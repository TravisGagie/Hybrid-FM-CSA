#!/bin/bash
# ropebwt3 (Li) on the same collection and queries as run_desktop.sh, for comparison.
#
#   bash run_ropebwt3.sh                 # quick: one round, 1 thread
#   MODE=full bash run_ropebwt3.sh       # 5 rounds, 1 and 2 threads
#   DATA=.../n32 bash run_ropebwt3.sh    # another size (default n16)
#
# ropebwt3's mem command finds SMEMs (= our MEMs) against both strands of the
# collection (it needs both).  Its default algorithm is BML (-l is the minimum
# length); --old-mem is Li's original forward-backward SMEM algorithm, and with
# -l 1 the default algorithm also reduces to forward-backward and reports every
# SMEM.  Query time = time from "loaded the BWT" to "processed", as ropebwt3
# logs them.  Run run_desktop.sh for the same size first (it cuts the collection).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=${DATA:-$HOME/shortcuts/chr19_scaling/hprc_chr19/n16}
MODE=${MODE:-quick}
SIZE=$(basename "$DATA")
OUT="$HERE/results/$SIZE/ropebwt3"
QUER=${QUER:-$DATA/../queries.txt}
COLL="$HERE/data/$SIZE.collection.txt"
[ -s "$COLL" ] || COLL="$DATA/../collection.txt"
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"

# ---- ropebwt3 itself (fixed version r285)
RB_DIR="$HERE/../ropebwt3"
if [ ! -x "$RB_DIR/ropebwt3" ]; then
  stamp "fetching and compiling ropebwt3"
  [ -d "$RB_DIR" ] || git clone -q https://github.com/lh3/ropebwt3 "$RB_DIR"
  (cd "$RB_DIR" && { [ ! -d .git ] || git checkout -q b408339; } && make -j4 > /dev/null 2>&1)
fi
RB="$RB_DIR/ropebwt3"
echo "ropebwt3: $("$RB" version 2>/dev/null || echo '?') ($RB)"
echo "collection: $COLL"
echo "queries: $QUER"

# ropebwt3 maps any other character to N, which breaks strand symmetry; check
BAD=$(LC_ALL=C tr -d 'ACGTN\n' < "$COLL" | wc -c)
[ "$BAD" = 0 ] || { echo "error: the collection has $BAD characters other than A, C, G, T, N"; exit 1; }

# ---- the index (both strands), with batches small enough for this machine
if [ ! -s "$OUT/index.fmd" ]; then
  stamp "building the ropebwt3 index (both strands)"
  "$RB" build -t "${BUILD_THREADS:-8}" -m "${BATCH:-1G}" -L -do "$OUT/index.fmd" "$COLL" 2> "$OUT/build.log"
  grep -E "Real time" "$OUT/build.log" || true
fi
echo "index: $(stat -c %s "$OUT/index.fmd") bytes = $(( $(stat -c %s "$OUT/index.fmd") / 1048576 )) MiB"
awk '{print ">q"NR-1; print}' "$QUER" > "$OUT/queries.fa"

PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}
PIN2=${PIN2:-$(lscpu -p=CPU,CORE 2>/dev/null | grep -v '^#' | awk -F, -v p="$PIN" '$1==p{c=$2} {cpu[$1]=$2} END{for (k in cpu) if (cpu[k]!=c) {print k; exit}}')}
CONFIGS=("bml100:-l 100" "bml1000:-l 1000" "fb1000:--old-mem -l 1000" "fball:-l 1")
THREADS=(1); ROUNDS=1
[ "$MODE" = full ] && { THREADS=(1 2); ROUNDS=5; }
echo "pinned to CPU $PIN (and $PIN2 for 2 threads); rounds: $ROUNDS"

# query time from ropebwt3's log: realtime at "processed" minus realtime at "loaded the BWT"
qtime() { awk -F'::' '/loaded the BWT/ {split($3,a,"*"); l=a[1]} /processed/ {split($3,b,"*"); p=b[1]} END {printf "%.3f\n", p-l}' "$1"; }
CSV="$OUT/times_$MODE.csv"
echo "config,threads,round,query_s,peak_rss" > "$CSV"
for r in $(seq 1 "$ROUNDS"); do
  for t in "${THREADS[@]}"; do
    CPUS=$PIN; [ "$t" = 2 ] && CPUS="$PIN,$PIN2"
    for c in "${CONFIGS[@]}"; do
      name=${c%%:*}; opts=${c#*:}
      out=/dev/null; [ "$r" = 1 ] && [ "$t" = 1 ] && out="$OUT/mems_$name.txt"
      taskset -c "$CPUS" "$RB" mem -t "$t" $opts "$OUT/index.fmd" "$OUT/queries.fa" > "$out" 2> "$OUT/log_${name}_t$t.txt"
      q=$(qtime "$OUT/log_${name}_t$t.txt"); rss=$(grep -o 'Peak RSS: [0-9.]* GB' "$OUT/log_${name}_t$t.txt" | awk '{print $3}')
      echo "$name,$t,$r,$q,$rss" >> "$CSV"
      stamp "round $r, $t thread(s), $name: $q s (peak RSS $rss GB)"
    done
  done
done
echo
echo "medians (s):"
python3 - "$CSV" <<'PY'
import csv, sys, statistics, collections
d = collections.defaultdict(list)
for row in csv.DictReader(open(sys.argv[1])): d[(row["config"], row["threads"])].append(float(row["query_s"]))
for (c, t), v in sorted(d.items()): print(f"  {c:8s} {t} thread(s): {statistics.median(v):9.2f}  (min {min(v):.2f}, max {max(v):.2f}, {len(v)} rounds)")
PY

# ---- exactness: compare with the MEMs of the two-level index (identical to plain search)
B="$HERE/build"; LV="$HERE/results/$SIZE"
if [ -x "$B/query_twolevel" ] && [ -s "$LV/v4.text.lvl" ]; then
  stamp "comparing with our MEMs"
  mkdir -p "$OUT/ours"
  "$B/query_twolevel" --index "$DATA/fwd_t3" --reversed-index "$DATA/rev_t3" --level "$LV/v4.text.lvl" --reversed-level "$LV/v4.rev.lvl" \
    --patterns "$QUER" --configs 2LJ --algorithm bml,fb --min-length 100,1000 --repeat 1 --dump "$OUT/ours" > /dev/null 2>&1
  cmpm() {  # $1: our files prefix, $2: ropebwt3 output, $3: min length
    (awk '{print "q0\t"$1"\t"$2}' "$OUT/ours/$1_q0.tsv"; awk '{print "q1\t"$1"\t"$2}' "$OUT/ours/$1_q1.tsv") | sort > "$OUT/a.txt"
    awk -v L="$3" '$3-$2>=L {print $1"\t"$2"\t"$3}' "$2" | sort > "$OUT/b.txt"
    echo "  $1 vs $(basename "$2"): ours $(wc -l < "$OUT/a.txt"), ropebwt3 $(wc -l < "$OUT/b.txt"), common $(comm -12 "$OUT/a.txt" "$OUT/b.txt" | wc -l)"
  }
  cmpm bml_L100 "$OUT/mems_bml100.txt" 100
  cmpm bml_L1000 "$OUT/mems_bml1000.txt" 1000
  cmpm fb_L1000 "$OUT/mems_fb1000.txt" 1000
  cmpm fb_L1000 "$OUT/mems_fball.txt" 1000
  echo "  (ropebwt3 also matches the reverse-complement strand, so it can find MEMs we do not, or longer ones)"
fi
stamp "done"
