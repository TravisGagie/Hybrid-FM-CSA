#!/bin/bash
# Movi 2 (Zakeri et al.) MEM finding on the same collection and queries, for comparison.
#
#   bash run_movi.sh                       # quick: regular index, one round
#   MODE=full bash run_movi.sh             # regular and blocked indexes, 5 rounds
#   DATA=.../n32 bash run_movi.sh          # another size (default n16)
#
# Movi 2 indexes both strands (like ropebwt3) and finds MEMs with BML (--min-mem-length L > 1)
# or all MEMs (L = 1, forward-backward style).  NB: use --min-mem-length; the short option -l
# is ignored.  MEM finding needs an ftab of k-mers (FTAB_K, default 10: 32 bytes per entry).
# Query time is what Movi 2 reports ("Time measured for processing the reads"), excluding loading.
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=${DATA:-$HOME/shortcuts/chr19_scaling/hprc_chr19/n16}
MODE=${MODE:-quick}
FTAB_K=${FTAB_K:-10}
SIZE=$(basename "$DATA")
OUT="$HERE/results/$SIZE/movi"
QUER=${QUER:-$DATA/../queries.txt}
COLL="$HERE/data/$SIZE.collection.txt"
[ -s "$COLL" ] || COLL="$DATA/../collection.txt"
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"

# ---- Movi 2 itself (fixed version 9bc0acc = 2.0.0); its build fetches sdsl-lite, pfp-thresholds, r-permute
MV_DIR="$HERE/../Movi"
if [ ! -x "$MV_DIR/build/movi" ]; then
  stamp "fetching and compiling Movi 2 (about 10 minutes)"
  [ -d "$MV_DIR" ] || git clone -q https://github.com/mohsenzakeri/Movi "$MV_DIR"
  (cd "$MV_DIR" && git checkout -q 9bc0acc && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 && make -j8 > make.log 2>&1) \
    || { echo "error: Movi 2 did not compile; see $MV_DIR/build/make.log"; exit 1; }
fi
MV="$MV_DIR/build/movi"
echo "collection: $COLL"; echo "queries: $QUER"
awk '{print ">q"NR-1; print}' "$QUER" > "$OUT/queries.fa"

TYPES=(regular); ROUNDS=1
[ "$MODE" = full ] && { TYPES=(regular blocked); ROUNDS=5; }
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}

for ty in "${TYPES[@]}"; do
  IDX="$OUT/index_$ty"
  if [ ! -s "$IDX/index.movi" ]; then
    [ -s "$OUT/collection.fa" ] || awk '{print ">s"NR; print}' "$COLL" > "$OUT/collection.fa"
    stamp "building the Movi 2 index ($ty, both strands)"
    TIME=(); [ -x /usr/bin/time ] && TIME=(/usr/bin/time -v)
    ${TIME[@]+"${TIME[@]}"} "$MV" build -f "$OUT/collection.fa" -i "$IDX" --type "$ty" > "$OUT/build_$ty.log" 2>&1 \
      || { tail -20 "$OUT/build_$ty.log"; exit 1; }
    grep -E "Elapsed|Maximum resident" "$OUT/build_$ty.log" || true
  fi
  if [ ! -s "$IDX/ftab.$FTAB_K.bin" ]; then
    "$MV" ftab -i "$IDX" --ftab-k "$FTAB_K" > "$OUT/ftab_$ty.log" 2>&1
  fi
  echo "index ($ty): $(( $(stat -c %s "$IDX/index.movi") / 1048576 )) MiB, ftab (k = $FTAB_K): $(( $(stat -c %s "$IDX/ftab.$FTAB_K.bin") / 1048576 )) MiB"
done
[ -s "$OUT/collection.fa" ] && rm -f "$OUT/collection.fa"

CSV="$OUT/times_$MODE.csv"
echo "type,L,round,query_s" > "$CSV"
for r in $(seq 1 "$ROUNDS"); do
  for ty in "${TYPES[@]}"; do
    for L in 100 1000 1; do
      out=/dev/null; [ "$r" = 1 ] && out="$OUT/mems_${ty}_L$L.txt"
      taskset -c "$PIN" "$MV" query -i "$OUT/index_$ty" -r "$OUT/queries.fa" --mem --min-mem-length "$L" --ftab-k "$FTAB_K" \
        --stdout -n -t 1 > "$out" 2> "$OUT/log_${ty}_L$L.txt"
      q=$(sed 's/\x1b\[[0-9;]*m//g' "$OUT/log_${ty}_L$L.txt" | grep -o 'processing the reads: [0-9.]*' | awk '{print $4}')
      echo "$ty,$L,$r,$q" >> "$CSV"
      stamp "round $r, $ty, L = $L: $q s"
    done
  done
done
echo; echo "medians (s):"
python3 - "$CSV" <<'PY'
import csv, sys, statistics, collections
d = collections.defaultdict(list)
for row in csv.DictReader(open(sys.argv[1])): d[(row["type"], int(row["L"]))].append(float(row["query_s"]))
for (t, L), v in sorted(d.items()): print(f"  {t:8s} L = {L:4d}: {statistics.median(v):9.2f}  (min {min(v):.2f}, max {max(v):.2f}, {len(v)} rounds)")
PY

# ---- compare with ropebwt3 (both strands, so they should agree)
RB="$HERE/results/$SIZE/ropebwt3"
for ty in "${TYPES[@]}"; do for x in "100 bml100" "1000 bml1000" "1000 fball"; do set -- $x
  m="$OUT/mems_${ty}_L$1.txt"; [ "$2" = fball ] && m="$OUT/mems_${ty}_L1.txt"
  [ -s "$m" ] && [ -s "$RB/mems_$2.txt" ] || continue
  awk -v L="$1" '$3-$2>=L {print $1"\t"$2"\t"$3}' "$m" | sort > "$OUT/a.txt"
  awk -v L="$1" '$3-$2>=L {print $1"\t"$2"\t"$3}' "$RB/mems_$2.txt" | sort > "$OUT/b.txt"
  echo "  Movi 2 ($ty, $(basename "$m")) vs ropebwt3 ($2), MEMs of length >= $1: $(wc -l < "$OUT/a.txt") vs $(wc -l < "$OUT/b.txt"), common $(comm -12 "$OUT/a.txt" "$OUT/b.txt" | wc -l)"
done; done
stamp "done"
