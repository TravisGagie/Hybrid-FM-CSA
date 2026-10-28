#!/bin/bash
# Only ropebwt3 (alone, both strands, 1 thread), at the sizes in SIZES (default n32 n64),
# after run_final.sh, which ran it only at 16 individuals.
#
#   nohup bash run_ropebwt3_final.sh > ropebwt3_final.log 2>&1 &
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
H=${H:-$HOME/shortcuts/chr19_scaling/hprc_chr19}
export RS_DIR=${RS_DIR:-$HOME/shortcuts_build/reportable-shortcuts}
SIZES=${SIZES:-"n32 n64"}
ROUNDS=${ROUNDS:-5}
OUT="$HERE/final"
stamp() { echo "[$(date '+%F %T')] $*"; }
SUM="$OUT/summary_ropebwt3.txt"
: > "$SUM"
note() { echo "$*" | tee -a "$SUM"; }
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}
B="$HERE/build"
[ -x "$B/count_strings" ] || { cmake -S "$HERE" -B "$B" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null; cmake --build "$B" -j4 --target count_strings > /dev/null; }
RB_DIR="$HERE/../ropebwt3"
[ -x "$RB_DIR/ropebwt3" ] || (cd "$RB_DIR" && make -j4 > /dev/null 2>&1)
RB="$RB_DIR/ropebwt3"
Q="$OUT/queries.txt"
[ -s "$Q" ] || { echo "run run_final.sh first"; exit 1; }
collection_for() {
  local n; n=$(od -An -t d8 -N8 "$1/fwd_t3/rlbwt.bin" | tr -d ' ')
  local c="$H/collection.txt"
  if [ "$(( $(stat -c %s "$c") + 1 ))" != "$n" ]; then c="$HERE/data/$(basename "$1").collection.txt"; fi
  echo "$c"
}
for SIZE in $SIZES; do
  D="$H/$SIZE"; R="$OUT/$SIZE"; COLL=$(collection_for "$D")
  note ""; note "== $SIZE"
  # ---- ropebwt3, alone and behind KeBaB
  RBD="$R/ropebwt3"; mkdir -p "$RBD"
  RBI="$HERE/results/$SIZE/ropebwt3/index.fmd"
  if [ ! -s "$RBI" ]; then
    stamp "$SIZE: ropebwt3 index"; mkdir -p "$(dirname "$RBI")"
    TIME=(); [ -x /usr/bin/time ] && TIME=(/usr/bin/time -v)
    ${TIME[@]+"${TIME[@]}"} "$RB" build -t 8 -m 1G -L -do "$RBI" "$COLL" 2> "$RBD/build.log" || { note "$SIZE: ropebwt3 index build failed (see $RBD/build.log)"; continue; }
    grep -h "Maximum resident\|Elapsed" "$RBD/build.log" | sed "s/^/  build: /" | tee -a "$SUM" || true
  fi
  note "ropebwt3 index (both strands): $(( $(stat -c %s "$RBI") / 1048576 )) MiB"
  qtime() { awk -F'::' '/loaded the BWT/ {split($3,a,"*"); l=a[1]} /processed/ {split($3,b,"*"); p=b[1]} END {printf "%.3f\n", p-l}' "$1"; }
  rbrun() {  # $1 name, $2 input fasta, $3.. options
    local name=$1 in=$2; shift 2
    local ts=()
    for r in $(seq 1 "$ROUNDS"); do
      local o=/dev/null; [ "$r" = 1 ] && o="$RBD/mems_$name.txt"
      taskset -c "$PIN" "$RB" mem -t 1 "$@" "$RBI" "$in" > "$o" 2> "$RBD/log_$name.txt"
      ts+=("$(qtime "$RBD/log_$name.txt")")
    done
    printf '%s\n' "${ts[@]}" | sort -g | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2}'
  }
  stamp "$SIZE: ropebwt3"
  note "ropebwt3 BML, L = 100: $(rbrun bml100 "$OUT/queries.fa" -l 100) s"
  note "ropebwt3 BML, L = 1000: $(rbrun bml1000 "$OUT/queries.fa" -l 1000) s"
  note "ropebwt3 forward-backward (--old-mem), L = 1000: $(rbrun fb1000 "$OUT/queries.fa" --old-mem -l 1000) s"
  note "ropebwt3 all SMEMs (-l 1): $(rbrun fball "$OUT/queries.fa" -l 1) s"
  ln -sfn "$R/ours" "$RBD/ours"
  note "check, ropebwt3 against ours (differences must be reverse-strand matches):"
  python3 "$HERE/check_strands.py" "$RBD" "$Q" "$D/fwd_t3/rlbwt.bin" "$B/count_strings" | tee -a "$SUM"

done
stamp "done"
