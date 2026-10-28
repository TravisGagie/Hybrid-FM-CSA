#!/bin/bash
# Smallest versus most frequent successor as the continuation of a phrase, on the
# desktop, after run_final.sh (it reuses that script's collections and levels).
#
#   nohup bash run_continuation.sh > continuation.log 2>&1 &
#
#   1. block statistics (four choices of continuation) at each size;
#   2. at SIZE_Q (default n16): the levels with the most frequent successor, then
#      two-level search with block jumps on both kinds of levels (BML, L = 100 and
#      1000; forward-backward, L = 1000), checked against plain search.
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
H=${H:-$HOME/shortcuts/chr19_scaling/hprc_chr19}
export RS_DIR=${RS_DIR:-$HOME/shortcuts_build/reportable-shortcuts}
SIZES=${SIZES:-"n16 n32 n64"}
SIZE_Q=${SIZE_Q:-n16}
ROUNDS=${ROUNDS:-3}
FINAL=${FINAL:-$HOME/twolevel_runs/twolevel}     # where run_final.sh ran
OUT="$HERE/continuation"
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}

stamp "compiling"
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null
cmake --build "$HERE/build" -j4 --target build_twolevel query_twolevel > /dev/null
B="$HERE/build"

collection_for() {  # $1 size: the collection run_final.sh used
  local s=$1 n; n=$(od -An -t d8 -N8 "$H/$s/fwd_t3/rlbwt.bin" | tr -d ' ')
  if [ "$(( $(stat -c %s "$H/collection.txt") + 1 ))" = "$n" ]; then echo "$H/collection.txt"; else echo "$FINAL/data/$s.collection.txt"; fi
}

for S in $SIZES; do
  C=$(collection_for "$S")
  stamp "$S: block statistics ($C)"
  TIME=(); [ -x /usr/bin/time ] && TIME=(/usr/bin/time -v)
  ${TIME[@]+"${TIME[@]}"} "$B/build_twolevel" --text "$C" --block-stats 2> "$OUT/stats_$S.log"
  grep -E "block stats|successor|Maximum resident" "$OUT/stats_$S.log" | sed "s/^/$S: /"
done

S=$SIZE_Q; D="$H/$S"; C=$(collection_for "$S")
LV="$OUT/$S.freq"
if [ ! -s "$LV.text.lvl" ]; then
  stamp "$S: levels with the most frequent successor"
  "$B/build_twolevel" --text "$C" --rlbwt "$D/fwd_t3/rlbwt.bin" --reversed-rlbwt "$D/rev_t3/rlbwt.bin" --out "$LV" \
    --continuation frequent 2> "$OUT/build_freq_$S.log"
fi
grep -E "^blocks|sizes" "$OUT/build_freq_$S.log" || true
Q="$FINAL/final/queries.txt"
for kind in smallest frequent; do
  if [ $kind = smallest ]; then L="$FINAL/results/$S/v5"; else L="$LV"; fi
  for alg in "bml 100,1000" "fb 1000"; do
    set -- $alg
    stamp "$S, $kind: $1"
    taskset -c "$PIN" "$B/query_twolevel" --index "$D/fwd_t3" --reversed-index "$D/rev_t3" \
      --level "$L.text.lvl" --reversed-level "$L.rev.lvl" --patterns "$Q" --configs plain,2L,2LJ \
      --algorithm "$1" --min-length "$2" --repeat "$ROUNDS" > "$OUT/q_${kind}_$1.txt" 2> "$OUT/q_${kind}_$1.err"
    sed "s/^/$kind: /" "$OUT/q_${kind}_$1.txt"
    grep -h "extra space" "$OUT/q_${kind}_$1.err" | sed "s/^/$kind: /" || true
    echo "$kind $1: mismatches against plain search: $(grep -c MISMATCH "$OUT/q_${kind}_$1.err" || true)"
  done
done
stamp "done"
