#!/bin/bash
# Two-level indexing against plain search and the shortcuts, on one collection of
# the chr19 scaling run (default: 16 individuals).  Self-contained: run it from the
# folder it was unpacked into.
#
#   MODE=probe bash run_desktop.sh               # show what data was found (no building or timing)
#   bash run_desktop.sh                          # build, self-test, build the levels, quick round (BML, L = 1000)
#   MODE=full bash run_desktop.sh                # all shortcut configurations found, BML L = 100 and 1000,
#                                                # forward-backward L = 1000, 5 rounds each (hours)
#   MODE=lean bash run_desktop.sh                # as full, but only plain, the fastest shortcuts and two-level
#   MODE=twolevel bash run_desktop.sh            # only the two two-level versions (minutes)
#   DATA=.../n32 MODE=full bash run_desktop.sh   # another collection size
#
# Nothing under ~/shortcuts is changed; levels and results go to ./results/<size>.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=${DATA:-$HOME/shortcuts/chr19_scaling/hprc_chr19/n16}
MODE=${MODE:-quick}
RS_DIR=${RS_DIR:-$HOME/shortcuts/reportable-shortcuts}
stamp() { echo "[$(date '+%F %T')] $*"; }

# ---- find the data
find_one() {  # first existing file among the arguments
  for f in "$@"; do [ -s "$f" ] && { echo "$f"; return; }; done; echo ""
}
COLL=${COLL:-$(find_one "$DATA/collection.txt" "$DATA/../collection.txt" $(ls "$DATA"/*collection*.txt 2>/dev/null | grep -v rev || true))}
QUER=${QUER:-$(find_one "$DATA/queries.txt" "$DATA/../queries.txt" "$HOME/shortcuts/chr19/queries.txt")}
FWD=${FWD:-$DATA/fwd_t3}
REV=${REV:-$DATA/rev_t3}
SCS=""; BEST=""
for t in t3 t6; do for e in 16 25 50; do
  f="$DATA/fwd_$t/ell${e}_theta0.shortcuts"; r="$DATA/rev_$t/ell${e}_theta0.shortcuts"
  if [ -s "$f" ] && [ -s "$r" ]; then
    SCS="$SCS,ell${e}_theta0@$DATA/fwd_$t@$DATA/rev_$t"
    if [ "$t$e" = t625 ]; then BEST="ell25_theta0@$DATA/fwd_t6@$DATA/rev_t6"; fi
  fi
done; done
SCS=${SCS#,}
[ -z "$BEST" ] && BEST=${SCS##*,}

echo "machine:     $(hostname), $(nproc) CPUs"
echo "data:        $DATA"
echo "collection:  ${COLL:-NOT FOUND}"
echo "queries:     ${QUER:-NOT FOUND}"
echo "RLBWTs:      $FWD/rlbwt.bin $( [ -s "$FWD/rlbwt.bin" ] && echo ok || echo MISSING ), $REV/rlbwt.bin $( [ -s "$REV/rlbwt.bin" ] && echo ok || echo MISSING )"
echo "shortcuts:   ${SCS:-NONE FOUND}"
echo "source:      $RS_DIR $( [ -s "$RS_DIR/include/rlbwt.hpp" ] && echo ok || echo MISSING )"
if [ "$MODE" = probe ] || [ -z "$COLL" ] || [ -z "$QUER" ] || [ ! -s "$FWD/rlbwt.bin" ] || [ ! -s "$REV/rlbwt.bin" ] || [ ! -s "$RS_DIR/include/rlbwt.hpp" ]; then
  echo; echo "contents of $DATA:"; find "$DATA" -maxdepth 2 -printf '%10s  %p\n' 2>/dev/null | head -60
  [ "$MODE" = probe ] && exit 0
  echo "error: something is missing (set COLL, QUER, FWD, REV or RS_DIR)"; exit 1
fi

# The collections of the scaling run are prefixes of one file: the collection of
# this size is the first n - 1 bytes, where n (the first 8 bytes of rlbwt.bin) is
# the length of T = S_1 # ... # S_k # $.
NROWS=$(od -An -t d8 -N8 "$FWD/rlbwt.bin" | tr -d ' ')
CSIZE=$(stat -c %s "$COLL")
if [ "$((CSIZE + 1))" != "$NROWS" ]; then
  SUB="$HERE/data/$(basename "$DATA").collection.txt"
  if [ ! -s "$SUB" ] || [ "$(( $(stat -c %s "$SUB") + 1 ))" != "$NROWS" ]; then
    [ "$CSIZE" -ge "$NROWS" ] || { echo "error: $COLL is shorter than the RLBWT ($NROWS rows)"; exit 1; }
    mkdir -p "$HERE/data"
    stamp "writing the first $((NROWS - 1)) bytes of $COLL to $SUB"
    head -c "$((NROWS - 1))" "$COLL" > "$SUB"
  fi
  [ "$(tail -c 1 "$SUB" | od -An -c | tr -d ' ')" = '\n' ] || { echo "error: $SUB does not end at a line break"; exit 1; }
  COLL=$SUB
fi
echo "collection used: $COLL ($(( NROWS - 1 )) bytes, $(wc -l < "$COLL") sequences)"

PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}
OUT=${OUT:-$HERE/results/$(basename "$DATA")}
mkdir -p "$OUT"

stamp "compiling"
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null
cmake --build "$HERE/build" -j4 --target build_twolevel query_twolevel build_index > /dev/null
B="$HERE/build"

stamp "self-test on a small synthetic collection"
T="$HERE/selftest"
rm -rf "$T"; mkdir -p "$T"
(cd "$T" && python3 "$HERE/gen_test.py" 3 200000 0.0005 > /dev/null)
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idx" --ell 25 --theta 0 > /dev/null 2>&1
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idxr" --ell 25 --theta 0 --reversed > /dev/null 2>&1
"$B/build_twolevel" --text "$T/collection.txt" --rlbwt "$T/idx/rlbwt.bin" --reversed-rlbwt "$T/idxr/rlbwt.bin" --out "$T/lv" 2> "$T/build.log"
"$B/query_twolevel" --index "$T/idx" --reversed-index "$T/idxr" --level "$T/lv.text.lvl" --reversed-level "$T/lv.rev.lvl" \
  --patterns "$T/queries.txt" --shortcuts ell25_theta0 --min-length 1,100,1000 --repeat 1 > "$T/query.txt" 2> "$T/query.log" \
  || { echo "self-test FAILED (see $T/query.log)"; exit 1; }
stamp "self-test passed: the same MEMs and intervals as plain search"

TIME=(); if [ -x /usr/bin/time ]; then TIME=(/usr/bin/time -v); fi
if [ ! -s "$OUT/v5.text.lvl" ] || [ ! -s "$OUT/v5.rev.lvl" ]; then
  stamp "building the levels of T and of the reversed collection"
  ${TIME[@]+"${TIME[@]}"} "$B/build_twolevel" --text "$COLL" --rlbwt "$FWD/rlbwt.bin" --reversed-rlbwt "$REV/rlbwt.bin" \
    --out "$OUT/v5" 2> "$OUT/build.log" || { grep -v '^\s' "$OUT/build.log" | tail -5; exit 1; }
  grep -E "^==|boundaries|distinct|B:|blocks|BWT|shared|sizes|Maximum resident|Elapsed" "$OUT/build.log" || true
fi

Q=(taskset -c "$PIN" "$B/query_twolevel" --index "$FWD" --reversed-index "$REV"
   --level "$OUT/v5.text.lvl" --reversed-level "$OUT/v5.rev.lvl" --patterns "$QUER")
stamp "queries pinned to CPU $PIN"
case "$MODE" in
  quick)
    stamp "quick run: BML, L = 1000, one round"
    "${Q[@]}" --shortcuts "$BEST" --algorithm bml --min-length 1000 --repeat 1 --csv "$OUT/quick.csv" 2>&1 | tee "$OUT/quick.txt" ;;
  full|lean)
    S=$SCS; [ "$MODE" = lean ] && S=$BEST
    stamp "$MODE run: BML, L = 100 and 1000"
    "${Q[@]}" --shortcuts "$S" --algorithm bml --min-length 100,1000 --repeat 5 --csv "$OUT/${MODE}_bml.csv" 2>&1 | tee "$OUT/${MODE}_bml.txt"
    stamp "$MODE run: forward-backward, L = 1000"
    "${Q[@]}" --shortcuts "$S" --algorithm fb --min-length 1000 --repeat 5 --csv "$OUT/${MODE}_fb.csv" 2>&1 | tee "$OUT/${MODE}_fb.txt" ;;
  twolevel)   # only the two-level index (5 rounds); plain search and the shortcuts are in the lean/full results
    stamp "two-level only: BML, L = 100 and 1000"
    "${Q[@]}" --configs 2L,2LJ --algorithm bml --min-length 100,1000 --repeat 5 --csv "$OUT/tl_bml.csv" 2>&1 | tee "$OUT/tl_bml.txt"
    stamp "two-level only: forward-backward, L = 1000"
    "${Q[@]}" --configs 2L,2LJ --algorithm fb --min-length 1000 --repeat 5 --csv "$OUT/tl_fb.csv" 2>&1 | tee "$OUT/tl_fb.txt" ;;
  *) echo "unknown MODE $MODE"; exit 1 ;;
esac
stamp "done"
