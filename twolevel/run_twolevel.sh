#!/bin/bash
# Two-level indexing on chr19, against plain search and the shortcuts, on this laptop.
#
#   bash run_twolevel.sh            # build, self-test, build the levels, quick run (about 15 minutes)
#   MODE=full bash run_twolevel.sh  # the full comparison (several hours): all six shortcut
#                                   # configurations, BML with L = 100 and 1000 and forward-backward
#                                   # with L = 1000, 5 rounds each, pinned to CPU $PIN (default 3)
#
# Close other programs before the timed part.  Nothing in ~/shortcuts is changed;
# the levels and results go to ./chr19 (or $OUT).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SC="$HOME/shortcuts"
C="$SC/chr19"
OUT=${OUT:-$HERE/chr19}
PIN=${PIN:-3}
MODE=${MODE:-quick}
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"

stamp "compiling"
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$HERE/build" -j4 --target build_twolevel query_twolevel build_index > /dev/null
B="$HERE/build"

stamp "self-test on a small synthetic collection"
T="$HERE/selftest"
rm -rf "$T"; mkdir -p "$T"
(cd "$T" && python3 "$HERE/gen_test.py" 3 200000 0.0005 > /dev/null)
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idx" --ell 25 --theta 0 > /dev/null 2>&1
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idxr" --ell 25 --theta 0 --reversed > /dev/null 2>&1
"$B/build_twolevel" --text "$T/collection.txt" --rlbwt "$T/idx/rlbwt.bin" --out "$T/t.lvl" 2> "$T/build.log"
"$B/build_twolevel" --text "$T/collection.txt" --rlbwt "$T/idxr/rlbwt.bin" --out "$T/r.lvl" --reversed 2>> "$T/build.log"
"$B/query_twolevel" --index "$T/idx" --reversed-index "$T/idxr" --level "$T/t.lvl" --reversed-level "$T/r.lvl" \
  --patterns "$T/queries.txt" --shortcuts ell25_theta0 --min-length 1,100,1000 --repeat 1 > "$T/query.txt" 2> "$T/query.log" \
  || { echo "self-test FAILED (see $T/query.log)"; exit 1; }
stamp "self-test passed: the same MEMs and intervals as plain search"

TIME=(); if [ -x /usr/bin/time ]; then TIME=(/usr/bin/time -v); fi
if [ ! -s "$OUT/text_v3.lvl" ]; then
  stamp "building the levels of T"
  ${TIME[@]+"${TIME[@]}"} "$B/build_twolevel" --text "$C/collection.txt" --rlbwt "$C/packed_t3/rlbwt.bin" --out "$OUT/text_v3.lvl" 2> "$OUT/build_text.log"
  grep -E "boundaries|distinct|B:|blocks|BWT|sizes|Maximum resident|Elapsed" "$OUT/build_text.log" || true
fi
if [ ! -s "$OUT/rev_v3.lvl" ]; then
  stamp "building the levels of the reversed collection"
  ${TIME[@]+"${TIME[@]}"} "$B/build_twolevel" --text "$C/collection.txt" --rlbwt "$C/packed_rev_t3/rlbwt.bin" --out "$OUT/rev_v3.lvl" --reversed 2> "$OUT/build_rev.log"
  grep -E "boundaries|distinct|B:|blocks|BWT|sizes|Maximum resident|Elapsed" "$OUT/build_rev.log" || true
fi

Q=(taskset -c "$PIN" "$B/query_twolevel" --index "$C/packed_t3" --reversed-index "$C/packed_rev_t3"
   --level "$OUT/text_v3.lvl" --reversed-level "$OUT/rev_v3.lvl" --patterns "$C/queries.txt")
if [ "$MODE" = quick ]; then
  stamp "quick run: BML, L = 1000, one round, best shortcuts (t6, ell 25)"
  "${Q[@]}" --shortcuts "ell25_theta0@$C/packed_t6@$C/packed_rev_t6" --algorithm bml --min-length 1000 --repeat 1 \
    --csv "$OUT/quick.csv" 2>&1 | tee "$OUT/quick.txt"
else
  SCS=""
  for t in t3 t6; do for e in 16 25 50; do
    SCS="$SCS,ell${e}_theta0@$C/packed_$t@$C/packed_rev_$t"
  done; done
  SCS=${SCS#,}
  stamp "full run: BML, L = 100 and 1000"
  "${Q[@]}" --shortcuts "$SCS" --algorithm bml --min-length 100,1000 --repeat 5 --csv "$OUT/full_bml.csv" 2>&1 | tee "$OUT/full_bml.txt"
  stamp "full run: forward-backward, L = 1000"
  "${Q[@]}" --shortcuts "$SCS" --algorithm fb --min-length 1000 --repeat 5 --csv "$OUT/full_fb.csv" 2>&1 | tee "$OUT/full_fb.txt"
fi
stamp "done"
