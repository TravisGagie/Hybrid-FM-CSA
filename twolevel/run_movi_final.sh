#!/bin/bash
# Only the Movi 2 part of run_final.sh (at 16 individuals), for when the Movi clone
# failed during the final run.  Needs run_final.sh's results in ./final.
#
#   nohup bash run_movi_final.sh > movi_final.log 2>&1 &
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
H=${H:-$HOME/shortcuts/chr19_scaling/hprc_chr19}
export RS_DIR=${RS_DIR:-$HOME/shortcuts_build/reportable-shortcuts}
ROUNDS=${ROUNDS:-5}
OUT="$HERE/final"
stamp() { echo "[$(date '+%F %T')] $*"; }
SUM="$OUT/summary_movi.txt"
: > "$SUM"
note() { echo "$*" | tee -a "$SUM"; }
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}
B="$HERE/build"
[ -x "$B/count_strings" ] || { cmake -S "$HERE" -B "$B" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null; cmake --build "$B" -j4 --target count_strings > /dev/null; }
Q="$OUT/queries.txt"
SIZE=${SIZE:-n16}
D="$H/$SIZE"; R="$OUT/$SIZE"
[ -s "$Q" ] && [ -d "$R/ours" ] || { echo "run run_final.sh first"; exit 1; }
collection_for() {
  local n; n=$(od -An -t d8 -N8 "$1/fwd_t3/rlbwt.bin" | tr -d ' ')
  local c="$H/collection.txt"
  if [ "$(( $(stat -c %s "$c") + 1 ))" != "$n" ]; then c="$HERE/data/$(basename "$1").collection.txt"; fi
  echo "$c"
}
COLL=$(collection_for "$D")
  # ---- Movi 2 (last)
MV_DIR="$HERE/../Movi"; MVD="$R/movi"; mkdir -p "$MVD"
[ -d "$MV_DIR" ] && [ ! -f "$MV_DIR/CMakeLists.txt" ] && rm -rf "$MV_DIR"   # left over from a failed clone
if [ ! -x "$MV_DIR/build/movi" ]; then
  stamp "fetching and compiling Movi 2"
  { [ -d "$MV_DIR" ] || git clone -q https://github.com/mohsenzakeri/Movi "$MV_DIR"; } && \
  (cd "$MV_DIR" && git checkout -q 9bc0acc && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 && make -j8 > make.log 2>&1) \
    || note "Movi 2 did not compile (see $MV_DIR/build/make.log)"
fi
MV="$MV_DIR/build/movi"; MVI="$HERE/results/$SIZE/movi_regular"
if [ -x "$MV" ] && [ ! -s "$MVI/index.movi" ]; then
  stamp "$SIZE: Movi 2 index (regular, both strands)"
  awk '{print ">s"NR; print}' "$COLL" > "$MVD/collection.fa"
  "$MV" build -f "$MVD/collection.fa" -i "$MVI" --type regular > "$MVD/build.log" 2>&1 || note "Movi 2 index build failed (see $MVD/build.log)"
  rm -f "$MVD/collection.fa"
fi
if [ -s "$MVI/index.movi" ]; then
  [ -s "$MVI/ftab.10.bin" ] || "$MV" ftab -i "$MVI" --ftab-k 10 > "$MVD/ftab.log" 2>&1
  note "Movi 2 index (regular, both strands): $(( $(stat -c %s "$MVI/index.movi") / 1048576 )) MiB, plus ftab (k = 10) $(( $(stat -c %s "$MVI/ftab.10.bin") / 1048576 )) MiB"
  for L in 100 1000 1; do
    ts=(); nr=$ROUNDS; [ "$L" = 1 ] && nr=1   # all MEMs (L = 1) can be much slower in Movi 2: one round
    for r in $(seq 1 "$nr"); do
      o=/dev/null; [ "$r" = 1 ] && o="$MVD/mems_movi$L.txt"
      taskset -c "$PIN" "$MV" query -i "$MVI" -r "$OUT/queries.fa" --mem --min-mem-length "$L" --ftab-k 10 --stdout -n -t 1 > "$o" 2> "$MVD/log_L$L.txt"
      ts+=("$(sed 's/\x1b\[[0-9;]*m//g' "$MVD/log_L$L.txt" | grep -o 'processing the reads: [0-9.]*' | awk '{print $4}')")
    done
    note "Movi 2 MEMs, L = $L: $(printf '%s\n' "${ts[@]}" | sort -g | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2}') s"
  done
  ln -sfn "$R/ours" "$MVD/ours"
  note "check, Movi 2 against ours (differences must be reverse-strand matches, or N in the collection read as A):"
  python3 "$HERE/check_strands.py" "$MVD" "$Q" "$D/fwd_t3/rlbwt.bin" "$B/count_strings" \
    "bml_L100:movi100:100,bml_L1000:movi1000:1000,fb_L1000:movi1:1000" | tee -a "$SUM"
fi
stamp "done"
