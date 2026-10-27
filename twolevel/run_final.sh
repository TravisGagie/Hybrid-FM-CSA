#!/bin/bash
# The final comparison for the paper, on the desktop.  MEMs are over {A,C,G,T}: the
# queries are split at their N runs (4 patterns), so every method answers the same
# question and the worst-case bound applies.
#
#   nohup bash run_final.sh > final.log 2>&1 &
#
# For each size in SIZES (default n16 n32 n64), with 5 rounds per configuration:
#   - plain search, the shortcuts (tau = 6, ell = 25), two-level without and with jumps
#     (BML with L = 100 and 1000, forward-backward with L = 1000);
# and, at 16 individuals only:
#   - KeBaB (k = 20) in front of plain search and two-level with jumps;
#   - ropebwt3, alone and behind KeBaB (both strands, 1 thread);
#   - Movi 2 (regular index, both strands, 1 thread), last, since its index build may not
#     fit in memory; if it fails, everything else is already done.
# Every method's MEMs are checked: two-level and shortcuts against plain search (inside
# query_twolevel), KeBaB pipelines against the full queries, ropebwt3 and Movi 2 against
# ours with the strand check.  Results go to ./final; a summary to final/summary.txt.
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
H=${H:-$HOME/shortcuts/chr19_scaling/hprc_chr19}
export RS_DIR=${RS_DIR:-$HOME/shortcuts_build/reportable-shortcuts}
SIZES=${SIZES:-"n16 n32 n64"}
ROUNDS=${ROUNDS:-5}
KMER=${KMER:-20}
OUT="$HERE/final"
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"
SUM="$OUT/summary.txt"
: > "$SUM"
note() { echo "$*" | tee -a "$SUM"; }
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}

# ---------------------------------------------------------------- tools
stamp "compiling"
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null
cmake --build "$HERE/build" -j4 --target build_twolevel query_twolevel count_strings build_index > /dev/null
B="$HERE/build"
RB_DIR="$HERE/../ropebwt3"
[ -x "$RB_DIR/ropebwt3" ] || (cd "$RB_DIR" && make -j4 > /dev/null 2>&1)
RB="$RB_DIR/ropebwt3"
KB_DIR="$HERE/../kebab"
if [ ! -x "$KB_DIR/kebab" ]; then
  stamp "fetching and compiling KeBaB"
  [ -d "$KB_DIR" ] || git clone -q https://github.com/drnatebrown/kebab "$KB_DIR"
  (cd "$KB_DIR" && git checkout -q 4bed76f && make -j4 > make.log 2>&1)
fi
KB="$KB_DIR/kebab"

# ---------------------------------------------------------------- queries over {A,C,G,T}
Q="$OUT/queries.txt"
tr -c 'ACGT\n' '\n' < "$H/queries.txt" | awk 'length($0) > 0' > "$Q"
awk '{print ">p"NR-1; print}' "$Q" > "$OUT/queries.fa"
note "queries: $(wc -l < "$Q") patterns over {A,C,G,T}, $(tr -d '\n' < "$Q" | wc -c) characters (from $H/queries.txt, split at N runs)"
note "machine: $(hostname), queries pinned to CPU $PIN, $ROUNDS rounds, medians"

# KeBaB fragments (FASTA, names pI:start-end, 1-based) -> one per line, plus offsets
frag2lines() {  # $1 fragments.fa, $2 out.txt, $3 out.off
  awk -v T="$2" -v O="$3" '/^>/ { split(substr($0, 2), a, ":"); split(a[2], b, "-"); p = substr(a[1], 2); s = b[1] - 1; next }
                           { print > T; print p "\t" s > O }' "$1"
}
# the MEMs of fragments, mapped back, against those of the full patterns: $1 frag dir, $2 offsets, $3 full dir, $4 name
cmpfrag() {
  python3 - "$1" "$2" "$3" "$4" "$(wc -l < "$Q")" <<'PY'
import os, sys
fd, offp, full, name, nq = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
offs = [l.split() for l in open(offp)]
got = set()
for i, (p, o) in enumerate(offs):
    f = os.path.join(fd, f"{name}_q{i}.tsv")
    for l in open(f):
        s, e, c = map(int, l.split()); got.add((int(p), s + int(o), e + int(o), c))
ref = set()
for q in range(nq):
    for l in open(os.path.join(full, f"{name}_q{q}.tsv")):
        s, e, c = map(int, l.split()); ref.add((q, s, e, c))
print(f"{len(got)} MEMs from the fragments, {len(ref)} from the full patterns: {'identical (positions and counts)' if got == ref else 'DIFFERENT'}")
PY
}
# collection of this size: the first n - 1 bytes of the 64-individual collection
collection_for() {  # $1 size dir
  local n; n=$(od -An -t d8 -N8 "$1/fwd_t3/rlbwt.bin" | tr -d ' ')
  local c="$H/collection.txt"
  if [ "$(( $(stat -c %s "$c") + 1 ))" != "$n" ]; then
    local sub="$HERE/data/$(basename "$1").collection.txt"
    if [ ! -s "$sub" ] || [ "$(( $(stat -c %s "$sub") + 1 ))" != "$n" ]; then
      mkdir -p "$HERE/data"; head -c "$((n - 1))" "$c" > "$sub"
    fi
    c=$sub
  fi
  echo "$c"
}

for SIZE in $SIZES; do
  D="$H/$SIZE"; R="$OUT/$SIZE"; mkdir -p "$R"
  stamp "===== $SIZE"
  COLL=$(collection_for "$D")
  LV="$HERE/results/$SIZE/v5"
  if [ ! -s "$LV.text.lvl" ] || [ ! -s "$LV.rev.lvl" ]; then
    stamp "building the two-level index"
    mkdir -p "$HERE/results/$SIZE"
    "$B/build_twolevel" --text "$COLL" --rlbwt "$D/fwd_t3/rlbwt.bin" --reversed-rlbwt "$D/rev_t3/rlbwt.bin" --out "$LV" 2> "$R/build.log"
  fi
  QT=(taskset -c "$PIN" "$B/query_twolevel" --index "$D/fwd_t3" --reversed-index "$D/rev_t3"
      --level "$LV.text.lvl" --reversed-level "$LV.rev.lvl")
  SC="ell25_theta0@$D/fwd_t6@$D/rev_t6"; SCN="ell25_theta0@$D/fwd_t6"

  # ---- plain, shortcuts, two-level
  mkdir -p "$R/ours"
  stamp "$SIZE: plain, shortcuts, two-level (BML)"
  "${QT[@]}" --patterns "$Q" --shortcuts "$SC" --configs "plain,$SCN,2L,2LJ" --algorithm bml --min-length 100,1000 \
    --repeat "$ROUNDS" --csv "$R/bml.csv" --dump "$R/ours" > "$R/bml.txt" 2> "$R/bml.err" || note "$SIZE BML: query_twolevel exited with $? (see $R/bml.err)"
  stamp "$SIZE: plain, shortcuts, two-level (forward-backward)"
  "${QT[@]}" --patterns "$Q" --shortcuts "$SC" --configs "plain,$SCN,2L,2LJ" --algorithm fb --min-length 1000 \
    --repeat "$ROUNDS" --csv "$R/fb.csv" --dump "$R/ours" > "$R/fb.txt" 2> "$R/fb.err" || note "$SIZE FB: query_twolevel exited with $? (see $R/fb.err)"
  note ""; note "== $SIZE ($(basename "$COLL"), $(wc -l < "$COLL") sequences)"
  grep -h "extra space" "$R/bml.err" | tee -a "$SUM" || true
  grep -h "MiB on disk" "$R/bml.err" | tee -a "$SUM" || true
  cat "$R/bml.txt" "$R/fb.txt" | tee -a "$SUM"
  note "mismatches against plain search: $(cat "$R/bml.err" "$R/fb.err" | grep -c MISMATCH || true)"

  [ "$SIZE" = "${EXTRA_SIZE:-n16}" ] || continue

  # ---- KeBaB in front of plain search and two-level with jumps
  K="$R/kebab"; mkdir -p "$K"
  if [ ! -s "$K/fwd.kbb" ]; then
    stamp "$SIZE: KeBaB index (forward strand, k = $KMER)"
    awk '{print ">s"NR; print}' "$COLL" > "$K/collection.fa"
    "$KB" build -k "$KMER" --kmer-mode forward -t 8 -o "$K/fwd" "$K/collection.fa" > "$K/build_fwd.log" 2>&1
    "$KB" build -k "$KMER" -t 8 -o "$K/can" "$K/collection.fa" > "$K/build_can.log" 2>&1   # canonical, for ropebwt3
    rm -f "$K/collection.fa"
  fi
  note "KeBaB filters (k = $KMER): forward $(( $(stat -c %s "$K/fwd.kbb") / 1048576 )) MiB, canonical $(( $(stat -c %s "$K/can.kbb") / 1048576 )) MiB"
  for L in 100 1000; do
    for idx in fwd can; do
      # scan time: median of ROUNDS, single thread
      ts=()
      for r in $(seq 1 "$ROUNDS"); do
        t0=$(date +%s.%N); taskset -c "$PIN" "$KB" scan -i "$K/$idx.kbb" -l "$L" -t 1 -o "$K/frag_${idx}_L$L.fa" "$OUT/queries.fa" > /dev/null 2>&1
        ts+=("$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN {printf "%.3f", b - a}')")
      done
      med=$(printf '%s\n' "${ts[@]}" | sort -g | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2}')
      note "KeBaB scan ($idx, L = $L): $med s; $(grep -c '>' "$K/frag_${idx}_L$L.fa") fragments, $(grep -v '>' "$K/frag_${idx}_L$L.fa" | tr -d '\n' | wc -c) characters kept"
    done
    frag2lines "$K/frag_fwd_L$L.fa" "$K/frag_L$L.txt" "$K/frag_L$L.off"
    mkdir -p "$K/ours_L$L"
    algs=bml; [ "$L" = 1000 ] && algs=bml,fb
    stamp "$SIZE: plain and two-level with jumps on KeBaB fragments, L = $L"
    "${QT[@]}" --patterns "$K/frag_L$L.txt" --configs plain,2LJ --algorithm "$algs" --min-length "$L" --repeat "$ROUNDS" \
      --dump "$K/ours_L$L" > "$K/kebab_L$L.txt" 2> "$K/kebab_L$L.err" || note "KeBaB L = $L: query_twolevel exited with $?"
    sed 's/^/KeBaB + /' "$K/kebab_L$L.txt" | tee -a "$SUM"
    for a in ${algs//,/ }; do note "  check, KeBaB + $a, L = $L: $(cmpfrag "$K/ours_L$L" "$K/frag_L$L.off" "$R/ours" "${a}_L$L")"; done
  done

  # ---- ropebwt3, alone and behind KeBaB
  RBD="$R/ropebwt3"; mkdir -p "$RBD"
  RBI="$HERE/results/$SIZE/ropebwt3/index.fmd"
  if [ ! -s "$RBI" ]; then
    stamp "$SIZE: ropebwt3 index"; mkdir -p "$(dirname "$RBI")"
    "$RB" build -t 8 -m 1G -L -do "$RBI" "$COLL" 2> "$RBD/build.log"
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
  for L in 100 1000; do
    note "KeBaB + ropebwt3 BML, L = $L: $(rbrun kebab_bml$L "$K/frag_can_L$L.fa" -l "$L") s (plus the canonical scan above)"
  done
  ln -sfn "$R/ours" "$RBD/ours"
  note "check, ropebwt3 against ours (differences must be reverse-strand matches):"
  python3 "$HERE/check_strands.py" "$RBD" "$Q" "$D/fwd_t3/rlbwt.bin" "$B/count_strings" | tee -a "$SUM"

  # ---- Movi 2 (last)
  MV_DIR="$HERE/../Movi"; MVD="$R/movi"; mkdir -p "$MVD"
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
done
stamp "done"
note ""; note "finished $(date '+%F %T')"
