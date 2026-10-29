#!/bin/bash
# Simulated long reads (10 kb, 1% and 5% errors, forward strand) from the two query
# haplotypes, after run_final.sh and run_extras_final.sh: plain search, shortcuts,
# two-level indexing with and without block jumps, KeBaB in front of plain search and
# block jumps, ropebwt3 alone and behind KeBaB, and Movi 2, at the sizes in SIZES.
# MEMs of length at least 30 and 100 (BML) and 30 (forward-backward).  No network needed.
#
#   nohup bash run_reads.sh > reads.log 2>&1 &
set -euo pipefail
export LC_ALL=C
HERE=$(cd "$(dirname "$0")" && pwd)
H=${H:-$HOME/shortcuts/chr19_scaling/hprc_chr19}
export RS_DIR=${RS_DIR:-$HOME/shortcuts_build/reportable-shortcuts}
SIZES=${SIZES:-"n16 n64"}
RATES=${RATES:-"1 5"}
NREADS=${NREADS:-2000}
RLEN=${RLEN:-10000}
ROUNDS=${ROUNDS:-5}
KMER=${KMER:-20}
LS=${LS:-"30 100"}
FINAL="$HERE/final"
OUT="$HERE/reads"
stamp() { echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"
SUM="$OUT/summary.txt"
: > "$SUM"
note() { echo "$*" | tee -a "$SUM"; }
PIN=${PIN:-$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)}
PIN=${PIN:-0}
B="$HERE/build"
[ -x "$B/query_twolevel" ] && [ -x "$B/count_strings" ] || { cmake -S "$HERE" -B "$B" -DCMAKE_BUILD_TYPE=Release -DRS_DIR="$RS_DIR" > /dev/null; cmake --build "$B" -j4 --target query_twolevel count_strings > /dev/null; }
RB="$HERE/../ropebwt3/ropebwt3"; KB="$HERE/../kebab/kebab"; MV="$HERE/../Movi/build/movi"
[ -s "$FINAL/queries.txt" ] || { echo "run run_final.sh first"; exit 1; }
# KeBaB fragments (FASTA, names pI:start-end, 1-based) -> one per line, plus offsets
frag2lines() {  # $1 fragments.fa, $2 out.txt, $3 out.off
  awk -v T="$2" -v O="$3" '/^>/ { split(substr($0, 2), a, ":"); split(a[2], b, "-"); p = substr(a[1], 2); s = b[1] - 1; next }
                           { print > T; print p "\t" s > O }' "$1"
}
# the MEMs of fragments, mapped back, against those of the full patterns: $1 frag dir, $2 offsets, $3 full dir, $4 name
cmpfrag() {
  python3 - "$1" "$2" "$3" "$4" "$NQ" <<'PY'
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

collection_for() {
  local n; n=$(od -An -t d8 -N8 "$1/fwd_t3/rlbwt.bin" | tr -d ' ')
  local c="$H/collection.txt"
  if [ "$(( $(stat -c %s "$c") + 1 ))" != "$n" ]; then c="$HERE/data/$(basename "$1").collection.txt"; fi
  echo "$c"
}
median() { sort -g | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2}'; }
qtime() { awk -F'::' '/loaded the BWT/ {split($3,a,"*"); l=a[1]} /processed/ {split($3,b,"*"); p=b[1]} END {printf "%.3f\n", p-l}' "$1"; }
LCSV=$(echo $LS | tr ' ' ',')
LMIN=$(echo $LS | tr ' ' '\n' | sort -n | head -1)

# ---- the reads
for E in $RATES; do
  R="$OUT/reads_e$E"
  [ -s "$R.txt" ] || python3 "$HERE/sim_reads.py" "$FINAL/queries.txt" "$R" "$NREADS" "$RLEN" "0.0$E" "$E"
  [ "$E" -ge 10 ] && { echo "RATES are percentages below 10"; exit 1; }
done
note "reads: $NREADS of $RLEN bases per error rate ($RATES %), from $FINAL/queries.txt; MEMs of length at least $LS"
note "machine: $(hostname), pinned to CPU $PIN, $ROUNDS rounds, medians"

for SIZE in $SIZES; do
  D="$H/$SIZE"; COLL=$(collection_for "$D")
  LV="$HERE/results/$SIZE/v5"
  QT=(taskset -c "$PIN" "$B/query_twolevel" --index "$D/fwd_t3" --reversed-index "$D/rev_t3"
      --level "$LV.text.lvl" --reversed-level "$LV.rev.lvl")
  SC="ell25_theta0@$D/fwd_t6@$D/rev_t6"; SCN="ell25_theta0@$D/fwd_t6"
  K="$FINAL/$SIZE/kebab"; mkdir -p "$K"
  if [ -x "$KB" ] && { [ ! -s "$K/fwd.kbb" ] || [ ! -s "$K/can.kbb" ]; }; then
    stamp "$SIZE: KeBaB filters"
    awk '{print ">s"NR; print}' "$COLL" > "$K/collection.fa"
    [ -s "$K/fwd.kbb" ] || "$KB" build -k "$KMER" --kmer-mode forward -t 8 -o "$K/fwd" "$K/collection.fa" > "$K/build_fwd.log" 2>&1
    [ -s "$K/can.kbb" ] || "$KB" build -k "$KMER" -t 8 -o "$K/can" "$K/collection.fa" > "$K/build_can.log" 2>&1
    rm -f "$K/collection.fa"
  fi
  RBI="$HERE/results/$SIZE/ropebwt3/index.fmd"; MVI="$HERE/results/$SIZE/movi_regular"
  for E in $RATES; do
    Q="$OUT/reads_e$E.txt"; QF="$OUT/reads_e$E.fa"; NQ=$(wc -l < "$Q")
    W="$OUT/$SIZE/e$E"; mkdir -p "$W/ours"
    note ""; note "== $SIZE, $E% errors"
    stamp "$SIZE, $E%: plain, shortcuts, two-level"
    "${QT[@]}" --patterns "$Q" --shortcuts "$SC" --configs "plain,$SCN,2L,2LJ" --algorithm bml --min-length "$LCSV" \
      --repeat "$ROUNDS" --dump "$W/ours" > "$W/bml.txt" 2> "$W/bml.err" || note "BML: query_twolevel exited with $?"
    "${QT[@]}" --patterns "$Q" --shortcuts "$SC" --configs "plain,$SCN,2L,2LJ" --algorithm fb --min-length "$LMIN" \
      --repeat "$ROUNDS" --dump "$W/ours" > "$W/fb.txt" 2> "$W/fb.err" || note "FB: query_twolevel exited with $?"
    cat "$W/bml.txt" "$W/fb.txt" | sed "s|@$D/fwd_t6||" | tee -a "$SUM"
    note "mismatches against plain search: $(cat "$W/bml.err" "$W/fb.err" | grep -c MISMATCH || true)"
    note "MEMs (BML): $(for L in $LS; do echo -n "L = $L: $(cat "$W"/ours/bml_L${L}_q*.tsv | wc -l); "; done)"

    # ---- KeBaB in front of plain search and block jumps
    if [ -x "$KB" ]; then
      for L in $LS; do
        ts=()
        for r in $(seq 1 "$ROUNDS"); do
          t0=$(date +%s.%N); taskset -c "$PIN" "$KB" scan -i "$K/fwd.kbb" -l "$L" -t 1 -o "$W/frag_fwd_L$L.fa" "$QF" > /dev/null 2>&1
          ts+=("$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN {printf "%.3f", b - a}')")
        done
        taskset -c "$PIN" "$KB" scan -i "$K/can.kbb" -l "$L" -t 1 -o "$W/frag_can_L$L.fa" "$QF" > /dev/null 2>&1
        note "KeBaB scan (fwd, L = $L): $(printf '%s\n' "${ts[@]}" | median) s; $(grep -c '>' "$W/frag_fwd_L$L.fa") fragments, $(grep -v '>' "$W/frag_fwd_L$L.fa" | tr -d '\n' | wc -c) of $(tr -d '\n' < "$Q" | wc -c) characters kept"
        frag2lines "$W/frag_fwd_L$L.fa" "$W/frag_L$L.txt" "$W/frag_L$L.off"
        mkdir -p "$W/kebab_L$L"
        algs=bml; [ "$L" = "$LMIN" ] && algs=bml,fb
        "${QT[@]}" --patterns "$W/frag_L$L.txt" --configs plain,2LJ --algorithm "$algs" --min-length "$L" --repeat "$ROUNDS" \
          --dump "$W/kebab_L$L" > "$W/kebab_L$L.txt" 2> "$W/kebab_L$L.err" || note "KeBaB L = $L: query_twolevel exited with $?"
        sed 's/^/KeBaB + /' "$W/kebab_L$L.txt" | tee -a "$SUM"
        for a in ${algs//,/ }; do note "  check, KeBaB + $a, L = $L: $(cmpfrag "$W/kebab_L$L" "$W/frag_L$L.off" "$W/ours" "${a}_L$L")"; done
      done
    fi

    # ---- ropebwt3, alone and behind KeBaB (canonical filter)
    if [ -x "$RB" ] && [ -s "$RBI" ]; then
      RBD="$W/ropebwt3"; mkdir -p "$RBD"; ln -sfn "$W/ours" "$RBD/ours"
      rbrun() {  # $1 name, $2 input fasta, $3.. options
        local name=$1 in=$2; shift 2
        local ts=()
        for r in $(seq 1 "$ROUNDS"); do
          local o=/dev/null; [ "$r" = 1 ] && o="$RBD/mems_$name.txt"
          taskset -c "$PIN" "$RB" mem -t 1 "$@" "$RBI" "$in" > "$o" 2> "$RBD/log_$name.txt"
          ts+=("$(qtime "$RBD/log_$name.txt")")
        done
        printf '%s\n' "${ts[@]}" | median
      }
      cases=""
      for L in $LS; do
        note "ropebwt3 BML, L = $L: $(rbrun bml$L "$QF" -l "$L") s"
        [ -x "$KB" ] && note "KeBaB + ropebwt3 BML, L = $L: $(rbrun kebab_bml$L "$W/frag_can_L$L.fa" -l "$L") s (plus the canonical scan)"
        cases="$cases,bml_L$L:bml$L:$L"
      done
      note "ropebwt3 forward-backward (--old-mem), L = $LMIN: $(rbrun fb$LMIN "$QF" --old-mem -l "$LMIN") s"
      cases="${cases#,},fb_L$LMIN:fb$LMIN:$LMIN"
      note "check, ropebwt3 against ours:"
      python3 "$HERE/check_strands.py" "$RBD" "$Q" "$D/fwd_t3/rlbwt.bin" "$B/count_strings" "$cases" | tee -a "$SUM"
    else
      note "ropebwt3 or its index ($RBI) not found"
    fi

    # ---- Movi 2
    if [ -x "$MV" ] && [ -s "$MVI/index.movi" ]; then
      MVD="$W/movi"; mkdir -p "$MVD"; ln -sfn "$W/ours" "$MVD/ours"
      cases=""
      for L in $LS; do
        ts=()
        for r in $(seq 1 "$ROUNDS"); do
          o=/dev/null; [ "$r" = 1 ] && o="$MVD/mems_movi$L.txt"
          taskset -c "$PIN" "$MV" query -i "$MVI" -r "$QF" --mem --min-mem-length "$L" --ftab-k 10 --stdout -n -t 1 > "$o" 2> "$MVD/log_L$L.txt"
          ts+=("$(sed 's/\x1b\[[0-9;]*m//g' "$MVD/log_L$L.txt" | grep -o 'processing the reads: [0-9.]*' | awk '{print $4}')")
        done
        note "Movi 2 MEMs, L = $L: $(printf '%s\n' "${ts[@]}" | median) s"
        cases="$cases,bml_L$L:movi$L:$L"
      done
      note "check, Movi 2 against ours:"
      python3 "$HERE/check_strands.py" "$MVD" "$Q" "$D/fwd_t3/rlbwt.bin" "$B/count_strings" "${cases#,}" | tee -a "$SUM"
    else
      note "Movi 2 or its index ($MVI) not found"
    fi
  done
done
stamp "done"
note ""; note "finished $(date '+%F %T')"
