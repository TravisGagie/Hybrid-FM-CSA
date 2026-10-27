#!/usr/bin/env python3
# Explains the MEMs that only one of ropebwt3 and our index reports.
#
#   python3 check_strands.py RESULTS_DIR/ropebwt3 QUERIES FWD_RLBWT COUNT_STRINGS
#
# ropebwt3 searches both strands, so:
#  - a MEM only ropebwt3 reports should not occur in T, while its reverse
#    complement should (it matches only the reverse strand);
#  - a MEM only we report should lie inside a longer ropebwt3 MEM that does not
#    occur in T (ropebwt3 extended it through the reverse strand).
import os, subprocess, sys, bisect
d, qpath, rlbwt, tool = sys.argv[1:5]
Q = [l.rstrip("\n\r") for l in open(qpath) if l.strip()]
tr = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(tr)[::-1]
cases = [("bml_L100", "bml100", 100), ("bml_L1000", "bml1000", 1000), ("fb_L1000", "fb1000", 1000), ("fb_L1000", "fball", 1000)]
for ours_name, rb_name, L in cases:
    ours = set()
    for q in (0, 1):
        for line in open(os.path.join(d, "ours", f"{ours_name}_q{q}.tsv")):
            s, e, _ = line.split(); ours.add((q, int(s), int(e)))
    rb = set()
    for line in open(os.path.join(d, f"mems_{rb_name}.txt")):
        f = line.split("\t"); q, s, e = int(f[0][1:]), int(f[1]), int(f[2])
        if e - s >= L: rb.add((q, s, e))
    only_ours, only_rb = sorted(ours - rb), sorted(rb - ours)
    # strings to count: each ropebwt3-only MEM and its reverse complement
    strs = []
    for q, s, e in only_rb:
        x = Q[q][s:e]; strs += [x, rc(x)]
    with open(os.path.join(d, "check.txt"), "w") as f:
        f.write("".join(x + "\n" for x in strs))
    out = subprocess.run([tool, rlbwt, os.path.join(d, "check.txt")], capture_output=True, text=True, check=True).stdout.split()
    cnt = list(map(int, out))
    rb_rev_only = sum(1 for i in range(len(only_rb)) if cnt[2 * i] == 0 and cnt[2 * i + 1] > 0)
    rb_bad = [only_rb[i] for i in range(len(only_rb)) if not (cnt[2 * i] == 0 and cnt[2 * i + 1] > 0)]
    # each of our MEMs missing from ropebwt3 should lie inside a ropebwt3-only MEM
    starts = {q: sorted((s, e) for qq, s, e in only_rb if qq == q) for q in (0, 1)}
    inside = 0; ours_bad = []
    for q, s, e in only_ours:
        lst = starts[q]; k = bisect.bisect_right(lst, (s, 1 << 62))
        ok = any(lst[j][0] <= s and lst[j][1] >= e for j in range(max(0, k - 50), k))
        inside += ok
        if not ok: ours_bad.append((q, s, e))
    print(f"{ours_name} vs {rb_name}: only ropebwt3 {len(only_rb)}, of which {rb_rev_only} occur only on the reverse strand; "
          f"only ours {len(only_ours)}, of which {inside} lie inside a ropebwt3-only MEM")
    for x in rb_bad[:5]: print("   unexplained ropebwt3 MEM", x, "length", x[2] - x[1])
    for x in ours_bad[:5]: print("   unexplained MEM of ours", x, "length", x[2] - x[1])
