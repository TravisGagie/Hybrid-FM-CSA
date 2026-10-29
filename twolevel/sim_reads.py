#!/usr/bin/env python3
# Simulated long reads from the query haplotypes (forward strand only, since our index
# finds only forward-strand MEMs): N reads of length LEN, starting uniformly at random,
# each base then substituted, deleted or followed by an inserted base, each with
# probability RATE / 3.  Writes PREFIX.txt (one read per line) and PREFIX.fa.
#
#   python3 sim_reads.py QUERIES.txt PREFIX N LEN RATE SEED
import random, sys
qpath, prefix, n, length, rate, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
rng = random.Random(seed)
Q = [l.strip() for l in open(qpath) if len(l.strip()) > length]
weights = [len(q) - length for q in Q]
total = sum(weights)
B = "ACGT"
with open(prefix + ".txt", "w") as t, open(prefix + ".fa", "w") as f:
    for i in range(n):
        x = rng.randrange(total)
        k = 0
        while x >= weights[k]: x -= weights[k]; k += 1
        src = Q[k][x:x + length]
        out = []
        for c in src:
            u = rng.random()
            if u < rate / 3:                      # substitution
                out.append(rng.choice([b for b in B if b != c]))
            elif u < 2 * rate / 3:                # deletion
                pass
            elif u < rate:                        # insertion after c
                out.append(c); out.append(rng.choice(B))
            else:
                out.append(c)
        r = "".join(out)
        t.write(r + "\n")
        f.write(">p%d\n%s\n" % (i, r))
