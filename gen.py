# Copyright (C) 2026 Travis Gagie. Licensed under the GNU GPL v3 or later; see LICENSE.
"""Generate a repetitive file of 4-byte little-endian ints for testing.

A base sequence is drawn from an order-1 Markov source over a large alphabet
(each value has a few successors with Zipf-like weights), then copied many
times with point mutations, like a collection of similar genomes/documents.
"""
import sys, numpy as np
out, sigma, base_len, copies, mut = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
rng = np.random.default_rng(1)
vals = rng.choice(2**32 - 1, size=sigma, replace=False).astype(np.uint32)
deg = 8
succ = rng.integers(0, sigma, size=(sigma, deg))
w = 1.0 / np.arange(1, deg + 1) ** 1.5; w /= w.sum()
choice = rng.choice(deg, size=base_len, p=w)
base = np.empty(base_len, dtype=np.int64); base[0] = 0
for i in range(1, base_len):
    base[i] = succ[base[i - 1], choice[i]]
parts = []
for c in range(copies):
    x = base.copy()
    k = rng.binomial(base_len, mut)
    pos = rng.integers(0, base_len, size=k)
    x[pos] = rng.integers(0, sigma, size=k)
    parts.append(x)
seq = vals[np.concatenate(parts)]
seq.astype('<u4').tofile(out)
print(f"{out}: {len(seq)} ints, sigma={sigma}")
