// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Run-length rank structure over an integer sequence.
//
// For each symbol c we keep the sorted list of the runs of c in the sequence as
// (start, cum) pairs, where cum = number of c's before that run, plus a sentinel
// (n, total).  rank_c(i) = cum + min(len, i - start) for the last run of c
// starting before i.  This is the run-length form of the CSA's per-character
// position lists (Psi restricted to one block of F).
//
// Symbols with many runs also get a small jump table indexed by the top bits of
// the position, so a rank query binary-searches only a short slice of the list.
// The table has about R/8 entries for a symbol with R runs (~6% overhead).
#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <stdexcept>

struct RunRank {
    struct Run { uint32_t start, cum; };

    uint64_t n = 0;
    uint32_t sigma = 0;
    std::vector<uint64_t> off;      // sigma+1: runs of c are runs[off[c] .. off[c+1]) incl. sentinel
    std::vector<Run> runs;
    std::vector<uint64_t> joff;     // sigma+1: jump table of c is jumps[joff[c] .. joff[c+1])
    std::vector<uint8_t> shift;     // per symbol
    std::vector<uint32_t> jumps;
    uint64_t nruns = 0;             // number of runs (without sentinels)

    static constexpr uint64_t JUMP_MIN_RUNS = 64;

    void build(const std::vector<uint32_t>& s, uint32_t sig) {
        n = s.size(); sigma = sig;
        if (n >= (1ull << 32)) throw std::runtime_error("sequence too long for 32-bit runs");
        std::vector<uint64_t> cnt(sigma + 1, 0);
        nruns = 0;
        for (uint64_t i = 0; i < n; ++i)
            if (i == 0 || s[i] != s[i - 1]) { cnt[s[i]]++; nruns++; }
        off.assign(sigma + 1, 0);
        for (uint32_t c = 0; c < sigma; ++c) off[c + 1] = off[c] + cnt[c] + 1;
        runs.assign(off[sigma], Run{0, 0});
        std::vector<uint64_t> pos(off.begin(), off.end() - 1);
        std::vector<uint32_t> seen(sigma, 0);
        for (uint64_t i = 0; i < n;) {
            uint64_t j = i + 1;
            while (j < n && s[j] == s[i]) ++j;
            uint32_t c = s[i];
            runs[pos[c]++] = Run{(uint32_t)i, seen[c]};
            seen[c] += (uint32_t)(j - i);
            i = j;
        }
        for (uint32_t c = 0; c < sigma; ++c) runs[pos[c]] = Run{(uint32_t)n, seen[c]};

        // jump tables
        joff.assign(sigma + 1, 0);
        shift.assign(sigma, 0);
        jumps.clear();
        for (uint32_t c = 0; c < sigma; ++c) {
            joff[c] = jumps.size();
            uint64_t R = off[c + 1] - off[c] - 1;
            if (R < JUMP_MIN_RUNS) continue;
            uint64_t target = R / 8;                  // desired number of buckets
            uint8_t sh = 0;
            while ((n >> sh) > target) ++sh;
            shift[c] = sh;
            uint64_t nb = (n >> sh) + 2;
            const Run* r = &runs[off[c]];
            uint64_t k = 0;
            for (uint64_t b = 0; b < nb; ++b) {
                uint64_t lim = b << sh;               // #runs with start < lim
                while (k < R && r[k].start < lim) ++k;
                jumps.push_back((uint32_t)k);
            }
        }
        joff[sigma] = jumps.size();
    }

    // number of occurrences of c in s[0, i)
    inline uint64_t rank(uint32_t c, uint64_t i) const {
        if (c >= sigma || i == 0) return 0;
        const Run* r = runs.data() + off[c];
        uint64_t R = off[c + 1] - off[c] - 1;
        if (R == 0) return 0;
        uint64_t lo = 0, hi = R;                      // find j = #runs with start < i, in [lo, hi]
        if (joff[c + 1] != joff[c]) {
            const uint32_t* J = jumps.data() + joff[c];
            uint64_t b = i >> shift[c];
            lo = J[b]; hi = J[b + 1];
            if (hi < R && r[hi].start < i) hi = R;    // (cannot happen, defensive)
        }
        while (lo < hi) {                             // first index with start >= i
            uint64_t mid = (lo + hi) >> 1;
            if (r[mid].start < i) lo = mid + 1; else hi = mid;
        }
        if (lo == 0) return 0;
        const Run& a = r[lo - 1];
        uint64_t len = r[lo].cum - a.cum;
        uint64_t d = i - a.start;
        return a.cum + (d < len ? d : len);
    }

    inline void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t& r1, uint64_t& r2) const {
        r1 = rank(c, i1); r2 = rank(c, i2);
    }

    // count of c in the whole sequence
    inline uint64_t total(uint32_t c) const {
        return c < sigma ? runs[off[c + 1] - 1].cum : 0;
    }

    uint64_t bytes() const {
        return off.size() * 8 + runs.size() * sizeof(Run) + joff.size() * 8 +
               shift.size() + jumps.size() * 4;
    }

    template <class V> static void wv(FILE* f, const V& v) {
        uint64_t k = v.size(); fwrite(&k, 8, 1, f);
        if (k) fwrite(v.data(), sizeof(v[0]), k, f);
    }
    template <class V> static void rv(FILE* f, V& v) {
        uint64_t k; if (fread(&k, 8, 1, f) != 1) throw std::runtime_error("read");
        v.resize(k);
        if (k && fread(v.data(), sizeof(v[0]), k, f) != k) throw std::runtime_error("read");
    }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nruns, 8, 1, f);
        wv(f, off); wv(f, runs); wv(f, joff); wv(f, shift); wv(f, jumps);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1)
            throw std::runtime_error("read");
        rv(f, off); rv(f, runs); rv(f, joff); rv(f, shift); rv(f, jumps);
    }
};
