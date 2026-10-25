// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// GapRank: rank over a sequence by gap-coded position lists, without run-lengths.
//
// For each symbol c, the sorted positions of c in the sequence are stored as
// Elias-delta codes of their gaps, with an absolute sample every SAMPLE elements.
// Applied to a BWT, the list of c is Psi restricted to c's block of F, so this is
// a Sadakane-style CSA with a gap-coded Psi; rank_c(i) is a binary search over c's
// samples followed by at most SAMPLE-1 decodes.  The gaps are small where c's
// occurrences cluster, so the space adapts to the local skew (high-order entropy)
// of the sequence, not to runs.  The number of occurrences of each symbol is kept
// in an Elias-Fano list, which also gives C[c], so lf2 needs no C array.
#pragma once
#include "efrank.hpp"
#include <vector>
#include <cstdint>
#include <stdexcept>

struct GapRank {
    static constexpr uint64_t SAMPLE = 64;
    uint64_t n = 0;
    uint32_t sigma = 0;
    std::vector<uint64_t> bits;        // delta codes of the gaps
    uint64_t nbits = 0;
    std::vector<uint32_t> sval;        // sampled values (positions)
    std::vector<uint64_t> soff;        // bit offset of the code after each sample
    EliasFano cnt;                     // sigma+1: number of occurrences of symbols < c (= C[c])
    EliasFano sfirst;                  // sigma+1: index of c's first sample

    void put(uint64_t x, int w) {      // append the w low bits of x, most significant first
        for (int b = w - 1; b >= 0; --b) {
            if ((x >> b) & 1) bits[nbits >> 6] |= 1ull << (63 - (nbits & 63));
            ++nbits;
            if ((nbits >> 6) >= bits.size()) bits.resize(bits.size() * 2 + 1, 0);
        }
    }
    static int lg(uint64_t x) { return 63 - __builtin_clzll(x); }
    void put_gamma(uint64_t x) { int l = lg(x); put(0, l); put(x, l + 1); }
    void put_delta(uint64_t x) { int l = lg(x); put_gamma(l + 1); put(x, l); }
    inline uint64_t peek64(uint64_t p) const {      // 64 bits starting at bit p
        uint64_t w = p >> 6, o = p & 63;
        uint64_t a = bits[w] << o;
        if (o && w + 1 < bits.size()) a |= bits[w + 1] >> (64 - o);
        return a;
    }
    inline uint64_t get(uint64_t& p, int w) const {
        if (w == 0) return 0;
        uint64_t x = peek64(p) >> (64 - w);
        p += w;
        return x;
    }
    inline uint64_t get_gamma(uint64_t& p) const {
        int z = __builtin_clzll(peek64(p) | 1);       // codes here are < 2^63
        p += z;
        return get(p, z + 1);
    }
    inline uint64_t get_delta(uint64_t& p) const {
        int l = (int)get_gamma(p) - 1;
        if (l == 0) return 1;
        return (1ull << l) | get(p, l);
    }

    void build(const std::vector<uint32_t>& s, uint32_t sig) {
        n = s.size(); sigma = sig;
        if (n >= (1ull << 32)) throw std::runtime_error("sequence too long for GapRank");
        std::vector<uint64_t> c(sigma + 1, 0);
        for (auto x : s) c[x + 1]++;
        for (uint32_t a = 0; a < sigma; ++a) c[a + 1] += c[a];
        std::vector<uint32_t> pos(n);
        {
            std::vector<uint64_t> fill(c.begin(), c.end() - 1);
            for (uint64_t i = 0; i < n; ++i) pos[fill[s[i]]++] = (uint32_t)i;
        }
        bits.assign(n / 16 + 2, 0); nbits = 0;
        std::vector<uint64_t> sf(sigma + 1, 0);
        sval.clear(); soff.clear();
        for (uint32_t a = 0; a < sigma; ++a) {
            sf[a] = sval.size();
            uint64_t prev = 0;
            for (uint64_t j = c[a]; j < c[a + 1]; ++j) {
                uint64_t k = j - c[a];
                if (k % SAMPLE == 0) { sval.push_back(pos[j]); soff.push_back(nbits); }
                else put_delta(pos[j] - prev);
                prev = pos[j];
            }
        }
        sf[sigma] = sval.size();
        bits.resize(nbits / 64 + 2);
        cnt.build(c, n + 1);
        sfirst.build(sf, sval.size() + 1);
    }
    // occurrences of a in s[0, i)
    inline uint64_t rank(uint32_t a, uint64_t i) const {
        if (a >= sigma || i == 0) return 0;
        uint64_t c0, c1, f0, f1;
        cnt.access2(a, c0, c1);
        if (c0 == c1) return 0;
        sfirst.access2(a, f0, f1);
        // last sample of a with value < i
        uint64_t lo = f0, hi = f1;                     // samples [f0, f1)
        if (sval[lo] >= i) return 0;
        while (hi - lo > 1) {
            uint64_t mid = (lo + hi) >> 1;
            if (sval[mid] < i) lo = mid; else hi = mid;
        }
        uint64_t k = (lo - f0) * SAMPLE;               // index of the sample within a's list
        uint64_t v = sval[lo], p = soff[lo], len = c1 - c0, r = k + 1;
        uint64_t stop = std::min(len, k + SAMPLE);
        while (r < stop) {
            v += get_delta(p);
            if (v >= i) break;
            ++r;
        }
        return r;
    }
    inline void rank2(uint32_t a, uint64_t i1, uint64_t i2, uint64_t& r1, uint64_t& r2) const {
        r1 = rank(a, i1); r2 = rank(a, i2);
    }
    inline void lf2(uint32_t a, uint64_t i1, uint64_t i2, uint64_t& o1, uint64_t& o2) const {
        if (a >= sigma) { o1 = o2 = n; return; }
        uint64_t c0 = cnt.access(a);
        o1 = c0 + rank(a, i1); o2 = c0 + rank(a, i2);
    }
    uint64_t bytes() const {
        return bits.size() * 8 + sval.size() * 4 + soff.size() * 8 + cnt.bytes() + sfirst.bytes();
    }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nbits, 8, 1, f);
        EliasFano::wv(f, bits); EliasFano::wv(f, sval); EliasFano::wv(f, soff); cnt.save(f); sfirst.save(f);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nbits, 8, 1, f) != 1)
            throw std::runtime_error("read");
        EliasFano::rv(f, bits); EliasFano::rv(f, sval); EliasFano::rv(f, soff); cnt.load(f); sfirst.load(f);
    }
};
