// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Compressed run-length rank structure, in the spirit of the RLCSA of
// Brown, Gagie, Manzini, Navarro and Sciortino ("Faster Run-Length Compressed
// Suffix Arrays"), simplified for practice.
//
// Runs of the sequence are listed in F'-order: grouped by symbol, and by start
// position within a symbol.  We keep two Elias-Fano sequences over them:
//
//   KS[j] = c * n + start_j   (globally increasing).  A predecessor query on
//                             KS replaces the binary search in the increasing
//                             interval of Psi' for c, and also locates that
//                             interval, so no per-symbol table is needed.
//   FP[j] = position in F of the first character of run j (the paper's B_F).
//
// rank_c(i): j = #{KS < c*n + i} - 1 is the last run of c starting before i;
// the answer is FP[j] - FP[first run of c] + min(FP[j+1] - FP[j], i - start_j).
// Because we search on positions directly, the rank on B_L (and the B_FL trick
// that avoids it) is not needed.
//
// Space is about r(log(sigma*n/r) + 2) + r(log(n/r) + 2) bits.
#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>
#include <stdexcept>
#ifdef __BMI2__
#include <immintrin.h>
#endif

static inline uint64_t select_in_word(uint64_t w, uint64_t k) {   // position of k-th (0-based) 1
#ifdef __BMI2__
    return __builtin_ctzll(_pdep_u64(1ull << k, w));
#else
    for (uint64_t i = 0;; ++i) if ((w >> i) & 1) { if (!k) return i; --k; }
#endif
}

struct EliasFano {
    uint64_t n = 0, u = 0;              // count, universe
    uint32_t l = 0;                     // low bits
    std::vector<uint64_t> low, high;
    uint64_t hbits = 0;
    std::vector<uint64_t> s0, s1;       // positions of every SAMPLE-th zero / one
    static constexpr uint64_t SAMPLE = 256;

    void build(const std::vector<uint64_t>& v, uint64_t universe) {
        n = v.size(); u = universe;
        l = 0;
        if (n && u / n > 1) while ((u / n) >> (l + 1)) ++l;
        low.assign((n * l + 63) / 64 + 1, 0);
        hbits = n + (u >> l) + 2;
        high.assign(hbits / 64 + 2, 0);
        for (uint64_t j = 0; j < n; ++j) {
            if (j && v[j] < v[j - 1]) throw std::runtime_error("EF: not sorted");
            uint64_t x = v[j];
            if (l) set_low(j, x & ((1ull << l) - 1));
            uint64_t p = (x >> l) + j;
            high[p >> 6] |= 1ull << (p & 63);
        }
        s0.clear(); s1.clear();
        uint64_t c0 = 0, c1 = 0;
        for (uint64_t p = 0; p < hbits; ++p) {
            if ((high[p >> 6] >> (p & 63)) & 1) { if (c1 % SAMPLE == 0) s1.push_back(p); ++c1; }
            else { if (c0 % SAMPLE == 0) s0.push_back(p); ++c0; }
        }
    }
    inline void set_low(uint64_t j, uint64_t x) {
        uint64_t b = j * l;
        low[b >> 6] |= x << (b & 63);
        if ((b & 63) + l > 64) low[(b >> 6) + 1] |= x >> (64 - (b & 63));
    }
    inline uint64_t get_low(uint64_t j) const {
        if (!l) return 0;
        uint64_t b = j * l, w = b >> 6, o = b & 63;
        uint64_t x = low[w] >> o;
        if (o + l > 64) x |= low[w + 1] << (64 - o);
        return x & ((1ull << l) - 1);
    }
    inline bool hbit(uint64_t p) const { return (high[p >> 6] >> (p & 63)) & 1; }

    inline uint64_t select1(uint64_t k) const {        // position of k-th one (0-based)
        uint64_t p = s1[k / SAMPLE];
        uint64_t rem = k % SAMPLE;
        uint64_t w = p >> 6;
        uint64_t word = high[w] & (~0ull << (p & 63));
        for (;;) {
            uint64_t c = __builtin_popcountll(word);
            if (rem < c) return (w << 6) + select_in_word(word, rem);
            rem -= c; word = high[++w];
        }
    }
    inline uint64_t select0(uint64_t k) const {        // position of k-th zero (0-based)
        uint64_t p = s0[k / SAMPLE];
        uint64_t rem = k % SAMPLE;
        uint64_t w = p >> 6;
        uint64_t word = ~high[w] & (~0ull << (p & 63));
        for (;;) {
            uint64_t c = __builtin_popcountll(word);
            if (rem < c) return (w << 6) + select_in_word(word, rem);
            rem -= c; word = ~high[++w];
        }
    }
    inline uint64_t access(uint64_t j) const {
        return ((select1(j) - j) << l) | get_low(j);
    }
    // access(j) and access(j+1) together (j+1 < n)
    inline void access2(uint64_t j, uint64_t& a, uint64_t& b) const {
        uint64_t p = select1(j);
        a = ((p - j) << l) | get_low(j);
        uint64_t q = p + 1;
        while (!hbit(q)) ++q;
        b = ((q - j - 1) << l) | get_low(j + 1);
    }
    // number of elements < x; also returns the high-bit position of the last one (if any)
    inline uint64_t count_less(uint64_t x) const {
        if (x >= u) return n;
        uint64_t h = x >> l, xl = x & ((1ull << l) - 1);
        uint64_t idx, pos;
        if (h == 0) { idx = 0; pos = 0; }
        else { uint64_t p = select0(h - 1); idx = p - (h - 1); pos = p + 1; }
        while (pos < hbits && hbit(pos) && get_low(idx) < xl) { ++idx; ++pos; }
        return idx;
    }
    uint64_t bytes() const { return (low.size() + high.size() + s0.size() + s1.size()) * 8 + 40; }

    template <class V> static void wv(FILE* f, const V& v) {
        uint64_t k = v.size(); fwrite(&k, 8, 1, f); if (k) fwrite(v.data(), sizeof(v[0]), k, f);
    }
    template <class V> static void rv(FILE* f, V& v) {
        uint64_t k; if (fread(&k, 8, 1, f) != 1) throw std::runtime_error("read");
        v.resize(k); if (k && fread(v.data(), sizeof(v[0]), k, f) != k) throw std::runtime_error("read");
    }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&u, 8, 1, f); fwrite(&l, 4, 1, f); fwrite(&hbits, 8, 1, f);
        wv(f, low); wv(f, high); wv(f, s0); wv(f, s1);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&u, 8, 1, f) != 1 || fread(&l, 4, 1, f) != 1 ||
            fread(&hbits, 8, 1, f) != 1) throw std::runtime_error("read");
        rv(f, low); rv(f, high); rv(f, s0); rv(f, s1);
    }
};

struct EFRank {
    uint64_t n = 0;
    uint32_t sigma = 0;
    uint64_t nruns = 0;
    EliasFano KS, FP;

    void build(const std::vector<uint32_t>& s, uint32_t sig) {
        n = s.size(); sigma = sig;
        if ((unsigned __int128)sigma * (n + 1) >= ((unsigned __int128)1 << 63))
            throw std::runtime_error("sigma*n too large");
        std::vector<uint64_t> cnt(sigma + 1, 0), runsof(sigma + 1, 0);
        nruns = 0;
        for (uint64_t i = 0; i < n; ++i) {
            cnt[s[i]]++;
            if (i == 0 || s[i] != s[i - 1]) { runsof[s[i]]++; nruns++; }
        }
        // F-start of each symbol, and slot of each symbol's first run
        std::vector<uint64_t> Cf(sigma + 1, 0), slot(sigma + 1, 0);
        for (uint32_t c = 0; c < sigma; ++c) { Cf[c + 1] = Cf[c] + cnt[c]; slot[c + 1] = slot[c] + runsof[c]; }
        std::vector<uint64_t> ks(nruns), fp(nruns);
        std::vector<uint64_t> seen(sigma, 0);
        for (uint64_t i = 0; i < n;) {
            uint64_t j = i + 1;
            while (j < n && s[j] == s[i]) ++j;
            uint32_t c = s[i];
            uint64_t k = slot[c]++;
            ks[k] = (uint64_t)c * n + i;
            fp[k] = Cf[c] + seen[c];
            seen[c] += j - i;
            i = j;
        }
        KS.build(ks, (uint64_t)sigma * n);
        FP.build(fp, n);
    }

    // occurrences of c in s[0, i)
    inline uint64_t rank(uint32_t c, uint64_t i) const {
        if (c >= sigma || i == 0) return 0;
        uint64_t base = (uint64_t)c * n;
        uint64_t b = KS.count_less(base);          // first run of c
        uint64_t a = KS.count_less(base + (i > n ? n : i));
        if (a == b) return 0;
        return rank_from(b, a, c, i);
    }
    // when the index b of c's first run is already known (backward search)
    inline uint64_t rank_from(uint64_t b, uint64_t a, uint32_t c, uint64_t i) const {
        uint64_t j = a - 1;
        uint64_t start = KS.access(j) - (uint64_t)c * n;
        uint64_t f, fn;
        if (j + 1 < nruns) FP.access2(j, f, fn); else { f = FP.access(j); fn = n; }
        uint64_t cc = FP.access(b);
        uint64_t len = fn - f, d = i - start;
        return f - cc + (d < len ? d : len);
    }
    // rank at two positions for the same symbol (one lookup of c's first run)
    inline void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t& r1, uint64_t& r2) const {
        if (c >= sigma) { r1 = r2 = 0; return; }
        uint64_t base = (uint64_t)c * n;
        uint64_t b = KS.count_less(base);
        uint64_t a1 = i1 ? KS.count_less(base + i1) : b;
        uint64_t a2 = i2 ? KS.count_less(base + i2) : b;
        r1 = a1 == b ? 0 : rank_from(b, a1, c, i1);
        r2 = a2 == b ? 0 : rank_from(b, a2, c, i2);
    }
    uint64_t bytes() const { return KS.bytes() + FP.bytes(); }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nruns, 8, 1, f); KS.save(f); FP.save(f);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1)
            throw std::runtime_error("read");
        KS.load(f); FP.load(f);
    }
};
