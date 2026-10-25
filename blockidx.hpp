// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Block index: backward search over chains of most frequent successors.
//
// Let s(x) be the most frequent successor of x in T (ties: smaller id).  Encode T as E:
// E[i] = # if T[i] = s(T[i-1]) (T[i-1] not a separator), else the real character T[i];
// # is smaller than every real character.  A block is a real character d followed by its
// j hashes, so it stands for d, s(d), s(s(d)), ... (j+1 characters).  The suffixes of E
// that start with a real character, in E's suffix order, are the suffixes of the block
// sequence, with block symbols ordered by (d ascending, j descending); the suffixes that
// start with # come before them (after the terminator), ordered by the number of leading
// hashes (more first) and then by the block suffix that follows.  We store
//   - the BWT of the block sequence (Elias-Fano run lists, C-free LF): one exact step per
//     block of the pattern;
//   - a bit per row of E's suffix order saying whether the row is preceded by a hash, so
//     that I(#Y) = 1 + [rank_#(lo), rank_#(hi)) (the terminator row is row 0);
//   - H: for the rows that start with a hash, the character that hash stands for, which
//     resolves a pattern's leading characters exactly (they may be hashes in the text);
//   - Ncnt[k]: the number of rows that start with at least k hashes;
//   - s, and the block symbols of each d (their j values), to encode patterns.
// Counting and longest-suffix search are exact.
#pragma once
#include "efrank.hpp"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include "libsais.h"

struct BlockIndex {
    uint32_t sigma = 0;                 // real characters 1..sigma (dense ids)
    uint32_t sep = 0;                   // separator id (never predicted, never predicts)
    std::vector<uint32_t> top;          // sigma+1: most frequent successor (0 = none)
    EliasFano dstart;                   // sigma+2: first block symbol (minus 1) of each d
    std::vector<uint16_t> jv;           // j value of each block symbol (descending within d)
    uint64_t nsym = 0;                  // block symbols (ids 1..nsym; 0 = terminator)
    uint64_t n = 0, nb = 0, nh = 0;     // characters, blocks, hash positions
    EFRank bbwt;                        // BWT of the block sequence
    EFRank hashbit;                     // over E's rows: 1 if preceded by a hash
    uint64_t hb_c1 = 0;                 // C[1] of hashbit (number of zeros)
    EFRank H;                           // over the hash rows: the character the hash stands for
    std::vector<uint64_t> ncnt;         // ncnt[k] = rows starting with >= k hashes (k >= 1)
    // statistics
    uint64_t steps = 0, hsteps = 0, hchecks = 0;

    static bool lessblk(uint64_t a, uint64_t b) { return a < b; }

    // block symbol id of (d, j), 0 if absent
    inline uint64_t sym(uint32_t d, uint32_t j) const {
        if (d == 0 || d > sigma) return 0;
        uint64_t a, b;
        dstart.access2(d, a, b);
        // jv[a..b) descending
        uint64_t lo = a, hi = b;
        while (lo < hi) { uint64_t m = (lo + hi) >> 1; if (jv[m] > j) lo = m + 1; else hi = m; }
        if (lo < b && jv[lo] == j) return lo + 1;
        return 0;
    }
    // block symbols (d, j') with j' >= j: ids [first, last] (empty if first > last)
    inline void syms_atleast(uint32_t d, uint32_t j, uint64_t& first, uint64_t& last) const {
        first = 1; last = 0;
        if (d == 0 || d > sigma) return;
        uint64_t a, b;
        dstart.access2(d, a, b);
        uint64_t lo = a, hi = b;            // first index with jv < j
        while (lo < hi) { uint64_t m = (lo + hi) >> 1; if (jv[m] >= j) lo = m + 1; else hi = m; }
        if (lo == a) return;
        first = a + 1; last = lo;           // ids a+1 .. lo
    }
    // rows of the block BWT whose suffix starts with a block symbol in [first, last]
    inline void sym_range(uint64_t first, uint64_t last, uint64_t& l, uint64_t& r) const {
        uint64_t x, y;
        bbwt.lf2((uint32_t)first, 0, 0, l, x);
        bbwt.lf2((uint32_t)(last + 1), 0, 0, r, y);
        if (last + 1 > nsym) r = nb + 1;
    }
    inline bool bstep(uint64_t c, uint64_t& l, uint64_t& r) {      // block-BWT step
        ++steps;
        uint64_t a, b;
        bbwt.lf2((uint32_t)c, l, r, a, b);
        l = a; r = b;
        return l < r;
    }
    // E-row interval of block rows [l, r)
    inline void to_erows(uint64_t& l, uint64_t& r) const { l += nh; r += nh; }
    inline bool hstep(uint64_t& l, uint64_t& r) {                  // I(#Y) from I(Y), in E-rows
        ++hsteps;
        uint64_t a, b;
        hashbit.lf2(1, l, r, a, b);
        l = a - hb_c1 + 1; r = b - hb_c1 + 1;
        return l < r;
    }
    inline uint64_t hcount(uint32_t c, uint64_t l, uint64_t r) {   // H = c among E-rows [l, r) (hash rows)
        ++hchecks;
        uint64_t a, b;
        H.rank2(c, l - 1, r - 1, a, b);
        return b - a;
    }

    // T: dense ids 1..sig; sep_ = separator id (0 if none).  Appends a separator (a new
    // symbol if there is none) so that the text ends with a real character.
    void build(std::vector<uint32_t> T, uint32_t sig, uint32_t sep_, bool verbose = false) {
        sigma = sig; sep = sep_;
        if (!sep) { sep = ++sigma; }
        T.push_back(sep);
        n = T.size();
        top.assign(sigma + 1, 0);
        {
            std::vector<uint64_t> pr; pr.reserve(n);
            for (uint64_t i = 0; i + 1 < n; ++i) if (T[i] != sep && T[i + 1] != sep) pr.push_back((uint64_t)T[i] << 32 | T[i + 1]);
            std::sort(pr.begin(), pr.end());
            std::vector<uint64_t> best(sigma + 1, 0);
            for (size_t i = 0; i < pr.size();) {
                size_t k = i; while (k < pr.size() && pr[k] == pr[i]) ++k;
                uint32_t a = pr[i] >> 32, b = (uint32_t)pr[i];
                if (k - i > best[a]) { best[a] = k - i; top[a] = b; }
                i = k;
            }
        }
        std::vector<uint8_t> hash(n, 0);
        nh = 0;
        for (uint64_t i = 1; i < n; ++i)
            if (T[i - 1] != sep && T[i] != sep && T[i] == top[T[i - 1]]) { hash[i] = 1; ++nh; }
        nb = n - nh;
        // block symbols
        std::vector<uint64_t> key; key.reserve(nb);
        uint64_t maxj = 0;
        for (uint64_t i = 0; i < n;) {
            uint64_t k = i + 1; while (k < n && hash[k]) ++k;
            uint64_t j = k - i - 1; maxj = std::max(maxj, j);
            key.push_back((uint64_t)T[i] << 32 | (0xFFFFFFFFull - j));
            i = k;
        }
        if (maxj >= 65535) throw std::runtime_error("chain of more than 65534 hashes");
        std::vector<uint64_t> u(key); std::sort(u.begin(), u.end()); u.erase(std::unique(u.begin(), u.end()), u.end());
        nsym = u.size();
        jv.resize(nsym);
        std::vector<uint64_t> ds(sigma + 2, 0);
        for (uint64_t k = 0; k < nsym; ++k) { jv[k] = (uint16_t)(0xFFFFFFFFull - (uint32_t)u[k]); ds[(u[k] >> 32) + 1]++; }
        for (uint32_t d = 0; d <= sigma; ++d) ds[d + 1] += ds[d];
        dstart.build(ds, nsym + 1);
        std::vector<int32_t> B(nb + 1);
        for (uint64_t k = 0; k < nb; ++k) B[k] = (int32_t)(std::lower_bound(u.begin(), u.end(), key[k]) - u.begin()) + 1;
        B[nb] = 0;
        std::vector<uint64_t>().swap(key); std::vector<uint64_t>().swap(u);
        {   // BWT of the block sequence
            uint64_t m = nb + 1, k = nsym + 1;
            std::vector<int32_t> SA(m + 6 * k);
            if (libsais_int(B.data(), SA.data(), (int32_t)m, (int32_t)k, (int32_t)(6 * k)) != 0) throw std::runtime_error("libsais");
            std::vector<uint32_t> L(m);
            for (uint64_t q = 0; q < m; ++q) L[q] = (uint32_t)B[SA[q] ? SA[q] - 1 : m - 1];
            bbwt.build(L, (uint32_t)k);
            if (verbose) fprintf(stderr, "  block BWT: %llu blocks, %llu symbols (%.1f MB)\n",
                                 (unsigned long long)nb, (unsigned long long)nsym, bbwt.bytes() / 1048576.0);
        }
        std::vector<int32_t>().swap(B);
        {   // E's suffix order: hashbit, H, ncnt
            uint64_t m = n + 1;
            std::vector<int32_t> X(m);
            for (uint64_t i = 0; i < n; ++i) X[i] = hash[i] ? 1 : (int32_t)T[i] + 1;
            X[n] = 0;
            uint64_t k = sigma + 2;
            std::vector<int32_t> SA(m + 6 * k);
            if (libsais_int(X.data(), SA.data(), (int32_t)m, (int32_t)k, (int32_t)(6 * k)) != 0) throw std::runtime_error("libsais");
            std::vector<uint32_t> hb(m), hv; hv.reserve(nh);
            uint64_t zeros = 0;
            for (uint64_t q = 0; q < m; ++q) {
                uint64_t p = SA[q];
                uint32_t prev = p ? (uint32_t)(X[p - 1] == 1) : 0;
                hb[q] = prev; zeros += !prev;
                if (p < n && hash[p]) hv.push_back(T[p]);
                if (q >= 1 && q <= nh && !(p < n && hash[p])) throw std::runtime_error("hash rows not contiguous");
            }
            hashbit.build(hb, 2);
            hb_c1 = zeros;
            H.build(hv, sigma + 1);
            // ncnt[k]: rows with >= k leading hashes = sum over m >= k of #(blocks with j >= m)
            std::vector<uint64_t> cj(maxj + 2, 0);
            for (uint64_t i = 0; i < n;) {
                uint64_t e = i + 1; while (e < n && hash[e]) ++e;
                cj[e - i - 1]++;           // a block with j hashes
                i = e;
            }
            // blocks with j >= m
            std::vector<uint64_t> ge(maxj + 2, 0);
            for (int64_t m2 = (int64_t)maxj; m2 >= 0; --m2) ge[m2] = cj[m2] + (m2 + 1 <= (int64_t)maxj ? ge[m2 + 1] : 0);
            ncnt.assign(maxj + 2, 0);
            for (int64_t m2 = (int64_t)maxj; m2 >= 1; --m2) ncnt[m2] = ge[m2] + ncnt[m2 + 1];
            if (ncnt[1] != nh) throw std::runtime_error("hash count mismatch");
            if (verbose) fprintf(stderr, "  E: %llu hash rows; hash bits %.1f MB, H %.1f MB\n",
                                 (unsigned long long)nh, hashbit.bytes() / 1048576.0, H.bytes() / 1048576.0);
        }
    }
    uint64_t rows_atleast(uint64_t k) const { return k < ncnt.size() ? ncnt[k] : 0; }

    // Number of occurrences of p (dense ids), exactly.
    uint64_t count(const std::vector<uint32_t>& p) {
        size_t m = p.size();
        if (m == 0) return n;
        for (auto x : p) if (x == 0 || x > sigma || x == sep) return 0;
        // block starts in the pattern
        std::vector<uint8_t> hs(m, 0);
        for (size_t t = 1; t < m; ++t) hs[t] = (p[t] == top[p[t - 1]]);
        // last block [a, m-1]
        size_t a = m - 1; while (a > 0 && hs[a]) --a;
        uint32_t len = (uint32_t)(m - 1 - a);      // hashes after p[a]
        if (a == 0) {                              // the whole pattern is one block
            uint64_t f, g, l, r;
            syms_atleast(p[0], len, f, g);
            uint64_t c = 0;
            if (f <= g) { sym_range(f, g, l, r); c += r - l; }
            uint64_t k = rows_atleast(len + 1);
            if (k) c += hcount(p[0], 1, 1 + k);
            return c;
        }
        uint64_t f, g, l, r;
        syms_atleast(p[a], len, f, g);
        if (f > g) return 0;
        sym_range(f, g, l, r);
        while (true) {
            size_t b = a - 1; while (b > 0 && hs[b]) --b;          // previous block [b, a-1]
            uint32_t j = (uint32_t)(a - 1 - b);
            if (b == 0) {                                          // the first block: real or hidden first character
                uint64_t c = 0;
                uint64_t s = sym(p[0], j);
                if (s) { uint64_t l2 = l, r2 = r; if (bstep(s, l2, r2)) c += r2 - l2; }
                uint64_t el = l, er = r; to_erows(el, er);
                bool ok = true;
                for (uint32_t t = 0; t <= j && ok; ++t) ok = hstep(el, er);
                if (ok) c += hcount(p[0], el, er);
                return c;
            }
            uint64_t s = sym(p[b], j);
            if (!s || !bstep(s, l, r)) return 0;
            a = b;
        }
    }

    // Longest-suffix search over a query of dense ids Q (0 = absent), with hs[t] = 1 if
    // Q[t] = s(Q[t-1]) within its segment and cut[t] = 1 if t starts a segment: the smallest
    // i >= floor such that Q[i..j] occurs (Q[j] present).  Exact.
    //   Q[x..j] with x inside the last block occurs iff a real Q[x] is followed by at least
    //   j-x hashes, or a hash row with at least j-x+1 leading hashes stands for Q[x].
    //   Extending left across whole blocks is exact with the block BWT.  When a block does
    //   not match whole, its longest matching suffix is found from the right: Q[x..j] with
    //   x = a-i occurs iff a real Q[x] followed by exactly i-1 hashes precedes the current
    //   rows, or the rows #^i Y (i hash steps) include one whose hash stands for Q[x].
    uint64_t longest_suffix(const std::vector<uint32_t>& Q, const std::vector<uint8_t>& hs,
                            const std::vector<uint8_t>& cut, uint64_t j, uint64_t floor) {
        uint64_t a = j; while (a > floor && hs[a]) --a;
        uint64_t l = 0, r = 0;
        bool real_ok = false;
        if (!hs[a]) {
            uint64_t f, g;
            syms_atleast(Q[a], (uint32_t)(j - a), f, g);
            if (f <= g) { sym_range(f, g, l, r); real_ok = r > l; }
        }
        if (!real_ok) return refine_partial(Q, j, a);
        if (a == floor || cut[a]) return a;
        while (true) {
            if (a == floor || cut[a] || !Q[a - 1]) return a;
            uint64_t b = a - 1; while (b > floor && hs[b]) --b;
            uint32_t jj = (uint32_t)(a - 1 - b);
            if (!hs[b]) {
                uint64_t sy = sym(Q[b], jj);
                uint64_t l2 = l, r2 = r;
                if (sy && bstep(sy, l2, r2)) { l = l2; r = r2; a = b; continue; }
            }
            // take as much of the block [b, a-1] as possible, from the right
            uint64_t best = a;
            uint64_t el = l, er = r; to_erows(el, er);
            for (uint64_t i = 1; i <= a - b; ++i) {
                uint64_t x = a - i;
                if (!Q[x]) break;
                bool ok = false;
                if (x != b || hs[b]) {              // a real Q[x] then i-1 hashes (for x = b with a real
                    uint64_t sy = sym(Q[x], (uint32_t)(i - 1));   //  start, that step just failed)
                    if (sy) { uint64_t l3 = l, r3 = r; ok = bstep(sy, l3, r3); }
                }
                bool hok = hstep(el, er);
                if (!ok && hok) ok = hcount(Q[x], el, er) > 0;
                if (!ok) break;
                best = x;
                if (!hok) break;                    // no longer hash rows to extend into
            }
            return best;
        }
    }
    // the smallest x in [a, j] such that Q[x..j] occurs, where x lies inside one block of the
    // query (no real start to extend from)
    uint64_t refine_partial(const std::vector<uint32_t>& Q, uint64_t j, uint64_t a) {
        uint64_t best = j;
        for (uint64_t x = j; x-- > a;) {
            if (!Q[x]) break;
            uint32_t len = (uint32_t)(j - x);
            bool ok = false;
            uint64_t f, g, l, r;
            syms_atleast(Q[x], len, f, g);
            if (f <= g) { sym_range(f, g, l, r); ok = r > l; }
            if (!ok) { uint64_t k = rows_atleast(len + 1); if (k) ok = hcount(Q[x], 1, 1 + k) > 0; }
            if (!ok) break;
            best = x;
        }
        return best;
    }

    uint64_t bytes() const {
        return top.size() * 4 + dstart.bytes() + jv.size() * 2 + bbwt.bytes() + hashbit.bytes() + H.bytes() + ncnt.size() * 8;
    }
    void save(FILE* f) const {
        fwrite(&sigma, 4, 1, f); fwrite(&sep, 4, 1, f); fwrite(&nsym, 8, 1, f); fwrite(&n, 8, 1, f);
        fwrite(&nb, 8, 1, f); fwrite(&nh, 8, 1, f); fwrite(&hb_c1, 8, 1, f);
        EliasFano::wv(f, top); dstart.save(f); EliasFano::wv(f, jv); bbwt.save(f); hashbit.save(f); H.save(f);
        EliasFano::wv(f, ncnt);
    }
    void load(FILE* f) {
        if (fread(&sigma, 4, 1, f) != 1 || fread(&sep, 4, 1, f) != 1 || fread(&nsym, 8, 1, f) != 1 ||
            fread(&n, 8, 1, f) != 1 || fread(&nb, 8, 1, f) != 1 || fread(&nh, 8, 1, f) != 1 ||
            fread(&hb_c1, 8, 1, f) != 1) throw std::runtime_error("read");
        EliasFano::rv(f, top); dstart.load(f); EliasFano::rv(f, jv); bbwt.load(f); hashbit.load(f); H.load(f);
        EliasFano::rv(f, ncnt);
    }
};
