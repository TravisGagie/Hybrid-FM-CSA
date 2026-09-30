// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Huffman-shaped wavelet trees.
//
// Backward search only needs rank_c, never order between symbols, so a Huffman
// shape (optimal expected depth) is used rather than Hu-Tucker.
//
//   HuffWT       : Huffman-shaped wavelet tree over the whole sequence.
//                  n(H0 + 1) bits plus ~12.5% for rank; rank_c costs |code(c)|
//                  bitvector ranks.
//   RLHuffRank   : RLFM-style. A Huffman wavelet tree over the r run heads L',
//                  Elias-Fano run starts (B_L), and Elias-Fano F-positions of the
//                  runs in F'-order (B_F).  rank_c(i) = one predecessor query on
//                  B_L + one descent of the tree (at two adjacent positions at
//                  once) + one or two accesses to B_F.
#pragma once
#include "efrank.hpp"
#include <queue>

// Bitvector with rank1; per 512-bit block: [cumulative ones][8 data words].
struct RankBV {
    std::vector<uint64_t> w;        // 9 words per block
    uint64_t n = 0;
    void build(const std::vector<uint64_t>& bits, uint64_t nbits) {
        n = nbits;
        uint64_t nb = nbits / 512 + 1;
        w.assign(nb * 9, 0);
        uint64_t cum = 0;
        for (uint64_t b = 0; b < nb; ++b) {
            w[b * 9] = cum;
            for (int k = 0; k < 8; ++k) {
                uint64_t idx = b * 8 + k;
                uint64_t x = idx < bits.size() ? bits[idx] : 0;
                w[b * 9 + 1 + k] = x;
                cum += __builtin_popcountll(x);
            }
        }
    }
    inline uint64_t rank1(uint64_t i) const {            // ones in [0, i)
        uint64_t b = i >> 9, k = (i >> 6) & 7, o = i & 63;
        const uint64_t* p = &w[b * 9];
        uint64_t r = p[0];
        for (uint64_t t = 0; t < k; ++t) r += __builtin_popcountll(p[1 + t]);
        if (o) r += __builtin_popcountll(p[1 + k] << (64 - o));
        return r;
    }
    uint64_t bytes() const { return w.size() * 8; }
    void save(FILE* f) const { fwrite(&n, 8, 1, f); EliasFano::wv(f, w); }
    void load(FILE* f) { if (fread(&n, 8, 1, f) != 1) throw std::runtime_error("read"); EliasFano::rv(f, w); }
};

struct HuffWT {
    uint64_t n = 0;
    uint32_t sigma = 0;
    std::vector<uint64_t> code;     // per symbol, bit d = branch at depth d
    std::vector<uint8_t> len;       // 255 = symbol absent
    std::vector<uint32_t> child;    // 2 per internal node
    std::vector<uint64_t> off;      // start of node's bits in the concatenation
    std::vector<uint64_t> ones0;    // ones before node start
    RankBV bv;
    uint32_t root = 0;
    uint64_t sum_depth = 0;         // total bits

    void build(const std::vector<uint32_t>& s, uint32_t sig) {
        n = s.size(); sigma = sig;
        std::vector<uint64_t> f(sigma, 0);
        for (auto c : s) f[c]++;
        // Huffman tree; leaves are ids 0..sigma-1, internal nodes sigma..
        using P = std::pair<uint64_t, uint32_t>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        for (uint32_t c = 0; c < sigma; ++c) if (f[c]) pq.push({f[c], c});
        std::vector<uint32_t> L, R;             // children of internal node sigma+k
        code.assign(sigma, 0); len.assign(sigma, 255);
        uint32_t nint = 0;
        if (pq.size() == 1) { len[pq.top().second] = 0; }
        while (pq.size() > 1) {
            P a = pq.top(); pq.pop(); P b = pq.top(); pq.pop();
            L.push_back(a.second); R.push_back(b.second);
            pq.push({a.first + b.first, sigma + nint++});
        }
        child.assign(2 * (size_t)nint, 0);
        if (nint) {
            root = 0;
            // renumber internal nodes top-down (BFS from root), assign codes
            std::vector<uint32_t> newid(nint, 0);
            std::vector<uint32_t> order{nint - 1};
            newid[nint - 1] = 0;
            std::vector<std::pair<uint64_t, uint8_t>> pc(nint, {0, 0});   // code prefix, depth
            for (size_t q = 0; q < order.size(); ++q) {
                uint32_t k = order[q];
                for (int b = 0; b < 2; ++b) {
                    uint32_t ch = b ? R[k] : L[k];
                    uint64_t cd = pc[k].first | ((uint64_t)b << pc[k].second);
                    uint8_t dp = pc[k].second + 1;
                    if (ch < sigma) { code[ch] = cd; len[ch] = dp; child[2 * newid[k] + b] = ch; }
                    else {
                        uint32_t kk = ch - sigma;
                        newid[kk] = order.size(); order.push_back(kk);
                        pc[kk] = {cd, dp};
                        child[2 * newid[k] + b] = sigma + newid[kk];
                    }
                }
            }
        }
        // node bit counts, offsets
        std::vector<uint64_t> cnt(nint, 0);
        sum_depth = 0;
        for (uint32_t c = 0; c < sigma; ++c) if (f[c] && len[c] != 255) {
            uint32_t v = 0;
            for (uint8_t d = 0; d < len[c]; ++d) {
                cnt[v] += f[c];
                uint32_t ch = child[2 * v + ((code[c] >> d) & 1)];
                v = ch - sigma;
            }
            sum_depth += f[c] * len[c];
        }
        off.assign(nint + 1, 0);
        for (uint32_t v = 0; v < nint; ++v) off[v + 1] = off[v] + cnt[v];
        std::vector<uint64_t> bits(off[nint] / 64 + 1, 0), pos(off.begin(), off.end());
        for (auto c : s) {
            uint32_t v = 0;
            for (uint8_t d = 0; d < len[c]; ++d) {
                uint64_t b = (code[c] >> d) & 1;
                uint64_t p = pos[v]++;
                if (b) bits[p >> 6] |= 1ull << (p & 63);
                v = child[2 * v + b] - sigma;
            }
        }
        bv.build(bits, off[nint]);
        ones0.assign(nint + 1, 0);
        for (uint32_t v = 0; v <= nint; ++v) ones0[v] = bv.rank1(off[v]);
    }

    inline uint64_t rank(uint32_t c, uint64_t i) const {
        if (c >= sigma || len[c] == 255) return 0;
        uint32_t v = 0;
        uint64_t cd = code[c];
        for (uint8_t d = 0; d < len[c]; ++d) {
            uint64_t o = bv.rank1(off[v] + i) - ones0[v];
            uint64_t b = cd & 1; cd >>= 1;
            i = b ? o : i - o;
            v = child[2 * v + b] - sigma;
        }
        return i;
    }
    // ranks at two positions in one descent
    inline void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t& r1, uint64_t& r2) const {
        if (c >= sigma || len[c] == 255) { r1 = r2 = 0; return; }
        uint32_t v = 0;
        uint64_t cd = code[c];
        for (uint8_t d = 0; d < len[c]; ++d) {
            uint64_t o1 = bv.rank1(off[v] + i1) - ones0[v];
            uint64_t o2 = bv.rank1(off[v] + i2) - ones0[v];
            uint64_t b = cd & 1; cd >>= 1;
            i1 = b ? o1 : i1 - o1;
            i2 = b ? o2 : i2 - o2;
            v = child[2 * v + b] - sigma;
        }
        r1 = i1; r2 = i2;
    }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&sum_depth, 8, 1, f);
        EliasFano::wv(f, code); EliasFano::wv(f, len); EliasFano::wv(f, child); EliasFano::wv(f, off);
        EliasFano::wv(f, ones0); bv.save(f);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&sum_depth, 8, 1, f) != 1)
            throw std::runtime_error("read");
        EliasFano::rv(f, code); EliasFano::rv(f, len); EliasFano::rv(f, child); EliasFano::rv(f, off);
        EliasFano::rv(f, ones0); bv.load(f); root = 0;
    }
    double avg_depth() const { return n ? (double)sum_depth / n : 0; }
    uint64_t bytes() const {
        return bv.bytes() + code.size() * 8 + len.size() + child.size() * 4 + off.size() * 8 + ones0.size() * 8;
    }
};

struct RLHuffRank {
    uint64_t n = 0;
    uint32_t sigma = 0;
    uint64_t nruns = 0;
    HuffWT heads;                   // L'
    EliasFano BL;                   // run start positions
    EliasFano FP;                   // F positions of runs in F'-order, plus sentinel n
    std::vector<uint32_t> first;    // sigma+1: index in F'-order of each symbol's first run

    void build(const std::vector<uint32_t>& s, uint32_t sig) {
        n = s.size(); sigma = sig;
        std::vector<uint32_t> Lp;
        std::vector<uint64_t> st, cnt(sigma + 1, 0), runsof(sigma + 1, 0);
        for (uint64_t i = 0; i < n; ++i) {
            cnt[s[i]]++;
            if (i == 0 || s[i] != s[i - 1]) { Lp.push_back(s[i]); st.push_back(i); runsof[s[i]]++; }
        }
        nruns = Lp.size();
        first.assign(sigma + 1, 0);
        std::vector<uint64_t> Cf(sigma + 1, 0);
        for (uint32_t c = 0; c < sigma; ++c) { first[c + 1] = first[c] + runsof[c]; Cf[c + 1] = Cf[c] + cnt[c]; }
        std::vector<uint64_t> fp(nruns + 1);
        std::vector<uint64_t> slot(first.begin(), first.end()), seen(sigma, 0);
        for (uint64_t k = 0; k < nruns; ++k) {
            uint32_t c = Lp[k];
            uint64_t len = (k + 1 < nruns ? st[k + 1] : n) - st[k];
            fp[slot[c]++] = Cf[c] + seen[c];
            seen[c] += len;
        }
        fp[nruns] = n;
        heads.build(Lp, sigma);
        BL.build(st, n);
        FP.build(fp, n + 1);
    }

    // occurrences of c in s[0, i)
    inline uint64_t rank(uint32_t c, uint64_t i) const {
        if (c >= sigma || i == 0 || first[c] == first[c + 1]) return 0;
        uint64_t j = BL.count_less(i) - 1;                 // run containing i-1
        uint64_t q0, q1;
        heads.rank2(c, j, j + 1, q0, q1);
        uint64_t f0 = FP.access(first[c]);
        uint64_t fq = FP.access(first[c] + q0);
        if (q1 > q0) return fq - f0 + (i - BL.access(j));  // i-1 lies in a run of c
        return fq - f0;
    }
    inline void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t& r1, uint64_t& r2) const {
        r1 = rank(c, i1); r2 = rank(c, i2);
    }
    uint64_t bytes() const { return heads.bytes() + BL.bytes() + FP.bytes() + first.size() * 4; }
    void save(FILE* f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nruns, 8, 1, f);
        heads.save(f); BL.save(f); FP.save(f); EliasFano::wv(f, first);
    }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1)
            throw std::runtime_error("read");
        heads.load(f); BL.load(f); FP.load(f); EliasFano::rv(f, first);
    }
};
