// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Two-level indexing over the RLBWT of the reportable-shortcuts code.
//
// The text X (T, or the reversed collection T') is parsed at the starts of
// closed syncmers (k = 51, s = 12), and at position 0.  Phrase i is
// X[b_i .. b_{i+1} + k), so consecutive phrases overlap by k characters, and
// the last phrase is X[b_last .. n), which ends with $.  No phrase has a
// syncmer start strictly inside it, so the dictionary is prefix-free, and if
// phrases are numbered in the byte order of their strings ($ = 0x00 <
// # = 0x01 < A < C < G < N < T, as in the RLBWT), parse suffixes are in the
// same order as the text suffixes at the boundaries.  Hence:
//
//   B     marks the rows of the character-level BWT whose suffixes start at
//         boundaries; the r-th 1 (0-based) is parse row r + 1 (row 0 of the
//         parse is the sentinel).  Every occurrence of a syncmer is a
//         boundary, so B is the union of the BWT intervals of the distinct
//         syncmers (plus the row of position 0), and is stored run-length.
//   pbwt  the BWT of the parse (EFRank, no C array).
//
// Blocks (the third level).  succ(x) is the smallest phrase that follows x in
// the parse.  Phrase i is a continuation if it equals succ of phrase i-1.  A
// block is a phrase that is not a continuation, followed by its j >= 0
// continuations; its symbol is (d, j), ordered by d ascending and j
// descending, so the suffixes starting at blocks are in the same order in the
// block sequence as in the parse.
//
//   BS    marks the parse rows that start blocks (and row 0);
//   bbwt  the BWT of the block sequence, over the rows BS marks.
//
// In a search, a block step prepends a whole block of the pattern.  It finds
// only the occurrences in which the block's first phrase starts a block of the
// text, which holds for all occurrences once the phrase before it is matched
// too (it is not a continuation of that phrase); so the last block taken is
// redone with parse steps from the interval before it, which is exact.
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "efrank.hpp"
#include "rlbwt.hpp"
#include "search.hpp"

namespace tl {

constexpr int K = 51, SM = 12, W = K - SM + 1;

inline uint64_t mix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27; x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    return x;
}

// Calls f(q) for each start q of a closed syncmer in x[0 .. len), in order: the
// k-mer x[q .. q+K) is over A/C/G/T, and its smallest s-mer, by g = a hash of
// the smaller of the 2-bit codes of the s-mer and of its reverse, is its first
// or its last (ties included).  Since g is the same for an s-mer and its
// reverse, the syncmers of the reverse of a string are the reverses of its
// syncmers, so a pattern and its reverse are scanned only once.  The window minima
// are computed by the van Herk / Gil-Werman method, in chunks.
template <class F>
void for_each_syncmer(const uint8_t* x, uint64_t len, F f) {
    static int8_t code[256];
    static bool init = false;
    if (!init) {
        for (int c = 0; c < 256; ++c) code[c] = -1;
        code['A'] = 0; code['C'] = 1; code['G'] = 2; code['T'] = 3;
        init = true;
    }
    if (len < (uint64_t)K) return;
    const uint64_t nq = len - K + 1;                 // k-mer starts 0 .. nq-1
    const uint64_t CH = (uint64_t)W * 1024;          // k-mer starts per chunk (a multiple of W)
    const uint64_t mask = (1ull << (2 * SM)) - 1;
    std::vector<uint64_t> g(CH + W), pre(CH + W), suf(CH + W);
    for (uint64_t c0 = 0; c0 < nq; c0 += CH) {
        const uint64_t cq = std::min(CH, nq - c0);   // k-mer starts c0 .. c0+cq-1
        const uint64_t ns = cq + W - 1;              // s-mer starts c0 .. c0+ns-1
        // g over the s-mers: 0 if the s-mer has a character other than A/C/G/T
        uint64_t smer = 0, rsmer = 0;                // codes of the s-mer and of its reverse
        int64_t bad = -1;                            // last bad character, relative to c0
        for (int64_t i = 0; i < SM - 1; ++i) {
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint64_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint64_t)cc << (2 * (SM - 1)));
        }
        for (uint64_t t = 0; t < ns; ++t) {
            int64_t i = (int64_t)t + SM - 1;
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint64_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint64_t)cc << (2 * (SM - 1)));
            uint64_t cm = smer < rsmer ? smer : rsmer;
            g[t] = bad >= (int64_t)t ? 0 : (mix64(cm) >> 1) + 1;
        }
        // block prefix and suffix minima, blocks of W s-mers
        for (uint64_t b0 = 0; b0 < ns; b0 += W) {
            uint64_t b1 = std::min(ns, b0 + W);
            uint64_t mn = ~0ull;
            for (uint64_t t = b0; t < b1; ++t) { mn = g[t] < mn ? g[t] : mn; pre[t] = mn; }
            mn = ~0ull;
            for (uint64_t t = b1; t-- > b0;) { mn = g[t] < mn ? g[t] : mn; suf[t] = mn; }
        }
        for (uint64_t q = 0; q < cq; ++q) {
            uint64_t a = suf[q], bb = pre[q + W - 1];
            uint64_t mn = a < bb ? a : bb;
            if (mn && (mn == g[q] || mn == g[q + W - 1])) f(c0 + q);
        }
    }
}

// Hash of a byte string (never 0).
inline uint64_t phrase_hash(const uint8_t* p, uint64_t len) {
    uint64_t hv = mix64(len ^ 0x9e3779b97f4a7c15ull);
    uint64_t i = 0;
    for (; i + 8 <= len; i += 8) { uint64_t w; std::memcpy(&w, p + i, 8); hv = mix64(hv ^ w) + i; }
    uint64_t w = 0;
    std::memcpy(&w, p + i, len - i);
    hv = mix64(hv ^ w ^ (len << 56));
    return hv ? hv : 1;
}

// A bitvector stored as runs of 1s: rank and select by Elias-Fano lists of
// the run starts and of the number of 1s before each run.
struct RLBits {
    uint64_t n = 0, ones = 0, nruns = 0;
    EliasFano RS, RC;
    void build(const std::vector<std::pair<uint64_t, uint64_t>>& runs, uint64_t size) {   // (start, length)
        n = size; nruns = runs.size();
        std::vector<uint64_t> s(nruns), c(nruns + 1);
        uint64_t acc = 0;
        for (uint64_t k = 0; k < nruns; ++k) { s[k] = runs[k].first; c[k] = acc; acc += runs[k].second; }
        c[nruns] = acc; ones = acc;
        RS.build(s, n + 1); RC.build(c, ones + 1);
    }
    inline uint64_t rank1(uint64_t pos) const {        // 1s in [0, pos)
        uint64_t k = RS.count_less(pos);
        if (k == 0) return 0;
        --k;
        uint64_t s = RS.access(k), c0, c1;
        RC.access2(k, c0, c1);
        uint64_t d = pos - s;
        return c0 + (d < c1 - c0 ? d : c1 - c0);
    }
    inline uint64_t select1(uint64_t r) const {        // position of the r-th 1 (0-based)
        uint64_t k = RC.count_less(r + 1) - 1;
        return RS.access(k) + (r - RC.access(k));
    }
    uint64_t bytes() const { return RS.bytes() + RC.bytes() + 24; }
    void save(FILE* f) const { fwrite(&n, 8, 1, f); fwrite(&ones, 8, 1, f); fwrite(&nruns, 8, 1, f); RS.save(f); RC.save(f); }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&ones, 8, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1) throw std::runtime_error("read");
        RS.load(f); RC.load(f);
    }
};

// Fixed-width packed integers.
struct PackedInts {
    uint64_t n = 0;
    uint32_t w = 1;
    std::vector<uint64_t> d;
    void build(const std::vector<uint64_t>& v) {
        n = v.size();
        uint64_t mx = 0;
        for (auto x : v) mx = std::max(mx, x);
        w = 1;
        while (w < 64 && (mx >> w)) ++w;
        d.assign((n * w + 63) / 64 + 1, 0);
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t bpos = i * w, x = v[i];
            d[bpos >> 6] |= x << (bpos & 63);
            if ((bpos & 63) + w > 64) d[(bpos >> 6) + 1] |= x >> (64 - (bpos & 63));
        }
    }
    inline uint64_t get(uint64_t i) const {
        uint64_t bpos = i * w, k = bpos >> 6, o = bpos & 63;
        uint64_t x = d[k] >> o;
        if (o + w > 64) x |= d[k + 1] << (64 - o);
        return w == 64 ? x : x & ((1ull << w) - 1);
    }
    uint64_t bytes() const { return d.size() * 8 + 12; }
    void save(FILE* f) const { fwrite(&n, 8, 1, f); fwrite(&w, 4, 1, f); EliasFano::wv(f, d); }
    void load(FILE* f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&w, 4, 1, f) != 1) throw std::runtime_error("read");
        EliasFano::rv(f, d);
    }
};

struct Level {
    uint64_t n = 0;          // rows of the character-level BWT
    uint64_t m = 0;          // phrases; the parse has m + 1 rows
    uint32_t sigma = 0;      // distinct phrases, ids 1..sigma
    uint64_t nblocks = 0;    // blocks; the block BWT has nblocks + 1 rows
    uint64_t nsym = 0;       // distinct block symbols, ids 1..nsym
    RLBits B;
    EFRank pbwt;
    EliasFano BS;
    EFRank bbwt;
    PackedInts succp;               // succ(x), x = 0..sigma (0: none)
    EliasFano dstart;               // sigma + 2 values: symbols of d are dstart[d] .. dstart[d+1]-1 (0-based)
    PackedInts jv;                  // j of each symbol, descending within each d
    std::vector<uint64_t> ccont;    // bit x: phrase x is a continuation somewhere in the parse
    // Phrase dictionary for queries: the 64-bit hashes sorted, in buckets by their
    // top hb bits; each slot stores the other 64 - hb bits (so the whole hash is
    // checked) and the phrase id.
    uint32_t hb = 1;
    EliasFano hbucket;              // 2^hb + 1 values: first slot of each bucket
    PackedInts hlow, hid;

    inline uint32_t succ(uint32_t x) const { return (uint32_t)succp.get(x); }
    inline bool can_continue(uint32_t x) const { return (ccont[x >> 6] >> (x & 63)) & 1; }
    inline uint32_t sym(uint32_t d, uint64_t j) const {
        uint64_t a0, b0;
        dstart.access2(d, a0, b0);
        uint64_t a = a0, b = b0;
        while (a < b) {                                 // jv descending
            uint64_t mid = (a + b) >> 1;
            if (jv.get(mid) > j) a = mid + 1; else b = mid;
        }
        return (a < b0 && jv.get(a) == j) ? (uint32_t)a + 1 : 0;
    }
    inline uint32_t lookup(uint64_t hv) const {
        uint64_t bk = hv >> (64 - hb), lo = hb == 64 ? 0 : hv & ((1ull << (64 - hb)) - 1);
        uint64_t a, b;
        hbucket.access2(bk, a, b);
        for (; a < b; ++a) {
            uint64_t x = hlow.get(a);
            if (x == lo) return (uint32_t)hid.get(a);
            if (x > lo) return 0;
        }
        return 0;
    }

    // From plain arrays, at build time.
    void set_succ(const std::vector<uint32_t>& s) { succp.build(std::vector<uint64_t>(s.begin(), s.end())); }
    void set_symbols(const std::vector<uint64_t>& ds, const std::vector<uint64_t>& j) {
        dstart.build(ds, nsym + 1); jv.build(j);
    }
    void set_dictionary(const std::vector<uint64_t>& phash) {      // phash[id], id = 1..sigma
        hb = 1;
        while (hb < 32 && (1ull << hb) < sigma) ++hb;
        std::vector<std::pair<uint64_t, uint32_t>> v;
        v.reserve(sigma);
        for (uint32_t id = 1; id <= sigma; ++id) v.push_back({phash[id], id});
        std::sort(v.begin(), v.end());
        for (size_t k = 1; k < v.size(); ++k)
            if (v[k].first == v[k - 1].first) throw std::runtime_error("two phrases with the same hash");
        uint64_t nb = 1ull << hb;
        std::vector<uint64_t> st(nb + 1), lw(sigma), ids(sigma);
        uint64_t k = 0;
        for (uint64_t bk = 0; bk < nb; ++bk) {
            st[bk] = k;
            while (k < v.size() && (v[k].first >> (64 - hb)) == bk) ++k;
        }
        st[nb] = k;
        for (size_t t = 0; t < v.size(); ++t) { lw[t] = v[t].first & ((1ull << (64 - hb)) - 1); ids[t] = v[t].second; }
        hbucket.build(st, sigma + 1); hlow.build(lw); hid.build(ids);
    }

    uint64_t table_bytes() const { return hbucket.bytes() + hlow.bytes() + hid.bytes(); }
    uint64_t succ_bytes() const { return succp.bytes() + ccont.size() * 8; }
    uint64_t sym_bytes() const { return dstart.bytes() + jv.bytes(); }

    void save(const std::string& path) const {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write " + path);
        fwrite("TWOLVL03", 1, 8, f);
        fwrite(&n, 8, 1, f); fwrite(&m, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nblocks, 8, 1, f); fwrite(&nsym, 8, 1, f);
        B.save(f); pbwt.save(f); BS.save(f); bbwt.save(f);
        succp.save(f); dstart.save(f); jv.save(f); EliasFano::wv(f, ccont);
        fwrite(&hb, 4, 1, f); hbucket.save(f); hlow.save(f); hid.save(f);
        fclose(f);
    }
    void load(const std::string& path) {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path);
        char mg[8];
        if (fread(mg, 1, 8, f) != 8 || std::memcmp(mg, "TWOLVL03", 8))
            throw std::runtime_error(path + ": not a two-level file of this version (rebuild it)");
        if (fread(&n, 8, 1, f) != 1 || fread(&m, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 ||
            fread(&nblocks, 8, 1, f) != 1 || fread(&nsym, 8, 1, f) != 1) throw std::runtime_error("read");
        B.load(f); pbwt.load(f); BS.load(f); bbwt.load(f);
        succp.load(f); dstart.load(f); jv.load(f); EliasFano::rv(f, ccont);
        if (fread(&hb, 4, 1, f) != 1) throw std::runtime_error("read");
        hbucket.load(f); hlow.load(f); hid.load(f);
        fclose(f);
    }
};

// The parse of a pattern against a level's dictionary.
struct PatInfo {
    std::vector<int64_t> S;         // syncmer starts
    std::vector<uint32_t> pid;      // pid[i]: id of P[S[i] .. S[i+1]+K), 0 if not a phrase (or i is last)
    std::vector<uint8_t> cont;      // cont[i]: phrase i is a continuation of phrase i-1
    std::vector<int32_t> bstart;    // bstart[i]: the block of phrase i starts at phrase bstart[i]
    // If mirror is given, it is the PatInfo of the reverse of P, whose syncmers
    // are reused.
    void build(const Level& L, const std::string& P, const PatInfo* mirror = nullptr) {
        S.clear();
        const uint8_t* x = reinterpret_cast<const uint8_t*>(P.data());
        if (mirror) {
            const int64_t m = (int64_t)P.size();
            S.resize(mirror->S.size());
            for (size_t i = 0, z = S.size(); i < z; ++i) S[i] = m - K - mirror->S[z - 1 - i];
        } else {
            for_each_syncmer(x, P.size(), [&](uint64_t q) { S.push_back((int64_t)q); });
        }
        size_t z = S.size();
        pid.assign(z, 0); cont.assign(z, 0); bstart.assign(z, 0);
        for (size_t i = 0; i + 1 < z; ++i) pid[i] = L.lookup(phrase_hash(x + S[i], S[i + 1] + K - S[i]));
        for (size_t i = 0; i < z; ++i) {
            cont[i] = i > 0 && pid[i] && pid[i - 1] && L.succ(pid[i - 1]) == pid[i];
            bstart[i] = cont[i] ? bstart[i - 1] : (int32_t)i;
        }
    }
};

struct Counters {
    int64_t searches = 0, char_steps = 0, parse_steps = 0, block_steps = 0;
};

// The longest suffix of P[begin .. end) that occurs in X, and its interval in
// the character-level BWT: the same result as rs::search_range without
// shortcuts.  With jumps = false, only the parse is used.
inline rs::SearchResult search2(const rs::RLBWT& rl, const Level& L, const PatInfo& pi, const std::string& P,
                                int64_t begin, int64_t end, bool jumps, Counters* cnt) {
    rs::SearchResult r;
    r.pattern_length = end - begin;
    int64_t lo = 0, hi = rl.size() - 1, j = end;
    int64_t csteps = 0, psteps = 0, bsteps = 0;
    auto char_steps_to = [&](int64_t stop) -> bool {      // false if a step fails
        while (j > stop) {
            int64_t nl = lo, nh = hi;
            rl.backward_step(static_cast<uint8_t>(P[(size_t)(j - 1)]), nl, nh);
            ++csteps;
            if (nh < nl) return false;
            lo = nl; hi = nh; --j;
        }
        return true;
    };
    const auto& S = pi.S;
    int64_t idx = (int64_t)(std::upper_bound(S.begin(), S.end(), end - K) - S.begin()) - 1;
    if (idx >= 0 && S[idx] >= begin && char_steps_to(S[idx])) {
        // parse level: rows [plo, phi] of the parse, for P[S[i] .. end)
        uint64_t plo = L.B.rank1((uint64_t)lo) + 1, phi = L.B.rank1((uint64_t)hi + 1);
        int64_t i = idx;
        int64_t noblock_above = idx - 1;      // block steps only from boundaries i <= this
        while (i > 0 && S[i - 1] >= begin && pi.pid[i - 1]) {
            if (jumps && i <= noblock_above && !pi.cont[i]) {
                uint64_t blo = L.BS.count_less(plo), bhi = L.BS.count_less(phi + 1);   // block rows [blo, bhi)
                int64_t cur = i, taken = 0;
                uint64_t pblo = 0, pbhi = 0; int64_t pcur = i;
                while (cur > 0 && blo < bhi) {
                    int32_t t0 = pi.bstart[cur - 1];
                    if (S[t0] < begin || !pi.pid[t0]) break;
                    uint32_t y = L.sym(pi.pid[t0], (uint64_t)(cur - 1 - t0));
                    if (!y) break;
                    uint64_t o1, o2;
                    L.bbwt.lf2(y, blo, bhi, o1, o2);
                    ++bsteps;
                    if (o1 >= o2) break;
                    pblo = blo; pbhi = bhi; pcur = cur;
                    blo = o1; bhi = o2; cur = t0; ++taken;
                }
                if (taken >= 1 && !L.can_continue(pi.pid[cur])) {
                    // the last block's first phrase is never a continuation, so every
                    // occurrence has it at the start of a block: the interval is exact
                    plo = L.BS.access(blo); phi = L.BS.access(bhi - 1);
                    i = cur;
                    noblock_above = i;
                    continue;
                }
                if (taken >= 2) {             // back up to the interval before the last block
                    plo = L.BS.access(pblo); phi = L.BS.access(pbhi - 1);
                    i = pcur;
                }
                // taken <= 1: stay at the entry interval, which is exact
                noblock_above = pi.bstart[i - 1];
            }
            uint64_t o1, o2;
            L.pbwt.lf2(pi.pid[i - 1], plo, phi + 1, o1, o2);
            ++psteps;
            if (o1 >= o2) break;
            plo = o1; phi = o2 - 1; --i;
        }
        lo = (int64_t)L.B.select1(plo - 1); hi = (int64_t)L.B.select1(phi - 1);
        j = S[i];
        char_steps_to(begin);
    } else if (idx < 0 || S[idx] < begin) {
        char_steps_to(begin);
    }
    r.lo = lo; r.hi = hi; r.matched = end - j;
    r.backward_steps = csteps;
    if (cnt) { ++cnt->searches; cnt->char_steps += csteps; cnt->parse_steps += psteps; cnt->block_steps += bsteps; }
    return r;
}

}  // namespace tl
