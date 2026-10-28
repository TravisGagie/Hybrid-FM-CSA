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
// Blocks (the third level).  succ(x) is the continuation of phrase x: by default
// the smallest phrase that follows x in the parse, or (built with --continuation
// frequent) the most frequent one.  Phrase i is a continuation if it equals succ
// of phrase i-1.  A block is a phrase d that is not a continuation, followed by
// its j >= 0 continuations; its symbol is (d, j, side), where side says whether
// the first phrase of the next block is larger (+) or smaller (-) than succ of
// the block's last phrase.  Symbols are ordered by d, then all - before all +,
// the -'s by j ascending and the +'s by j descending, so the suffixes starting
// at blocks are in the same order in the block sequence as in the parse.  With
// the smallest successor every side is + (the next block's first phrase is
// larger than succ), and the sides are not stored.
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
    const uint64_t CH = (uint64_t)W * 256;           // k-mer starts per chunk (a multiple of W)
    const uint32_t mask = (1u << (2 * SM)) - 1;
    std::vector<uint32_t> g(CH + W), pre(CH + W), suf(CH + W);
    for (uint64_t c0 = 0; c0 < nq; c0 += CH) {
        const uint64_t cq = std::min(CH, nq - c0);   // k-mer starts c0 .. c0+cq-1
        const uint64_t ns = cq + W - 1;              // s-mer starts c0 .. c0+ns-1
        // g over the s-mers: 0 if the s-mer has a character other than A/C/G/T
        uint32_t smer = 0, rsmer = 0;                // codes of the s-mer and of its reverse
        int64_t bad = -1;                            // last bad character, relative to c0
        for (int64_t i = 0; i < SM - 1; ++i) {
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint32_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint32_t)cc << (2 * (SM - 1)));
        }
        for (uint64_t t = 0; t < ns; ++t) {
            int64_t i = (int64_t)t + SM - 1;
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint32_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint32_t)cc << (2 * (SM - 1)));
            uint32_t h = smer < rsmer ? smer : rsmer;   // murmur3 fmix32 (a bijection) of the smaller code
            h ^= h >> 16; h *= 0x85ebca6bu; h ^= h >> 13; h *= 0xc2b2ae35u; h ^= h >> 16;
            g[t] = bad >= (int64_t)t ? 0u : (h >> 1) + 1u;
        }
        // block prefix and suffix minima, blocks of W s-mers
        for (uint64_t b0 = 0; b0 < ns; b0 += W) {
            uint64_t b1 = std::min(ns, b0 + W);
            uint32_t mn = ~0u;
            for (uint64_t t = b0; t < b1; ++t) { mn = g[t] < mn ? g[t] : mn; pre[t] = mn; }
            mn = ~0u;
            for (uint64_t t = b1; t-- > b0;) { mn = g[t] < mn ? g[t] : mn; suf[t] = mn; }
        }
        for (uint64_t q = 0; q < cq; ++q) {
            uint32_t a = suf[q], bb = pre[q + W - 1];
            uint32_t mn = a < bb ? a : bb;
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

// The phrase dictionary for queries, shared by the two directions: the 64-bit
// hashes of T's phrases, sorted and split into buckets by their top hb bits
// (about 8 per bucket); each slot stores the other 64 - hb bits (so the whole
// hash is checked), the phrase's id in T, and the id of its reverse in the
// reversed collection.
struct PhraseDict {
    uint32_t hb = 1;
    uint64_t sigma = 0;
    std::vector<uint32_t> bstart;   // 2^hb + 1: first slot of each bucket
    PackedInts low, idf, idr;
    void build(const std::vector<uint64_t>& hf, const std::vector<uint64_t>& revid, uint64_t sig) {   // indexed by id 1..sig
        sigma = sig;
        hb = 1;
        while (hb < 32 && (8ull << hb) < sigma) ++hb;
        std::vector<std::pair<uint64_t, uint32_t>> v;
        v.reserve(sigma);
        for (uint32_t id = 1; id <= sigma; ++id) v.push_back({hf[id], id});
        std::sort(v.begin(), v.end());
        for (size_t k = 1; k < v.size(); ++k)
            if (v[k].first == v[k - 1].first) throw std::runtime_error("two phrases with the same hash");
        uint64_t nb = 1ull << hb;
        bstart.assign(nb + 1, 0);
        std::vector<uint64_t> lw(sigma), f(sigma), r(sigma);
        uint64_t k = 0;
        for (uint64_t bk = 0; bk < nb; ++bk) {
            bstart[bk] = (uint32_t)k;
            while (k < v.size() && (v[k].first >> (64 - hb)) == bk) ++k;
        }
        bstart[nb] = (uint32_t)k;
        for (size_t t = 0; t < v.size(); ++t) {
            lw[t] = v[t].first & ((1ull << (64 - hb)) - 1); f[t] = v[t].second; r[t] = revid[v[t].second];
        }
        low.build(lw); idf.build(f); idr.build(r);
    }
    inline uint64_t lowbits(uint64_t hv) const { return hv & ((1ull << (64 - hb)) - 1); }
    inline void prefetch(uint64_t hv) const {
        uint64_t a = bstart[hv >> (64 - hb)];
        __builtin_prefetch(&low.d[(a * low.w) >> 6]);
    }
    // slot of hv, or ~0
    inline uint64_t slot(uint64_t hv) const {
        uint64_t bk = hv >> (64 - hb), lo = lowbits(hv);
        for (uint64_t a = bstart[bk], b = bstart[bk + 1]; a < b; ++a) {
            uint64_t x = low.get(a);
            if (x == lo) return a;
            if (x > lo) break;
        }
        return ~0ull;
    }
    uint64_t bytes() const { return bstart.size() * 4 + low.bytes() + idf.bytes() + idr.bytes(); }
    void save(FILE* f) const { fwrite(&hb, 4, 1, f); fwrite(&sigma, 8, 1, f); EliasFano::wv(f, bstart); low.save(f); idf.save(f); idr.save(f); }
    void load(FILE* f) {
        if (fread(&hb, 4, 1, f) != 1 || fread(&sigma, 8, 1, f) != 1) throw std::runtime_error("read");
        EliasFano::rv(f, bstart); low.load(f); idf.load(f); idr.load(f);
    }
};

// The block-start rows: Elias-Fano or run-length coded, whichever is smaller.
struct BlockStarts {
    uint8_t rl = 0;
    uint64_t alt_bytes = 0;         // size of the encoding not chosen (at build time)
    EliasFano ef;
    RLBits r;
    void build(const std::vector<uint64_t>& rows, uint64_t universe) {
        std::vector<std::pair<uint64_t, uint64_t>> runs;
        for (uint64_t x : rows) {
            if (!runs.empty() && runs.back().first + runs.back().second == x) ++runs.back().second;
            else runs.push_back({x, 1});
        }
        r.build(runs, universe);
        ef.build(rows, universe);
        rl = r.bytes() < ef.bytes();
        alt_bytes = rl ? ef.bytes() : r.bytes();
        if (rl) ef = EliasFano(); else r = RLBits();
    }
    inline uint64_t count_less(uint64_t x) const { return rl ? r.rank1(x) : ef.count_less(x); }
    inline uint64_t access(uint64_t k) const { return rl ? r.select1(k) : ef.access(k); }
    uint64_t bytes() const { return rl ? r.bytes() : ef.bytes(); }
    void save(FILE* f) const { fwrite(&rl, 1, 1, f); if (rl) r.save(f); else ef.save(f); }
    void load(FILE* f) {
        if (fread(&rl, 1, 1, f) != 1) throw std::runtime_error("read");
        if (rl) r.load(f); else ef.load(f);
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
    BlockStarts BS;
    EFRank bbwt;
    PackedInts succp;               // succ(x), x = 0..sigma (0: none)
    EliasFano dstart;               // sigma + 2 values: symbols of d are dstart[d] .. dstart[d+1]-1 (0-based)
    PackedInts jv;                  // j of each symbol, in symbol order within each d
    uint8_t sided = 0;              // 1: symbols have sides (stored in sidev); 0: all +
    std::vector<uint64_t> sidev;    // bit k: symbol k + 1 has side +
    std::vector<uint64_t> ccont;    // bit x: phrase x is a continuation somewhere in the parse
    PhraseDict dict;                // only in the level of T (shared by both directions)

    inline uint32_t succ(uint32_t x) const { return (uint32_t)succp.get(x); }
    inline bool can_continue(uint32_t x) const { return (ccont[x >> 6] >> (x & 63)) & 1; }
    inline bool sym_plus(uint64_t k) const { return !sided || ((sidev[k >> 6] >> (k & 63)) & 1); }
    // the order of (j, side) within a d: - by j ascending, then + by j descending
    static inline uint64_t sym_key(uint64_t j, bool plus) { return plus ? (1ull << 32) | (0xffffffffull - j) : j; }
    // the symbol (d, j, side) (1-based), or 0 if it does not occur
    inline uint32_t sym(uint32_t d, uint64_t j, bool plus) const {
        if (!sided && !plus) return 0;
        uint64_t a0, b0;
        dstart.access2(d, a0, b0);
        const uint64_t want = sym_key(j, plus);
        uint64_t a = a0, b = b0;
        while (a < b) {
            uint64_t mid = (a + b) >> 1;
            if (sym_key(jv.get(mid), sym_plus(mid)) < want) a = mid + 1; else b = mid;
        }
        return (a < b0 && sym_key(jv.get(a), sym_plus(a)) == want) ? (uint32_t)a + 1 : 0;
    }
    // From plain arrays, at build time.
    void set_succ(const std::vector<uint32_t>& s) { succp.build(std::vector<uint64_t>(s.begin(), s.end())); }
    void set_symbols(const std::vector<uint64_t>& ds, const std::vector<uint64_t>& j,
                     const std::vector<uint8_t>& plus, bool with_sides) {
        dstart.build(ds, nsym + 1); jv.build(j);
        sided = with_sides;
        sidev.assign(with_sides ? j.size() / 64 + 1 : 0, 0);
        if (with_sides) for (size_t k = 0; k < j.size(); ++k) if (plus[k]) sidev[k >> 6] |= 1ull << (k & 63);
    }
    uint64_t table_bytes() const { return dict.bytes(); }
    uint64_t succ_bytes() const { return succp.bytes() + ccont.size() * 8; }
    uint64_t sym_bytes() const { return dstart.bytes() + jv.bytes() + sidev.size() * 8; }

    void save(const std::string& path) const {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write " + path);
        fwrite(sided ? "TWOLVL06" : "TWOLVL05", 1, 8, f);
        fwrite(&n, 8, 1, f); fwrite(&m, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nblocks, 8, 1, f); fwrite(&nsym, 8, 1, f);
        B.save(f); pbwt.save(f); BS.save(f); bbwt.save(f);
        succp.save(f); dstart.save(f); jv.save(f); EliasFano::wv(f, ccont);
        if (sided) EliasFano::wv(f, sidev);
        dict.save(f);
        fclose(f);
    }
    void load(const std::string& path) {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path);
        char mg[8];
        if (fread(mg, 1, 8, f) != 8 || (std::memcmp(mg, "TWOLVL05", 8) && std::memcmp(mg, "TWOLVL06", 8)))
            throw std::runtime_error(path + ": not a two-level file of this version (rebuild it)");
        if (fread(&n, 8, 1, f) != 1 || fread(&m, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 ||
            fread(&nblocks, 8, 1, f) != 1 || fread(&nsym, 8, 1, f) != 1) throw std::runtime_error("read");
        B.load(f); pbwt.load(f); BS.load(f); bbwt.load(f);
        succp.load(f); dstart.load(f); jv.load(f); EliasFano::rv(f, ccont);
        sided = mg[7] == '6';
        sidev.clear();
        if (sided) EliasFano::rv(f, sidev);
        dict.load(f);
        fclose(f);
    }
};

// The parse of a pattern against a level's dictionary.
struct PatInfo {
    std::vector<int64_t> S;         // syncmer starts
    std::vector<uint32_t> pid;      // pid[i]: id of P[S[i] .. S[i+1]+K), 0 if not a phrase (or i is last)
    std::vector<uint8_t> cont;      // cont[i]: phrase i is a continuation of phrase i-1
    std::vector<int32_t> bstart;    // bstart[i]: the block of phrase i starts at phrase bstart[i]
    void finish(const Level& L) {
        size_t z = S.size();
        cont.assign(z, 0); bstart.assign(z, 0);
        for (size_t i = 0; i < z; ++i) {
            cont[i] = i > 0 && pid[i] && pid[i - 1] && L.succ(pid[i - 1]) == pid[i];
            bstart[i] = cont[i] ? bstart[i - 1] : (int32_t)i;
        }
    }
};

// Parses P against T's level (into pi) and its reverse Q against the reversed
// collection's (into qi), with one scan for syncmers and one dictionary lookup per
// phrase: the phrases of Q are the reverses of P's, in reverse order.
inline void parse_pattern(const Level& L, const Level& R, const std::string& P, PatInfo& pi, PatInfo& qi) {
    const PhraseDict& D = L.dict;
    const uint8_t* x = reinterpret_cast<const uint8_t*>(P.data());
    pi.S.clear();
    for_each_syncmer(x, P.size(), [&](uint64_t q) { pi.S.push_back((int64_t)q); });
    const size_t z = pi.S.size(), np = z ? z - 1 : 0;
    const int64_t m = (int64_t)P.size();
    qi.S.resize(z);
    for (size_t i = 0; i < z; ++i) qi.S[i] = m - K - pi.S[z - 1 - i];
    std::vector<uint64_t> h(np), sl(np);
    for (size_t i = 0; i < np; ++i) h[i] = phrase_hash(x + pi.S[i], pi.S[i + 1] + K - pi.S[i]);
    constexpr size_t AHEAD = 16;
    for (size_t i = 0; i < np && i < AHEAD; ++i) D.prefetch(h[i]);
    for (size_t i = 0; i < np; ++i) {
        if (i + AHEAD < np) D.prefetch(h[i + AHEAD]);
        sl[i] = D.slot(h[i]);
    }
    pi.pid.assign(z, 0); qi.pid.assign(z, 0);
    for (size_t i = 0; i < np; ++i) {
        if (i + AHEAD < np && sl[i + AHEAD] != ~0ull) {
            __builtin_prefetch(&D.idf.d[(sl[i + AHEAD] * D.idf.w) >> 6]);
            __builtin_prefetch(&D.idr.d[(sl[i + AHEAD] * D.idr.w) >> 6]);
        }
        if (sl[i] == ~0ull) continue;
        pi.pid[i] = (uint32_t)D.idf.get(sl[i]);
        qi.pid[np - 1 - i] = (uint32_t)D.idr.get(sl[i]);
    }
    pi.finish(L); qi.finish(R);
}

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
                    // side: the phrase after the block against succ of the block's last phrase
                    const bool plus = pi.pid[cur] > L.succ(pi.pid[cur - 1]);
                    uint32_t y = L.sym(pi.pid[t0], (uint64_t)(cur - 1 - t0), plus);
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
