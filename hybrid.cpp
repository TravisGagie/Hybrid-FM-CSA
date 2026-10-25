// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Psi_E index: counting over an order-1 frequency-rank encoding of the text.
//
//   T[0..n0-1] : input integers, remapped to dense ids 1..sigma; $ = 0.
//   S[k]       : rank (1 = most frequent) of T[k+1] among the successors of T[k],
//                for k = 0..n0-2; S[n0-1] = 0 is the terminator.  |S| = n0.
//                Suffix S[k..] encodes T[k+1..], and is preceded in T by T[k].
//   BWT(S)     : run-length encoded, used for backward search on enc(P)[2..m].
//   Psi_E[i]   : T[SA_S[i]], the unencoded character preceding the i-th encoded
//                suffix.
//
//   count(P) = rank_{P[0]}(Psi_E, r) - rank_{P[0]}(Psi_E, l), where [l, r) is the
//   SA_S interval of enc(P)[2..m].
//
// Both BWT(S) and Psi_E are stored with a run-length rank structure, either
//   RunRank (explicit per-symbol run lists, 64 bits per run, fast), or
//   EFRank  (two Elias-Fano sequences over the runs in F'-order, compact).
//
// Modes:
//   hybrid build <input> <index> [width] [explicit]   width in bytes: 1, 2, 4 (default), 8
//   hybrid count <index> <patterns.txt>               one pattern per line, decimal values
//   hybrid mems <index> [--no-psi] [-L len] <query.qenc> ...   pseudo-MEMs of queries parsed by syncparse
//            (--no-psi: from encoded matches only, without checks on Psi_E; still safe, but looser;
//             -L: only the pseudo-MEMs of at least len characters, found with Boyer-Moore-Li skipping;
//             an index built as bihuff also has the reversed encoding, for forward-backward search;
//             --one-way: do not use it)
//   hybrid bench <input> [width] [queries] [m ...]
//            builds the Psi_E index and a baseline RLFM-index of T, each with both
//            rank structures, checks they agree, and times counting.
#include "runrank.hpp"
#include "efrank.hpp"
#include "huffwt.hpp"
#include "gaprank.hpp"
using EscRank = EscRankT<RLHuffRank, EFRank>;     // run-length: for repetitive inputs
using EscRankP = EscRankT<HuffWT, GapRank>;       // plain: for inputs with few long runs
#include <functional>
#include "libsais.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <memory>

using namespace std;
using clk = chrono::steady_clock;
static double secs(clk::time_point a) { return chrono::duration<double>(clk::now() - a).count(); }

// ---------------------------------------------------------------- input
static vector<uint64_t> read_ints(const string& fn, int width) {
    FILE* f = fopen(fn.c_str(), "rb");
    if (!f) { perror(fn.c_str()); exit(1); }
    fseek(f, 0, SEEK_END); long long sz = ftell(f); fseek(f, 0, SEEK_SET);
    vector<uint8_t> buf(sz);
    if (sz && fread(buf.data(), 1, sz, f) != (size_t)sz) { perror("read"); exit(1); }
    fclose(f);
    if (sz % width) fprintf(stderr, "warning: file size not a multiple of %d; ignoring tail\n", width);
    size_t n = sz / width;
    vector<uint64_t> v(n);
    for (size_t i = 0; i < n; ++i) {
        uint64_t x = 0;
        memcpy(&x, &buf[i * width], width);           // little-endian
        v[i] = x;
    }
    return v;
}

static double H0(const vector<uint32_t>& s, uint32_t sig) {
    vector<uint64_t> c(sig, 0);
    for (auto x : s) c[x]++;
    double h = 0, n = s.size();
    for (auto x : c) if (x) h -= x / n * log2(x / n);
    return h;
}

static vector<uint32_t> densify(const vector<uint64_t>& raw, vector<uint64_t>& vals) {
    vals = raw;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    vector<uint32_t> T(raw.size());
    for (size_t i = 0; i < raw.size(); ++i)
        T[i] = (lower_bound(vals.begin(), vals.end(), raw[i]) - vals.begin()) + 1;
    return T;
}

// ---------------------------------------------------------------- encoding
// A bit array with reads and writes of up to 64 bits at any bit position.
struct BitArray {
    vector<uint64_t> a;
    void init(uint64_t nbits) { a.assign(nbits / 64 + 2, 0); }
    inline void set(uint64_t pos, uint32_t w, uint64_t x) {
        if (!w) return;
        uint64_t q = pos >> 6, o = pos & 63;
        a[q] |= x << o;
        if (o + w > 64) a[q + 1] |= x >> (64 - o);
    }
    inline uint64_t get(uint64_t pos, uint32_t w) const {
        if (!w) return 0;
        uint64_t q = pos >> 6, o = pos & 63;
        uint64_t x = a[q] >> o;
        if (o + w > 64) x |= a[q + 1] << (64 - o);
        return w == 64 ? x : x & ((1ULL << w) - 1);
    }
};
static inline uint32_t bits_for(uint64_t x) { return x ? 64 - __builtin_clzll(x) : 0; }

struct Encoding {
    // Map from input values to dense ids 1..sigma.  If the input values are
    // consecutive (as phrase ids from syncparse are), only the first is stored.
    vector<uint64_t> vals;          // vals[id-1] = raw value of dense id (empty if consecutive)
    bool consecutive = false;
    uint64_t vals_base = 0;
    uint32_t sigma = 0;             // number of distinct input values
    uint32_t kE = 0;                // alphabet size of S = max rank + 1
    uint64_t n0 = 0;                // |S|
    vector<uint64_t> C;             // kE+1

    // Rank lists while building: for a, pairs (b, rank) sorted by b in soff[a]..soff[a+1].
    vector<uint64_t> soff;
    vector<uint32_t> sb, srank;
    // Compact rank lists, used after building:
    //   offs: Elias-Fano list of the index of each phrase's first successor (about
    //     2 + lg(bigrams/sigma) bits per phrase);
    //   csucc: the successors, ceil(lg(sigma+1)) bits each.  For a phrase with
    //     d <= SMALL successors they are stored in rank order, so the rank of b
    //     is its position + 1 (found by a short linear scan, with no rank stored).
    //     For d > SMALL they are sorted by id (found by binary search), and
    //     their ranks are stored in bigranks, wbig bits each; bigA lists those
    //     phrases and bigstart says where their ranks start.
    // A lookup costs about one cache miss for small degrees.
    static constexpr uint64_t SMALL = 16;
    bool compact = false;
    EliasFano offs;
    uint32_t wsucc = 0, wbig = 0;
    BitArray csucc, bigranks;
    vector<uint32_t> bigA;
    vector<uint64_t> bigstart;
    uint64_t nbigrams = 0;

    void compress_vals() {
        consecutive = true;
        for (size_t i = 0; i < vals.size(); ++i) if (vals[i] != vals[0] + i) { consecutive = false; break; }
        if (consecutive) { vals_base = vals.empty() ? 0 : vals[0]; vector<uint64_t>().swap(vals); }
    }
    int64_t id_of(uint64_t raw) const {
        if (consecutive) return (raw >= vals_base && raw - vals_base < sigma) ? (int64_t)(raw - vals_base) + 1 : -1;
        auto it = lower_bound(vals.begin(), vals.end(), raw);
        if (it == vals.end() || *it != raw) return -1;
        return (it - vals.begin()) + 1;
    }
    // Capped rank lists (cap > 0): only each phrase's cap most frequent successors are
    // ranked.  A pair (a, b) whose b is not among them is encoded as the escape cap + b,
    // which names b itself, so E[i] still determines T[i] given T[i-1]; the rank lists
    // hold at most sigma * cap entries.  cap = 0: no cap.
    uint32_t cap = 0;
    // the symbol of E for the pair (a, b): its rank, its escape, or 0 if (a, b) is known not to occur
    inline uint32_t encode(uint32_t a, uint32_t b) const {
        uint32_t r = succ_rank(a, b);
        if (r || !cap || !b || b > sigma) return r;
        return cap + b;
    }
    uint32_t succ_rank(uint32_t a, uint32_t b) const {
        if (!compact) {
            const uint32_t* lo = sb.data() + soff[a];
            const uint32_t* hi = sb.data() + soff[a + 1];
            const uint32_t* it = lower_bound(lo, hi, b);
            if (it == hi || *it != b) return 0;
            return srank[it - sb.data()];
        }
        if (a > sigma) return 0;
        uint64_t lo, hi;
        offs.access2(a, lo, hi);
        uint64_t d = hi - lo;
        if (d == 0) return 0;
        if (d <= SMALL) {
            for (uint64_t j = 0; j < d; ++j) if (csucc.get((lo + j) * wsucc, wsucc) == b) return (uint32_t)j + 1;
            return 0;
        }
        uint64_t l = lo, h = hi;
        while (l < h) {
            uint64_t mid = (l + h) >> 1;
            if (csucc.get(mid * wsucc, wsucc) < b) l = mid + 1; else h = mid;
        }
        if (l == hi || csucc.get(l * wsucc, wsucc) != b) return 0;
        uint64_t k = lower_bound(bigA.begin(), bigA.end(), a) - bigA.begin();
        return (uint32_t)bigranks.get((bigstart[k] + (l - lo)) * wbig, wbig);
    }
    void compress_ranks() {
        nbigrams = sb.size();
        if (nbigrams >= (1ull << 32)) { fprintf(stderr, "too many distinct pairs\n"); exit(1); }
        offs.build(soff, nbigrams + 1);
        wsucc = bits_for(sigma);
        csucc.init(nbigrams * wsucc);
        uint64_t maxd = 0, nbig = 0;
        for (uint32_t a = 0; a <= sigma; ++a) {
            uint64_t d = soff[a + 1] - soff[a];
            maxd = max(maxd, d);
            if (d > SMALL) { bigA.push_back(a); bigstart.push_back(nbig); nbig += d; }
        }
        wbig = bits_for(maxd);
        bigranks.init(nbig * wbig);
        uint64_t k = 0;
        for (uint32_t a = 0; a <= sigma; ++a) {
            uint64_t lo = soff[a], d = soff[a + 1] - lo;
            if (d <= SMALL) {
                for (uint64_t j = 0; j < d; ++j) csucc.set((lo + srank[lo + j] - 1) * wsucc, wsucc, sb[lo + j]);
            } else {
                for (uint64_t j = 0; j < d; ++j) {
                    csucc.set((lo + j) * wsucc, wsucc, sb[lo + j]);
                    bigranks.set((bigstart[k] + j) * wbig, wbig, srank[lo + j]);
                }
                ++k;
            }
        }
        compact = true;
        vector<uint64_t>().swap(soff); vector<uint32_t>().swap(sb); vector<uint32_t>().swap(srank);
    }

    // From T (dense ids, no terminator) compute the rankings, S, SA(S), BWT(S) and Psi_E.
    void prepare(const vector<uint32_t>& T, vector<uint32_t>& L, vector<uint32_t>& P, bool verbose,
                 vector<uint32_t>* Lrev = nullptr) {
        n0 = T.size();
        if (n0 >= (1ull << 31) - 16) { fprintf(stderr, "input too long for 32-bit libsais\n"); exit(1); }
        auto t0 = clk::now();
        vector<uint64_t> pr;
        pr.reserve(n0 ? n0 - 1 : 0);
        for (uint64_t k = 0; k + 1 < n0; ++k) pr.push_back((uint64_t)T[k] << 32 | T[k + 1]);
        sort(pr.begin(), pr.end());
        soff.assign(sigma + 2, 0);
        sb.clear(); srank.clear();
        vector<pair<uint64_t, uint32_t>> grp;         // (count, b)
        size_t i = 0;
        kE = 1;
        for (uint32_t a = 1; a <= sigma; ++a) {
            soff[a] = sb.size();
            grp.clear();
            while (i < pr.size() && (pr[i] >> 32) == a) {
                size_t j = i;
                while (j < pr.size() && pr[j] == pr[i]) ++j;
                grp.push_back({j - i, (uint32_t)pr[i]});
                i = j;
            }
            vector<size_t> ord(grp.size());
            for (size_t q = 0; q < ord.size(); ++q) ord[q] = q;
            sort(ord.begin(), ord.end(), [&](size_t x, size_t y) {
                return grp[x].first != grp[y].first ? grp[x].first > grp[y].first
                                                    : grp[x].second < grp[y].second; });
            vector<uint32_t> rk(grp.size());
            for (size_t q = 0; q < ord.size(); ++q) rk[ord[q]] = q + 1;
            for (size_t q = 0; q < grp.size(); ++q) {
                if (cap && rk[q] > cap) continue;          // escaped: not in the lists
                sb.push_back(grp[q].second); srank.push_back(rk[q]);
            }
            kE = max<uint32_t>(kE, (cap ? min<uint32_t>(grp.size(), cap) : grp.size()) + 1);
        }
        soff[sigma + 1] = sb.size();
        vector<uint64_t>().swap(pr);
        if (cap) kE = cap + sigma + 1;              // ranks 1..cap, escapes cap+1..cap+sigma

        vector<int32_t> S(n0);
        uint64_t nesc = 0;
        for (uint64_t k = 0; k + 1 < n0; ++k) { S[k] = encode(T[k], T[k + 1]); nesc += (uint32_t)S[k] > cap && cap; }
        if (n0) S[n0 - 1] = 0;
        if (verbose) fprintf(stderr, "  encoded (%.2fs): %zu ranked bigrams, E alphabet %u%s\n",
                             secs(t0), sb.size(), kE, cap ? (", cap " + to_string(cap) + ", " + to_string(nesc) + " escapes").c_str() : "");
        compress_ranks();
        t0 = clk::now();
        vector<int32_t> SA(n0 + 6 * (size_t)kE);
        if (libsais_int(S.data(), SA.data(), n0, kE, 6 * kE) != 0) { fprintf(stderr, "libsais failed\n"); exit(1); }
        SA.resize(n0);
        if (verbose) fprintf(stderr, "  SA(S) (%.2fs)\n", secs(t0));

        L.assign(n0, 0); P.assign(n0, 0);
        C.assign(kE + 1, 0);
        for (uint64_t k = 0; k < n0; ++k) C[S[k] + 1]++;
        for (uint32_t c = 0; c < kE; ++c) C[c + 1] += C[c];
        for (uint64_t r = 0; r < n0; ++r) {
            uint64_t p = SA[r];
            L[r] = S[p ? p - 1 : n0 - 1];
            P[r] = T[p];
        }
        if (verbose) {
            vector<uint32_t> Su(S.begin(), S.end());
            fprintf(stderr, "  H0(S) = %.3f bits\n", H0(Su, kE));
        }
        if (Lrev) {                     // BWT of S reversed (terminator kept at the end)
            vector<int32_t>().swap(SA);
            vector<int32_t> R(n0);
            for (uint64_t k = 0; k + 1 < n0; ++k) R[k] = S[n0 - 2 - k];
            if (n0) R[n0 - 1] = 0;
            vector<int32_t> SR(n0 + 6 * (size_t)kE);
            if (libsais_int(R.data(), SR.data(), n0, kE, 6 * kE) != 0) { fprintf(stderr, "libsais failed\n"); exit(1); }
            Lrev->assign(n0, 0);
            for (uint64_t r = 0; r < n0; ++r) { uint64_t p = SR[r]; (*Lrev)[r] = R[p ? p - 1 : n0 - 1]; }
            if (verbose) fprintf(stderr, "  BWT of the reversed encoding built\n");
        }
    }
    // (C is needed only by rank structures that cannot supply C[c] themselves)
    uint64_t bytes(bool withC = true) const {
        return vals.size() * 8 + (withC ? C.size() * 8 : 0) + offs.bytes() + csucc.a.size() * 8 + bigranks.a.size() * 8 +
               bigA.size() * 4 + bigstart.size() * 8 + 64;
    }
    void save(FILE* f) const {
        fwrite(&sigma, 4, 1, f); fwrite(&kE, 4, 1, f); fwrite(&n0, 8, 1, f);
        uint8_t c = consecutive; fwrite(&c, 1, 1, f); fwrite(&vals_base, 8, 1, f);
        fwrite(&wsucc, 4, 1, f); fwrite(&wbig, 4, 1, f); fwrite(&nbigrams, 8, 1, f); fwrite(&cap, 4, 1, f);
        RunRank::wv(f, vals); RunRank::wv(f, C); offs.save(f);
        RunRank::wv(f, csucc.a); RunRank::wv(f, bigranks.a); RunRank::wv(f, bigA); RunRank::wv(f, bigstart);
    }
    void load(FILE* f, int version = 4) {
        uint8_t c;
        if (fread(&sigma, 4, 1, f) != 1 || fread(&kE, 4, 1, f) != 1 || fread(&n0, 8, 1, f) != 1 ||
            fread(&c, 1, 1, f) != 1 || fread(&vals_base, 8, 1, f) != 1 || fread(&wsucc, 4, 1, f) != 1 ||
            fread(&wbig, 4, 1, f) != 1 || fread(&nbigrams, 8, 1, f) != 1) exit(1);
        cap = 0;
        if (version >= 4 && fread(&cap, 4, 1, f) != 1) exit(1);
        consecutive = c;
        RunRank::rv(f, vals); RunRank::rv(f, C);
        if (version >= 4) offs.load(f);
        else { vector<uint32_t> o; RunRank::rv(f, o); offs.build(vector<uint64_t>(o.begin(), o.end()), nbigrams + 1); }
        RunRank::rv(f, csucc.a); RunRank::rv(f, bigranks.a); RunRank::rv(f, bigA); RunRank::rv(f, bigstart);
        compact = true;
    }
};

// LF step [l, r) -> [C[c] + rank_c(l), C[c] + rank_c(r)).  Rank structures that keep
// F positions of their runs (Elias-Fano run lists, run-length Huffman WTs) do this
// without a C array (method lf2); for the others C is passed in.
template <class R, class = void> struct has_lf2 : std::false_type {};
template <class R> struct has_lf2<R, std::void_t<decltype(std::declval<const R&>().lf2(0u, 0ull, 0ull,
    std::declval<uint64_t&>(), std::declval<uint64_t&>()))>> : std::true_type {};
template <class R, class CF>
static inline void lf_step(const R& bwt, CF Cof, uint32_t c, uint64_t& l, uint64_t& r) {
    if constexpr (has_lf2<R>::value) { uint64_t a, b; bwt.lf2(c, l, r, a, b); l = a; r = b; }
    else { uint64_t rl, rr, cc = Cof(c); bwt.rank2(c, l, r, rl, rr); l = cc + rl; r = cc + rr; }
}

// ---------------------------------------------------------------- index
template <class RB, class RP = RB>
struct PsiEIndex {
    const Encoding* enc = nullptr;
    RB bwt;
    RP psi;

    void build(const Encoding& e, const vector<uint32_t>& L, const vector<uint32_t>& P) {
        enc = &e;
        if constexpr (std::is_same_v<RB, EscRank> || std::is_same_v<RB, EscRankP>) bwt.cap = e.cap;
        bwt.build(L, e.kE);
        psi.build(P, e.sigma + 1);
    }
    // pattern as dense ids (0 = value absent from T)
    uint64_t count_ids(const vector<uint32_t>& p) const {
        size_t m = p.size();
        if (m == 0) return enc->n0;
        for (auto x : p) if (x == 0) return 0;
        uint64_t l = 0, r = enc->n0;
        for (size_t j = m - 1; j >= 1; --j) {          // backward search enc(P)[2..m]
            uint32_t c = enc->encode(p[j - 1], p[j]);
            if (c == 0) return 0;
            lf_step(bwt, [this](uint32_t x) { return enc->C[x]; }, c, l, r);
            if (l >= r) return 0;
        }
        uint64_t rl, rr;
        psi.rank2(p[0], l, r, rl, rr);
        return rr - rl;
    }
    uint64_t count_raw(const vector<uint64_t>& p) const {
        vector<uint32_t> q(p.size());
        for (size_t j = 0; j < p.size(); ++j) { int64_t id = enc->id_of(p[j]); q[j] = id < 0 ? 0 : id; }
        return count_ids(q);
    }
    uint64_t bytes() const { return enc->bytes(!has_lf2<RB>::value) + bwt.bytes() + psi.bytes(); }
};

// Baseline: run-length FM-index of T$ with the same rank structure (with Elias-Fano
// run lists this is a compressed RLCSA of T).  A C array, Elias-Fano coded, is kept
// only when the rank structure cannot supply C[c] itself.
template <class R>
struct RLFM {
    uint64_t n = 0;
    EliasFano C;
    R bwt;
    void build(const vector<uint32_t>& L, const vector<uint64_t>& C_, uint32_t sigma) {
        n = L.size();
        if (!has_lf2<R>::value) C.build(C_, n + 1);
        bwt.build(L, sigma + 1);
    }
    inline void step(uint32_t c, uint64_t& l, uint64_t& r) const {
        lf_step(bwt, [this](uint32_t x) { return C.access(x); }, c, l, r);
    }
    // the interval of the single symbol c
    inline void range(uint32_t c, uint64_t& l, uint64_t& r) const { l = 0; r = n; step(c, l, r); }
    uint64_t count_ids(const vector<uint32_t>& p) const {
        uint64_t l = 0, r = n;
        for (size_t j = p.size(); j-- > 0;) {
            uint32_t c = p[j];
            if (c == 0) return 0;
            step(c, l, r);
            if (l >= r) return 0;
        }
        return r - l;
    }
    uint64_t bytes() const { return (has_lf2<R>::value ? 0 : C.bytes()) + bwt.bytes(); }
    void save(FILE* f) const { fwrite(&n, 8, 1, f); bwt.save(f); if (!has_lf2<R>::value) C.save(f); }
    void load(FILE* f) { if (fread(&n, 8, 1, f) != 1) exit(1); bwt.load(f); if (!has_lf2<R>::value) C.load(f); }
};

static void bwt_of_T(const vector<uint32_t>& T, uint32_t sigma, vector<uint32_t>& L, vector<uint64_t>& C) {
    uint64_t n = T.size() + 1;
    vector<int32_t> X(n), SA(n + 6 * (size_t)(sigma + 1));
    for (size_t k = 0; k < T.size(); ++k) X[k] = T[k];
    X[n - 1] = 0;
    if (libsais_int(X.data(), SA.data(), n, sigma + 1, 6 * (sigma + 1)) != 0) { fprintf(stderr, "libsais failed\n"); exit(1); }
    L.assign(n, 0);
    C.assign(sigma + 2, 0);
    for (uint64_t k = 0; k < n; ++k) C[X[k] + 1]++;
    for (uint32_t c = 0; c <= sigma; ++c) C[c + 1] += C[c];
    for (uint64_t r = 0; r < n; ++r) L[r] = X[SA[r] ? SA[r] - 1 : n - 1];
}

static string mb(uint64_t b) { char s[64]; snprintf(s, 64, "%.2f MB", b / 1048576.0); return s; }
static uint64_t count_runs(const vector<uint32_t>& s) {
    uint64_t r = 0; for (size_t i = 0; i < s.size(); ++i) r += (i == 0 || s[i] != s[i - 1]); return r;
}

// ---------------------------------------------------------------- modes
static int do_build(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: hybrid build <input> <index> [width] [huff|ef|explicit|bihuff|rlcsa|birlcsa] [cap]\n"); return 1; }
    int width = argc > 4 ? atoi(argv[4]) : 4;
    string kind = argc > 5 ? argv[5] : "huff";
    uint32_t cap = argc > 6 ? (uint32_t)atoi(argv[6]) : 0;
    if (cap && kind != "huff" && kind != "ef") { fprintf(stderr, "a cap needs kind huff or ef\n"); return 1; }
    if (kind != "explicit" && kind != "ef" && kind != "huff" && kind != "bihuff" && kind != "rlcsa" && kind != "birlcsa") {
        fprintf(stderr, "kind must be explicit, ef, huff, bihuff, rlcsa or birlcsa\n"); return 1;
    }
    auto t0 = clk::now();
    auto raw = read_ints(argv[2], width);
    Encoding e;
    e.cap = cap;
    auto T = densify(raw, e.vals);
    e.sigma = e.vals.size();
    e.compress_vals();
    vector<uint64_t>().swap(raw);
    if (kind == "rlcsa" || kind == "birlcsa") {   // RLCSA (Elias-Fano run lists) of T, and of T reversed with birlcsa
        FILE* f = fopen(argv[3], "wb");
        if (!f) { perror(argv[3]); return 1; }
        fwrite("RLCSA002", 1, 8, f);
        uint8_t c = e.consecutive, hr = kind == "birlcsa";
        fwrite(&e.sigma, 4, 1, f); fwrite(&c, 1, 1, f); fwrite(&e.vals_base, 8, 1, f); fwrite(&hr, 1, 1, f);
        RunRank::wv(f, e.vals);
        uint64_t sz = e.vals.size() * 8, rf, rr = 0;
        {
            vector<uint32_t> LT; vector<uint64_t> CT;
            bwt_of_T(T, e.sigma, LT, CT);
            rf = count_runs(LT);
            RLFM<EFRank> fm; fm.build(LT, CT, e.sigma); fm.save(f); sz += fm.bytes();
            printf("RLCSA of T: %llu runs, %s\n", (unsigned long long)rf, mb(fm.bytes()).c_str());
        }
        if (hr) {
        reverse(T.begin(), T.end());
            vector<uint32_t> LT; vector<uint64_t> CT;
            bwt_of_T(T, e.sigma, LT, CT);
            rr = count_runs(LT);
            RLFM<EFRank> fm; fm.build(LT, CT, e.sigma); fm.save(f); sz += fm.bytes();
            printf("RLCSA of T reversed: %llu runs, %s\n", (unsigned long long)rr, mb(fm.bytes()).c_str());
        }
        (void)rr;
        fclose(f);
        printf("n = %zu, sigma = %u; rlcsa index: %s, built in %.2fs\n", T.size() + 1, e.sigma, mb(sz).c_str(), secs(t0));
        return 0;
    }
    vector<uint32_t> L, P, Lrev;
    e.prepare(T, L, P, true, kind == "bihuff" ? &Lrev : nullptr);
    FILE* f = fopen(argv[3], "wb");
    if (!f) { perror(argv[3]); return 1; }
    bool esc = cap && kind == "huff";
    fwrite(kind == "explicit" ? "PSIEX004" : kind == "ef" ? "PSIEF004" : esc ? "PSIEC004" : "PSIEH004", 1, 8, f);
    e.save(f);
    uint64_t sz;
    if (kind == "explicit") { PsiEIndex<RunRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else if (kind == "ef")  { PsiEIndex<EFRank> ix;  ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else if (esc) {
        PsiEIndex<EscRank, EFRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes();
        printf("BWT(E) %s (ranks and ESC %s, escapes %s), Psi_E %s, rank lists %s\n", mb(ix.bwt.bytes()).c_str(),
               mb(ix.bwt.top.bytes()).c_str(), mb(ix.bwt.has_sub ? ix.bwt.sub.bytes() : 0).c_str(),
               mb(ix.psi.bytes()).c_str(), mb(e.bytes(false)).c_str());
    }
    else {
        PsiEIndex<RLHuffRank, EFRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes();
        printf("BWT(E) %s, Psi_E %s, rank lists %s\n", mb(ix.bwt.bytes()).c_str(), mb(ix.psi.bytes()).c_str(), mb(e.bytes(false)).c_str());
    }
    if (kind == "bihuff") {             // index of the reversed encoding, for forward search
        RLHuffRank rev; rev.build(Lrev, e.kE);
        fwrite("REVBWT01", 1, 8, f); rev.save(f);
        printf("reversed encoding: %llu runs, %s\n", (unsigned long long)count_runs(Lrev), mb(rev.bytes()).c_str());
        sz += rev.bytes();
    }
    fclose(f);
    printf("n = %llu, sigma = %u, E alphabet = %u, runs: BWT(E) = %llu, Psi_E = %llu\n",
           (unsigned long long)e.n0, e.sigma, e.kE,
           (unsigned long long)count_runs(L), (unsigned long long)count_runs(P));
    printf("%s index: %s, built in %.2fs\n", kind.c_str(), mb(sz).c_str(), secs(t0));
    return 0;
}

template <class RB, class RP>
static int count_with(FILE* f, const Encoding& e, const char* pfile) {
    PsiEIndex<RB, RP> ix;
    ix.enc = &e;
    ix.bwt.load(f); ix.psi.load(f);
    ifstream in(pfile);
    string line;
    auto t0 = clk::now();
    uint64_t q = 0;
    while (getline(in, line)) {
        istringstream ss(line);
        vector<uint64_t> p; uint64_t x;
        while (ss >> x) p.push_back(x);
        if (p.empty()) continue;
        printf("%llu\n", (unsigned long long)ix.count_raw(p));
        ++q;
    }
    fprintf(stderr, "%llu patterns in %.3fs\n", (unsigned long long)q, secs(t0));
    return 0;
}

static int do_count(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: hybrid count <index> <patterns.txt>\n"); return 1; }
    FILE* f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    char mg[8];
    if (fread(mg, 1, 8, f) != 8) { fprintf(stderr, "bad index file\n"); return 1; }
    Encoding e;
    if (memcmp(mg, "PSIE", 4)) { fprintf(stderr, "bad index file\n"); return 1; }
    e.load(f, mg[7] - '0');
    char k = mg[4];
    if (k == 'X') return count_with<RunRank, RunRank>(f, e, argv[3]);
    if (k == 'F') return count_with<EFRank, EFRank>(f, e, argv[3]);
    if (k == 'H') return count_with<RLHuffRank, EFRank>(f, e, argv[3]);
    if (k == 'C') return count_with<EscRank, EFRank>(f, e, argv[3]);
    fprintf(stderr, "bad index file\n");
    return 1;
}

// ---------------------------------------------------------------- MEMs and pseudo-MEMs
// Reads a query parse written by `syncparse query` (a .qenc file), finds the MEMs of
// the query's parse with respect to the indexed parse using only backward search on
// the hybrid index, and reports pseudo-MEMs.
//
// Let the query's phrases be Q[0..t-1], phrase a covering positions [s_a, s_{a+1})
// of the query (non-overlapping phrases; the last phrase ends at the query's end).
// For each MEM Q[i..j] of the query's parse, the pseudo-MEM runs from the start of
// phrase i-1 to k characters past the end of phrase j+1 (clipped to the query).
// For each pair of consecutive phrases a, a+1 that are both absent from the
// dataset, the pseudo-MEM runs from s_a to k characters past the end of phrase a+1.
// Every MEM of the query with respect to the dataset lies inside a pseudo-MEM.
// The lower bound of a pseudo-MEM from a parse MEM Q[i..j] is s_{j+1} - s_i: the
// phrases i..j occur in the dataset's parse, so their concatenation occurs in the
// dataset.  (With overlapping phrases, from `syncparse build --overlap`, the
// pseudo-MEMs extend just one phrase to the right, as in prefix-free parsing.)
// N runs: syncparse splits the query at maximal runs of characters other than
// A, C, G, T and parses each segment separately.  N counts as a mismatch, so no
// MEM crosses a segment boundary: parse MEMs stop at boundaries and pseudo-MEMs are
// clipped to their segment (at its ends there is no next phrase to extend into).
struct QueryParse {
    int k = 0, overlap = 0;
    uint64_t L = 0;
    vector<uint64_t> raw, start;        // phrase ids as written by syncparse (0xFFFFFFFF = absent)
    vector<uint64_t> segend, segbeg;    // per phrase: end and start of its segment
    vector<uint64_t> rcraw;             // overlap = 2 (strand-symmetric parse): id of each phrase's reverse complement
    vector<uint8_t> first, last;        // first / last phrase of its segment
    uint64_t nruns = 0, nbases = 0;     // N runs and their total length
};
static QueryParse read_qenc(const char* fn) {
    QueryParse q;
    FILE* f = fopen(fn, "rb");
    if (!f) { perror(fn); exit(1); }
    char mg[8]; uint64_t np;
    bool v3 = false;
    if (fread(mg, 1, 8, f) != 8 || (memcmp(mg, "QENC0002", 8) && !(v3 = !memcmp(mg, "QENC0003", 8))) || fread(&q.k, 4, 1, f) != 1 ||
        fread(&q.overlap, 4, 1, f) != 1 || fread(&q.L, 8, 1, f) != 1 || fread(&np, 8, 1, f) != 1) {
        fprintf(stderr, "%s is not a .qenc file from this version of syncparse\n", fn); exit(1);
    }
    vector<uint32_t> rec(4 * np);
    if (np && fread(rec.data(), 4, rec.size(), f) != rec.size()) { fprintf(stderr, "%s is truncated\n", fn); exit(1); }
    vector<uint64_t> nr;
    if (v3) {
        if (fread(&q.nruns, 8, 1, f) != 1) { fprintf(stderr, "%s is truncated\n", fn); exit(1); }
        nr.resize(2 * q.nruns);
        if (q.nruns && fread(nr.data(), 8, nr.size(), f) != nr.size()) { fprintf(stderr, "%s is truncated\n", fn); exit(1); }
        for (uint64_t r = 0; r < q.nruns; ++r) q.nbases += nr[2 * r + 1] - nr[2 * r];
        if (q.overlap == 2) {
            vector<uint32_t> rc(np);
            if (np && fread(rc.data(), 4, np, f) != np) { fprintf(stderr, "%s is truncated\n", fn); exit(1); }
            q.rcraw.assign(rc.begin(), rc.end());
        }
    }
    fclose(f);
    q.raw.resize(np); q.start.resize(np); q.segend.resize(np); q.segbeg.resize(np); q.first.resize(np); q.last.resize(np);
    size_t r = 0;
    for (uint64_t t = 0; t < np; ++t) {
        q.raw[t] = rec[4 * t]; q.start[t] = rec[4 * t + 3];
        uint32_t fl = rec[4 * t + 2];
        q.first[t] = t == 0 || (v3 && (fl & 1));
        q.last[t] = t + 1 == np || (v3 && (fl & 2));
        while (r < q.nruns && nr[2 * r] <= q.start[t]) ++r;
        q.segend[t] = r < q.nruns ? nr[2 * r] : q.L;
        q.segbeg[t] = r > 0 ? nr[2 * r - 1] : 0;
    }
    return q;
}

// Search backends for MemFinder.
//   EncBackend: the hybrid.  The interval of state a (match Q[a..j]) is that of
//     R[a+1..j] in BWT(E), starting from the whole range; a backward step from a
//     to a-1 is on R[a]; the first phrase Q[a] is checked against Psi_E.
//   PlainBackend: RLCSAs (RLFMs with Elias-Fano run lists) of T and of T
//     reversed.  The interval of state a is that of Q[a..j] itself, starting from
//     the range of Q[j]; a backward step from a to a-1 is on Q[a-1].  Matches are
//     exact, so there is nothing to check.
template <class RB, class RP>
struct EncBackend {
    const PsiEIndex<RB, RP>& ix;
    const Encoding& e;
    const RB* rev = nullptr;
    static constexpr bool exact = false;
    EncBackend(const PsiEIndex<RB, RP>& ix_, const RB* rev_) : ix(ix_), e(*ix_.enc), rev(rev_) {}
    int64_t id_of(uint64_t raw) const { return e.id_of(raw); }
    string tag() const { return e.cap ? ".cap" + to_string(e.cap) : string(); }
    uint32_t pair_sym(uint32_t a, uint32_t b) const { return e.encode(a, b); }   // R for (a, b), 0 if known absent
    void init(uint32_t, uint64_t& l, uint64_t& r) const { l = 0; r = e.n0; }
    void finit(uint32_t, uint64_t& l, uint64_t& r) const { l = 0; r = e.n0; }
    inline bool step(uint32_t c, uint64_t& l, uint64_t& r) const {
        lf_step(ix.bwt, [this](uint32_t x) { return e.C[x]; }, c, l, r); return l < r;
    }
    inline bool fstep(uint32_t c, uint64_t& l, uint64_t& r) const {   // same symbol counts as E
        lf_step(*rev, [this](uint32_t x) { return e.C[x]; }, c, l, r); return l < r;
    }
    inline bool precedes(uint32_t a, uint64_t l, uint64_t r) const {
        uint64_t rl, rr;
        ix.psi.rank2(a, l, r, rl, rr);
        return rr > rl;
    }
};
struct PlainIndex {
    Encoding map;                       // only the map from input values to dense ids
    RLFM<EFRank> fwd, rev;
    bool have_rev = false;
    uint64_t bytes() const { return map.vals.size() * 8 + fwd.bytes() + (have_rev ? rev.bytes() : 0); }
};
struct PlainBackend {
    const PlainIndex& ix;
    const RLFM<EFRank>* rev = nullptr;
    static constexpr bool exact = true;
    PlainBackend(const PlainIndex& ix_, const RLFM<EFRank>* rev_) : ix(ix_), rev(rev_) {}
    int64_t id_of(uint64_t raw) const { return ix.map.id_of(raw); }
    string tag() const { return ".rlcsa"; }
    uint32_t pair_sym(uint32_t a, uint32_t) const { return a; }
    void init(uint32_t q, uint64_t& l, uint64_t& r) const { ix.fwd.range(q, l, r); }
    void finit(uint32_t q, uint64_t& l, uint64_t& r) const { rev->range(q, l, r); }
    inline bool step(uint32_t c, uint64_t& l, uint64_t& r) const { ix.fwd.step(c, l, r); return l < r; }
    inline bool fstep(uint32_t c, uint64_t& l, uint64_t& r) const { rev->step(c, l, r); return l < r; }
    inline bool precedes(uint32_t, uint64_t, uint64_t) const { return true; }
};

template <class BE>
struct MemFinder {
    const BE& be;
    bool rev = false;                   // forward search available
    vector<uint32_t> Q;                 // dense ids, 0 = absent
    // R[t]: backward step from state t to t-1 (0 = impossible: a phrase absent, the pair
    // never occurs, or t starts a segment); F[t]: forward step from t-1 to t.
    // Hybrid: both are the rank of Q[t] after Q[t-1].  Plain: R[t] = Q[t-1], F[t] = Q[t].
    vector<uint32_t> R, F;
    // Strand-symmetric mode: the index holds the dataset and its reverse complements,
    // so extending Q[p..j] to the right is a backward step for rc(Q[p..j]) in the same
    // index.  rcQ[t] = dense id of rc(Q[t]); F[t] then steps from rc(Q[p..t-1]) to
    // rc(Q[p..t]), and a forward match can be checked with Psi_E like a backward one.
    bool rcmode = false;
    vector<uint32_t> rcQ;
    uint64_t steps = 0, fsteps = 0, checks = 0;   // backward steps, forward steps, rank queries on Psi_E
    static constexpr uint64_t CHECK_EVERY = 16;
    bool use_psi = true;                // false: encoded matches count as matches (no Psi_E)
    uint64_t last_l = 0, last_r = 0;    // interval of the state last returned by an encoded longest_suffix

    // "e-occurs": Q[a..b] e-occurs if Q[a] is present and R[a+1..b] occurs in the
    // encoding.  Every occurrence is an e-occurrence; Q[a..b] occurs if moreover some
    // e-occurrence is preceded by Q[a] (one rank query on Psi_E).  With the plain
    // backend e-occurrence is occurrence.
    // (Without Psi_E, matching is relaxed differently in the two directions -- pairs of
    // consecutive phrases forward, pairs of their reverse complements backward -- and
    // the two relaxations cannot be mixed; so strand-symmetric forward search is used
    // only with Psi_E, or with the exact plain backend.)
    MemFinder(const BE& be_, const QueryParse& q, bool fwd_ok, bool have_rev, bool psi = true)
        : be(be_), rev(have_rev && fwd_ok) {
        use_psi = psi;
        uint64_t np = q.raw.size();
        Q.assign(np, 0); R.assign(np, 0); F.assign(np, 0);
        for (uint64_t t = 0; t < np; ++t) {
            if (q.raw[t] == 0xFFFFFFFFu) continue;
            int64_t id = be.id_of(q.raw[t]);
            Q[t] = id < 0 ? 0 : (uint32_t)id;
        }
        rcmode = q.overlap == 2 && fwd_ok && (BE::exact || use_psi);
        if (rcmode) {
            rcQ.assign(np, 0);
            for (uint64_t t = 0; t < np; ++t) if (Q[t] && q.rcraw[t] != 0xFFFFFFFFu) {
                int64_t id = be.id_of(q.rcraw[t]);
                rcQ[t] = id < 0 ? 0 : (uint32_t)id;
            }
            for (uint64_t t = 0; t < np; ++t) if (!rcQ[t]) Q[t] = 0;   // (cannot happen with a symmetric parse)
        }
        for (uint64_t t = 1; t < np; ++t) if (Q[t - 1] && Q[t] && !q.first[t]) {
            if (BE::exact) R[t] = Q[t - 1];
            else R[t] = be.pair_sym(Q[t - 1], Q[t]);
            if (rcmode) F[t] = BE::exact ? rcQ[t] : be.pair_sym(rcQ[t], rcQ[t - 1]);
            else F[t] = BE::exact ? Q[t] : R[t];
        }
        if (rcmode) rev = true;             // forward search is always available
        if (BE::exact) use_psi = false;
    }
    inline bool step(uint32_t c, uint64_t& l, uint64_t& r) { ++steps; return be.step(c, l, r); }
    inline bool fstep(uint32_t c, uint64_t& l, uint64_t& r) { ++fsteps; return be.fstep(c, l, r); }
    inline bool precedes(uint32_t a, uint64_t l, uint64_t r, bool truth) {
        if (!truth) return true;
        ++checks;
        return be.precedes(a, l, r);
    }
    bool occurs(uint64_t a, uint64_t b, bool truth) {
        if (!Q[a] || !Q[b]) return false;
        uint64_t l, r;
        be.init(Q[b], l, r);
        for (uint64_t x = b; x > a; --x) if (!R[x] || !step(R[x], l, r)) return false;
        return a == b || precedes(Q[a], l, r, truth);
    }
    // smallest i >= floor such that Q[i..j] occurs (truth) or e-occurs; Q[j] must be present
    uint64_t longest_suffix(uint64_t j, uint64_t floor, bool truth) {
        uint64_t l, r;
        be.init(Q[j], l, r);
        vector<pair<uint64_t, uint64_t>> ivl{{l, r}};      // ivl[d]: interval of state j-d
        uint64_t good = j, a = j;
        bool bad = false;
        while (a > floor && Q[a - 1] && R[a]) {
            if (!step(R[a], l, r)) break;
            --a;
            ivl.push_back({l, r});
            if (truth && (j - a) % CHECK_EVERY == 0) {
                if (precedes(Q[a], l, r, true)) good = a; else { bad = true; break; }
            }
        }
        if (!truth) { last_l = ivl.back().first; last_r = ivl.back().second; return a; }
        if (!bad && a < good) {
            auto [ll, rr] = ivl[j - a];
            if (precedes(Q[a], ll, rr, true)) return a;
        }
        uint64_t lo = a + 1, hi = good;
        while (lo < hi) {
            uint64_t mid = (lo + hi) >> 1;
            auto [ll, rr] = ivl[j - mid];
            if (precedes(Q[mid], ll, rr, true)) hi = mid; else lo = mid + 1;
        }
        return hi;
    }
    // largest j <= cap such that Q[p..j] occurs (truth) or e-occurs, by forward search
    // (Q[p] present).  In strand-symmetric mode the intervals are those of rc(Q[p..j])
    // in the same index, and Q[j] is checked with Psi_E (as rc(Q[j]) is the first phrase
    // of rc(Q[p..j])), every CHECK_EVERY steps and then by binary search.
    uint64_t forward_extend(uint64_t p, uint64_t cap, bool truth = false) {
        uint64_t j = p, l, r;
        if (!rcmode) {
            be.finit(Q[p], l, r);
            while (j < cap && F[j + 1] && fstep(F[j + 1], l, r)) ++j;
            return j;
        }
        be.init(rcQ[p], l, r);
        vector<pair<uint64_t, uint64_t>> ivl{{l, r}};      // ivl[d]: interval of rc(Q[p..p+d])
        uint64_t good = p;
        bool bad = false;
        while (j < cap && F[j + 1]) {
            if (!fstep_same(F[j + 1], l, r)) break;
            ++j;
            ivl.push_back({l, r});
            if (truth && (j - p) % CHECK_EVERY == 0) {
                if (precedes(rcQ[j], l, r, true)) good = j; else { bad = true; break; }
            }
        }
        if (!truth) return j;
        if (!bad && j > good) { auto [ll, rr] = ivl[j - p]; if (precedes(rcQ[j], ll, rr, true)) return j; }
        uint64_t lo = good, hi = bad ? j : j;   // largest x in [good, hi) that occurs; good occurs
        while (hi - lo > 1) {
            uint64_t mid = (lo + hi) >> 1;
            auto [ll, rr] = ivl[mid - p];
            if (precedes(rcQ[mid], ll, rr, true)) lo = mid; else hi = mid;
        }
        return lo;
    }
    inline bool fstep_same(uint32_t c, uint64_t& l, uint64_t& r) { ++fsteps; return be.step(c, l, r); }
    // largest j' in [p, cap] such that Q[p..j'] occurs (truth) or e-occurs
    uint64_t right_end(uint64_t p, uint64_t cap, bool truth) {
        if (rcmode) return forward_extend(p, cap, truth);
        if (rev && !truth) return forward_extend(p, cap);
        uint64_t lo = p, hi = cap + 1;          // Q[p..lo] occurs; Q[p..hi] does not (or hi = cap + 1)
        uint64_t d = 1;
        while (lo + d < hi && occurs(p, lo + d, truth)) { lo += d; d *= 2; }
        if (lo + d < hi) hi = lo + d;
        while (hi - lo > 1) {
            uint64_t mid = (lo + hi) >> 1;
            if (occurs(p, mid, truth)) lo = mid; else hi = mid;
        }
        return lo;
    }
    // MEMs (truth) or e-MEMs of Q[lo..hi], right to left.  For e-MEMs, if etrue is given,
    // also records whether each e-MEM occurs (one rank query on Psi_E, using the interval
    // left by the backward search).
    void mems_range(uint64_t lo, uint64_t hi, bool truth, vector<pair<uint64_t, uint64_t>>& out,
                    vector<bool>* etrue = nullptr) {
        int64_t j = (int64_t)hi;
        while (j >= (int64_t)lo && !Q[j]) --j;
        while (j >= (int64_t)lo) {
            uint64_t i = longest_suffix(j, lo, truth);
            out.push_back({i, (uint64_t)j});
            if (etrue) etrue->push_back(i == (uint64_t)j || precedes(Q[i], last_l, last_r, true));
            if (i == lo) break;
            uint64_t p = i - 1;
            if (!Q[p]) {
                j = (int64_t)p - 1;
                while (j >= (int64_t)lo && !Q[j]) --j;
                continue;
            }
            j = (int64_t)right_end(p, (uint64_t)j - 1, truth);
        }
    }
    // the MEMs inside an e-MEM [a, b] (all of it, if it occurs)
    void refine(uint64_t a, uint64_t b, vector<pair<uint64_t, uint64_t>>& out) {
        if (occurs(a, b, true)) { out.push_back({a, b}); return; }
        mems_range(a, b, true, out);
    }
    // sort, and drop duplicates and intervals contained in others.  (Searching inside
    // each e-MEM separately can report a truncated piece of a MEM that sticks out of
    // that e-MEM; the whole MEM is found from the e-MEM that contains it.)
    static void tidy(vector<pair<uint64_t, uint64_t>>& v) {
        sort(v.begin(), v.end(), [](auto& x, auto& y) { return x.first != y.first ? x.first < y.first : x.second > y.second; });
        vector<pair<uint64_t, uint64_t>> w;
        for (auto& x : v) if (w.empty() || x.second > w.back().second) w.push_back(x);
        v.swap(w);
    }
    vector<pair<uint64_t, uint64_t>> mems() {
        vector<pair<uint64_t, uint64_t>> out;
        if (Q.empty()) return out;
        if (!rev || rcmode) { mems_range(0, Q.size() - 1, use_psi, out); tidy(out); return out; }
        vector<pair<uint64_t, uint64_t>> em;
        vector<bool> et;
        mems_range(0, Q.size() - 1, false, em, use_psi ? &et : nullptr);   // e-MEMs by forward-backward
        if (!use_psi) { tidy(em); return em; }
        for (size_t t = 0; t < em.size(); ++t) {
            if (et[t]) out.push_back(em[t]);        // the e-MEM occurs, so it is a MEM
            else mems_range(em[t].first, em[t].second, true, out);
        }
        tidy(out);
        return out;
    }
    // Boyer-Moore-Li over phrases: the MEMs Q[i..j] with ext(i, j) >= L.  For each start a
    // the window is Q[a..b(a)], b(a) the smallest b with ext(a, b) >= L; every such MEM
    // starting at i contains the window for a = i, and b(a) is nondecreasing in a.  If
    // Q[s..b(a)] does not (e-)occur, no such MEM starts in [a, s].  With the reversed
    // index, windows are tested for e-occurrence, a hit is extended by forward search to
    // its e-MEM, and (with Psi_E) the MEMs inside that e-MEM are found.
    template <class Ext>
    vector<pair<uint64_t, uint64_t>> long_mems(uint64_t L, Ext ext) {
        vector<pair<uint64_t, uint64_t>> out;
        uint64_t np = Q.size();
        bool wtruth = (!rev || rcmode) && use_psi;   // what windows are tested for
        uint64_t a = 0, b = 0;
        while (a < np) {
            if (b < a) b = a;
            while (b < np && ext(a, b) < L) ++b;
            if (b == np) break;
            if (!Q[b]) { a = b + 1; continue; }
            uint64_t s = longest_suffix(b, a, wtruth);
            if (s > a) { a = s; continue; }
            uint64_t j = right_end(a, np - 1, wtruth);
            uint64_t i = longest_suffix(j, 0, wtruth);
            if (rev && use_psi && !rcmode) {
                vector<pair<uint64_t, uint64_t>> in;
                if (i == j || precedes(Q[i], last_l, last_r, true)) in.push_back({i, j});
                else mems_range(i, j, true, in);
                for (auto x : in) if (ext(x.first, x.second) >= L) out.push_back(x);
            } else out.push_back({i, j});
            if (j + 1 >= np) break;
            if (!Q[j + 1]) { a = j + 2; continue; }
            a = max<uint64_t>(i + 1, longest_suffix(j + 1, i + 1, wtruth));
        }
        tidy(out);
        return out;
    }
};

template <class BE>
static int mems_queries(const BE& be, bool have_rev, int argc, char** argv) {
    bool use_psi = true, one_way = false;
    uint64_t minlen = 0;                // -L: only pseudo-MEMs of at least this many characters
    vector<string> qfiles;
    for (int a = 3; a < argc; ++a) {
        string x = argv[a];
        if (x == "--no-psi") use_psi = false;
        else if (x == "--one-way") one_way = true;
        else if (x == "-L" && a + 1 < argc) minlen = strtoull(argv[++a], nullptr, 10);
        else qfiles.push_back(x);
    }
    for (auto& qfn : qfiles) {
        QueryParse q = read_qenc(qfn.c_str());
        auto t0 = clk::now();
        MemFinder<BE> mf(be, q, !one_way, have_rev, BE::exact || use_psi);
        bool exact = BE::exact || use_psi;      // reported parse MEMs are true parse MEMs
        uint64_t np = q.raw.size();
        // Pseudo-MEM bounds.  overlap = 0: a boundary depends on the k characters from it
        // on, so a parse MEM Q[i..j] is extended by one phrase to the left and by one phrase
        // plus k characters to the right.  overlap = 2 (strand-symmetric): phrases run from
        // one syncmer centre to the next, inclusive, and a boundary depends on h = (k-1)/2
        // characters on each side, so the extension is one phrase plus h characters on
        // each side.  (overlap = 1: one phrase, plus the k characters shared, to the right.)
        bool sym = q.overlap == 2;
        uint64_t h = sym ? (uint64_t)(q.k - 1) / 2 : 0;
        auto pstart = [&](uint64_t t) { return q.start[t]; };
        auto pend = [&](uint64_t t) -> uint64_t {      // end (exclusive) of phrase t
            if (q.last[t]) return q.segend[t];
            return q.start[t + 1] + (sym ? 1 : 0);
        };
        uint64_t ext = (uint64_t)q.k;    // extra characters past the end of the extension phrase
        auto left_of = [&](uint64_t t) -> uint64_t {   // start of phrase t minus the margin
            return q.first[t] ? q.start[t] : max<uint64_t>(q.segbeg[t], q.start[t] >= h ? q.start[t] - h : 0);
        };
        // the pseudo-MEM of Q[i..j], clipped to the segment
        auto pm_begin = [&](uint64_t i) -> uint64_t {
            if (q.first[i]) return pstart(i);
            return sym ? left_of(i - 1) : pstart(i - 1);
        };
        auto pm_end = [&](uint64_t j) -> uint64_t {
            if (q.last[j]) return q.segend[j];
            if (sym) return min<uint64_t>(q.segend[j], pend(j + 1) + h);
            if (q.last[j + 1]) return q.segend[j];
            if (q.overlap) return min<uint64_t>(q.segend[j], q.start[j + 2] + ext);
            return min<uint64_t>(q.segend[j], pend(j + 1) + ext);
        };
        auto extlen = [&](uint64_t i, uint64_t j) { return pm_end(j) - pm_begin(i); };
        auto ms = minlen ? mf.long_mems(minlen, extlen) : mf.mems();
        double tm = secs(t0);
        struct PM { uint64_t b, e, lb, i, j; bool pair; };
        vector<PM> pm;
        for (auto [i, j] : ms) {
            uint64_t b = pm_begin(i), en = pm_end(j);
            uint64_t lb = (q.overlap == 1 ? (q.last[j] ? q.segend[j] : min<uint64_t>(q.segend[j], q.start[j + 1] + ext)) : pend(j)) - pstart(i);
            if (!exact) {               // only single phrases and pairs of consecutive phrases are known to occur
                lb = pend(i) - pstart(i);
                for (uint64_t t = i + 1; t <= j; ++t) lb = max<uint64_t>(lb, pend(t) - pstart(t - 1));
            }
            pm.push_back({b, en, lb, i, j, false});
        }
        for (uint64_t t = 0; t < np; ++t) {
            if (mf.Q[t]) continue;
            uint64_t en;
            if (q.last[t]) {            // a segment that is one absent phrase
                if (!q.first[t]) continue;
                en = q.segend[t];
            } else if (!mf.Q[t + 1]) {
                if (sym) en = q.last[t + 1] ? q.segend[t] : min<uint64_t>(q.segend[t], pend(t + 1) + h);
                else en = q.last[t + 1] ? q.segend[t] : min<uint64_t>(q.segend[t], q.overlap ? q.start[t + 2] + ext : pend(t + 1) + ext);
            } else continue;
            uint64_t bg = sym ? left_of(t) : pstart(t);
            if (en - bg >= minlen) pm.push_back({bg, en, 0, t, q.last[t] ? t : t + 1, true});
        }
        sort(pm.begin(), pm.end(), [](const PM& x, const PM& y) { return x.b != y.b ? x.b < y.b : x.e < y.e; });
        string ofn = qfn + be.tag() + (BE::exact || use_psi ? "" : ".nopsi") + (minlen ? ".L" + to_string(minlen) : string()) +
                     (one_way && (have_rev || sym) ? ".oneway" : "") + ".pmems.tsv";
        FILE* fo = fopen(ofn.c_str(), "w");
        if (!fo) { perror(ofn.c_str()); return 1; }
        fprintf(fo, "start\tend\tlength\tlower_bound\tfirst_phrase\tlast_phrase\tkind\n");
        uint64_t covered = 0, cur_b = 0, cur_e = 0, total = 0, maxlb = 0, maxlen = 0, mem_phrases = 0;
        for (auto& x : pm) {
            fprintf(fo, "%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%s\n", (unsigned long long)x.b, (unsigned long long)x.e,
                    (unsigned long long)(x.e - x.b), (unsigned long long)x.lb, (unsigned long long)x.i,
                    (unsigned long long)x.j, x.pair ? "absent-pair" : "parse-MEM");
            total += x.e - x.b; maxlb = max(maxlb, x.lb); maxlen = max(maxlen, x.e - x.b);
            if (!x.pair) mem_phrases += x.j - x.i + 1;
            if (x.b > cur_e) { covered += cur_e - cur_b; cur_b = x.b; cur_e = x.e; }
            else cur_e = max(cur_e, x.e);
        }
        covered += cur_e - cur_b;
        fclose(fo);
        uint64_t nabs = 0; for (auto x : mf.Q) nabs += (x == 0);
        printf("%s: %llu bp, %llu phrases (%llu absent from the dataset); %llu N runs (%llu bp)\n", qfn.c_str(),
               (unsigned long long)q.L, (unsigned long long)np, (unsigned long long)nabs,
               (unsigned long long)q.nruns, (unsigned long long)q.nbases);
        printf("  parse MEMs: %zu, total %llu phrases; absent pairs: %zu\n", ms.size(),
               (unsigned long long)mem_phrases, pm.size() - ms.size());
        printf("  pseudo-MEMs: %zu, total length %llu bp (%.2f x the query), covering %.2f%% of the query; "
               "longest %llu bp, largest lower bound %llu bp\n", pm.size(), (unsigned long long)total,
               (double)total / max<uint64_t>(q.L, 1), 100.0 * covered / max<uint64_t>(q.L, 1),
               (unsigned long long)maxlen, (unsigned long long)maxlb);
        printf("  %s: %llu backward steps, %llu forward steps, %llu Psi_E checks, %.3fs -> %s\n",
               mf.rev ? "forward-backward" : "backward only", (unsigned long long)mf.steps,
               (unsigned long long)mf.fsteps, (unsigned long long)mf.checks, tm, ofn.c_str());
    }
    return 0;
}

template <class RB, class RP>
static int mems_with(FILE* f, const Encoding& e, int argc, char** argv) {
    PsiEIndex<RB, RP> ix;
    ix.enc = &e;
    ix.bwt.load(f); ix.psi.load(f);
    RB rev;
    bool have_rev = false;
    {
        char mg[8];
        if (fread(mg, 1, 8, f) == 8 && !memcmp(mg, "REVBWT01", 8)) { rev.load(f); have_rev = true; }
    }
    fclose(f);
    printf("hybrid index: %s (%s)\n", mb(ix.bytes() + (have_rev ? rev.bytes() : 0)).c_str(),
           have_rev ? "with BWT of the reversed encoding" : "no reversed BWT; forward search via reverse complements with --rc parses");
    EncBackend<RB, RP> be(ix, have_rev ? &rev : nullptr);
    return mems_queries(be, have_rev, argc, argv);
}
static int mems_plain(FILE* f, int argc, char** argv) {
    PlainIndex ix;
    uint8_t c;
    uint8_t hr;
    if (fread(&ix.map.sigma, 4, 1, f) != 1 || fread(&c, 1, 1, f) != 1 || fread(&ix.map.vals_base, 8, 1, f) != 1 ||
        fread(&hr, 1, 1, f) != 1) exit(1);
    ix.map.consecutive = c;
    RunRank::rv(f, ix.map.vals);
    ix.fwd.load(f);
    if (hr) { ix.rev.load(f); ix.have_rev = true; }
    fclose(f);
    if (hr) printf("RLCSA index: %s (forward %s, reversed %s)\n", mb(ix.bytes()).c_str(), mb(ix.fwd.bytes()).c_str(), mb(ix.rev.bytes()).c_str());
    else printf("RLCSA index: %s\n", mb(ix.bytes()).c_str());
    PlainBackend be(ix, hr ? &ix.rev : nullptr);
    return mems_queries(be, hr, argc, argv);
}

static int do_mems(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: hybrid mems <index> [--no-psi] [--one-way] [-L min_length] <query.qenc> ...\n"); return 1; }
    FILE* f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    char mg[8];
    if (fread(mg, 1, 8, f) != 8) { fprintf(stderr, "bad index file\n"); return 1; }
    if (!memcmp(mg, "RLCSA002", 8)) return mems_plain(f, argc, argv);
    Encoding e;
    if (memcmp(mg, "PSIE", 4)) { fprintf(stderr, "bad index file\n"); return 1; }
    e.load(f, mg[7] - '0');
    char k = mg[4];
    if (k == 'X') return mems_with<RunRank, RunRank>(f, e, argc, argv);
    if (k == 'F') return mems_with<EFRank, EFRank>(f, e, argc, argv);
    if (k == 'H') return mems_with<RLHuffRank, EFRank>(f, e, argc, argv);
    if (k == 'C') return mems_with<EscRank, EFRank>(f, e, argc, argv);
    fprintf(stderr, "bad index file\n");
    return 1;
}

template <class I>
static double time_counts(const I& ix, const vector<vector<uint32_t>>& pats, vector<uint64_t>& out) {
    out.assign(pats.size(), 0);
    auto t0 = clk::now();
    for (size_t q = 0; q < pats.size(); ++q) out[q] = ix.count_ids(pats[q]);
    return secs(t0) / pats.size() * 1e6;
}

static int do_bench(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: hybrid bench <input> [width] [queries] [m ...]\n"); return 1; }
    int width = argc > 3 ? atoi(argv[3]) : 4;
    size_t Q = argc > 4 ? atoll(argv[4]) : 100000;
    vector<size_t> ms;
    vector<uint32_t> caps;              // cap=N arguments: also benchmark hybrids with capped rank lists
    for (int a = 5; a < argc; ++a) {
        if (!strncmp(argv[a], "cap=", 4)) caps.push_back((uint32_t)atoi(argv[a] + 4));
        else ms.push_back(atoll(argv[a]));
    }
    if (ms.empty()) ms = {1, 2, 4, 8, 16, 32, 64};

    auto raw = read_ints(argv[2], width);
    Encoding e;
    auto T = densify(raw, e.vals);
    e.sigma = e.vals.size();
    e.compress_vals();
    vector<uint64_t>().swap(raw);
    size_t n0 = T.size();

    fprintf(stderr, "building Psi_E indexes\n");
    vector<uint32_t> L, P;
    e.prepare(T, L, P, true);
    uint64_t rE = count_runs(L), rP = count_runs(P);
    struct Row { string name; uint64_t bytes; function<uint64_t(const vector<uint32_t>&)> count; };
    vector<Row> rows;
    auto t0 = clk::now();
    auto pX = make_shared<PsiEIndex<RunRank>>();              pX->build(e, L, P);
    auto pF = make_shared<PsiEIndex<EFRank>>();               pF->build(e, L, P);
    auto pH = make_shared<PsiEIndex<HuffWT, EFRank>>();       pH->build(e, L, P);
    auto pR = make_shared<PsiEIndex<RLHuffRank, EFRank>>();   pR->build(e, L, P);
    auto pG = make_shared<PsiEIndex<HuffWT, GapRank>>();      pG->build(e, L, P);
    fprintf(stderr, "  rank structures (%.2fs)\n", secs(t0));
    double dE = pH->bwt.avg_depth(), dEr = pR->bwt.heads.avg_depth();
    rows.push_back({"PsiE: BWT(E) explicit", pX->bytes(), [pX](auto& p) { return pX->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) EF", pF->bytes(), [pF](auto& p) { return pF->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) HuffWT", pH->bytes(), [pH](auto& p) { return pH->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) RL-HuffWT", pR->bytes(), [pR](auto& p) { return pR->count_ids(p); }});
    rows.push_back({"PsiE: HuffWT, gap Psi_E", pG->bytes(), [pG](auto& p) { return pG->count_ids(p); }});
    uint64_t psiGap = pG->psi.bytes();
    uint64_t psiEF = pF->psi.bytes(), bwtEF = pF->bwt.bytes(), bwtH = pH->bwt.bytes(), bwtR = pR->bwt.bytes();
    vector<uint32_t>().swap(L); vector<uint32_t>().swap(P);

    vector<shared_ptr<Encoding>> encs;
    vector<string> capinfo;
    for (uint32_t cap : caps) {
        auto ec = make_shared<Encoding>();
        ec->vals = e.vals; ec->consecutive = e.consecutive; ec->vals_base = e.vals_base; ec->sigma = e.sigma;
        ec->cap = cap;
        fprintf(stderr, "building capped hybrids, cap = %u\n", cap);
        vector<uint32_t> Lc, Pc;
        ec->prepare(T, Lc, Pc, true);
        auto cR = make_shared<PsiEIndex<EscRank, EFRank>>(); cR->build(*ec, Lc, Pc);
        auto cF = make_shared<PsiEIndex<EFRank>>();          cF->build(*ec, Lc, Pc);
        auto cP = make_shared<PsiEIndex<EscRankP, GapRank>>(); cP->build(*ec, Lc, Pc);
        uint64_t nesc = 0; for (auto x : Lc) nesc += x > cap;
        char buf[512];
        snprintf(buf, sizeof buf, "cap %u: runs BWT(E) = %llu, Psi_E = %llu; escapes %.2f%%; BWT(E) RL-HuffWT+escapes %s "
                 "(escapes %s), EF %s; Psi_E %s; rank lists %s",
                 cap, (unsigned long long)count_runs(Lc), (unsigned long long)count_runs(Pc), 100.0 * nesc / max<size_t>(n0, 1),
                 mb(cR->bwt.bytes()).c_str(), mb(cR->bwt.has_sub ? cR->bwt.sub.bytes() : 0).c_str(), mb(cF->bwt.bytes()).c_str(),
                 mb(cR->psi.bytes()).c_str(), mb(ec->bytes(false)).c_str());
        capinfo.push_back(string(buf));
        snprintf(buf, sizeof buf, "cap %u, plain: BWT(E) HuffWT+escapes %s (escapes %s), Psi_E gap-coded %s",
                 cap, mb(cP->bwt.bytes()).c_str(), mb(cP->bwt.has_sub ? cP->bwt.sub.bytes() : 0).c_str(), mb(cP->psi.bytes()).c_str());
        capinfo.push_back(buf);
        rows.push_back({"PsiE cap " + to_string(cap) + ": RL-HuffWT+esc", cR->bytes(), [cR](auto& p) { return cR->count_ids(p); }});
        rows.push_back({"PsiE cap " + to_string(cap) + ": EF", cF->bytes(), [cF](auto& p) { return cF->count_ids(p); }});
        rows.push_back({"PsiE cap " + to_string(cap) + ": HuffWT+esc, gap", cP->bytes(), [cP](auto& p) { return cP->count_ids(p); }});
        encs.push_back(ec);
    }

    fprintf(stderr, "building baselines\n");
    vector<uint32_t> LT; vector<uint64_t> CT;
    bwt_of_T(T, e.sigma, LT, CT);
    uint64_t rT = count_runs(LT);
    auto fX = make_shared<RLFM<RunRank>>();    fX->build(LT, CT, e.sigma);
    auto fF = make_shared<RLFM<EFRank>>();     fF->build(LT, CT, e.sigma);
    auto fH = make_shared<RLFM<HuffWT>>();     fH->build(LT, CT, e.sigma);
    auto fR = make_shared<RLFM<RLHuffRank>>(); fR->build(LT, CT, e.sigma);
    auto fG = make_shared<RLFM<GapRank>>();    fG->build(LT, CT, e.sigma);
    double dT = fH->bwt.avg_depth(), dTr = fR->bwt.heads.avg_depth();
    rows.push_back({"FM(T): explicit", fX->bytes(), [fX](auto& p) { return fX->count_ids(p); }});
    rows.push_back({"FM(T): EF", fF->bytes(), [fF](auto& p) { return fF->count_ids(p); }});
    rows.push_back({"FM(T): HuffWT", fH->bytes(), [fH](auto& p) { return fH->count_ids(p); }});
    rows.push_back({"FM(T): RL-HuffWT", fR->bytes(), [fR](auto& p) { return fR->count_ids(p); }});
    rows.push_back({"CSA(T): gap Psi", fG->bytes(), [fG](auto& p) { return fG->count_ids(p); }});
    vector<uint32_t>().swap(LT);

    printf("n = %zu, sigma = %u, H0(T) = %.2f bits, E alphabet = %u, distinct bigrams = %zu\n",
           n0, e.sigma, H0(T, e.sigma + 1), e.kE, (size_t)e.nbigrams);
    printf("runs: BWT(T) = %llu; BWT(E) = %llu, Psi_E = %llu\n",
           (unsigned long long)rT, (unsigned long long)rE, (unsigned long long)rP);
    printf("Huffman WT depth (bits/symbol): BWT(E) %.2f, run heads of BWT(E) %.2f; BWT(T) %.2f, run heads of BWT(T) %.2f\n",
           dE, dEr, dT, dTr);
    printf("BWT(E) alone: explicit %s, EF %s, HuffWT %s, RL-HuffWT %s; Psi_E (EF) %s; rankings %s\n",
           mb(pX->bwt.bytes()).c_str(), mb(bwtEF).c_str(), mb(bwtH).c_str(), mb(bwtR).c_str(),
           mb(psiEF).c_str(), mb(e.bytes(false)).c_str());
    printf("Psi_E gap-coded (plain, no runs) %s\n", mb(psiGap).c_str());
    for (auto& x : capinfo) printf("%s\n", x.c_str());

    mt19937_64 rng(12345);
    if (n0 <= 5000000) {
        size_t bad = 0;
        for (int t = 0; t < 30; ++t) {
            size_t m = 1 + rng() % 8;
            if (m > n0) break;
            size_t s = rng() % (n0 - m + 1);
            vector<uint32_t> p(T.begin() + s, T.begin() + s + m);
            if (t % 3 == 0) p[rng() % m] = 1 + rng() % e.sigma;
            uint64_t bf = 0;
            for (size_t k = 0; k + m <= n0; ++k) bf += equal(p.begin(), p.end(), T.begin() + k);
            for (auto& r : rows) if (r.count(p) != bf) ++bad;
        }
        printf("brute-force check: %s\n", bad ? "FAILED" : "ok");
    }

    vector<vector<vector<uint32_t>>> pats;
    for (size_t m : ms) {
        if (m > n0) continue;
        vector<vector<uint32_t>> ps(Q);
        for (auto& p : ps) { size_t s = rng() % (n0 - m + 1); p.assign(T.begin() + s, T.begin() + s + m); }
        pats.push_back(move(ps));
    }
    printf("\ncount time, microseconds per query (%zu random substrings of T per m)\n", Q);
    printf("%-24s %10s", "index", "size");
    for (auto& ps : pats) printf(" %8s", ("m=" + to_string(ps[0].size())).c_str());
    printf("\n");
    vector<vector<uint64_t>> ref(pats.size());
    bool mism = false;
    for (auto& r : rows) {
        printf("%-24s %10s", r.name.c_str(), mb(r.bytes).c_str());
        for (size_t k = 0; k < pats.size(); ++k) {
            vector<uint64_t> out(Q);
            auto t1 = clk::now();
            for (size_t q = 0; q < Q; ++q) out[q] = r.count(pats[k][q]);
            printf(" %8.3f", secs(t1) / Q * 1e6);
            if (ref[k].empty()) ref[k] = out; else if (out != ref[k]) mism = true;
        }
        printf("\n");
        fflush(stdout);
    }
    if (mism) printf("MISMATCH between indexes\n");
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: hybrid build|count|bench|mems ...\n"); return 1; }
    string mode = argv[1];
    if (mode == "build") return do_build(argc, argv);
    if (mode == "count") return do_count(argc, argv);
    if (mode == "bench") return do_bench(argc, argv);
    if (mode == "mems") return do_mems(argc, argv);
    fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 1;
}
