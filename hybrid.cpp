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
    //   off[a]: index of a's first successor (32 bits per phrase);
    //   csucc: the successors, ceil(lg(sigma+1)) bits each.  For a phrase with
    //     d <= SMALL successors they are stored in rank order, so the rank of b
    //     is its position + 1 (found by a short linear scan, with no rank stored).
    //     For d > SMALL they are sorted by id (found by binary search), and
    //     their ranks are stored in bigranks, wbig bits each; bigA lists those
    //     phrases and bigstart says where their ranks start.
    // A lookup costs about one cache miss for small degrees.
    static constexpr uint64_t SMALL = 16;
    bool compact = false;
    vector<uint32_t> off;
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
    uint32_t succ_rank(uint32_t a, uint32_t b) const {
        if (!compact) {
            const uint32_t* lo = sb.data() + soff[a];
            const uint32_t* hi = sb.data() + soff[a + 1];
            const uint32_t* it = lower_bound(lo, hi, b);
            if (it == hi || *it != b) return 0;
            return srank[it - sb.data()];
        }
        if (a > sigma) return 0;
        uint64_t lo = off[a], hi = off[a + 1], d = hi - lo;
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
        off.assign(soff.begin(), soff.end());
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
            for (size_t q = 0; q < grp.size(); ++q) { sb.push_back(grp[q].second); srank.push_back(rk[q]); }
            kE = max<uint32_t>(kE, grp.size() + 1);
        }
        soff[sigma + 1] = sb.size();
        vector<uint64_t>().swap(pr);

        vector<int32_t> S(n0);
        for (uint64_t k = 0; k + 1 < n0; ++k) S[k] = succ_rank(T[k], T[k + 1]);
        if (n0) S[n0 - 1] = 0;
        if (verbose) fprintf(stderr, "  encoded (%.2fs): %zu distinct bigrams, E alphabet %u\n",
                             secs(t0), sb.size(), kE);
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
    uint64_t bytes() const {
        return vals.size() * 8 + C.size() * 8 + off.size() * 4 + csucc.a.size() * 8 + bigranks.a.size() * 8 +
               bigA.size() * 4 + bigstart.size() * 8 + 64;
    }
    void save(FILE* f) const {
        fwrite(&sigma, 4, 1, f); fwrite(&kE, 4, 1, f); fwrite(&n0, 8, 1, f);
        uint8_t c = consecutive; fwrite(&c, 1, 1, f); fwrite(&vals_base, 8, 1, f);
        fwrite(&wsucc, 4, 1, f); fwrite(&wbig, 4, 1, f); fwrite(&nbigrams, 8, 1, f);
        RunRank::wv(f, vals); RunRank::wv(f, C); RunRank::wv(f, off);
        RunRank::wv(f, csucc.a); RunRank::wv(f, bigranks.a); RunRank::wv(f, bigA); RunRank::wv(f, bigstart);
    }
    void load(FILE* f) {
        uint8_t c;
        if (fread(&sigma, 4, 1, f) != 1 || fread(&kE, 4, 1, f) != 1 || fread(&n0, 8, 1, f) != 1 ||
            fread(&c, 1, 1, f) != 1 || fread(&vals_base, 8, 1, f) != 1 || fread(&wsucc, 4, 1, f) != 1 ||
            fread(&wbig, 4, 1, f) != 1 || fread(&nbigrams, 8, 1, f) != 1) exit(1);
        consecutive = c;
        RunRank::rv(f, vals); RunRank::rv(f, C); RunRank::rv(f, off);
        RunRank::rv(f, csucc.a); RunRank::rv(f, bigranks.a); RunRank::rv(f, bigA); RunRank::rv(f, bigstart);
        compact = true;
    }
};

// ---------------------------------------------------------------- index
template <class RB, class RP = RB>
struct PsiEIndex {
    const Encoding* enc = nullptr;
    RB bwt;
    RP psi;

    void build(const Encoding& e, const vector<uint32_t>& L, const vector<uint32_t>& P) {
        enc = &e;
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
            uint32_t c = enc->succ_rank(p[j - 1], p[j]);
            if (c == 0) return 0;
            uint64_t rl, rr;
            bwt.rank2(c, l, r, rl, rr);
            l = enc->C[c] + rl; r = enc->C[c] + rr;
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
    uint64_t bytes() const { return enc->bytes() + bwt.bytes() + psi.bytes(); }
};

// Baseline: run-length FM-index of T$ with the same rank structure.
template <class R>
struct RLFM {
    uint64_t n = 0;
    vector<uint64_t> C;
    R bwt;
    void build(const vector<uint32_t>& L, const vector<uint64_t>& C_, uint32_t sigma) {
        n = L.size(); C = C_;
        bwt.build(L, sigma + 1);
    }
    uint64_t count_ids(const vector<uint32_t>& p) const {
        uint64_t l = 0, r = n;
        for (size_t j = p.size(); j-- > 0;) {
            uint32_t c = p[j];
            if (c == 0) return 0;
            uint64_t rl, rr;
            bwt.rank2(c, l, r, rl, rr);
            l = C[c] + rl; r = C[c] + rr;
            if (l >= r) return 0;
        }
        return r - l;
    }
    uint64_t bytes() const { return C.size() * 8 + bwt.bytes(); }
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
    if (argc < 4) { fprintf(stderr, "usage: hybrid build <input> <index> [width] [huff|ef|explicit|bihuff]\n"); return 1; }
    int width = argc > 4 ? atoi(argv[4]) : 4;
    string kind = argc > 5 ? argv[5] : "huff";
    if (kind != "explicit" && kind != "ef" && kind != "huff" && kind != "bihuff") {
        fprintf(stderr, "kind must be explicit, ef, huff or bihuff\n"); return 1;
    }
    auto t0 = clk::now();
    auto raw = read_ints(argv[2], width);
    Encoding e;
    auto T = densify(raw, e.vals);
    e.sigma = e.vals.size();
    e.compress_vals();
    vector<uint64_t>().swap(raw);
    vector<uint32_t> L, P, Lrev;
    e.prepare(T, L, P, true, kind == "bihuff" ? &Lrev : nullptr);
    FILE* f = fopen(argv[3], "wb");
    if (!f) { perror(argv[3]); return 1; }
    fwrite(kind == "explicit" ? "PSIEX003" : kind == "ef" ? "PSIEF003" : "PSIEH003", 1, 8, f);
    e.save(f);
    uint64_t sz;
    if (kind == "explicit") { PsiEIndex<RunRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else if (kind == "ef")  { PsiEIndex<EFRank> ix;  ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else { PsiEIndex<RLHuffRank, EFRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
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
    e.load(f);
    if (!memcmp(mg, "PSIEX003", 8)) return count_with<RunRank, RunRank>(f, e, argv[3]);
    if (!memcmp(mg, "PSIEF003", 8)) return count_with<EFRank, EFRank>(f, e, argv[3]);
    if (!memcmp(mg, "PSIEH003", 8)) return count_with<RLHuffRank, EFRank>(f, e, argv[3]);
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
    vector<uint64_t> segend;            // per phrase: end of its segment
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
    }
    fclose(f);
    q.raw.resize(np); q.start.resize(np); q.segend.resize(np); q.first.resize(np); q.last.resize(np);
    size_t r = 0;
    for (uint64_t t = 0; t < np; ++t) {
        q.raw[t] = rec[4 * t]; q.start[t] = rec[4 * t + 3];
        uint32_t fl = rec[4 * t + 2];
        q.first[t] = t == 0 || (v3 && (fl & 1));
        q.last[t] = t + 1 == np || (v3 && (fl & 2));
        while (r < q.nruns && nr[2 * r] <= q.start[t]) ++r;
        q.segend[t] = r < q.nruns ? nr[2 * r] : q.L;
    }
    return q;
}

template <class RB, class RP>
struct MemFinder {
    const PsiEIndex<RB, RP>& ix;
    const Encoding& e;
    const RB* rev = nullptr;            // index of the reversed encoding (forward search), if any
    vector<uint32_t> Q;                 // dense ids, 0 = absent
    vector<uint32_t> R;                 // R[a] = rank of Q[a] after Q[a-1], 0 if the pair never occurs
    uint64_t steps = 0, fsteps = 0, checks = 0;   // backward steps, forward steps, rank queries on Psi_E
    static constexpr uint64_t CHECK_EVERY = 16;
    bool use_psi = true;                // false: encoded matches count as matches (no Psi_E)
    uint64_t last_l = 0, last_r = 0;    // interval of R[a+1..j] for the start a last returned by an encoded longest_suffix

    // "e-occurs": Q[a..b] e-occurs if Q[a] is present and R[a+1..b] occurs in the
    // encoding.  Every occurrence is an e-occurrence; Q[a..b] occurs if moreover some
    // e-occurrence is preceded by Q[a] (one rank query on Psi_E).
    MemFinder(const PsiEIndex<RB, RP>& ix_, const QueryParse& q) : ix(ix_), e(*ix_.enc) {
        uint64_t np = q.raw.size();
        Q.assign(np, 0); R.assign(np, 0);
        for (uint64_t t = 0; t < np; ++t) {
            if (q.raw[t] == 0xFFFFFFFFu) continue;
            int64_t id = e.id_of(q.raw[t]);
            Q[t] = id < 0 ? 0 : (uint32_t)id;
        }
        for (uint64_t t = 1; t < np; ++t) if (Q[t - 1] && Q[t] && !q.first[t]) R[t] = e.succ_rank(Q[t - 1], Q[t]);
    }
    inline bool step(uint32_t c, uint64_t& l, uint64_t& r) {        // backward step in BWT(E)
        ++steps;
        uint64_t rl, rr;
        ix.bwt.rank2(c, l, r, rl, rr);
        l = e.C[c] + rl; r = e.C[c] + rr;
        return l < r;
    }
    inline bool fstep(uint32_t c, uint64_t& l, uint64_t& r) {       // backward step in BWT(E reversed)
        ++fsteps;
        uint64_t rl, rr;
        rev->rank2(c, l, r, rl, rr);
        l = e.C[c] + rl; r = e.C[c] + rr;       // same symbol counts as E
        return l < r;
    }
    inline bool precedes(uint32_t a, uint64_t l, uint64_t r, bool truth) {
        if (!truth) return true;
        ++checks;
        uint64_t rl, rr;
        ix.psi.rank2(a, l, r, rl, rr);
        return rr > rl;
    }
    bool occurs(uint64_t a, uint64_t b, bool truth) {
        if (!Q[a]) return false;
        uint64_t l = 0, r = e.n0;
        for (uint64_t x = b; x > a; --x) if (!R[x] || !step(R[x], l, r)) return false;
        return a == b || precedes(Q[a], l, r, truth);
    }
    // smallest i >= floor such that Q[i..j] occurs (truth) or e-occurs; Q[j] must be present
    uint64_t longest_suffix(uint64_t j, uint64_t floor, bool truth) {
        vector<pair<uint64_t, uint64_t>> ivl{{0, e.n0}};   // ivl[d]: interval of R[j-d+1..j]
        uint64_t good = j, a = j, l = 0, r = e.n0;
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
    // largest j <= cap such that Q[p..j] e-occurs, by forward search (Q[p] present)
    uint64_t forward_extend(uint64_t p, uint64_t cap) {
        uint64_t j = p, l = 0, r = e.n0;
        while (j < cap && R[j + 1] && fstep(R[j + 1], l, r)) ++j;
        return j;
    }
    // largest j' in [p, cap] such that Q[p..j'] occurs (truth) or e-occurs
    uint64_t right_end(uint64_t p, uint64_t cap, bool truth) {
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
        if (!rev) { mems_range(0, Q.size() - 1, use_psi, out); tidy(out); return out; }
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
        bool wtruth = !rev && use_psi;          // what windows are tested for
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
            if (rev && use_psi) {
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
        MemFinder<RB, RP> mf(ix, q);
        mf.use_psi = use_psi;
        if (have_rev && !one_way) mf.rev = &rev;
        uint64_t np = q.raw.size();
        auto pstart = [&](uint64_t t) { return q.start[t]; };
        auto pend = [&](uint64_t t) { return q.last[t] ? q.segend[t] : q.start[t + 1]; };   // non-overlapping end
        uint64_t ext = (uint64_t)q.k;    // extra characters past the end of the extension phrase
        // the pseudo-MEM of Q[i..j], clipped to the segment
        auto pm_begin = [&](uint64_t i) { return pstart(q.first[i] ? i : i - 1); };
        auto pm_end = [&](uint64_t j) -> uint64_t {
            if (q.last[j] || q.last[j + 1]) return q.segend[j];
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
            uint64_t lb = (q.overlap ? (q.last[j] ? q.segend[j] : min<uint64_t>(q.segend[j], q.start[j + 1] + ext)) : pend(j)) - pstart(i);
            if (!use_psi) {             // only single phrases and pairs of consecutive phrases are known to occur
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
                en = q.last[t + 1] ? q.segend[t] : min<uint64_t>(q.segend[t], q.overlap ? q.start[t + 2] + ext : pend(t + 1) + ext);
            } else continue;
            if (en - pstart(t) >= minlen) pm.push_back({pstart(t), en, 0, t, q.last[t] ? t : t + 1, true});
        }
        sort(pm.begin(), pm.end(), [](const PM& x, const PM& y) { return x.b != y.b ? x.b < y.b : x.e < y.e; });
        string ofn = qfn + (use_psi ? "" : ".nopsi") + (minlen ? ".L" + to_string(minlen) : string()) +
                     (have_rev && one_way ? ".oneway" : "") + ".pmems.tsv";
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

static int do_mems(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: hybrid mems <index> [--no-psi] [--one-way] [-L min_length] <query.qenc> ...\n"); return 1; }
    FILE* f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    char mg[8];
    if (fread(mg, 1, 8, f) != 8) { fprintf(stderr, "bad index file\n"); return 1; }
    Encoding e;
    e.load(f);
    if (!memcmp(mg, "PSIEX003", 8)) return mems_with<RunRank, RunRank>(f, e, argc, argv);
    if (!memcmp(mg, "PSIEF003", 8)) return mems_with<EFRank, EFRank>(f, e, argc, argv);
    if (!memcmp(mg, "PSIEH003", 8)) return mems_with<RLHuffRank, EFRank>(f, e, argc, argv);
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
    for (int a = 5; a < argc; ++a) ms.push_back(atoll(argv[a]));
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
    fprintf(stderr, "  rank structures (%.2fs)\n", secs(t0));
    double dE = pH->bwt.avg_depth(), dEr = pR->bwt.heads.avg_depth();
    rows.push_back({"PsiE: BWT(E) explicit", pX->bytes(), [pX](auto& p) { return pX->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) EF", pF->bytes(), [pF](auto& p) { return pF->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) HuffWT", pH->bytes(), [pH](auto& p) { return pH->count_ids(p); }});
    rows.push_back({"PsiE: BWT(E) RL-HuffWT", pR->bytes(), [pR](auto& p) { return pR->count_ids(p); }});
    uint64_t psiEF = pF->psi.bytes(), bwtEF = pF->bwt.bytes(), bwtH = pH->bwt.bytes(), bwtR = pR->bwt.bytes();
    vector<uint32_t>().swap(L); vector<uint32_t>().swap(P);

    fprintf(stderr, "building baselines\n");
    vector<uint32_t> LT; vector<uint64_t> CT;
    bwt_of_T(T, e.sigma, LT, CT);
    uint64_t rT = count_runs(LT);
    auto fX = make_shared<RLFM<RunRank>>();    fX->build(LT, CT, e.sigma);
    auto fF = make_shared<RLFM<EFRank>>();     fF->build(LT, CT, e.sigma);
    auto fH = make_shared<RLFM<HuffWT>>();     fH->build(LT, CT, e.sigma);
    auto fR = make_shared<RLFM<RLHuffRank>>(); fR->build(LT, CT, e.sigma);
    double dT = fH->bwt.avg_depth(), dTr = fR->bwt.heads.avg_depth();
    rows.push_back({"FM(T): explicit", fX->bytes(), [fX](auto& p) { return fX->count_ids(p); }});
    rows.push_back({"FM(T): EF", fF->bytes(), [fF](auto& p) { return fF->count_ids(p); }});
    rows.push_back({"FM(T): HuffWT", fH->bytes(), [fH](auto& p) { return fH->count_ids(p); }});
    rows.push_back({"FM(T): RL-HuffWT", fR->bytes(), [fR](auto& p) { return fR->count_ids(p); }});
    vector<uint32_t>().swap(LT);

    printf("n = %zu, sigma = %u, H0(T) = %.2f bits, E alphabet = %u, distinct bigrams = %zu\n",
           n0, e.sigma, H0(T, e.sigma + 1), e.kE, (size_t)e.nbigrams);
    printf("runs: BWT(T) = %llu; BWT(E) = %llu, Psi_E = %llu\n",
           (unsigned long long)rT, (unsigned long long)rE, (unsigned long long)rP);
    printf("Huffman WT depth (bits/symbol): BWT(E) %.2f, run heads of BWT(E) %.2f; BWT(T) %.2f, run heads of BWT(T) %.2f\n",
           dE, dEr, dT, dTr);
    printf("BWT(E) alone: explicit %s, EF %s, HuffWT %s, RL-HuffWT %s; Psi_E (EF) %s; rankings %s\n",
           mb(pX->bwt.bytes()).c_str(), mb(bwtEF).c_str(), mb(bwtH).c_str(), mb(bwtR).c_str(),
           mb(psiEF).c_str(), mb(e.bytes()).c_str());

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
