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
struct Encoding {
    vector<uint64_t> vals;          // vals[id-1] = raw value of dense id
    uint32_t sigma = 0;             // number of distinct input values
    vector<uint64_t> soff;          // sigma+2; for a: pairs (b, rank) sorted by b
    vector<uint32_t> sb, srank;
    uint32_t kE = 0;                // alphabet size of S = max rank + 1
    uint64_t n0 = 0;                // |S|
    vector<uint64_t> C;             // kE+1

    int64_t id_of(uint64_t raw) const {
        auto it = lower_bound(vals.begin(), vals.end(), raw);
        if (it == vals.end() || *it != raw) return -1;
        return (it - vals.begin()) + 1;
    }
    uint32_t succ_rank(uint32_t a, uint32_t b) const {
        const uint32_t* lo = sb.data() + soff[a];
        const uint32_t* hi = sb.data() + soff[a + 1];
        const uint32_t* it = lower_bound(lo, hi, b);
        if (it == hi || *it != b) return 0;
        return srank[it - sb.data()];
    }

    // From T (dense ids, no terminator) compute the rankings, S, SA(S), BWT(S) and Psi_E.
    void prepare(const vector<uint32_t>& T, vector<uint32_t>& L, vector<uint32_t>& P, bool verbose) {
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
    }
    uint64_t bytes() const { return vals.size() * 8 + soff.size() * 8 + sb.size() * 4 + srank.size() * 4 + C.size() * 8; }
    void save(FILE* f) const {
        fwrite(&sigma, 4, 1, f); fwrite(&kE, 4, 1, f); fwrite(&n0, 8, 1, f);
        RunRank::wv(f, vals); RunRank::wv(f, soff); RunRank::wv(f, sb); RunRank::wv(f, srank); RunRank::wv(f, C);
    }
    void load(FILE* f) {
        if (fread(&sigma, 4, 1, f) != 1 || fread(&kE, 4, 1, f) != 1 || fread(&n0, 8, 1, f) != 1) exit(1);
        RunRank::rv(f, vals); RunRank::rv(f, soff); RunRank::rv(f, sb); RunRank::rv(f, srank); RunRank::rv(f, C);
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
    if (argc < 4) { fprintf(stderr, "usage: hybrid build <input> <index> [width] [huff|ef|explicit]\n"); return 1; }
    int width = argc > 4 ? atoi(argv[4]) : 4;
    string kind = argc > 5 ? argv[5] : "huff";
    if (kind != "explicit" && kind != "ef" && kind != "huff") { fprintf(stderr, "kind must be explicit, ef or huff\n"); return 1; }
    auto t0 = clk::now();
    auto raw = read_ints(argv[2], width);
    Encoding e;
    auto T = densify(raw, e.vals);
    e.sigma = e.vals.size();
    vector<uint64_t>().swap(raw);
    vector<uint32_t> L, P;
    e.prepare(T, L, P, true);
    FILE* f = fopen(argv[3], "wb");
    if (!f) { perror(argv[3]); return 1; }
    fwrite(kind == "explicit" ? "PSIEX002" : kind == "ef" ? "PSIEF002" : "PSIEH002", 1, 8, f);
    e.save(f);
    uint64_t sz;
    if (kind == "explicit") { PsiEIndex<RunRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else if (kind == "ef")  { PsiEIndex<EFRank> ix;  ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
    else { PsiEIndex<RLHuffRank, EFRank> ix; ix.build(e, L, P); ix.bwt.save(f); ix.psi.save(f); sz = ix.bytes(); }
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
    if (!memcmp(mg, "PSIEX002", 8)) return count_with<RunRank, RunRank>(f, e, argv[3]);
    if (!memcmp(mg, "PSIEF002", 8)) return count_with<EFRank, EFRank>(f, e, argv[3]);
    if (!memcmp(mg, "PSIEH002", 8)) return count_with<RLHuffRank, EFRank>(f, e, argv[3]);
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
           n0, e.sigma, H0(T, e.sigma + 1), e.kE, e.sb.size());
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
    if (argc < 2) { fprintf(stderr, "usage: hybrid build|count|bench ...\n"); return 1; }
    string mode = argv[1];
    if (mode == "build") return do_build(argc, argv);
    if (mode == "count") return do_count(argc, argv);
    if (mode == "bench") return do_bench(argc, argv);
    fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 1;
}
