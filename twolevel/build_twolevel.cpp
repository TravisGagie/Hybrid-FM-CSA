// Copyright (C) 2026 Travis Gagie  (GPL v3 or later; see LICENSE)
//
// Builds the parse, block and mapping levels over an RLBWT of the
// reportable-shortcuts code.
//
//   build_twolevel --text COLLECTION --rlbwt DIR/rlbwt.bin --reversed-rlbwt REVDIR/rlbwt.bin --out PREFIX
//
// COLLECTION is the collection file the RLBWTs were built from (one sequence per
// line); the second RLBWT is that of the reversed collection
// rev(S_k) # ... # rev(S_1) # $.  Writes PREFIX.text.lvl and PREFIX.rev.lvl; the
// first holds the phrase dictionary of both, which give each phrase of T and
// its reverse in the reversed collection.
#include <chrono>
#include <iostream>
#include <unordered_map>

#include "libsais.h"
#include "text.hpp"
#include "twolevel.hpp"

using namespace tl;
using Clock = std::chrono::steady_clock;

static double secs(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// The text T = S_1 # ... # S_k # $ of a collection file with one sequence per
// line, as rs::read_text builds it (empty lines skipped, CR dropped), but read
// in one piece so that the memory is not doubled while it grows.
static std::vector<uint8_t> read_collection(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    fseek(f, 0, SEEK_END);
    long long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    std::vector<uint8_t> X((size_t)sz + 2);
    if (sz > 0 && fread(X.data(), 1, (size_t)sz, f) != (size_t)sz) throw std::runtime_error("cannot read " + path);
    fclose(f);
    if (sz > 0 && X[0] == '>') throw std::runtime_error(path + ": FASTA is not supported; use one sequence per line");
    size_t w = 0;
    bool in_line = false;
    for (size_t i = 0; i < (size_t)sz; ++i) {
        uint8_t c = X[i];
        if (c == '\r') continue;
        if (c == '\n') { if (in_line) X[w++] = rs::kSeparator; in_line = false; continue; }
        if (c == rs::kSeparator || c == rs::kTerminator) throw std::runtime_error(path + " contains a byte 0x00 or 0x01");
        X[w++] = c; in_line = true;
    }
    if (in_line) X[w++] = rs::kSeparator;
    X[w++] = rs::kTerminator;
    X.resize(w);
    X.shrink_to_fit();
    return X;
}

// Builds the levels of one direction.  hf[id] and hr[id] are the hashes of
// phrase id and of its reverse, for joining the two directions' dictionaries.
static void build_level(const std::string& text_path, bool reversed, const std::string& rlbwt_path, Level& L,
                        std::vector<uint64_t>& hf, std::vector<uint64_t>& hr, Clock::time_point t0) {
    std::vector<uint8_t> X = read_collection(text_path);
    if (reversed) {   // S_1 # ... # S_k # $  ->  rev(S_k) # ... # rev(S_1) # $
        std::reverse(X.begin(), X.end() - 1);
        std::rotate(X.begin(), X.begin() + 1, X.end() - 1);
    }
    const uint64_t n = X.size();
    std::cerr << "text: " << n << " characters (" << secs(t0) << " s)\n";

    // boundaries
    std::vector<uint64_t> b;
    bool zero_is_syncmer = false;
    for_each_syncmer(X.data(), n, [&](uint64_t q) { if (q == 0) zero_is_syncmer = true; b.push_back(q); });
    if (!zero_is_syncmer) b.insert(b.begin(), 0);
    const uint64_t m = b.size();
    if (m + 1 >= (1ull << 31)) throw std::runtime_error("too many phrases for a 32-bit parse suffix array");
    auto plen = [&](uint64_t i) -> uint64_t { return i + 1 < m ? b[i + 1] + K - b[i] : n - b[i]; };
    std::cerr << "boundaries: " << m << " (" << secs(t0) << " s)\n";

    // distinct phrases
    std::vector<uint32_t> pz(m);
    std::vector<uint32_t> rep;        // first occurrence (phrase index) of each distinct phrase
    std::vector<uint64_t> rhash;
    {
        std::unordered_map<uint64_t, uint32_t> seen;
        seen.reserve(m / 8);
        for (uint64_t i = 0; i < m; ++i) {
            uint64_t len = plen(i);
            uint64_t hv = phrase_hash(X.data() + b[i], len);
            auto it = seen.find(hv);
            if (it == seen.end()) {
                seen.emplace(hv, (uint32_t)rep.size());
                pz[i] = (uint32_t)rep.size();
                rep.push_back((uint32_t)i); rhash.push_back(hv);
            } else {
                uint32_t r0 = rep[it->second];
                if (plen(r0) != len || std::memcmp(X.data() + b[r0], X.data() + b[i], len))
                    throw std::runtime_error("phrase hash collision; change phrase_hash");
                pz[i] = it->second;
            }
        }
    }
    const uint32_t sigma = (uint32_t)rep.size();
    std::cerr << "distinct phrases: " << sigma << " (" << secs(t0) << " s)\n";

    // lexicographic ids 1..sigma (byte order; the dictionary is prefix-free)
    std::vector<uint32_t> order(sigma);
    for (uint32_t a = 0; a < sigma; ++a) order[a] = a;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t c) {
        uint64_t la = plen(rep[a]), lc = plen(rep[c]);
        int r = std::memcmp(X.data() + b[rep[a]], X.data() + b[rep[c]], std::min(la, lc));
        if (r == 0 && la != lc) {
            // a proper prefix: cannot happen if the dictionary is prefix-free
            throw std::runtime_error("dictionary is not prefix-free");
        }
        return r < 0;
    });
    std::vector<uint32_t> idof(sigma);
    for (uint32_t r = 0; r < sigma; ++r) idof[order[r]] = r + 1;
    for (auto& x : pz) x = idof[x];
    L.n = n; L.m = m; L.sigma = sigma;
    hf.assign(sigma + 1, 0); hr.assign(sigma + 1, 0);
    {
        std::vector<uint8_t> buf;
        for (uint32_t a = 0; a < sigma; ++a) {
            uint64_t len = plen(rep[a]);
            const uint8_t* p = X.data() + b[rep[a]];
            buf.assign(p, p + len);
            std::reverse(buf.begin(), buf.end());
            hf[idof[a]] = rhash[a];
            hr[idof[a]] = phrase_hash(buf.data(), len);
        }
    }
    std::cerr << "sorted (" << secs(t0) << " s)\n";

    // distinct syncmers: phrases with equal first K characters are adjacent in id order
    std::vector<uint64_t> synrep;     // a text position of each distinct syncmer
    {
        uint32_t skip = zero_is_syncmer ? 0 : pz[0];
        int64_t prev = -1;  // position of the previous distinct syncmer
        for (uint32_t r = 0; r < sigma; ++r) {
            uint32_t a = order[r];
            if (r + 1 == skip) continue;
            uint64_t pos = b[rep[a]];
            if (prev < 0 || std::memcmp(X.data() + prev, X.data() + pos, K)) { synrep.push_back(pos); prev = pos; }
        }
    }
    const size_t nsyn = synrep.size();
    std::vector<uint8_t> syn(nsyn * K);
    for (size_t k = 0; k < nsyn; ++k) std::memcpy(syn.data() + k * K, X.data() + synrep[k], K);
    std::vector<uint64_t>().swap(synrep);
    std::vector<uint32_t>().swap(rep); std::vector<uint64_t>().swap(rhash);
    std::vector<uint32_t>().swap(order); std::vector<uint32_t>().swap(idof);
    std::vector<uint8_t>().swap(X);
    std::vector<uint64_t>().swap(b);
    std::cerr << "distinct syncmers: " << nsyn << "\n";

    // B: the union of the syncmers' intervals, and the row of position 0
    {
        rs::RLBWT rl = rs::RLBWT::load(rlbwt_path);
        if ((uint64_t)rl.size() != n) throw std::runtime_error("the RLBWT does not have n rows; wrong --text, or the two RLBWTs swapped?");
        std::vector<std::pair<uint64_t, uint64_t>> iv;
        iv.reserve(nsyn + 1);
        for (size_t k = 0; k < nsyn; ++k) {
            const uint8_t* s = syn.data() + k * K;
            int64_t lo = 0, hi = rl.size() - 1;
            for (int c = K - 1; c >= 0 && lo <= hi; --c) rl.backward_step(s[c], lo, hi);
            if (hi < lo) throw std::runtime_error("a syncmer does not occur in the RLBWT");
            iv.push_back({(uint64_t)lo, (uint64_t)(hi - lo + 1)});
        }
        if (!zero_is_syncmer) {     // the row whose BWT character is $
            int64_t row = -1;
            for (int64_t run = 0; run < rl.num_runs(); ++run)
                if (rl.run_char(run) == rs::kTerminator) { row = rl.run_start(run); break; }
            if (row < 0) throw std::runtime_error("no $ in the RLBWT");
            iv.push_back({(uint64_t)row, 1});
        }
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<uint64_t, uint64_t>> runs;
        uint64_t total = 0;
        for (auto& x : iv) {
            if (!runs.empty() && runs.back().first + runs.back().second > x.first) throw std::runtime_error("overlapping syncmer intervals");
            if (!runs.empty() && runs.back().first + runs.back().second == x.first) runs.back().second += x.second;
            else runs.push_back(x);
            total += x.second;
        }
        if (total != m) throw std::runtime_error("B has " + std::to_string(total) + " ones, not m = " + std::to_string(m) +
                                                 "; is the RLBWT that of this text?");
        L.B.build(runs, n);
        std::cerr << "B: " << runs.size() << " runs, " << L.B.bytes() / 1e6 << " MB (" << secs(t0) << " s)\n";
    }
    std::vector<uint8_t>().swap(syn);

    // successors and blocks
    std::vector<uint32_t> succ(sigma + 1, 0);
    for (uint64_t i = 0; i + 1 < m; ++i) {
        uint32_t& s = succ[pz[i]];
        if (!s || pz[i + 1] < s) s = pz[i + 1];
    }
    std::vector<uint8_t> cont(m + 1, 0);
    L.ccont.assign(sigma / 64 + 1, 0);
    for (uint64_t i = 1; i < m; ++i) {
        cont[i] = succ[pz[i - 1]] == pz[i];
        if (cont[i]) L.ccont[pz[i] >> 6] |= 1ull << (pz[i] & 63);
    }
    // symbol of the block ending just before each block start (and m)
    std::vector<uint64_t> pairs;      // d << 32 | (~j), sorted gives d ascending, j descending
    {
        uint64_t start = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || !cont[i]) {
                uint64_t j = i - 1 - start;
                if (j >= (1ull << 32) - 1) throw std::runtime_error("block too long");
                pairs.push_back(((uint64_t)pz[start] << 32) | (uint32_t)~(uint32_t)j);
                start = i;
            }
    }
    L.nblocks = pairs.size();
    std::vector<uint64_t> syms(pairs);
    std::sort(syms.begin(), syms.end());
    syms.erase(std::unique(syms.begin(), syms.end()), syms.end());
    L.nsym = syms.size();
    {
        std::vector<uint64_t> ds(sigma + 2, 0), jvv(syms.size());
        for (size_t k = 0; k < syms.size(); ++k) { ds[(syms[k] >> 32) + 1]++; jvv[k] = ~(uint32_t)syms[k]; }
        for (uint32_t d = 0; d <= sigma; ++d) ds[d + 1] += ds[d];
        L.set_symbols(ds, jvv);
    }
    L.set_succ(succ);
    std::vector<uint32_t>().swap(succ);
    std::cerr << "blocks: " << L.nblocks << " (" << (double)m / L.nblocks << " phrases each), " << L.nsym << " symbols\n";

    // the parse BWT
    std::vector<int32_t> SA(m + 1);
    {
        std::vector<int32_t> Tz(m + 1);
        for (uint64_t i = 0; i < m; ++i) Tz[i] = (int32_t)pz[i];
        Tz[m] = 0;
        if (libsais_int(Tz.data(), SA.data(), (int32_t)(m + 1), (int32_t)sigma + 1, 0) != 0) throw std::runtime_error("libsais_int failed");
    }
    {
        std::vector<uint32_t> bw(m + 1);
        for (uint64_t r = 0; r <= m; ++r) bw[r] = SA[r] == 0 ? 0 : pz[SA[r] - 1];
        L.pbwt.build(bw, sigma + 1);
    }
    std::cerr << "parse BWT: " << L.pbwt.nruns << " runs, " << L.pbwt.bytes() / 1e6 << " MB (" << secs(t0) << " s)\n";

    // block rows and the block BWT
    {
        std::vector<uint32_t> prevsym(m + 1, 0);   // at block starts p > 0 (and m): symbol of the block ending at p-1
        uint64_t k = 0, start = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || !cont[i]) {
                uint64_t key = pairs[k++];
                prevsym[i] = (uint32_t)(std::lower_bound(syms.begin(), syms.end(), key) - syms.begin()) + 1;
                start = i;
            }
        (void)start;
        std::vector<uint64_t> rows;
        std::vector<uint32_t> bw;
        rows.reserve(L.nblocks + 1); bw.reserve(L.nblocks + 1);
        for (uint64_t r = 0; r <= m; ++r) {
            uint64_t p = (uint64_t)SA[r];
            if (p == m || !cont[p]) { rows.push_back(r); bw.push_back(p == 0 ? 0 : prevsym[p]); }
        }
        if (rows.size() != L.nblocks + 1) throw std::runtime_error("block row count");
        L.BS.build(rows, m + 2);
        L.bbwt.build(bw, (uint32_t)L.nsym + 1);
    }
    std::cerr << "block BWT: " << L.bbwt.nruns << " runs, " << L.bbwt.bytes() / 1e6 << " MB; BS " << L.BS.bytes() / 1e6
              << " MB (" << secs(t0) << " s)\n";
}

int main(int argc, char** argv) {
    std::string text_path, rlbwt_path, rev_path, out_prefix;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (i + 1 < argc && a == "--text") text_path = argv[++i];
        else if (i + 1 < argc && a == "--rlbwt") rlbwt_path = argv[++i];
        else if (i + 1 < argc && a == "--reversed-rlbwt") rev_path = argv[++i];
        else if (i + 1 < argc && a == "--out") out_prefix = argv[++i];
        else { std::cerr << "unknown argument " << a << "\n"; return 1; }
    }
    if (text_path.empty() || rlbwt_path.empty() || rev_path.empty() || out_prefix.empty()) {
        std::cerr << "usage: build_twolevel --text COLLECTION --rlbwt DIR/rlbwt.bin --reversed-rlbwt REVDIR/rlbwt.bin --out PREFIX\n";
        return 1;
    }
    try {
        auto t0 = Clock::now();
        Level L, R;
        std::vector<uint64_t> hf, hr, rhf, rhr;
        std::cerr << "== T\n";
        build_level(text_path, false, rlbwt_path, L, hf, hr, t0);
        std::cerr << "== the reversed collection\n";
        build_level(text_path, true, rev_path, R, rhf, rhr, t0);
        // The phrases of the reversed collection are the reverses of T's (the syncmers
        // are reversal-symmetric), except at the ends of the text, which contain $ or
        // start at position 0 and so never match a pattern's phrase.
        {
            std::vector<std::pair<uint64_t, uint32_t>> rv;
            rv.reserve(R.sigma);
            for (uint32_t y = 1; y <= R.sigma; ++y) rv.push_back({rhf[y], y});
            std::sort(rv.begin(), rv.end());
            std::vector<uint64_t> revid(L.sigma + 1, 0);
            uint64_t missing = 0;
            for (uint32_t x = 1; x <= L.sigma; ++x) {
                auto it = std::lower_bound(rv.begin(), rv.end(), std::make_pair(hr[x], (uint32_t)0));
                if (it != rv.end() && it->first == hr[x]) revid[x] = it->second; else ++missing;
            }
            std::cerr << "shared dictionary: " << L.sigma << " phrases, " << missing << " without a reverse in the reversed collection\n";
            if (missing > 4) throw std::runtime_error("too many phrases without reverses");
            L.dict.build(hf, revid, L.sigma);
        }
        L.save(out_prefix + ".text.lvl");
        R.save(out_prefix + ".rev.lvl");
        for (int d = 0; d < 2; ++d) {
            const Level& x = d ? R : L;
            std::cerr << (d ? "reversed" : "text") << " sizes (MB): B " << x.B.bytes() / 1e6 << ", parse BWT " << x.pbwt.bytes() / 1e6
                      << ", phrase table " << x.dict.bytes() / 1e6 << ", succ " << x.succ_bytes() / 1e6 << ", BS " << x.BS.bytes() / 1e6
                      << ", block BWT " << x.bbwt.bytes() / 1e6 << ", symbols " << x.sym_bytes() / 1e6 << "\n";
        }
        std::cerr << "wrote " << out_prefix << ".text.lvl and .rev.lvl (" << secs(t0) << " s)\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
