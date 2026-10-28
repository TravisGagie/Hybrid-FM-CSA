// Copyright (C) 2026 Travis Gagie  (GPL v3 or later; see LICENSE)
//
// Builds the parse, block and mapping levels over an RLBWT of the
// reportable-shortcuts code.
//
//   build_twolevel --text COLLECTION --rlbwt DIR/rlbwt.bin --reversed-rlbwt REVDIR/rlbwt.bin --out PREFIX
//   build_twolevel --text COLLECTION --block-stats     (block statistics only)
//   ... --continuation frequent   (continuations: the most frequent successor, not the smallest)
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

#ifndef TL_CHUNK
#define TL_CHUNK (1 << 24)
#endif

using namespace tl;
using Clock = std::chrono::steady_clock;
static bool g_stats_only = false;
static bool g_frequent = false;      // --continuation frequent: the most frequent successor

static double secs(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// The characters of T = S_1 # ... # S_k # $ (or, reversed, of
// T' = rev(S_k) # ... # rev(S_1) # $) for a collection file with one sequence per
// line, as rs::read_text builds them (empty lines skipped, CR dropped), in
// pieces, so that the text is never held in memory.
struct TextStream {
    FILE* f = nullptr;
    bool rev = false, done = false, in_line = false, finished = false;
    long long pos = 0, size = 0;              // reversed: bytes [0, pos) not yet read
    std::vector<uint8_t> raw;
    TextStream(const std::string& path, bool reversed) : rev(reversed) {
        f = fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path);
        fseeko(f, 0, SEEK_END); size = ftello(f); fseeko(f, 0, SEEK_SET);
        int c = fgetc(f);
        if (c == '>') throw std::runtime_error(path + ": FASTA is not supported; use one sequence per line");
        fseeko(f, 0, SEEK_SET);
        pos = rev ? size : 0;
        raw.resize(TL_CHUNK);
    }
    ~TextStream() { if (f) fclose(f); }
    void emit(uint8_t c, std::vector<uint8_t>& out) {
        if (c == '\r') return;
        if (c == '\n') { if (in_line) out.push_back(rs::kSeparator); in_line = false; return; }
        if (c == rs::kSeparator || c == rs::kTerminator) throw std::runtime_error("the collection contains a byte 0x00 or 0x01");
        out.push_back(c); in_line = true;
    }
    // Appends at least one character to out, or returns false at the end.
    bool read(std::vector<uint8_t>& out) {
        if (finished) return false;
        size_t before = out.size();
        while (out.size() == before) {
            size_t got;
            if (!rev) {
                got = fread(raw.data(), 1, raw.size(), f);
                for (size_t i = 0; i < got; ++i) emit(raw[i], out);
            } else {
                long long want = std::min<long long>((long long)raw.size(), pos);
                pos -= want;
                fseeko(f, pos, SEEK_SET);
                got = want ? fread(raw.data(), 1, (size_t)want, f) : 0;
                if ((long long)got != want) throw std::runtime_error("read error");
                // reading backwards, a newline ends the line that follows it
                for (size_t i = got; i-- > 0;) emit(raw[i], out);
            }
            if (got == 0) {
                if (in_line) out.push_back(rs::kSeparator);
                out.push_back(rs::kTerminator);
                finished = true;
                break;
            }
        }
        return true;
    }
};

// Block statistics for other choices of the continuation c(x) of each phrase x
// (--block-stats): the smallest successor (as built), the largest, whichever of
// those two is more frequent after x, and the most frequent successor.  Any
// choice keeps the block BWT in parse order if a block's symbol also records on
// which side of c(last phrase) the next block's first phrase lies; the smallest
// (largest) successor needs no side bit, since that phrase is always larger
// (smaller).  The most frequent successor maximizes the number of continuations.
static void block_stats(const std::vector<uint32_t>& pz, uint32_t sigma) {
    const uint64_t m = pz.size();
    std::vector<uint32_t> cmin(sigma + 1, 0), cmax(sigma + 1, 0), cfrq(sigma + 1, 0);
    std::vector<uint64_t> nmin(sigma + 1, 0), nmax(sigma + 1, 0);
    {
        std::vector<uint64_t> pr;
        pr.reserve(m ? m - 1 : 0);
        for (uint64_t i = 0; i + 1 < m; ++i) pr.push_back((uint64_t)pz[i] << 32 | pz[i + 1]);
        std::sort(pr.begin(), pr.end());
        uint64_t distinct = 0, single = 0;
        for (uint64_t a = 0; a < pr.size();) {
            uint32_t x = (uint32_t)(pr[a] >> 32);
            uint64_t b = a, best = 0, ny = 0;
            while (b < pr.size() && (uint32_t)(pr[b] >> 32) == x) {
                uint64_t c = b;
                while (c < pr.size() && pr[c] == pr[b]) ++c;
                uint32_t y = (uint32_t)pr[b];
                if (!cmin[x]) { cmin[x] = y; nmin[x] = c - b; }
                cmax[x] = y; nmax[x] = c - b;
                if (c - b > best) { best = c - b; cfrq[x] = y; }
                ++ny; b = c;
            }
            distinct += ny; if (ny == 1) ++single;
            a = b;
        }
        std::cerr << "block stats: " << m << " phrases, " << sigma << " distinct, " << distinct
                  << " distinct adjacent pairs; " << single << " phrases with a single successor\n";
    }
    std::vector<uint32_t> cmm(sigma + 1);
    for (uint32_t x = 0; x <= sigma; ++x) cmm[x] = nmax[x] > nmin[x] ? cmax[x] : cmin[x];
    const char* names[4] = {"smallest successor", "largest successor", "smaller or larger, more frequent", "most frequent successor"};
    const std::vector<uint32_t>* cs[4] = {&cmin, &cmax, &cmm, &cfrq};
    for (int r = 0; r < 4; ++r) {
        const auto& c = *cs[r];
        std::vector<uint64_t> keys;       // d, j, side
        uint64_t start = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || c[pz[i - 1]] != pz[i]) {
                uint64_t j = i - 1 - start;
                int side = (i < m && pz[i] > c[pz[i - 1]]) ? 1 : 0;
                if (r != 3) side = 0;  // the side is implied by the choice
                keys.push_back((uint64_t)pz[start] << 33 | (j & 0xffffffffull) << 1 | (uint64_t)side);
                start = i;
            }
        uint64_t nb = keys.size();
        std::sort(keys.begin(), keys.end());
        uint64_t ns = std::unique(keys.begin(), keys.end()) - keys.begin();
        std::cerr << "  " << names[r] << ": " << nb << " blocks (" << (double)m / nb << " phrases each), " << ns << " symbols\n";
    }
}

// Builds the levels of one direction.  hf[id] and hr[id] are the hashes of
// phrase id and of its reverse, for joining the two directions' dictionaries.
// The text is streamed: only the distinct phrases and the parse are kept.
static void build_level(const std::string& text_path, bool reversed, const std::string& rlbwt_path, Level& L,
                        std::vector<uint64_t>& hf, std::vector<uint64_t>& hr, Clock::time_point t0) {
    std::vector<uint32_t> pz;             // the parse, as indices of distinct phrases (later, ids)
    std::vector<uint8_t> dict;            // the distinct phrases, concatenated
    std::vector<uint64_t> doff;           // start of each in dict (one extra entry at the end)
    std::vector<uint64_t> rhash;
    std::unordered_map<uint64_t, uint32_t> seen;
    bool zero_is_syncmer = false;
    uint64_t n = 0;
    {
        TextStream ts(text_path, reversed);
        std::vector<uint8_t> buf;         // X[base .. base + buf.size())
        uint64_t base = 0, open = 0, scanned = 0;   // open: start of the unfinished phrase; scanned: k-mer starts evaluated
        auto add_phrase = [&](uint64_t s, uint64_t e) {
            const uint8_t* p = buf.data() + (s - base);
            uint64_t len = e - s, hv = phrase_hash(p, len);
            auto it = seen.find(hv);
            if (it == seen.end()) {
                uint32_t idx = (uint32_t)rhash.size();
                seen.emplace(hv, idx);
                rhash.push_back(hv);
                doff.push_back(dict.size());
                dict.insert(dict.end(), p, p + len);
                pz.push_back(idx);
            } else {
                uint32_t idx = it->second;
                uint64_t a = doff[idx], b2 = idx + 1 < doff.size() ? doff[idx + 1] : dict.size();
                if (b2 - a != len || std::memcmp(dict.data() + a, p, len)) throw std::runtime_error("phrase hash collision; change phrase_hash");
                pz.push_back(idx);
            }
        };
        bool more = true;
        while (more) {
            more = ts.read(buf);
            const uint64_t end = base + buf.size();
            if (end >= (uint64_t)K && end - K + 1 > scanned) {
                const uint64_t from = scanned;
                for_each_syncmer(buf.data() + (from - base), end - from, [&](uint64_t q) {
                    uint64_t a = from + q;
                    if (a == 0) { zero_is_syncmer = true; return; }
                    add_phrase(open, a + K);
                    open = a;
                });
                scanned = end - K + 1;
            }
            if (!more) { n = end; add_phrase(open, n); break; }
            uint64_t keep = std::min(open, scanned);
            if (keep > base) { buf.erase(buf.begin(), buf.begin() + (keep - base)); base = keep; }
        }
    }
    doff.push_back(dict.size());
    seen = std::unordered_map<uint64_t, uint32_t>();
    const uint64_t m = pz.size();
    if (m + 1 >= (1ull << 31)) throw std::runtime_error("too many phrases for a 32-bit parse suffix array");
    const uint32_t sigma = (uint32_t)rhash.size();
    auto ptr = [&](uint32_t a) { return dict.data() + doff[a]; };
    auto dlen = [&](uint32_t a) { return doff[a + 1] - doff[a]; };
    std::cerr << "text: " << n << " characters; " << m << " phrases, " << sigma << " distinct (" << dict.size() / 1e6
              << " MB) (" << secs(t0) << " s)\n";

    // lexicographic ids 1..sigma (byte order; the dictionary is prefix-free)
    std::vector<uint32_t> order(sigma);
    for (uint32_t a = 0; a < sigma; ++a) order[a] = a;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t c) {
        uint64_t la = dlen(a), lc = dlen(c);
        int r = std::memcmp(ptr(a), ptr(c), std::min(la, lc));
        if (r == 0 && la != lc) throw std::runtime_error("dictionary is not prefix-free");
        return r < 0;
    });
    std::vector<uint32_t> idof(sigma);
    for (uint32_t r = 0; r < sigma; ++r) idof[order[r]] = r + 1;
    for (auto& x : pz) x = idof[x];
    if (g_stats_only) { block_stats(pz, (uint32_t)rhash.size()); return; }
    L.n = n; L.m = m; L.sigma = sigma;
    hf.assign(sigma + 1, 0); hr.assign(sigma + 1, 0);
    {
        std::vector<uint8_t> tmp;
        for (uint32_t a = 0; a < sigma; ++a) {
            tmp.assign(ptr(a), ptr(a) + dlen(a));
            std::reverse(tmp.begin(), tmp.end());
            hf[idof[a]] = rhash[a];
            hr[idof[a]] = phrase_hash(tmp.data(), tmp.size());
        }
    }
    std::cerr << "sorted (" << secs(t0) << " s)\n";

    // distinct syncmers: phrases with equal first K characters are adjacent in id order
    std::vector<uint8_t> syn;
    {
        uint32_t skip = zero_is_syncmer ? 0 : pz[0];
        const uint8_t* prev = nullptr;
        for (uint32_t r = 0; r < sigma; ++r) {
            if (r + 1 == skip) continue;
            const uint8_t* p = ptr(order[r]);
            if (!prev || std::memcmp(prev, p, K)) { syn.insert(syn.end(), p, p + K); prev = p; }
        }
    }
    const size_t nsyn = syn.size() / K;
    std::vector<uint32_t>().swap(order); std::vector<uint32_t>().swap(idof);
    std::vector<uint64_t>().swap(rhash); std::vector<uint8_t>().swap(dict); std::vector<uint64_t>().swap(doff);
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
    if (!g_frequent) {
        for (uint64_t i = 0; i + 1 < m; ++i) {
            uint32_t& s = succ[pz[i]];
            if (!s || pz[i + 1] < s) s = pz[i + 1];
        }
    } else {      // the most frequent successor (ties: the smallest)
        std::vector<uint64_t> pr;
        pr.reserve(m ? m - 1 : 0);
        for (uint64_t i = 0; i + 1 < m; ++i) pr.push_back((uint64_t)pz[i] << 32 | pz[i + 1]);
        std::sort(pr.begin(), pr.end());
        std::vector<uint64_t> best(sigma + 1, 0);
        for (uint64_t a = 0; a < pr.size();) {
            uint64_t b = a;
            while (b < pr.size() && pr[b] == pr[a]) ++b;
            uint32_t x = (uint32_t)(pr[a] >> 32);
            if (b - a > best[x]) { best[x] = b - a; succ[x] = (uint32_t)pr[a]; }
            a = b;
        }
    }
    std::vector<uint8_t> cont(m + 1, 0);
    L.ccont.assign(sigma / 64 + 1, 0);
    for (uint64_t i = 1; i < m; ++i) {
        cont[i] = succ[pz[i - 1]] == pz[i];
        if (cont[i]) L.ccont[pz[i] >> 6] |= 1ull << (pz[i] & 63);
    }
    // symbol of the block ending just before each block start (and m)
    // keys d << 33 | Level::sym_key(j, side): sorted, the symbol order of Level::sym.  With
    // the smallest successor the next block's first phrase is always larger than succ
    // (+); the last block, followed by the sentinel, is given + too, which is safe since
    // its last phrase occurs only once.
    std::vector<uint64_t> pairs;
    {
        uint64_t start = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || !cont[i]) {
                uint64_t j = i - 1 - start;
                if (j >= (1ull << 32) - 1) throw std::runtime_error("block too long");
                bool plus = !g_frequent || (i < m && pz[i] > succ[pz[i - 1]]);
                pairs.push_back(((uint64_t)pz[start] << 33) | Level::sym_key(j, plus));
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
        std::vector<uint8_t> plus(syms.size());
        for (size_t k = 0; k < syms.size(); ++k) {
            ds[(syms[k] >> 33) + 1]++;
            plus[k] = (syms[k] >> 32) & 1;
            uint64_t low = syms[k] & 0xffffffffull;
            jvv[k] = plus[k] ? 0xffffffffull - low : low;
        }
        for (uint32_t d = 0; d <= sigma; ++d) ds[d + 1] += ds[d];
        L.set_symbols(ds, jvv, plus, g_frequent);
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
    std::cerr << "block BWT: " << L.bbwt.nruns << " runs, " << L.bbwt.bytes() / 1e6 << " MB; BS " << L.BS.bytes() / 1e6 << (L.BS.rl ? " (run-length, " + std::to_string(L.BS.r.nruns) + " runs" : std::string(" (Elias-Fano")) << "; the other encoding " << L.BS.alt_bytes / 1e6 << " MB)"
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
        else if (a == "--block-stats") g_stats_only = true;
        else if (i + 1 < argc && a == "--continuation") {
            std::string v = argv[++i];
            if (v == "frequent") g_frequent = true;
            else if (v != "smallest") { std::cerr << "--continuation smallest|frequent\n"; return 1; }
        }
        else { std::cerr << "unknown argument " << a << "\n"; return 1; }
    }
    if (g_stats_only && !text_path.empty()) {
        // statistics only: the parse of T, no RLBWT needed
        try {
            Level L; std::vector<uint64_t> hf, hr;
            build_level(text_path, false, "", L, hf, hr, Clock::now());
        } catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 1; }
        return 0;
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
