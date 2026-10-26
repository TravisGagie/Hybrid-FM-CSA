// Copyright (C) 2026 Travis Gagie  (GPL v3 or later; see LICENSE)
//
// Finds the MEMs of patterns by forward-backward search and Boyer-Moore-Li,
// as rs::mems_forward_backward and rs::mems_bml do, with each of
//   plain       plain backward search over the two RLBWTs;
//   NAME        the shortcuts NAME (DIR/NAME.shortcuts and REVDIR/NAME.shortcuts);
//   2L          two-level search over the parse, without block jumps;
//   2LJ         two-level search with block jumps,
// in rounds (plain, then each other configuration in turn, K times), checks
// that all find the same MEMs (with the same BWT intervals) as plain search,
// and reports the median time of each.  The two-level times include parsing
// the pattern and its reverse.
//
//   query_twolevel --index DIR --reversed-index DIR --level FILE --reversed-level FILE --patterns FILE
//                  [--shortcuts NAME,...] [--algorithm fb,bml] [--min-length 100,1000] [--repeat 5]
//                  [--configs plain,2L,2LJ] [--csv FILE]
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "mems.hpp"
#include "shortcuts.hpp"
#include "twolevel.hpp"

using Clock = std::chrono::steady_clock;

namespace {

struct Cnt {
    int64_t searches = 0, char_steps = 0, parse_steps = 0, block_steps = 0, shortcuts = 0;
    double prep_ns = 0;
};

// One configuration: how a search of P[begin .. end) on T, or of the reversed
// pattern on T', is done.
struct Engine {
    std::string name;
    const rs::RLBWT* text = nullptr;
    const rs::RLBWT* rev = nullptr;
    const rs::Shortcuts* sc = nullptr;
    const rs::Shortcuts* sc_rev = nullptr;
    const tl::Level* L = nullptr;
    const tl::Level* Lrev = nullptr;
    bool twolevel = false, jumps = false;
    tl::PatInfo pi, qi;

    void prepare(const std::string& P, const std::string& Q) {
        if (twolevel) tl::parse_pattern(*L, *Lrev, P, pi, qi);
    }
    rs::SearchResult search(bool reversed, const std::string& X, int64_t begin, int64_t end, Cnt& c) {
        if (twolevel) {
            tl::Counters t;
            rs::SearchResult r = tl::search2(reversed ? *rev : *text, reversed ? *Lrev : *L, reversed ? qi : pi, X, begin, end,
                                             jumps, &t);
            c.searches += t.searches; c.char_steps += t.char_steps; c.parse_steps += t.parse_steps; c.block_steps += t.block_steps;
            return r;
        }
        rs::SearchResult r = rs::search_range(reversed ? *rev : *text, reversed ? sc_rev : sc, X, begin, end);
        ++c.searches; c.char_steps += r.backward_steps; c.shortcuts += r.shortcuts;
        return r;
    }
};

// The MEM finder of mems.cpp, over an Engine.
struct Finder {
    Engine& E;
    const std::string& P;
    std::string Q;
    int64_t m;
    Cnt& c;
    Finder(Engine& e, const std::string& p, Cnt& cn) : E(e), P(p), Q(p.rbegin(), p.rend()), m((int64_t)p.size()), c(cn) {
        auto t0 = Clock::now();
        E.prepare(P, Q);
        c.prep_ns += std::chrono::duration<double, std::nano>(Clock::now() - t0).count();
    }
    int64_t forward(int64_t x) { return x + E.search(true, Q, 0, m - x, c).matched; }
    rs::SearchResult backward(int64_t b, int64_t e) { return E.search(false, P, b, e, c); }
    int64_t report_and_advance(int64_t x, int64_t e, int64_t L, std::vector<rs::Mem>& out) {
        rs::SearchResult r = backward(0, e);
        const int64_t s = e - r.matched;
        if (e - s >= L) out.push_back({s, e, r.lo, r.hi});
        if (e == m) return m;
        const int64_t s1 = e + 1 - backward(0, e + 1).matched;
        return std::max(s1, x + 1);
    }
    std::vector<rs::Mem> fb(int64_t L) {
        std::vector<rs::Mem> out;
        for (int64_t x = 0; x < m;) {
            const int64_t e = forward(x);
            if (e == x) { ++x; continue; }
            x = report_and_advance(x, e, L, out);
        }
        return out;
    }
    std::vector<rs::Mem> bml(int64_t L) {
        std::vector<rs::Mem> out;
        for (int64_t x = 0; x + L <= m;) {
            const int64_t we = x + L;
            const rs::SearchResult r = backward(x, we);
            if (r.matched < L) { x = we - r.matched; continue; }
            x = report_and_advance(x, forward(x), L, out);
        }
        return out;
    }
};

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    size_t k = v.size() / 2;
    return v.size() % 2 ? v[k] : (v[k - 1] + v[k]) / 2;
}

std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> v;
    std::stringstream ss(s);
    for (std::string x; std::getline(ss, x, ',');) if (!x.empty()) v.push_back(x);
    return v;
}

}  // namespace

int main(int argc, char** argv) {
    std::string idx, ridx, lev, rlev, pat, csvp;
    std::vector<std::string> scn, algs = {"fb", "bml"}, configs = {"plain", "2L", "2LJ"};
    std::vector<int64_t> Ls = {100, 1000};
    int64_t repeat = 5;
    bool configs_given = false;
    for (int i = 1; i + 1 < argc; i += 2) {
        std::string f = argv[i], v = argv[i + 1];
        if (f == "--index") idx = v;
        else if (f == "--reversed-index") ridx = v;
        else if (f == "--level") lev = v;
        else if (f == "--reversed-level") rlev = v;
        else if (f == "--patterns") pat = v;
        else if (f == "--shortcuts") scn = split(v);
        else if (f == "--algorithm") algs = split(v);
        else if (f == "--min-length") { Ls.clear(); for (auto& x : split(v)) Ls.push_back(std::stoll(x)); }
        else if (f == "--repeat") repeat = std::max<int64_t>(1, std::stoll(v));
        else if (f == "--configs") { configs = split(v); configs_given = true; }
        else if (f == "--csv") csvp = v;
        else { std::cerr << "unknown flag " << f << "\n"; return 1; }
    }
    if (argc % 2 == 0 || idx.empty() || ridx.empty() || pat.empty()) {
        std::cerr << "usage: query_twolevel --index DIR --reversed-index DIR --level FILE --reversed-level FILE --patterns FILE\n"
                     "       [--shortcuts NAME,...] [--algorithm fb,bml] [--min-length 100,1000] [--repeat 5]\n"
                     "       [--configs plain,2L,2LJ] [--csv FILE]\n";
        return 1;
    }
    try {
        auto t0 = Clock::now();
        const rs::RLBWT text = rs::RLBWT::load(idx + "/rlbwt.bin");
        const rs::RLBWT rev = rs::RLBWT::load(ridx + "/rlbwt.bin");
        std::vector<std::unique_ptr<rs::Shortcuts>> sc, scr;
        // NAME (in --index and --reversed-index) or NAME@DIR@REVDIR; the configuration is called NAME or NAME@DIR
        for (auto& spec : scn) {
            std::string nm = spec, d = idx, rd = ridx;
            size_t a = spec.find('@');
            if (a != std::string::npos) {
                size_t b2 = spec.find('@', a + 1);
                if (b2 == std::string::npos) throw std::runtime_error("--shortcuts: expected NAME@DIR@REVDIR in " + spec);
                nm = spec.substr(0, a); d = spec.substr(a + 1, b2 - a - 1); rd = spec.substr(b2 + 1);
                spec = nm + "@" + d;
            }
            sc.emplace_back(new rs::Shortcuts(rs::Shortcuts::load(d + "/" + nm + ".shortcuts")));
            scr.emplace_back(new rs::Shortcuts(rs::Shortcuts::load(rd + "/" + nm + ".shortcuts")));
            std::cerr << "shortcuts " << spec << ": " << (std::filesystem::file_size(d + "/" + nm + ".shortcuts") +
                                                           std::filesystem::file_size(rd + "/" + nm + ".shortcuts")) / 1048576.0
                      << " MiB on disk (both directions)\n";
        }
        tl::Level L, Lr;
        bool need2 = false;
        for (auto& c : configs) if (c == "2L" || c == "2LJ") need2 = true;
        if (need2) {
            if (lev.empty() || rlev.empty()) throw std::runtime_error("--level and --reversed-level are needed for 2L and 2LJ");
            L.load(lev); Lr.load(rlev);
            auto mb = [](const tl::Level& x, bool j) {
                double s = x.B.bytes() + x.pbwt.bytes() + x.table_bytes();
                if (j) s += x.succ_bytes() + x.BS.bytes() + x.bbwt.bytes() + x.sym_bytes();
                return s / 1048576.0;
            };
            for (int d = 0; d < 2; ++d) {
                const tl::Level& x = d ? Lr : L;
                std::cerr << (d ? "reversed" : "text") << " level: m = " << x.m << ", sigma = " << x.sigma << ", blocks = " << x.nblocks
                          << "; MiB: B " << x.B.bytes() / 1048576.0 << ", parse BWT " << x.pbwt.bytes() / 1048576.0
                          << ", phrase table " << x.table_bytes() / 1048576.0 << ", succ " << x.succ_bytes() / 1048576.0
                          << ", BS " << x.BS.bytes() / 1048576.0 << ", block BWT " << x.bbwt.bytes() / 1048576.0 << ", symbols "
                          << x.sym_bytes() / 1048576.0 << "\n";
            }
            std::cerr << "extra space over the two RLBWTs: 2L " << mb(L, false) + mb(Lr, false) << " MiB, 2LJ "
                      << mb(L, true) + mb(Lr, true) << " MiB\n";
        }
        std::cerr << "loaded in " << std::chrono::duration<double>(Clock::now() - t0).count() << " s\n";

        std::vector<Engine> E;
        if (!configs_given && !scn.empty()) configs.insert(configs.begin() + 1, scn.begin(), scn.end());
        for (auto& c : configs) {
            Engine e; e.name = c; e.text = &text; e.rev = &rev;
            if (c == "plain") {}
            else if (c == "2L" || c == "2LJ") { e.twolevel = true; e.jumps = c == "2LJ"; e.L = &L; e.Lrev = &Lr; }
            else {
                auto it = std::find(scn.begin(), scn.end(), c);
                if (it == scn.end()) throw std::runtime_error("unknown configuration " + c);
                e.sc = sc[it - scn.begin()].get(); e.sc_rev = scr[it - scn.begin()].get();
            }
            E.push_back(std::move(e));
        }
        if (E.empty() || E[0].name != "plain") throw std::runtime_error("the first configuration must be plain");

        std::vector<std::string> pats;
        {
            std::ifstream in(pat);
            if (!in) throw std::runtime_error("cannot open " + pat);
            for (std::string P; std::getline(in, P);) {
                if (!P.empty() && P.back() == '\r') P.pop_back();
                if (!P.empty()) pats.push_back(P);
            }
        }
        std::ofstream csv;
        if (!csvp.empty()) {
            csv.open(csvp);
            csv << "config,algorithm,min_length,pattern,length,mems,mem_chars,median_ns,min_ns,max_ns,searches,char_steps,parse_steps,"
                   "block_steps,shortcuts\n";
        }
        int64_t mism = 0;
        for (auto& alg : algs) {
            for (int64_t Lm : Ls) {
                std::vector<double> total(E.size(), 0);
                std::vector<Cnt> ctot(E.size());
                for (size_t q = 0; q < pats.size(); ++q) {
                    const std::string& P = pats[q];
                    std::vector<std::vector<double>> times(E.size());
                    std::vector<rs::Mem> ref;
                    std::vector<Cnt> cs(E.size());
                    for (int64_t k = 0; k < repeat; ++k) {
                        for (size_t c = 0; c < E.size(); ++c) {
                            Cnt cn;
                            auto t1 = Clock::now();
                            std::vector<rs::Mem> res;
                            {
                                Finder f(E[c], P, cn);
                                res = alg == "fb" ? f.fb(Lm) : f.bml(Lm);
                            }
                            times[c].push_back(std::chrono::duration<double, std::nano>(Clock::now() - t1).count());
                            cs[c] = cn;
                            if (c == 0) ref = std::move(res);
                            else if (res != ref) {
                                ++mism;
                                std::cerr << "MISMATCH: " << E[c].name << ", " << alg << ", L = " << Lm << ", pattern " << q << " ("
                                          << res.size() << " MEMs vs " << ref.size() << ")\n";
                            }
                        }
                    }
                    int64_t mc = 0;
                    for (auto& x : ref) mc += x.end - x.start;
                    for (size_t c = 0; c < E.size(); ++c) {
                        double md = median(times[c]);
                        total[c] += md;
                        ctot[c].searches += cs[c].searches; ctot[c].char_steps += cs[c].char_steps;
                        ctot[c].parse_steps += cs[c].parse_steps; ctot[c].block_steps += cs[c].block_steps;
                        ctot[c].shortcuts += cs[c].shortcuts; ctot[c].prep_ns += cs[c].prep_ns;
                        if (csv.is_open())
                            csv << E[c].name << ',' << alg << ',' << Lm << ',' << q << ',' << P.size() << ',' << ref.size() << ',' << mc
                                << ',' << md << ',' << *std::min_element(times[c].begin(), times[c].end()) << ','
                                << *std::max_element(times[c].begin(), times[c].end()) << ',' << cs[c].searches << ','
                                << cs[c].char_steps << ',' << cs[c].parse_steps << ',' << cs[c].block_steps << ',' << cs[c].shortcuts
                                << '\n';
                    }
                    if (csv.is_open()) csv.flush();
                }
                for (size_t c = 0; c < E.size(); ++c)
                    std::cout << alg << ", L = " << Lm << ", " << E[c].name << ": " << total[c] / 1e9 << " s (speedup "
                              << total[0] / total[c] << "); steps: char " << ctot[c].char_steps << ", parse " << ctot[c].parse_steps
                              << ", block " << ctot[c].block_steps << ", shortcuts " << ctot[c].shortcuts
                              << "; parsing the patterns " << ctot[c].prep_ns / 1e9 << " s" << std::endl;
            }
        }
        if (mism) { std::cerr << mism << " mismatches\n"; return 2; }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
