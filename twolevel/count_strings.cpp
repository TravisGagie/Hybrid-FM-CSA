// Copyright (C) 2026 Travis Gagie  (GPL v3 or later; see LICENSE)
//
// Counts the occurrences in T of each line of a file, by backward search on
// T's RLBWT:   count_strings DIR/rlbwt.bin STRINGS > COUNTS
#include <fstream>
#include <iostream>
#include "rlbwt.hpp"
#include "search.hpp"
int main(int argc, char** argv) {
    if (argc != 3) { std::cerr << "usage: count_strings DIR/rlbwt.bin STRINGS\n"; return 1; }
    const rs::RLBWT rl = rs::RLBWT::load(argv[1]);
    std::ifstream in(argv[2]);
    for (std::string s; std::getline(in, s);) {
        rs::SearchResult r = rs::search_plain(rl, s);
        std::cout << (r.found() ? r.hi - r.lo + 1 : 0) << '\n';
    }
    return 0;
}
