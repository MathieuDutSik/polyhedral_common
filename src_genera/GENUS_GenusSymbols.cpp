// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GenusSymbol.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  The genus symbols (GenusSymbol.h): the genera of positive definite
  lattices of given rank and determinant, and the genera given by labels.

    GENUS_GenusSymbols count n dmin dmax
        a line "n d number" per d, the even genera of rank n and
        determinant d (EnumerateGenusSymbols), as genera() of SageMath
        counts them; countodd: the odd genera
    GENUS_GenusSymbols list n d [build [gram]]
        the symbols of the even genera (listodd: the odd ones); with
        build, a lattice of each genus is constructed
        (GenusSymbolRepresentative, which checks its genus) and its symbol
        compared, printed with gram
    GENUS_GenusSymbols labels File
        File has lines "label n g_11 g_12 ... g_1n g_22 ... g_nn" (a
        label and the upper triangle of a Gram matrix of the genus, as
        DefLatDB gives them); for each, the symbol of the Gram matrix is
        compared with the parsed label and a lattice of the genus is built
        from the label alone: a line "OK label time" or "FAIL label stage
        time", then the counts
 */

using T = mpq_class;
using Tint = mpz_class;

void PrintUsage() {
  std::cerr << "GENUS_GenusSymbols count|countodd n dmin dmax\n";
  std::cerr << "GENUS_GenusSymbols list|listodd n d [build [gram]]\n";
  std::cerr << "GENUS_GenusSymbols labels File\n";
}

void ProcessLabels(std::string const &file) {
  std::ifstream is(file);
  std::string line;
  int n_ok = 0, n_fail = 0;
  while (std::getline(is, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream ls(line);
    std::string label;
    int n;
    ls >> label >> n;
    MyMatrix<T> G(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        long x;
        ls >> x;
        G(i, j) = T(x);
        G(j, i) = T(x);
      }
    }
    auto start = std::chrono::steady_clock::now();
    std::string stage = "parse";
    auto elapsed = [&]() -> double {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                           start)
          .count();
    };
    try {
      GenusSymbol s = ParseGenusSymbol<Tint>(label);
      stage = "symbol";
      GenusSymbol s_G = GenusSymbolOfGram<T, Tint>(G);
      if (s_G != s) {
        std::cerr << "symbol of the Gram matrix: " << GenusSymbolString(s_G)
                  << "\nparsed label:              " << GenusSymbolString(s)
                  << "\n";
        throw TerminalException{1};
      }
      stage = "representative";
      MyMatrix<T> R = GenusSymbolRepresentative<T, Tint>(s, std::cerr);
      std::cout << "OK " << label << " " << elapsed() << "\n";
      n_ok++;
    } catch (TerminalException const &e) {
      std::cout << "FAIL " << label << " " << stage << " " << elapsed()
                << "\n";
      n_fail++;
    }
    std::cout.flush();
  }
  std::cout << "n_ok=" << n_ok << " n_fail=" << n_fail << "\n";
}

int main(int argc, char *argv[]) {
  try {
    if (argc < 3) {
      PrintUsage();
      return -1;
    }
    std::string mode = argv[1];
    if (mode == "labels") {
      ProcessLabels(argv[2]);
      return 0;
    }
    bool even = true;
    if (mode == "countodd" || mode == "listodd") {
      even = false;
      mode = mode.substr(0, mode.size() - 3);
    }
    if (mode == "count" && argc == 5) {
      int n = std::stoi(argv[2]);
      int64_t dmin = std::stoll(argv[3]);
      int64_t dmax = std::stoll(argv[4]);
      for (int64_t d = dmin; d <= dmax; d++) {
        std::vector<GenusSymbol> l =
            EnumerateGenusSymbols<Tint>(n, Tint(static_cast<long>(d)), even);
        std::cout << n << " " << d << " " << l.size() << "\n";
      }
      return 0;
    }
    if (mode == "list" && argc >= 4) {
      int n = std::stoi(argv[2]);
      int64_t d = std::stoll(argv[3]);
      bool build = (argc > 4 && std::string(argv[4]) == "build");
      bool gram = (argc > 5 && std::string(argv[5]) == "gram");
      std::vector<GenusSymbol> l =
          EnumerateGenusSymbols<Tint>(n, Tint(static_cast<long>(d)), even);
      for (auto const &g : l) {
        std::cout << GenusSymbolString(g);
        if (build) {
          MyMatrix<T> G = GenusSymbolRepresentative<T, Tint>(g, std::cerr);
          GenusSymbol g2 = GenusSymbolOfGram<T, Tint>(G);
          std::cout << (g2 == g ? " built" : " MISMATCH");
          if (gram) {
            std::cout << "\n" << G.rows() << " " << G.cols() << "\n";
            for (int i = 0; i < G.rows(); i++) {
              for (int j = 0; j < G.cols(); j++) {
                std::cout << (j > 0 ? " " : "") << G(i, j);
              }
              std::cout << "\n";
            }
          }
        }
        std::cout << "\n";
      }
      return 0;
    }
    PrintUsage();
    return -1;
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_GenusSymbols\n";
    exit(e.eVal);
  }
}
