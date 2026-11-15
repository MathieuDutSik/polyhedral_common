// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GenusSymbol.h"
#include <string>
#include <vector>
// clang-format on

/*
  The genera of even positive definite lattices of rank n and determinant
  d (EnumerateEvenGenusSymbols), for comparison with genera() of SageMath:
    count n dmin dmax : a line "n d number" per d (countodd, listodd:
                        the odd lattices)
    list n d [build [gram]] : the symbols; with build, a lattice of each
                        genus is constructed (GenusSymbolRepresentative
                        checks its genus) and its symbol compared, printed
                        with gram
 */

using T = mpq_class;
using Tint = mpz_class;

int main(int argc, char *argv[]) {
  try {
    if (argc < 4) {
      std::cerr << "Test_GenusSymbolEnumeration count n dmin dmax\n";
      std::cerr << "Test_GenusSymbolEnumeration list n d [build]\n";
      return -1;
    }
    std::string mode = argv[1];
    // countodd / listodd: the odd lattices
    bool even = true;
    if (mode == "countodd" || mode == "listodd") {
      even = false;
      mode = mode.substr(0, mode.size() - 3);
    }
    int n = std::stoi(argv[2]);
    if (mode == "count") {
      int64_t dmin = std::stoll(argv[3]);
      int64_t dmax = std::stoll(argv[4]);
      for (int64_t d = dmin; d <= dmax; d++) {
        std::vector<GenusSymbol> l =
            EnumerateGenusSymbols<Tint>(n, Tint(static_cast<long>(d)), even);
        std::cout << n << " " << d << " " << l.size() << "\n";
      }
      return 0;
    }
    int64_t d = std::stoll(argv[3]);
    bool build = (argc > 4 && std::string(argv[4]) == "build");
    std::vector<GenusSymbol> l =
        EnumerateGenusSymbols<Tint>(n, Tint(static_cast<long>(d)), even);
    for (auto const &g : l) {
      std::cout << GenusSymbolString(g);
      if (build) {
        MyMatrix<T> G = GenusSymbolRepresentative<T, Tint>(g, std::cerr);
        GenusSymbol g2 = GenusSymbolOfGram<T, Tint>(G);
        std::cout << (g2 == g ? " built" : " MISMATCH");
        if (argc > 5 && std::string(argv[5]) == "gram") {
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
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_GenusSymbolEnumeration\n";
    exit(e.eVal);
  }
}
