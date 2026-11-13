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
  The genus symbols against a list of genera given by a label and one
  lattice, as DefLatDB gives them: per line
      label n g_11 g_12 ... g_1n g_22 ... g_nn
  (the upper triangle of a Gram matrix of the genus). For every line:
  - the symbol of the Gram matrix is the parsed label;
  - GenusSymbolRepresentative builds a lattice of that genus from the label
    alone (it checks its genus itself).
  One line of output per genus: OK or FAIL, the stage, the time.
 */

using T = mpq_class;
using Tint = mpz_class;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Test_GenusSymbol [file]\n";
    std::cerr << "with lines: label n upper-triangle-of-a-Gram-matrix\n";
    return -1;
  }
  std::ifstream is(argv[1]);
  std::string line;
  int n_ok = 0, n_fail = 0;
  while (std::getline(is, line)) {
    if (line.empty()) {
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
      double sec = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count();
      std::cout << "OK " << label << " " << sec << "\n";
      n_ok++;
    } catch (TerminalException const &e) {
      double sec = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count();
      std::cout << "FAIL " << label << " " << stage << " " << sec << "\n";
      n_fail++;
    }
    std::cout.flush();
  }
  std::cout << "n_ok=" << n_ok << " n_fail=" << n_fail << "\n";
  return 0;
}
