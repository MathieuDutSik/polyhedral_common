// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GenusSpinor.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  The spinor genera of a genus of positive definite lattices of rank >= 3
  (GenusSpinor.h): their number, proper and improper, and for the small
  odd primes l prime to det the number of connected components of the
  l-neighbour graph of the genus, 1 when a single seed suffices for the
  enumeration at l.
 */

using T = mpq_class;
using Tint = mpz_class;

static void PrintUsage() {
  std::cerr << R"(
GENUS_SpinorGenera symbol [GenusSymbol]
GENUS_SpinorGenera gram [FileGram]
GENUS_SpinorGenera list [FileList]

The spinor genera of a genus of positive definite lattices of rank >= 3.

  symbol   the genus as a symbol in the notation of the DefLatDB labels
           (see GenusSymbol.h), e.g. "II{1}^{8}_{II}{25}^{-1}"; a lattice
           of the genus is built for the improper spinor genera
  gram     a Gram matrix of a lattice of the genus
  list     per line: label n g_11 g_12 ... g_nn (the upper triangle of a
           Gram matrix), one output line per genus:
           label proper improper
)";
}

static void PrintReport(GenusSymbol const &s, MyMatrix<T> const &G) {
  std::cout << "genus " << GenusSymbolString(s) << "\n";
  std::cout << "rank " << s.rank << ", determinant "
            << GenusSymbolDeterminant<Tint>(s) << "\n";
  if (s.rank < 3) {
    std::cout << "rank below 3: the spinor genera do not control the "
                 "neighbour graphs\n";
    return;
  }
  SpinorGenusData data = GenusSpinorData<T, Tint>(s, G);
  std::cout << "proper spinor genera " << data.n_proper() << "\n";
  std::cout << "improper spinor genera " << *data.n_improper() << "\n";
  Tint d = GenusSymbolDeterminant<Tint>(s);
  std::cout << "components of the l-neighbour graph, odd l prime to det:";
  int n_shown = 0;
  for (int64_t l = 3; n_shown < 12; l += 2) {
    bool is_prime = true;
    for (int64_t f = 3; f * f <= l; f += 2) {
      if (l % f == 0) {
        is_prime = false;
      }
    }
    if (!is_prime || ResInt(d, Tint(static_cast<long>(l))) == 0) {
      continue;
    }
    std::cout << " " << l << ":" << data.neighbor_components<Tint>(l);
    n_shown++;
  }
  std::cout << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      PrintUsage();
      throw TerminalException{1};
    }
    std::string mode = argv[1];
    std::string arg = argv[2];
    if (mode == "symbol") {
      GenusSymbol s = ParseGenusSymbol<Tint>(arg);
      MyMatrix<T> G = GenusSymbolRepresentative<T, Tint>(s, std::cerr);
      PrintReport(s, G);
    } else if (mode == "gram") {
      MyMatrix<T> G = ReadMatrixFile<T>(arg);
      if (!IsSymmetricMatrix(G) || !IsIntegralMatrix(G)) {
        std::cerr << "GENUS_SpinorGenera: the Gram matrix must be symmetric "
                     "and integral\n";
        throw TerminalException{1};
      }
      PrintReport(GenusSymbolOfGram<T, Tint>(G), G);
    } else if (mode == "list") {
      std::ifstream is(arg);
      std::string line;
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
        if (n < 3) {
          std::cout << label << " rank_below_3\n";
          continue;
        }
        SpinorGenusData data =
            GenusSpinorData<T, Tint>(GenusSymbolOfGram<T, Tint>(G), G);
        std::cout << label << " " << data.n_proper() << " "
                  << *data.n_improper() << "\n";
        std::cout.flush();
      }
    } else {
      PrintUsage();
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of GENUS_SpinorGenera\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_SpinorGenera\n";
    exit(e.eVal);
  }
  runtime(time);
}
