// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GenusSymbol.h"
#include <string>
// clang-format on

/*
  The genus of a positive definite integral lattice, as the symbol that
  GENUS_Enumerate takes (GenusSymbol in its DATA block): the Conway-Sloane
  symbol in the notation of the DefLatDB labels (GenusSymbol.h).
 */

using T = mpq_class;
using Tint = mpz_class;

static void PrintUsage() {
  std::cerr << R"(
GENUS_GetGenusFromLattice [FileGram] [FileO]
or
GENUS_GetGenusFromLattice [FileGram]

Prints the genus symbol of the Gram matrix in FileGram, e.g.
II{1}^{8}_{II}{25}^{-1}, the input of GENUS_Enumerate. The Gram matrix must
be integral and positive definite.

 ------- FileO --------

  stdout       write to std::cout
  stderr       write to std::cerr
  anything     else is taken as a file name to write to
               Default: stdout
)";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 3) {
      PrintUsage();
      throw TerminalException{1};
    }
    std::string FileGram = argv[1];
    std::string FileO = (argc == 3) ? argv[2] : "stdout";
    MyMatrix<T> G = ReadMatrixFile<T>(FileGram);
    if (!IsSymmetricMatrix(G) || !IsIntegralMatrix(G) ||
        !IsPositiveDefinite(G, std::cerr)) {
      std::cerr << "GENUS_GetGenusFromLattice: the Gram matrix must be "
                   "symmetric, integral and positive definite\n";
      throw TerminalException{1};
    }
    GenusSymbol s = GenusSymbolOfGram<T, Tint>(G);
    std::cerr << "rank " << s.rank << ", determinant "
              << GenusSymbolDeterminant<Tint>(s) << ", "
              << (s.even ? "even" : "odd") << "\n";
    auto f = [&](std::ostream &os) -> void {
      os << GenusSymbolString(s) << "\n";
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of GENUS_GetGenusFromLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_GetGenusFromLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}
