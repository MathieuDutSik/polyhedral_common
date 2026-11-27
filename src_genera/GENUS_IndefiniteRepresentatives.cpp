// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GenusIndefinite.h"
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_IndefiniteRepresentatives: representatives of the classes of the
  genus of an indefinite lattice over Z of rank >= 3, one per improper
  spinor genus (GenusIndefinite.h).
    GENUS_IndefiniteRepresentatives [GramFile] [OutFormat] [OutFile]
  GramFile holds the Gram matrix of a lattice of the genus. OutFormat is
  GAP (a record with the number of classes, the signature, the Gram
  matrices and the primes of the neighbors used) or CPP (the number of
  classes, then the Gram matrices).
 */

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    if (argc != 2 && argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_IndefiniteRepresentatives [GramFile] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "with OutFormat GAP or CPP, OutFile a file or stdout / "
                   "stderr\n";
      return -1;
    }
    std::string GramFile = argv[1];
    std::string OutFormat = (argc == 4) ? argv[2] : "GAP";
    std::string OutFile = (argc == 4) ? argv[3] : "stdout";
    if (OutFormat != "GAP" && OutFormat != "CPP") {
      std::cerr << "GENUS_IndefiniteRepresentatives: OutFormat is GAP or CPP\n";
      throw TerminalException{1};
    }
    MyMatrix<T> G_T = ReadMatrixFile<T>(GramFile);
    int n = G_T.rows();
    if (G_T.cols() != n || !IsSymmetricMatrix(G_T) || !IsIntegralMatrix(G_T)) {
      std::cerr << "GENUS_IndefiniteRepresentatives: the Gram matrix must be "
                   "square, symmetric and integral\n";
      throw TerminalException{1};
    }
    if (DeterminantMat(G_T) == 0) {
      std::cerr << "GENUS_IndefiniteRepresentatives: the form is degenerate\n";
      throw TerminalException{1};
    }
    std::pair<int, int> sig = GenusIndefiniteSignature(G_T);
    if (sig.first == 0 || sig.second == 0) {
      std::cerr << "GENUS_IndefiniteRepresentatives: the form is definite "
                   "(GENUS_Enumerate enumerates the definite genera)\n";
      throw TerminalException{1};
    }
    if (n < 3) {
      std::cerr << "GENUS_IndefiniteRepresentatives: binary forms are not "
                   "treated (their classes are not the spinor genera)\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> G = UniversalMatrixConversion<Tint, T>(G_T);
    GenusIndefiniteResult<T, Tint> res =
        GenusIndefiniteRepresentatives<T, Tint>(G, std::cerr);
    std::cerr << "GENUS_IndefiniteRepresentatives: signature (" << sig.first
              << "," << sig.second << "), " << res.n_class << " classes\n";
    auto f = [&](std::ostream &os) -> void {
      if (OutFormat == "CPP") {
        os << res.l_gram.size() << "\n";
        for (auto const &M : res.l_gram) {
          WriteMatrix(os, M);
        }
        return;
      }
      os << "return rec(n_class:=" << res.n_class << ", signature:=["
         << sig.first << "," << sig.second << "],\nListGram:=[";
      for (size_t i = 0; i < res.l_gram.size(); i++) {
        os << (i > 0 ? ",\n" : "");
        WriteMatrixGAP(os, res.l_gram[i]);
      }
      os << "],\nListPrime:=[";
      for (size_t i = 0; i < res.l_prime.size(); i++) {
        os << (i > 0 ? ", " : "") << res.l_prime[i];
      }
      os << "]);\n";
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of GENUS_IndefiniteRepresentatives\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_IndefiniteRepresentatives\n";
    exit(e.eVal);
  }
  runtime(time);
}
