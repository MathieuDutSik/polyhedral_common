// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "HermitianIndefinite.h"
#include <string>
#include <utility>
// clang-format on

/*
  GENUS_HermitianIndefinite: representatives of the classes of the genus of
  an indefinite Hermitian lattice over an imaginary quadratic field, one
  per special genus (HermitianIndefinite.h).

    GENUS_HermitianIndefinite [description] [OutFormat] [OutFile]

  The description is that of GENUS_HermitianGenusOfLattice: the ring, the
  form, and the optional generators of the lattice. OutFormat is GAP (a
  record with the number of classes, the signature, the Gram matrices of
  the trace form Tr h on the Z-bases, the genus symbols and the primes of
  the neighbors used) or CPP (the number of classes, then the Gram
  matrices).
 */

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    if (argc != 2 && argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_HermitianIndefinite [description] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "with OutFormat GAP or CPP, OutFile a file or stdout / "
                   "stderr\n";
      return -1;
    }
    std::string OutFormat = (argc == 4) ? argv[2] : "GAP";
    std::string OutFile = (argc == 4) ? argv[3] : "stdout";
    if (OutFormat != "GAP" && OutFormat != "CPP") {
      std::cerr << "GENUS_HermitianIndefinite: OutFormat is GAP or CPP\n";
      throw TerminalException{1};
    }
    DescriptionTokenReader<Tint> tr{read_description_tokens(argv[1])};
    auto [ring, opt_ord] = read_cm_realified_ring<T, Tint>(tr, std::cerr);
    if (!opt_ord) {
      std::cerr << "GENUS_HermitianIndefinite: the ring must be given as a "
                   "NumberField\n";
      throw TerminalException{1};
    }
    CMOrder<T, Tint> const &ord = *opt_ord;
    if (ord.deg != 2) {
      std::cerr << "GENUS_HermitianIndefinite: only the imaginary quadratic "
                   "fields (K_0 = Q) are treated\n";
      throw TerminalException{1};
    }
    MyMatrix<T> A = read_herm_lattice_basis<T, Tint>(tr, ord, ring);
    MyMatrix<T> G = A * UniversalMatrixConversion<T, Tint>(ring.B) *
                    A.transpose();
    if (DeterminantMat(G) == 0) {
      std::cerr << "GENUS_HermitianIndefinite: the form is degenerate\n";
      throw TerminalException{1};
    }
    // the signature of Tr h is twice that of h
    DiagSymMat<T> dsm = DiagonalizeSymmetricMatrix(G, std::cerr);
    std::pair<int, int> sig{dsm.nbPlus / 2, dsm.nbMinus / 2};
    if (sig.first == 0 || sig.second == 0) {
      std::cerr << "GENUS_HermitianIndefinite: the form is definite "
                   "(GENUS_HermitianEnumerate enumerates the definite "
                   "genera)\n";
      throw TerminalException{1};
    }
    int64_t max_prime = 100000;
    HermIndefiniteResult<T, Tint> res =
        HermIndefiniteRepresentatives<T, Tint>(ord, ring, A, max_prime,
                                               std::cerr);
    if (!res.done) {
      std::cerr << "GENUS_HermitianIndefinite: not done: " << res.reason
                << "\n";
      throw TerminalException{1};
    }
    std::cerr << "GENUS_HermitianIndefinite: signature (" << sig.first << ","
              << sig.second << "), " << res.n_class << " classes\n";
    auto f = [&](std::ostream &os) -> void {
      std::vector<MyMatrix<T>> l_gram;
      for (auto const &B : res.l_basis) {
        l_gram.push_back(B * UniversalMatrixConversion<T, Tint>(ring.B) *
                         B.transpose());
      }
      if (OutFormat == "CPP") {
        os << l_gram.size() << "\n";
        for (auto const &M : l_gram) {
          WriteMatrix(os, M);
        }
        return;
      }
      os << "return rec(n_class:=" << res.n_class << ", signature:=["
         << sig.first << "," << sig.second << "],\nListGram:=[";
      for (size_t i = 0; i < l_gram.size(); i++) {
        os << (i > 0 ? ",\n" : "");
        WriteMatrixGAP(os, l_gram[i]);
      }
      os << "],\nListSymbol:=[";
      for (size_t i = 0; i < res.l_basis.size(); i++) {
        HermLatticeData<T, Tint> L = HermLatticeOfBasis(ord, ring, res.l_basis[i]);
        os << (i > 0 ? ",\n" : "") << "\""
           << HermGlobalSymbolString(HermGenusSymbolOfLattice(L)) << "\"";
      }
      os << "],\nListPrime:=[";
      for (size_t i = 0; i < res.l_prime.size(); i++) {
        os << (i > 0 ? ", " : "") << res.l_prime[i];
      }
      os << "]);\n";
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of GENUS_HermitianIndefinite\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianIndefinite\n";
    exit(e.eVal);
  }
  runtime(time);
}
