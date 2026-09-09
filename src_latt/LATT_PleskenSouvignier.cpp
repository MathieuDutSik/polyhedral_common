// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "LatticePleskenSouvignier.h"
#include "Positivity.h"
// clang-format on

template <typename T, typename Tint>
void process_automorphism(std::string const &FileGram,
                          std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileGram);
  if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
    std::cerr << "LATT_PleskenSouvignier: The input Gram matrix in "
              << FileGram << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListMat{eMat};
  PleskenSouvignierAutomResult<Tint> result =
      PleskenSouvignierLatticeAutomorphism<T, Tint>(ListMat, std::cerr);
  mpz_class order = PleskenSouvignierGroupOrder<mpz_class>(
      result.ListOrbitSize);
  if (OutFormat == "GAP") {
    os << "return rec(ListGen:=[";
    for (size_t i = 0; i < result.ListGen.size(); i++) {
      if (i > 0) {
        os << ",\n";
      }
      WriteMatrixGAP(os, result.ListGen[i]);
    }
    os << "], order:=" << order << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << order << "\n";
    os << result.ListGen.size() << "\n";
    for (auto &eGen : result.ListGen) {
      WriteMatrix(os, eGen);
    }
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
void process_isometry(std::string const &FileGram1,
                      std::string const &FileGram2,
                      std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> eMat1 = ReadMatrixFile<T>(FileGram1);
  MyMatrix<T> eMat2 = ReadMatrixFile<T>(FileGram2);
  for (auto &eMat : {eMat1, eMat2}) {
    if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
      std::cerr << "LATT_PleskenSouvignier: an input Gram matrix is not "
                << "symmetric positive definite\n";
      throw TerminalException{1};
    }
  }
  std::vector<MyMatrix<T>> ListMat1{eMat1};
  std::vector<MyMatrix<T>> ListMat2{eMat2};
  std::vector<MyMatrix<Tint>> ListGenAut2;
  std::optional<MyMatrix<Tint>> opt =
      PleskenSouvignierLatticeIsometry<T, Tint>(ListMat1, ListMat2,
                                                ListGenAut2, std::cerr);
  if (OutFormat == "GAP") {
    if (!opt) {
      os << "return fail;\n";
    } else {
      os << "return ";
      WriteMatrixGAP(os, *opt);
      os << ";\n";
    }
    return;
  }
  if (OutFormat == "CPP") {
    if (!opt) {
      os << "fail\n";
    } else {
      WriteMatrix(os, *opt);
    }
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6 && argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_PleskenSouvignier [arith] aut [GramFile] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_PleskenSouvignier [arith] iso [GramFile1] "
                << "[GramFile2] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: gmp for T=mpq, Tint=mpz\n";
      std::cerr << "aut: the automorphism group generators and order\n";
      std::cerr << "iso: an isometry P with P * M1 * P^T = M2, or fail\n";
      std::cerr << "OutFormat: GAP (default) or CPP\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string mode = argv[2];
    if (mode != "aut" && mode != "iso") {
      std::cerr << "The mode should be aut or iso\n";
      throw TerminalException{1};
    }
    int n_file = (mode == "aut") ? 1 : 2;
    std::string FileGram1 = argv[3];
    std::string FileGram2 = (n_file == 2) ? argv[4] : "unset";
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 3 + n_file + 2) {
      OutFormat = argv[3 + n_file];
      OutFile = argv[3 + n_file + 1];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        if (mode == "aut") {
          return process_automorphism<T, Tint>(FileGram1, OutFormat, os);
        }
        return process_isometry<T, Tint>(FileGram1, FileGram2, OutFormat, os);
      }
      if (arith == "safe") {
        // Machine integers with overflow detection: an overflow throws
        // instead of corrupting, so a caller can fall back to gmp.
        using T = Rational<SafeInt64>;
        using Tint = SafeInt64;
        if (mode == "aut") {
          return process_automorphism<T, Tint>(FileGram1, OutFormat, os);
        }
        return process_isometry<T, Tint>(FileGram1, FileGram2, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_PleskenSouvignier\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_PleskenSouvignier\n";
    exit(e.eVal);
  }
  runtime(time);
}
