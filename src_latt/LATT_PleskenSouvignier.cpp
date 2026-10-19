// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticePleskenSouvignier.h"
#include "Positivity.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process_automorphism(std::string const &FileGram,
                          std::string const &OutFormat, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
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

template <typename T>
void process_isometry(std::string const &FileGram1,
                      std::string const &FileGram2,
                      std::string const &OutFormat, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
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
    std::string mode, FileGram1;
    std::string FileGram2 = "unset";
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(
        dispatch, [&]<typename T>() -> void {
          if (mode == "aut") {
            return process_automorphism<T>(FileGram1, OutFormat, *os_ptr);
          }
          return process_isometry<T>(FileGram1, FileGram2, OutFormat,
                                     *os_ptr);
        });
    if (argc != 4 && argc != 6 && argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_PleskenSouvignier [arith] aut [GramFile] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_PleskenSouvignier [arith] iso [GramFile1] "
                << "[GramFile2] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: " << dispatch.choices() << "\n";
      std::cerr << "aut: the automorphism group generators and order\n";
      std::cerr << "iso: an isometry P with P * M1 * P^T = M2, or fail\n";
      std::cerr << "OutFormat: GAP (default) or CPP\n";
      return -1;
    }
    std::string arith = argv[1];
    mode = argv[2];
    if (mode != "aut" && mode != "iso") {
      std::cerr << "The mode should be aut or iso\n";
      throw TerminalException{1};
    }
    int n_file = (mode == "aut") ? 1 : 2;
    FileGram1 = argv[3];
    if (n_file == 2) {
      FileGram2 = argv[4];
    }
    if (argc == 3 + n_file + 2) {
      OutFormat = argv[3 + n_file];
      OutFile = argv[3 + n_file + 1];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_PleskenSouvignier\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_PleskenSouvignier\n";
    exit(e.eVal);
  }
  runtime(time);
}
