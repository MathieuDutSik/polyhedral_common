// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticeEichlerDecomposition.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

static void PrintUsage() {
  std::cerr << R"(
LATT_EichlerDecomposition [arith] [FileI] [OutFormat] [FileO]
or
LATT_EichlerDecomposition [arith] [FileI]

Decomposes a positive definite Gram matrix into its indecomposable orthogonal
summands (Eichler's theorem: the decomposition is unique). See
LatticeEichlerDecomposition.h.

 ------- arith --------

  gmp          mpq_class over mpz_class
)";
#ifdef ENABLE_BOOST_TYPES
  std::cerr << R"(  gmp_boost    boost mpq_rational over mpz_int
  multi_boost  boost cpp_rational over cpp_int
)";
#endif
#ifdef ENABLE_FLINT_SUPPORT
  std::cerr << R"(  flint        fmpq_class over fmpz_class
)";
#endif
  std::cerr << R"(
 ------- OutFormat --------

  GAP          the list of the summands, each with its basis (rows, in the
               coordinates of the input) and its Gram matrix, GAP readable
  CPP          the number of summands, then for each its basis and its Gram
               matrix
               Default: GAP

 ------- FileO --------

  stdout       write to std::cout
  stderr       write to std::cerr
  anything     else is taken as a file name to write to
               Default: stderr
)";
}

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileI);
  if (!IsSymmetricMatrix(GramMat) ||
      !IsPositiveDefinite(GramMat, std::cerr)) {
    std::cerr << "LATT_EichlerDecomposition: the Gram matrix must be "
                 "symmetric positive definite\n";
    throw TerminalException{1};
  }
  std::vector<EichlerComponent<T, Tint>> ListComp =
      EichlerDecomposition<T, Tint>(GramMat, std::cerr);
  std::cerr << "dimension " << GramMat.rows() << ", " << ListComp.size()
            << " indecomposable summands of ranks";
  for (auto const &comp : ListComp) {
    std::cerr << " " << comp.Gram.rows();
  }
  std::cerr << "\n";
  if (OutFormat == "GAP") {
    os << "return [";
    for (size_t i = 0; i < ListComp.size(); i++) {
      if (i > 0) {
        os << ",\n";
      }
      os << "rec(Basis:=";
      WriteMatrixGAP(os, ListComp[i].Basis);
      os << ", Gram:=";
      WriteMatrixGAP(os, ListComp[i].Gram);
      os << ")";
    }
    os << "];\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << ListComp.size() << "\n";
    for (auto const &comp : ListComp) {
      WriteMatrix(os, comp.Basis);
      WriteMatrix(os, comp.Gram);
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat=" << OutFormat
            << ". Allowed: GAP, CPP\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      PrintUsage();
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            process<T>(FileI, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_EichlerDecomposition\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_EichlerDecomposition\n";
    exit(e.eVal);
  }
  runtime(time);
}
