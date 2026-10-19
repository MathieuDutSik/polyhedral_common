// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeCanonicalMultiple(std::string const &FileI,
                              std::string const &OutFormat, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMatrix = ReadListMatrixFile<T>(FileI);
  MyMatrix<Tint> B =
      ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMatrix, std::cerr);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> Mat_red = B_T * ListMatrix[0] * B_T.transpose();
  if (OutFormat == "CPP") {
    WriteMatrix(os, Mat_red);
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, B);
    os << ", eG:=";
    WriteMatrixGAP(os, Mat_red);
    os << ");\n";
    return;
  }
  std::cerr << "OutFormat does not match anything\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_CanonicalizeMultiple [arith] [FileMatrices] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_CanonicalizeMultiple [arith] [FileMatrices]\n";
      std::cerr << "\n";
      std::cerr << "FileMatrices (input) : The list of matrices used on "
                   "input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
      std::cerr << "OutFormat values:\n";
      std::cerr << "  CPP : only the reduced matrix is in output\n";
      std::cerr << "  GAP : the basis and the reduced matrix in GAP format\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    //
    // A field is required: the canonicalization of a family that does not
    // span Z^n goes through LinPolytopeIntegral_Canonicalization_Subspaces,
    // which divides. Integral input is unaffected, a rational type prints an
    // integer as an integer.
    auto f = [&](std::ostream &os_out) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeCanonicalMultiple<T>(FileI, OutFormat, os_out);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_CanonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CanonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
