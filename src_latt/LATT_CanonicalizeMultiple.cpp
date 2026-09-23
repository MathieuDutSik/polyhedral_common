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
// clang-format on

template <typename T, typename Tint>
void ComputeCanonicalMultiple(std::string const &FileI,
                              std::string const &OutFormat, std::ostream &os) {
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
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
      std::cerr << "  multi_boost : the boost multiprecision types\n";
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeCanonicalMultiple<T, Tint>(FileI, OutFormat, os_out);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeCanonicalMultiple<T, Tint>(FileI, OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeCanonicalMultiple<T, Tint>(FileI, OutFormat, os_out);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ComputeCanonicalMultiple<T, Tint>(FileI, OutFormat, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Available possibilities: gmp, gmp_boost, "
                << "multi_boost, flint\n";
#else
      std::cerr << "Available possibilities: gmp, gmp_boost, "
                << "multi_boost (build with ENABLE_FLINT_SUPPORT=1 "
                << "for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_CanonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CanonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
