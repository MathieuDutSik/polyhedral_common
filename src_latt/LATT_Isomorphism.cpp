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
#include "SignatureSymmetric.h"
// clang-format on

template <typename T, typename Tint>
void ComputeIsomorphism(std::string const &FileListMat1,
                        std::string const &FileListMat2,
                        std::string const &OutFormat, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMat1 = ReadListMatrixFile<T>(FileListMat1);
  std::vector<MyMatrix<T>> ListMat2 = ReadListMatrixFile<T>(FileListMat2);
  if (ListMat1.empty()) {
    std::cerr << "LATT_Isomorphism: The input matrix list in " << FileListMat1
              << " is empty\n";
    throw TerminalException{1};
  }
  if (ListMat2.empty()) {
    std::cerr << "LATT_Isomorphism: The input matrix list in " << FileListMat2
              << " is empty\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat1[0]) ||
      !IsPositiveDefinite(ListMat1[0], std::cerr)) {
    std::cerr << "LATT_Isomorphism: The first input Gram matrix in "
              << FileListMat1 << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat2[0]) ||
      !IsPositiveDefinite(ListMat2[0], std::cerr)) {
    std::cerr << "LATT_Isomorphism: The first input Gram matrix in "
              << FileListMat2 << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  std::optional<MyMatrix<Tint>> equiv =
      ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                     std::cerr);
  if (OutFormat == "GAP") {
    if (equiv) {
      os << "return ";
      WriteMatrixGAP(os, *equiv);
      os << ";\n";
    } else {
      os << "return false;\n";
    }
    return;
  }
  if (OutFormat == "Oscar") {
    if (equiv) {
      WriteMatrix(os, *equiv);
    } else {
      os << "0 0\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Isomorphism [arith] [ListMat1] [ListMat2] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [arith] [ListMat1] [ListMat2]\n";
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
      std::cerr << "  GAP   : the equivalence matrix or false, GAP readable\n";
      std::cerr << "  Oscar : the equivalence matrix in the Oscar format\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileListMat1 = argv[2];
    std::string FileListMat2 = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeIsomorphism<T, Tint>(FileListMat1, FileListMat2,
                                           OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeIsomorphism<T, Tint>(FileListMat1, FileListMat2,
                                           OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeIsomorphism<T, Tint>(FileListMat1, FileListMat2,
                                           OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ComputeIsomorphism<T, Tint>(FileListMat1, FileListMat2,
                                           OutFormat, os);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of LATT_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
