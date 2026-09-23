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
#include "InvariantVectorFamily.h"
#include "SignatureSymmetric.h"
// clang-format on

template <typename T, typename Tint>
void ComputeAutomorphism(std::string const &FileListMat,
                         std::string const &OutFormat, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);
  if (ListMat.empty()) {
    std::cerr << "LATT_Automorphism: The input matrix list in " << FileListMat
              << " is empty\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat[0]) ||
      !IsPositiveDefinite(ListMat[0], std::cerr)) {
    std::cerr << "LATT_Automorphism: The first input Gram matrix in "
              << FileListMat << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  if (OutFormat == "GAP_order") {
    // Lightweight path: compute a full-rank invariant short-vector
    // family, take the permutation action of the lattice automorphism
    // group on it, and report just the group order. Matches what
    // ArithmeticAutomorphismGroupMultiple_inner computes internally
    // before lifting permutations back to matrices, but skips the
    // (sometimes expensive) integral matrix recovery step.
    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat[0], std::cerr);
    MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
    int n_row = SHV_T.rows();
    std::vector<T> Vdiag(n_row, T(0));
    std::vector<std::vector<Tidx>> ListPerm =
        GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(SHV_T, ListMat,
                                                           Vdiag, std::cerr);
    std::vector<Telt> ListPermGens;
    for (auto &eList : ListPerm) {
      ListPermGens.push_back(Telt(eList));
    }
    Tgroup grp(ListPermGens, n_row);
    os << "return " << grp.size() << ";\n";
    return;
  }
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, std::cerr);
  if (OutFormat == "GAP") {
    os << "return ";
    WriteListMatrixGAP(os, ListGen);
    os << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    os << ListGen.size() << "\n";
    for (auto &eMat : ListGen) {
      WriteMatrix(os, eMat);
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
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Automorphism [arith] [ListMat] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Automorphism [arith] [ListMat]\n";
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
      std::cerr << "  GAP       : ListGen returned as a GAP-readable list of\n";
      std::cerr << "              integral matrix generators (default)\n";
      std::cerr << "  Oscar     : ListGen returned in the Oscar format\n";
      std::cerr << "  GAP_order : the order |Aut(GramMat)| only, computed "
                << "via\n";
      std::cerr << "              the permutation action on a full-rank\n";
      std::cerr << "              invariant vector family (skips the matrix\n";
      std::cerr << "              lift, much cheaper)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileListMat = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeAutomorphism<T, Tint>(FileListMat, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeAutomorphism<T, Tint>(FileListMat, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeAutomorphism<T, Tint>(FileListMat, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ComputeAutomorphism<T, Tint>(FileListMat, OutFormat, os);
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
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
