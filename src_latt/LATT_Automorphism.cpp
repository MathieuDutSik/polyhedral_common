// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "InvariantVectorFamily.h"
#include "SignatureSymmetric.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Automorphism [ListMat] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Automorphism [ListMat]\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  GAP       : ListGen returned as a GAP-readable list of\n";
      std::cerr << "              integral matrix generators (default)\n";
      std::cerr << "  Oscar     : ListGen returned in the Oscar format\n";
      std::cerr << "  GAP_order : the order |Aut(GramMat)| only, computed via\n";
      std::cerr << "              the permutation action on a full-rank\n";
      std::cerr << "              invariant vector family (skips the matrix\n";
      std::cerr << "              lift, much cheaper)\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::string FileListMat = argv[1];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      OutFile = argv[3];
    }
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

    auto prt = [&](std::ostream &os) -> void {
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
            GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(
                SHV_T, ListMat, Vdiag, std::cerr);
        std::vector<Telt> ListPermGens;
        for (auto &eList : ListPerm) {
          ListPermGens.push_back(Telt(eList));
        }
        Tgroup grp(ListPermGens, n_row);
        os << "return " << grp.size() << ";\n";
        return;
      }
      std::vector<MyMatrix<Tint>> ListGen =
          ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat,
                                                                std::cerr);
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
    };
    print_stderr_stdout_file(OutFile, prt);
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
