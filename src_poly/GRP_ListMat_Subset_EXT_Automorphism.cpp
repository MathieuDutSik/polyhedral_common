// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Automorphism [FileI] [OutFormat] [FileO]\n";
      std::cerr << "        or\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Automorphism [FileI]\n";
      std::cerr << "\n";
      std::cerr << "FileI     : The file containing the group\n";
      std::cerr << "OutFormat : The output format that can be GAP or Oscar\n";
      std::cerr << "FileO     : The file containg the output\n";
      return -1;
    }
    //    using T = mpz_class;
    using T = long;
    using Tidx = int32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
#else
    using Tint = mpz_class;
#endif
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
    //
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    std::ifstream is(FileI);
    int nbMat, len;
    is >> nbMat;
    std::vector<MyMatrix<T>> ListMat;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat.push_back(eMat);
    }
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    for (auto &eMat : ListMat) {
      if (!IsSymmetricMatrix(eMat)) {
        std::cerr << "The matrix eMat should be symmetric\n";
        throw TerminalException{1};
      }
      if (eMat.cols() != EXT.cols()) {
        std::cerr << "|eMat|=" << eMat.cols() << " |EXT|=" << EXT.cols()
                  << "\n";
        throw TerminalException{1};
      }
    }
    int n_rows = EXT.rows();
    std::cerr << "n_rows=" << n_rows << "\n";
    is >> len;
    if (len != n_rows) {
      std::cerr << "We have n_rows=" << n_rows << " but len=" << len << "\n";
      throw TerminalException{1};
    }
    std::vector<T> Vdiag(n_rows);
    for (int i = 0; i < n_rows; i++) {
      T val;
      is >> val;
      Vdiag[i] = val;
    }
    //
    const bool use_scheme = true;
    std::vector<std::vector<Tidx>> ListGen =
        GetListGenAutomorphism_ListMat_Vdiag<T, Tidx, use_scheme>(EXT, ListMat,
                                                                  Vdiag);
    //
    std::vector<Telt> LGen;
    for (auto &eList : ListGen)
      LGen.push_back(Telt(eList));
    Tgroup GRP(LGen, n_rows);
    std::cerr << "|GRP|=" << GRP.size() << "\n";
    auto f_print=[&](std::ostream & os) -> void {
      if (OutFormat == "GAP") {
        os << "return " << GRP.GapString() << ";\n";
        return;
      }
      if (OutFormat == "Oscar") {
        WriteGroup(os, GRP);
        return;
      }
      std::cerr << "Failed to find a matching OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      f_print(std::cerr);
    } else {
      if (FileO == "stdout") {
        f_print(std::cout);
      } else {
        std::ofstream os(FileO);
        f_print(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_ListMat_Subset_EXT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
