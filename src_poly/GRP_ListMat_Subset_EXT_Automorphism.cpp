// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheoryRational.h"
#endif
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Automorphism [INfile] [OUTfile]\n";
      std::cerr << "        or\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Automorphism [INfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
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
    std::ifstream is(argv[1]);
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
    if (argc == 3) {
      std::ofstream os(argv[2]);
      os << "return " << GRP.GapString() << ";\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_ListMat_Subset_EXT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
