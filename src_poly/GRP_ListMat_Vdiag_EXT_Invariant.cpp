// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_ListMat_Vdiag_EXT_Invariant [INfile] [OUTfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpz_int;
    using Tfield = boost::multiprecision::mpq_rational;
#else
    using T = mpz_class;
    using Tfield = mpq_class;
#endif
    //
    std::cerr << "GRP_ListMat_Vdiag_EXT_Invariant : Reading input\n";
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
    size_t e_hash = GetInvariant_ListMat_Vdiag<T, Tfield>(EXT, ListMat, Vdiag, std::cerr);
    //
    std::ofstream os(argv[2]);
    os << "return " << e_hash << ";\n";
    std::cerr << "Normal termination of GRP_ListMat_Vdiag_EXT_Invariant\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_ListMat_Vdiag_EXT_Invariant\n";
    exit(e.eVal);
  }
  runtime(time1);
}
