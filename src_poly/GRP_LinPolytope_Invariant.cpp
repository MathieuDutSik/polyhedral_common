// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Invariant [EXTIN] [OutCan]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      std::cerr << "OutCan : The canonicalization file\n";
      return -1;
    }
    //
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tfield = T;
#else
    using T = mpq_class;
    using Tfield = T;
#endif
    std::string FileExt = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileExt);
    MyMatrix<T> EXTred = RowReduction(EXT);
    int n_rows = EXT.rows();
    //
    MyMatrix<T> Qinv = GetQmatrix(EXTred);
    std::vector<MyMatrix<T>> ListMat = { Qinv };
    std::vector<T> Vdiag(n_rows, 0);
    //
    size_t e_hash = GetInvariant_ListMat_Vdiag<T,Tfield>(EXTred, ListMat, Vdiag);
    //
    std::string FileOut = argv[2];
    std::ofstream os(FileOut);
    os << e_hash << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time);
}
