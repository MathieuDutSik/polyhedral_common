// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Canonic [EXTIN] [OutCan]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      std::cerr << "OutCan : The canonicalization file\n";
      return -1;
    }
    //
    using T = mpq_class;
    const bool use_scheme = true;
    using Tidx = int16_t;
    std::string FileExt = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileExt);
    size_t nbCol = EXT.cols();
    size_t nbRow = EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    std::vector<Tidx> CanonicOrd =
        LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT);
    //
    if (argc == 3) {
      std::ofstream os(argv[2]);
      for (size_t iRow = 0; iRow < nbRow; iRow++)
        os << " " << CanonicOrd[iRow];
      os << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
