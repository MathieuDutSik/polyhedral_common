// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Canonic [EXTIN] [OutCan]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      std::cerr << "OutCan : The canonicalization file\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tidx = int16_t;
    std::string FileExt = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileExt);
    size_t nbCol = EXT.cols();
    size_t nbRow = EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    size_t threshold = THRESHOLD_USE_SUBSET_SCHEME_CANONIC;
    std::vector<Tidx> CanonicOrd =
      LinPolytope_CanonicOrdering<T, Tidx>(EXT, threshold, std::cerr);
    //
    if (argc == 3) {
      std::ofstream os(argv[2]);
      for (size_t iRow = 0; iRow < nbRow; iRow++)
        os << " " << CanonicOrd[iRow];
      os << "\n";
    }
    std::cerr << "Normal termination of GRP_LinPolytope_Canonic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
