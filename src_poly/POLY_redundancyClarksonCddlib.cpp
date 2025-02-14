// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "POLY_c_cddlib.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancy [DATAIN] [DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATAOUT : The list of irredundant facets\n";
      return -1;
    }
    //
    using T = mpq_class;
    std::string eFile = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
    //
#ifdef USE_CDDLIB
    std::vector<int> ListIrred = cbased_cdd::RedundancyReductionClarkson(EXT);
    //
    std::ofstream os(argv[2]);
    os << "return [";
    int nbIrred = ListIrred.size();
    for (int i = 0; i < nbIrred; i++) {
      if (i > 0)
        os << ",";
      int eVal = ListIrred[i] + 1;
      os << eVal;
    }
    os << "];\n";
#else
    std::cerr << "You need to compile this program with USE_CDDLIB\n";
    throw TerminalException{1};
#endif
    std::cerr << "Normal termination of POLY_redundancyClarksonCddlib\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancyClarksonCddlib\n";
    exit(e.eVal);
  }
  runtime(time1);
}
