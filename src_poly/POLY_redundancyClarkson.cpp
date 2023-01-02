// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_cddlib.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancyClarskson [DATAIN] [DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      std::cerr << "DATAOUT : The list of irredundant facets\n";
      return -1;
    }
    //
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    //
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(EXT);
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
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
