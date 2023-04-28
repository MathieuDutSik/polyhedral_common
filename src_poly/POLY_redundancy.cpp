// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "POLY_RedundancyElimination.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancy [DATAIN] [DATAOUT]\n";
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
    std::vector<int> ListIrred = EliminationByRedundance_HitAndRun(EXT);
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
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancy\n";
    exit(e.eVal);
  }
  runtime(time1);
}
