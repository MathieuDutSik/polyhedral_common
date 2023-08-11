// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinearDetermineByInequalities [DATAFAC]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT: The input data of the polytope vertices\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    using T = mpq_class;
    std::string eFile = argv[1];
    MyMatrix<T> FAC = ReadMatrixFile<T>(eFile);
    //
    MyMatrix<T> LinSpace = LinearDeterminedByInequalities(FAC);
    std::cerr << "LinSpace=\n";
    WriteMatrix(std::cerr, LinSpace);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_LinearDetermineByInequalities\n";
    exit(e.eVal);
  }
  runtime(time1);
}
