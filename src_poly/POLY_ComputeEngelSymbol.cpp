// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_PolytopeFct.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ComputeEngelSymbol [eFileI] [eFileO]\n";
      std::cerr << "\n";
      std::cerr << "eFileI: the input file\n";
      std::cerr << "eFileO: the output file\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    std::ifstream is(argv[1]);
    MyMatrix<mpq_class> EXT = ReadMatrix<mpq_class>(is);
    MyMatrix<mpq_class> FAC = ReadMatrix<mpq_class>(is);
    std::cerr << "After read matrix\n";
    //
    std::string eFileO = argv[2];
    ComputeEngelPolyhedralSubordinationFile(eFileO, EXT, FAC);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
