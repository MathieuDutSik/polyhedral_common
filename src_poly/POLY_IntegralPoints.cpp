// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "POLY_PolytopeInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    using Tint = int;
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IntegralPoints [FAC] [EXT] [OUT]\n";
      std::cerr << "\n";
      std::cerr << "FAC: list of facets (IN)\n";
      std::cerr << "EXT: list of vertices (IN)\n";
      std::cerr << "OUT: list of integral points (OUT)\n";
      return -1;
    }
    //
    //  std::cerr << "Reading input\n";
    //
    MyMatrix<T> FAC = ReadMatrixFile<T>(argv[1]);
    //
    MyMatrix<T> EXT = ReadMatrixFile<T>(argv[2]);
    //
    std::vector<MyVector<Tint>> ListPoint =
        GetListIntegralPoint<T, Tint>(FAC, EXT);
    //
    std::ofstream os(argv[3]);
    os << "return [";
    bool IsFirst = true;
    for (auto &ePt : ListPoint) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      WriteVectorGAP(os, ePt);
    }
    os << "];\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_IntegralPoints\n";
    exit(e.eVal);
  }
  runtime(time1);
}
