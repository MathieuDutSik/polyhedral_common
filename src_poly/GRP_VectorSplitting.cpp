// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeDefinitions.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_VectorSplitting [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "FileIn    : The input file containing n, the modulo, and "
                   "the list of matrices\n";
      std::cerr << "FileOut   : The list of computed orbits\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    std::ifstream is(argv[1]);
    //
    int n;
    is >> n;
    //
    int pPrime;
    is >> pPrime;
    //
    int nbMat;
    is >> nbMat;
    std::vector<MyMatrix<int>> ListMat(nbMat);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<int> eMat(n, n);
      for (int iRow = 0; iRow < n; iRow++)
        for (int iCol = 0; iCol < n; iCol++) {
          int eVal;
          is >> eVal;
          int res = eVal % pPrime;
          if (res < 0)
            res += pPrime;
          eMat(iRow, iCol) = res;
        }
      ListMat[iMat] = eMat;
    }
    //
    std::ofstream os(argv[2]);
    os << "return [";
    bool IsFirst = true;
    std::function<void(std::vector<int> const &, int const &)> FCT =
        [&](std::vector<int> const &V, int const &OrbSize) -> void {
      if (IsFirst == false)
        os << ",\n";
      IsFirst = false;
      //
      os << "rec(OrbSize := " << OrbSize << ", eVect := [";
      for (int i = 0; i < n; i++) {
        if (i > 0)
          os << ",";
        os << V[i];
      }
      os << "])";
    };
    EnumerateOrbitPrimitiveVector<int>(ListMat, pPrime, n, FCT);
    os << "];\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_VectorSplitting\n";
    exit(e.eVal);
  }
  runtime(time1);
}
