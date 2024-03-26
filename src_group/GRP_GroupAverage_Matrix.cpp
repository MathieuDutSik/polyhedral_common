// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "MatrixGroupAverage.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_GroupAverage_Matrix [MatrixFile] [GroupFile]\n";
      return -1;
    }
    using T = mpq_class;
    //
    std::string MatrixFile = argv[1];
    std::string GroupFile = argv[1];
    MyMatrix<T> eMatrix = ReadMatrixFile<T>(MatrixFile);
    std::vector<MyMatrix<T>> LGen = ReadListMatrixFile<T>(GroupFile);
    MyMatrix<T> eMatAvg = OrbitBarycenterSymmetricMatrix(eMatrix, LGen);
    std::cerr << "eMatAvg=\n";
    WriteMatrix(std::cerr, eMatAvg);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_PermMatrTest_PolyStabilizer\n";
    exit(e.eVal);
  }
  runtime(time1);
}
