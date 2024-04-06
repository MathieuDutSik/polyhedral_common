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
      std::cerr << "GRP_GroupAverage_Matrix [VectorFile] [GroupFile]\n";
      return -1;
    }
    using T = mpq_class;
    //
    std::string VectorFile = argv[1];
    std::string GroupFile = argv[2];
    MyVector<T> eVector = ReadVectorFile<T>(VectorFile);
    std::vector<MyMatrix<T>> LGen = ReadListMatrixFile<T>(GroupFile);
    MyVector<T> eVecAvg = OrbitBarycenter(eVector, LGen);
    std::cerr << "eVecAvg=\n";
    WriteVector(std::cerr, eVecAvg);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_GroupAverage_Vector\n";
    exit(e.eVal);
  }
  runtime(time1);
}
