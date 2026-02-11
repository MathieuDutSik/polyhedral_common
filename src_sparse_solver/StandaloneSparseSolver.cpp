// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GampMatlab.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr,
              "StandaloneSparseSolver [inputMat] [inputVect] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MySparseMatrix<double> SpMat = ReadSparseMatrix<double>(INmat);
    std::cerr << "Read SpMat rows/cols = " << SpMat.rows() << " / "
              << SpMat.cols() << "\n";
    // reading the vector
    std::ifstream INvect(argv[2]);
    MyVector<double> eVect = ReadVector<double>(INvect);
    // creating the RecSparse operators.
    MyVector<double> x = AMP_SolutionSparseSystem(SpMat, eVect, std::cerr);
    //
    std::ofstream OUTvect(argv[3]);
    WriteVector(OUTvect, x);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
