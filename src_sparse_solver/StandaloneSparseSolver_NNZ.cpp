// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GampMatlab.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr,
              "StandaloneSparseSolver_NNZ [inputMat] [inputVect] [output]\n");
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
    std::ofstream os(argv[3]);
    os << "return [";
    double tol = 1.0e-5;
    int siz = x.size();
    int eNNZ = 0;
    for (int i = 0; i < siz; i++) {
      double eVal = x(i);
      int eVal_i = 0;
      if (fabs(eVal) > tol) {
        eVal_i = 1;
        eNNZ++;
      }
      if (i > 0)
        os << ",";
      os << eVal_i;
    }
    os << "];\n";
    std::cerr << "eNNZ=" << eNNZ << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
