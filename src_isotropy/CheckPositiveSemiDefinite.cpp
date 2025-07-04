// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Positivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "CheckPositiveSemiDefinite [inputMat]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    // computing the kernel
    DiagSymMat<T> eDiag = DiagonalizeSymmetricMatrix(TheMat, std::cerr);
    if (eDiag.nbMinus > 0) {
      std::cout << "The matrix is NOT positive semidefinite\n";
    } else {
      if (eDiag.nbZero > 0)
        std::cout << "The matrix is positive semidefinite\n";
      else
        std::cout << "The matrix is positive definite\n";
    }
    std::cerr << "Normal termination of CheckPositiveSemiDefinite\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CheckPositiveSemiDefinite\n";
    exit(e.eVal);
  }
  runtime(time);
}
