#include "Temp_Positivity.h"
int main(int argc, char *argv[])
{
  using T=mpq_class;
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "CheckPositiveSemiDefinite [inputMat]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat=ReadMatrix<T>(INmat);
    // computing the kernel
    DiagSymMat<T> eDiag = DiagonalizeNonDegenerateSymmetricMatrix(TheMat);
    if (eDiag.nbMinus > 0)
      std::cout << "The matrix is NOT positive semidefinite\n";
    else {
      if (eDiag.nbZero > 0)
        std::cout << "The matrix is positive semidefinite\n";
      else
        std::cout << "The matrix is positive definite\n";
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Normal termination of the program\n";
}
