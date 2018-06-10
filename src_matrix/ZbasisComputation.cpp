#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  using T=mpq_class;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "ZbasisComputation [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat=ReadMatrix<T>(INmat);
    // computing the kernel
    MyMatrix<T> TheBasis = GetZbasis(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, TheBasis);
    os << ";\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
