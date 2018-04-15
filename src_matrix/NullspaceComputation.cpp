#include "MAT_Matrix.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  using T=mpq_class;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputation [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat=ReadMatrix<T>(INmat);
    // computing the kernel
    Eigen::FullPivLU<MyMatrix<T>> lu(TheMat);
    MyMatrix<T> TheKer = lu.kernel();

    
    //    MyMatrix<T> TheKer=NullspaceMat(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, TheKer);
    os << ";\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
