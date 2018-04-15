#include "Temp_MatrixInt.h"
using namespace std;
int main(int argc, char *argv[])
{
  std::ifstream EXTfs;
  if (argc != 2) {
    fprintf(stderr, "Number of argument is = %d\n", argc);
    fprintf(stderr, "This program is used as\n");
    fprintf(stderr, "NullspaceIntMat [DATAEXT]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DATAEXT: The input data of the matrix\n");
    fprintf(stderr, "It returns the integral nullspace of the matrix\n");
    return -1;
  }
  //
  fprintf(stderr, "Reading input\n");
  EXTfs.open(argv[1]);
  MyMatrix<int> EXT=ReadMatrix<int>(EXTfs);
  EXTfs.close();
  MyMatrix<int> NSP=NullspaceIntMat(EXT);
  WriteMatrix(std::cout, NSP);
  //
  fprintf(stderr, "Completion of the program\n");
}
