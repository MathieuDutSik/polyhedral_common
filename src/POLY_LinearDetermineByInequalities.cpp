#include "POLY_LinearProgramming.h"

int main(int argc, char *argv[])
{
  try {
    std::vector<std::vector<int> > TheReturn;
    VectVectInt TheOutput;
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LinearDeterminedByInequalities [DATAFAC]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT: The input data of the polytope vertices\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream Ifs(argv[1]);
    MyMatrix<mpq_class> FAC=ReadMatrix<mpq_class>(Ifs);
    //
    MyMatrix<mpq_class> LinSpace=LinearDeterminedByInequalities(FAC);
    std::cerr << "LinSpace=\n";
    WriteMatrix(std::cerr, LinSpace);
    //
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
