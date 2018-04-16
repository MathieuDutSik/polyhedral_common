#include "NumberTheory.h"
#include "LatticeDefinitions.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_lll [DATAIN] [DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN  : The Gram matrix on input\n";
      std::cerr << "DATAOUT : The Gram reduced matrix + the transformation matrix (GAP readable)\n";
      return -1;
    }
    //
    std::ifstream is(argv[1]);
    MyMatrix<mpq_class> GramMat=ReadMatrix<mpq_class>(is);
    //
    std::ofstream os(argv[2]);
    LLLreduction<mpq_class,mpz_class> recLLL = LLLreducedBasis<mpq_class,mpz_class>(GramMat);
    os << "return rec(GramMat:=";
    WriteMatrixGAP(os, recLLL.GramMatRed);
    os << ", Pmat:=";
    WriteMatrixGAP(os, recLLL.Pmat);
    os << ");\n";
    std::cerr << "recLLL.GramMatRed=\n";
    WriteMatrix(std::cerr, recLLL.GramMatRed);
    std::cerr << "recLLL.Pmat=\n";
    WriteMatrix(std::cerr, recLLL.Pmat);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
