#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_lp_cdd [DATAIN\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The inequalities and function to be minimized\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    MyMatrix<mpq_class> ListIneq=ReadMatrixLrsCdd<mpq_class>(is);
    std::string eChoice;
    is >> eChoice;
    if (eChoice != "minimize") {
      std::cerr << "Only the minimize is available here\n";
      throw TerminalException{1};
    }
    int nbCol=ListIneq.cols();
    MyVector<mpq_class> ToBeMinimized(nbCol);
    for (int iCol=0; iCol<nbCol; iCol++) {
      mpq_class eVal;
      is >> eVal;
      ToBeMinimized(iCol)=eVal;
    }
    LpSolution<mpq_class> eSol=CDD_LinearProgramming(ListIneq, ToBeMinimized);
    //
    std::cerr << "PrimalDefined = " << eSol.PrimalDefined << "\n";
    std::cerr << "  DualDefined = " << eSol.DualDefined << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
