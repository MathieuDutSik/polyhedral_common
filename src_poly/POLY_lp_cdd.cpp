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
    using T=mpq_class;
    MyMatrix<T> ListIneq=ReadMatrixLrsCdd<T>(is);
    std::string eChoice;
    is >> eChoice;
    if (eChoice != "minimize") {
      std::cerr << "Only the minimize is available here\n";
      throw TerminalException{1};
    }
    int nbCol=ListIneq.cols();
    MyVector<T> ToBeMinimized(nbCol);
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal;
      is >> eVal;
      ToBeMinimized(iCol)=eVal;
    }
    LpSolution<T> eSol=CDD_LinearProgramming(ListIneq, ToBeMinimized);
    //
    std::cerr << "PrimalDefined = " << eSol.PrimalDefined << "\n";
    std::cerr << "  DualDefined = " << eSol.DualDefined << "\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
