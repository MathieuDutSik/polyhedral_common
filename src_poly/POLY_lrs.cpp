// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_lrslib.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_lrs [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      return -1;
    }
    //
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
    int nbCol = EXT.cols();
    //
    std::cout << "V-representation\n";
    std::cout << "begin\n";
    std::cout << "****** " << nbCol << " rational\n";
    long nVertices = 0;
    std::function<void(T *)> fPrint = [&](T *out) -> void {
      for (int iCol = 0; iCol < nbCol; iCol++)
        std::cout << " " << out[iCol];
      std::cout << "\n";
      nVertices++;
    };
    lrs::Kernel_DualDescription(EXT, fPrint);
    std::cout << "end\n";
    std::cout << "*Total: nvertices=" << nVertices << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
