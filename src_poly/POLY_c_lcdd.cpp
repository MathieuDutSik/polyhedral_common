// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_PolytopeFct.h"
#include "POLY_c_cddlib.h"

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_c_lcdd [DATAEXT]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT (in) : The polytope vertices\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
    //
    vectface ListIncd = cbased_cdd::DualDescription_incd(EXT);
    //
    int nbFace = ListIncd.size();
    std::ofstream os(argv[2]);
    std::cerr << "nbFace=" << nbFace << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_c_lcdd\n";
    exit(e.eVal);
  }
  runtime(time1);
}
