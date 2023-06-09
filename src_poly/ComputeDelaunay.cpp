// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeDelaunay.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    Eigen::initParallel();
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_DELAUNAY();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ComputeDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    using Tgroup = TheGroupFormat<mpz_class>;
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    TreatDelaunayEntry<T, Tint, Tgroup>(eFull);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ComputeDelaunay\n";
    exit(e.eVal);
  }
  runtime(time1);
}
