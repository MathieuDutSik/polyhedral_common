// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "POLY_ThreadDualDescription.h"
int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    Eigen::initParallel();
    FullNamelist eFull = NAMELIST_GetStandard_TEMP_THREADED_ADM();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ThreadedADM [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    using T = mpq_class;
    using Tgroup = TheGroupFormat<mpz_class>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionComputeDualDesc<T, Tgroup>(eFull);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_ThreadedADM\n";
    exit(e.eVal);
  }
  runtime(time1);
}
