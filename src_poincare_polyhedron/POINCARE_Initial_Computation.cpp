// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "poincare_polyhedron.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetPoincareInput();
    if (argc != 2) {
      std::cerr << "POINCARE_Initial_Computation [FileNML]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      throw TerminalException{1};
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    RecOption rec_option = ReadInitialData(eFull);
    Process_rec_option(rec_option);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IndefiniteReduction\n";
    exit(e.eVal);
  }
  runtime(time1);
}
