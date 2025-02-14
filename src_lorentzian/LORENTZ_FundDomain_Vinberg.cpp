// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "LatticeStabEquiCan.h"
#include "PolytopeEquiStab.h"
#include "vinberg.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_VINBERG();
    if (argc != 2) {
      std::cerr << "LORENTZ_FundDomain_Vinberg [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      throw TerminalException{1};
    }
    std::string eFileName = argv[1];
    using T = mpq_class;
    using Tint = mpz_class;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionVinberg<T, Tint>(eFull, std::cerr);
    std::cerr << "Normal termination of LORENTZ_FundDomain_Vinberg\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_FundDomain_Vinberg\n";
    exit(e.eVal);
  }
  runtime(time);
}
