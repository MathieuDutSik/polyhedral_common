// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_General.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    FullNamelist eFull = NAMELIST_GetOneTSPACE();
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ConvertTspace [TspaceNamelistFile] [FILEOUT]\n";
      std::cerr << "\n";
      std::cerr << "TspaceNamelistFile : The namelist file containing the "
                   "description of the T-space\n";
      std::cerr << "FILEOUT            : The output file of the T-space\n";
      return -1;
    }
    //
    std::string TspaceNamelistFile = argv[1];
    std::string FILEOUT = argv[2];
    //
    NAMELIST_ReadNamelistFile(TspaceNamelistFile, eFull);
    //
    SingleBlock BlockTSPACE = eFull.ListBlock.at("TSPACE");
    LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
    WriteLinSpaceFile(FILEOUT, LinSpa);
    std::cerr << "Normal termination of LATT_ConvertTspace\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ConvertTspace\n";
    exit(e.eVal);
  }
  runtime(time);
}
