// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
#include "Temp_PerfectForm.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    Eigen::initParallel();
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_PERFECT();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ComputePerfect [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    using Tidx = uint32_t;
    using T = mpq_class;
    using Tint = mpz_class;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    TreatPerfectLatticesEntry<T, Tint, Tgroup>(eFull);
    std::cerr << "Completion of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
