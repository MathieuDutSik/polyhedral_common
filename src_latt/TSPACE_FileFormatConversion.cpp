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
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ConvertTspace [TspaceNamelistFile] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_ConvertTspace [TspaceNamelistFile]\n";
      std::cerr << "\n";
      std::cerr << "TspaceNamelistFile : The namelist file containing the "
                   "description of the T-space\n";
      std::cerr << "OutFormat          : CPP or GAP\n";
      std::cerr << "OutFile            : The output file of the T-space\n";
      return -1;
    }
    //
    std::string TspaceNamelistFile = argv[1];
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      OutFile = argv[3];
    }
    //
    NAMELIST_ReadNamelistFile(TspaceNamelistFile, eFull);
    //
    SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
    LinSpaceMatrix<T> LinSpa =
        ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
    auto f = [&](std::ostream &os) -> void {
      if (OutFormat == "CPP") {
        WriteLinSpace(os, LinSpa);
        return;
      }
      if (OutFormat == "GAP") {
        os << "return ";
        WriteLinSpaceGAP(os, LinSpa);
        os << ";\n";
        return;
      }
      std::cerr << "Failed to find a matching entry for OutFormat\n";
      std::cerr << "Allowed choices: CPP, GAP\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of LATT_ConvertTspace\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ConvertTspace\n";
    exit(e.eVal);
  }
  runtime(time);
}
