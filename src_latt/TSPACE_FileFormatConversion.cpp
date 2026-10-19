// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Tspace_General.h"
#include "Group.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ConvertTspace(std::string const &TspaceNamelistFile,
                   std::string const &OutFormat, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  FullNamelist eFull = NAMELIST_GetOneTSPACE();
  NAMELIST_ReadNamelistFile(TspaceNamelistFile, eFull);
  //
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
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
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_FileFormatConversion [arith] "
                   "[TspaceNamelistFile] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "TSPACE_FileFormatConversion [arith] "
                   "[TspaceNamelistFile]\n";
      std::cerr << "\n";
      std::cerr << "arith              : The chosen arithmetic (see below)\n";
      std::cerr << "TspaceNamelistFile : The namelist file containing the "
                   "description of the T-space\n";
      std::cerr << "OutFormat          : CPP or GAP\n";
      std::cerr << "OutFile            : The output file of the T-space\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string TspaceNamelistFile = argv[2];
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ConvertTspace<T>(TspaceNamelistFile, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of TSPACE_FileFormatConversion\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_FileFormatConversion\n";
    exit(e.eVal);
  }
  runtime(time);
}
