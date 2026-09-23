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
// clang-format on

template <typename T, typename Tint>
void ConvertTspace(std::string const &TspaceNamelistFile,
                   std::string const &OutFormat, std::ostream &os) {
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ConvertTspace<T, Tint>(TspaceNamelistFile, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ConvertTspace<T, Tint>(TspaceNamelistFile, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ConvertTspace<T, Tint>(TspaceNamelistFile, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ConvertTspace<T, Tint>(TspaceNamelistFile, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of TSPACE_FileFormatConversion\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_FileFormatConversion\n";
    exit(e.eVal);
  }
  runtime(time);
}
