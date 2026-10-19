// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "SHORT_Realizability.h"
#include "rational.h"
#include "Arithmetic_dispatch_types.h"

template <typename T>
void automorphism_group(std::string const& FileSHV, std::string const& OutFormat, std::string const& OutFile) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  MyMatrix<Tint> SHV = ReadMatrixFile<Tint>(FileSHV);
  if (RankMat(SHV) != SHV.cols()) {
    std::cerr << "SHORT_AutomorphismGroup: The input matrix in "
              << FileSHV << " is not of full rank\n";
    throw TerminalException{1};
  }

  std::vector<MyMatrix<Tint>> list_mat = SHORT_GetStabilizer<T,Tint,Tgroup>(SHV, std::cerr);

  auto f_out=[&](std::ostream& osf) -> void {
    if (OutFormat == "GAP") {
      osf << "return ";
      WriteListMatrixGAP(osf, list_mat);
      osf << ";\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    std::cerr << "Allowed OutFormat are GAP, TXT\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_out);
}



int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileSHV;
    std::string OutFormat;
    std::string OutFile;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      automorphism_group<T>(FileSHV, OutFormat, OutFile);
    });
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_AutomorphismGroup [arith] [FileSHV]\n";
      std::cerr << "or\n";
      std::cerr << "SHORT_AutomorphismGroup [arith] [FileSHV] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "[arith]        : The arithmetic being used:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "[FileIn]       : The input file of the system\n";
      std::cerr << "[OutFormat]    : The output format, TXT or GAP\n";
      std::cerr << "[OutFile]      : The output file of the program\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    FileSHV = argv[2];
    OutFormat = argv[3];
    OutFile = argv[4];
    //
    dispatch.run(arith);
    //
    //
    std::cerr << "Normal termination of SHORT_AutomorphismGroup\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_AutomorphismGroup\n";
    exit(e.eVal);
  }
  runtime(time1);
}
