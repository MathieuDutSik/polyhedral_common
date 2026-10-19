// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Permutation.h"
#include "SHORT_Realizability.h"
#include "rational.h"
#include "Arithmetic_dispatch_types.h"

template <typename T>
void test_realizability(std::string const& FileSHV, std::string const& OutFormat, std::string const& OutFile) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  MyMatrix<Tint> SHV = ReadMatrixFile<Tint>(FileSHV);
  if (RankMat(SHV) != SHV.cols()) {
    std::cerr << "SHORT_TestRealizability: The input matrix in "
              << FileSHV << " is not of full rank\n";
    throw TerminalException{1};
  }

  ReplyRealizability<T, Tint> eRes = SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(SHV, std::cerr);

  auto f_out=[&](std::ostream& osf) -> void {
    if (OutFormat == "TXT") {
      if (eRes.reply) {
        osf << "Realizable with following matrix\n";
        WriteMatrix(osf, eRes.eMat);
      } else {
        osf << "Non Realizable\n";
      }
      return;
    }
    if (OutFormat == "GAP") {
      osf << "return rec(realizable:=";
      if (eRes.reply) {
        osf << "true, matrix:=" << StringMatrixGAP(eRes.eMat);
      } else {
        osf << "false";
      }
      osf << ");\n";
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
    std::string OutFormat = "GAP";
    std::string OutFile = "stdout";
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      test_realizability<T>(FileSHV, OutFormat, OutFile);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_TestRealizability [arith] [FileIn]\n";
      std::cerr << "or\n";
      std::cerr << "SHORT_TestRealizability [arith] [FileIn] [OutFormat] "
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
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    dispatch.run(arith);
    //
    //
    std::cerr << "Normal termination of SHORT_TestRealizability\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_TestRealizability\n";
    exit(e.eVal);
  }
  runtime(time1);
}
