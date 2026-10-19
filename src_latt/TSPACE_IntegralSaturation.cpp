// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Tspace_ListMatSaturation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void IntegralSaturation(std::string const &FileI, std::string const &OutFormat,
                        std::ostream &os) {
  std::vector<MyMatrix<T>> ListMatA = ReadListMatrixFile<T>(FileI);
  std::vector<MyMatrix<T>> ListMatB = IntegralSaturationSpace(ListMatA, std::cerr);
  if (OutFormat == "GAP") {
    os << "return ";
    WriteListMatrixGAP(os, ListMatB);
    os << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_IntegralSaturation [arith] [ListGram]\n";
      std::cerr << "    or\n";
      std::cerr << "TSPACE_IntegralSaturation [arith] [ListGram] [OutFormat] "
                << "[OutFile]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics(
          dispatch, [&]<typename T>() -> void {
            IntegralSaturation<T>(FileI, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of TSPACE_IntegralSaturation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_IntegralSaturation\n";
    exit(e.eVal);
  }
  runtime(time);
}
