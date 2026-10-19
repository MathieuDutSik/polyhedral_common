// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "FiniteMatrixGroupTest.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const& FileListMat, std::string const& OutFormat, std::string const& OutFile) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::vector<MyMatrix<Tint>> ListMat = ReadListMatrixFile<Tint>(FileListMat);

  bool is_finite = test_finiteness_group<T,Tint>(ListMat, std::cerr);
  auto f_print=[&](std::ostream& os) -> void {
    if (OutFormat == "Raw") {
      os << "is_finite=" << is_finite << "\n";
      return;
    }
    if (OutFormat == "GAP") {
      os << "return rec(is_finite:=" << GAP_logical(is_finite) << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_TestFiniteness [arith] [ListMat] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "GRP_TestFiniteness [arith] [ListMat]\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileListMat = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f=[&]() -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            process<T>(FileListMat, OutFormat, OutFile);
          });
      dispatch.run(arith);
    };
    f();
    std::cerr << "Normal termination of GRP_TestFiniteness\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_TestFiniteness\n";
    exit(e.eVal);
  }
  runtime(time);
}
