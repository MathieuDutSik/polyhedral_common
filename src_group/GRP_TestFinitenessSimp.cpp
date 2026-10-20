// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "FiniteMatrixGroupFund.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template<typename Tint>
void process(std::string const& FileListMat, std::string const& OutFormat, std::string const& OutFile) {
  std::vector<MyMatrix<Tint>> ListMat = ReadListMatrixFile<Tint>(FileListMat);


  std::vector<MyMatrix<Tint>> ListMatRed = ExhaustiveReductionComplexityGroupMatrix(ListMat, std::cerr);
  FiniteMatrixGroupTest<Tint> fmgt(ListMatRed, std::cerr);
  bool is_finite = fmgt.full_resolution();
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
    std::string FileListMat;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process<Tint>(FileListMat, OutFormat, OutFile);
        });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_TestFinitenessSimp [arith] [ListMat] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "GRP_TestFinitenessSimp [arith] [ListMat]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileListMat = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    dispatch.run(arith);
    std::cerr << "Normal termination of GRP_TestFinitenessSimp\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_TestFinitenessSimp\n";
    exit(e.eVal);
  }
  runtime(time);
}
