// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_Functions.h"
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
  print_stderr_stdout_file(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_TestFiniteness [arith] [ListMat] [OutFormat] [OutFile]\n";
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint>(FileListMat, OutFormat, OutFile);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of GRP_TestFiniteness\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_TestFiniteness\n";
    exit(e.eVal);
  }
  runtime(time);
}
