// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "PolytopeEquiStab.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template<typename T>
void process(std::string const& FileExt, std::string const& OutFormat, std::string const& OutFile) {
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileExt);
    size_t threshold = THRESHOLD_USE_SUBSET_SCHEME_CANONIC;

    MyMatrix<T> EXT_can = LinPolytope_CanonicForm<T>(EXT, threshold, std::cerr);
    auto f_print=[&](std::ostream& osf) -> void {
      if (OutFormat == "GAP") {
        osf << "return ";
        WriteMatrixGAP(osf, EXT_can);
        osf << ";\n";
        return;
      }
      if (OutFormat == "CPP") {
        WriteMatrix(osf, EXT_can);
        return;
      }
      std::cerr << "GRP_LinPolytope_Canonic: No matching entry for OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f_print);
}



int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileExt;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    add_rational_arithmetics(dispatch, [&]<typename T>() -> void {
      process<T>(FileExt, OutFormat, OutFile);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Canonic [arith] [EXTIN] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Canonic [arith] [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      std::cerr << "OutCan : The canonicalization file\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileExt = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    dispatch.run(arith);
    std::cerr << "Normal termination of GRP_LinPolytope_Canonic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
