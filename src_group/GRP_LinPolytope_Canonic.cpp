// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PolytopeEquiStab.h"
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
    print_stderr_stdout_file(OutFile, f_print);
}



int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Canonic [arith] [EXTIN] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Canonic [arith] [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      std::cerr << "OutCan : The canonicalization file\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExt = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f_work=[&]() -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileExt, OutFormat, OutFile);
      }
      std::cerr << "GRP_LinPolytope_Canonic: No matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f_work();
    std::cerr << "Normal termination of GRP_LinPolytope_Canonic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
