// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "generalized_polytopes.h"
// clang-format on


template <typename T>
void process(std::string const &GenPolyFile,
             std::string const &OutFormat,
             std::string const &OutFile) {
  std::vector<MyMatrix<T>> list_ext = ReadListMatrixFile<T>(GenPolyFile);
  int dim = list_ext[0].cols();
  GeneralizedPolytope<T> gp = list_ext_to_generalizedpolytope(dim, list_ext);
  MyMatrix<T> EXT = get_vertices_gp(gp, std::cerr);

  auto f_print=[&](std::ostream& os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return " << StringMatrixGAP(EXT) << ";\n";
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
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PolyGen_vertices arith [ListEXT] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "PolyGen_vertices arith [ListEXT]\n";
      std::cerr << "\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string PolyFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f=[&]() -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(PolyFile, OutFormat, OutFile);
      }
      std::cerr << "Error for the template parameter arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination PolyGen_vertices\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PolyGen_difference\n";
    exit(e.eVal);
  }
  runtime(time);
}
