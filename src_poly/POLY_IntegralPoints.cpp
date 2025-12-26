// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "POLY_PolytopeInt.h"
// clang-format on

template<typename T, typename Tint>
void process(std::string const& FileFac, std::string const& method, std::string const& OutFormat, std::string const& OutFile) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(FileFac);
  int n_col = FAC.cols();

  auto get_vert=[&]() -> std::vector<MyVector<Tint>> {
    if (method == "LP_no_LLL") {
      return GetListIntegralPoint_LP<T,Tint>(FAC, std::cerr);
    }
    if (method == "ITER_no_LLL") {
      return GetListIntegralPoint_ITER<T,Tint>(FAC, std::cerr);
    }
    std::cerr << "POLYINT: Failed to find a matching entry for get_vert\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<Tint>> ListVert = get_vert();
  MyMatrix<Tint> M_vert = MatrixFromVectorFamilyDim(n_col, ListVert);
  auto f_print=[&](std::ostream& osf) -> void {
    if (OutFormat == "GAP") {
      osf << "return ";
      WriteMatrixGAP(osf, M_vert);
      osf << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      WriteMatrix(osf, M_vert);
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}




int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IntegralPoints arith method [FAC] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_IntegralPoints arith method [FAC]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string method = argv[2];
    std::string FileFac = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&]() -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint>(FileFac, method, OutFormat, OutFile);
      }
      std::cerr << "Failed to find a matching arithmetic for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of POLY_IntegralPoints\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_IntegralPoints\n";
    exit(e.eVal);
  }
  runtime(time1);
}
