// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_PolytopeInt.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const& FileFac, std::string const& method, std::string const& OutFormat, std::string const& OutFile) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> FAC = ReadMatrixFile<T>(FileFac);
  int n_col = FAC.cols();

  auto get_vert=[&]() -> std::vector<MyVector<Tint>> {
    if (method == "LP_no_LLL") {
      return GetListIntegralPoint_LP_no_LLL<T,Tint>(FAC, std::cerr);
    }
    if (method == "ITER_no_LLL") {
      return GetListIntegralPoint_ITER_no_LLL<T,Tint>(FAC, std::cerr);
    }
    if (method == "ITER") {
      return GetListIntegralPoint_ITER<T,Tint>(FAC, std::cerr);
    }
    if (method == "LP") {
      return GetListIntegralPoint_LP<T,Tint>(FAC, std::cerr);
    }
    std::cerr << "POLYINT: Failed to find a matching entry for get_vert\n";
    throw TerminalException{1};
  };
  std::vector<MyVector<Tint>> ListVert = get_vert();
  MyMatrix<Tint> M_vert = MatrixFromVectorFamilyDim(n_col-1, ListVert);
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
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}




int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string method, FileFac;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileFac, method, OutFormat, OutFile);
    };
    add_rational_arithmetics_by_library(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IntegralPoints arith method [FAC] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_IntegralPoints arith method [FAC]\n";
      std::cerr << "\n";
      std::cerr << "arith: " << dispatch.choices() << "\n";
      return -1;
    }
    std::string arith = argv[1];
    method = argv[2];
    FileFac = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    dispatch.run(arith);
    std::cerr << "Normal termination of POLY_IntegralPoints\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_IntegralPoints\n";
    exit(e.eVal);
  }
  runtime(time1);
}
