// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticeDelaunay.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &FileM, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileM);
  HumanTime time_total;
  CVPSolver<T, Tint> solver(GramMat, std::cerr);
  MyMatrix<Tint> EXT = FindDelaunayPolytope<T, Tint>(solver, std::cerr);
  std::cerr << "|FindDelaunayPolytope|=" << time_total << "\n";
  if (OutFormat == "GAP") {
    os << "return ";
    WriteMatrix(os, EXT);
    os << ";\n";
    return;
  }
  std::cerr << "No type available for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_FindOneVertex [arith] [FileM] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_FindOneVertex [arith] [FileM]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileM = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            process<T>(FileM, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_FindOneVertex\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindOneVertex\n";
    exit(e.eVal);
  }
  runtime(time1);
}
