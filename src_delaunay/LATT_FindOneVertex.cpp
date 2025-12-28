// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeDelaunay.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &FileM, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileM);
  HumanTime time_total;
  CVPSolver<T, Tint> solver(GramMat, std::cerr);
  MyMatrix<Tint> EXT =
      FindDelaunayPolytope<T, Tint>(GramMat, solver, std::cerr);
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileM, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of LATT_FindOneVertex\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindOneVertex\n";
    exit(e.eVal);
  }
  runtime(time1);
}
