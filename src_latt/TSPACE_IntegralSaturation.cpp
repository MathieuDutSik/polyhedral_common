// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_ListMatSaturation.h"
// clang-format on

template <typename T>
void IntegralSaturation(std::string const &FileI, std::string const &OutFormat,
                        std::ostream &os) {
  std::vector<MyMatrix<T>> ListMatA = ReadListMatrixFile<T>(FileI);
  std::vector<MyMatrix<T>> ListMatB = IntegralSaturationSpace(ListMatA, std::cerr);
  if (OutFormat == "GAP") {
    os << "return ";
    WriteListMatrixGAP(os, ListMatB);
    os << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_IntegralSaturation [arith] [ListGram]\n";
      std::cerr << "    or\n";
      std::cerr << "TSPACE_IntegralSaturation [arith] [ListGram] [OutFormat] [OutFile]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return IntegralSaturation<T>(FileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    //
    std::cerr << "Normal termination of TSPACE_IntegralSaturation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_IntegralSaturation\n";
    exit(e.eVal);
  }
  runtime(time);
}
