// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeStabEquiCan.h"
// clang-format on

template<typename T, typename Tint>
void ComputeCanonical(std::string const& FileI, std::string const& OutFormat, std::ostream& os) {
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
  Canonic_PosDef<T, Tint> RetF =
    ComputeCanonicalForm<T, Tint>(eMat, std::cerr);
  if (OutFormat == "CPP") {
    WriteMatrix(os, RetF.Mat);
    return;
  }
  if (OutFormat == "PYTHON") {
    os << "{\"Basis\":" << StringMatrixPYTHON(RetF.Basis) << ", \"SHV\":" << StringMatrixPYTHON(TransposedMat(RetF.SHV)) << ", \"eG\":" << StringMatrixPYTHON(RetF.Mat) << "}\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return ";
    WriteMatrixGAP(os, RetF.Mat);
    os << ";\n";
    return;
  }
  if (OutFormat == "GAP_full") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, RetF.Basis);
    os << ", SHV:=";
    WriteMatrixGAP(os, TransposedMat(RetF.SHV));
    os << ", eG:=";
    WriteMatrixGAP(os, RetF.Mat);
    os << ");\n";
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
      std::cerr << "LATT_Canonicalize [arith] [GramI]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f=[&](std::ostream& os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeCanonical<T,Tint>(FileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    //
    std::cerr << "Normal termination of LATT_Canonicalize\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Canonicalize\n";
    exit(e.eVal);
  }
  runtime(time);
}
