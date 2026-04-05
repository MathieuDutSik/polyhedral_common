// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeStabEquiCan.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_canonicalizeSymplectic [GramI] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_canonicalizeSymplectic [GramI]\n";
      std::cerr << "\n";
      std::cerr << "If opt=1 then only the matrix is in output.\n";
      std::cerr << "If opt=2 then the basis, list of vectors and matrix is in "
                   "GAP formatted output\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    //
    MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
    MyMatrix<Tint> B = ComputeCanonicalFormSymplectic<T, Tint>(eMat, std::cerr);
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
    MyMatrix<T> eMat_red = B_T * eMat * B_T.transpose();
    //
    auto f=[&](std::ostream& os_out) -> void {
      if (OutFormat == "CPP") {
        WriteMatrix(os_out, eMat_red);
        return;
      }
      if (OutFormat == "GAP") {
        os_out << "return rec(Basis:=";
        WriteMatrixGAP(os_out, B);
        os_out << ", eG:=";
        WriteMatrixGAP(os_out, eMat_red);
        os_out << ");\n";
        return;
      }
      std::cerr << "LATT_CanonicalizeSymplectic: No matching OutFormat\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, f);
    std::cerr << "Normal termination of LATT_canonicalizeSymplectic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_canonicalizeSymplectic\n";
    exit(e.eVal);
  }
  runtime(time);
}
