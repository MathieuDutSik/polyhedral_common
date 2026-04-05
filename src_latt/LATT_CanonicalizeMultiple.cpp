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
      std::cerr << "LATT_canonicalizeMultiple [FileMatrices] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_canonicalizeMultiple [FileMatrices]\n";
      std::cerr << "\n";
      std::cerr << "If opt=1 then only the matrix is in output.\n";
      std::cerr << "If opt=2 then the basis, list of vectors and matrix is in "
                   "GAP formatted output\n";
      std::cerr << "FileGrams (input) : The list of matrices used on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    using T = mpz_class;
    using Tint = mpz_class;
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    //
    std::vector<MyMatrix<T>> ListMatrix = ReadListMatrixFile<T>(FileI);
    MyMatrix<Tint> B = ComputeCanonicalFormMultiple<T, Tint>(ListMatrix, std::cerr);
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
    MyMatrix<T> Mat_red = B_T * ListMatrix[0] * B_T.transpose();
    //
    auto f=[&](std::ostream& os_out) -> void {
      if (OutFormat == "CPP") {
        WriteMatrix(os_out, Mat_red);
        return;
      }
      if (OutFormat == "GAP") {
        os_out << "return rec(Basis:=";
        WriteMatrixGAP(os_out, B);
        os_out << ", eG:=";
        WriteMatrixGAP(os_out, Mat_red);
        os_out << ");\n";
        return;
      }
      std::cerr << "OutFormat does not match anything\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, f);
    std::cerr << "Normal termination of LATT_canonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_canonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
