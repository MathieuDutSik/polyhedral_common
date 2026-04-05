// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveListCone [DATASYMM] [OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "CP_CopositiveListCone [DATASYMM]\n";
      std::cerr << "\n";
      std::cerr
          << "DATASYMM: The input data of the copositive symmetric matrix\n";
      std::cerr << "It returns the list of cones that is used to test that the "
                   "matrix is copositive\n";
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
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
    //
    int n = eSymmMat.rows();
    MyMatrix<Tint> TheBasis = IdentityMat<Tint>(n);
    ResultListCone<Tint> res =
        EnumerateListConeCopositive<T, Tint>(eSymmMat, TheBasis, std::cerr);
    auto f_process=[&](std::ostream& os_out) -> void {
      size_t nbCone = res.ListCone.size();
      if (OutFormat == "CPP") {
        if (res.test == false) {
          os_out << "Matrix is not Copositive\n";
        } else {
          for (size_t iCone = 0; iCone < nbCone; iCone++) {
            os_out << "iCone=" << iCone << "/" << nbCone << " Basis\n";
            MyMatrix<Tint> const &eMat = res.ListCone[iCone];
            WriteMatrix(os_out, eMat);
          }
        }
      }
      if (OutFormat == "GAP") {
        os_out << "return rec(";
        if (res.test == false) {
          os_out << "copositive:=false";
        } else {
          os_out << "copositive:=true";
          os_out << ", ListCone:=";
          WriteListMatrixGAP(os_out, res.ListCone);
        }
        os_out << ");\n";
      }
    };
    print_stderr_stdout_file(FileO, f_process);
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_CopositiveListCone\n";
    exit(e.eVal);
  }
  runtime(time1);
}
