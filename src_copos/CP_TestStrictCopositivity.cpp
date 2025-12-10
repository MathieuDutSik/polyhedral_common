// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

template <typename T, typename Tint>
void compute(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os_out) {
  std::cerr << "Reading input\n";
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  std::cerr << "eSymmMat=\n";
  WriteMatrix(std::cerr, eSymmMat);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  //
  CopositivityTestResult<Tint> eResult =
      TestStrictCopositivity<T, Tint>(eSymmMat, InitialBasis, std::cerr);
  //
  WriteCopositivityTestResult(os_out, OutFormat, eResult);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestStrictCopositivity [arith] [DATASYMM]\n";
      std::cerr << "or\n";
      std::cerr << "CP_TestStrictCopositivity [arith] [DATASYMM] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: The chosen arithmetic\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "OutFormat: classic or GAP. Default is classic\n";
      std::cerr << "OutFile: File to the utput. If absent then it goes to "
                   "std::cerr\n";
      std::cerr << "\n";
      std::cerr
          << "It returns true if the matrix is copositive. If not it returns a "
             "non-negative vector V with A[V] < 0\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f_print = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute<T, Tint>(FileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    //
    print_stderr_stdout_file(OutFile, f_print);
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestStrictCopositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
