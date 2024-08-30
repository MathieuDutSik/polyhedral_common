// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc < 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCopositivity [DATASYMM]\n";
      std::cerr << "CP_TestCopositivity [DATASYMM] [OutFormat]\n";
      std::cerr << "CP_TestCopositivity [DATASYMM] [OutFormat] [OutFile]\n";
      std::cerr << "CP_TestCopositivity [DATASYMM] [OutFormat] [OutFile] "
                   "[InitialBasis]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "OutFormat: classic or GAP. Default is classic\n";
      std::cerr << "OutFile: File to the utput. If absent then it goes to "
                   "std::cerr\n";
      std::cerr << "InitialBasis: If missing, the standard basis, otherwise "
                   "the basis from the file\n";
      std::cerr << "\n";
      std::cerr
          << "It returns true if the matrix is copositive. If not it returns a "
             "non-negative vector V with A[V] < 0\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T = mpq_class;
    using Tint = int;
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    std::string OutFormat = "classic";
    if (argc >= 3) {
      OutFormat = argv[2];
    }
    //
    auto process = [&](std::ostream &os) -> void {
      MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
      if (argc >= 5)
        InitialBasis = ReadMatrixFile<Tint>(argv[4]);
      //
      std::pair<SingleTestResult<Tint>, size_t> eResult =
          TestCopositivity<T, Tint>(eSymmMat, InitialBasis, std::cerr);
      //
      WriteSingleTestResult(os, OutFormat, eResult);
    };
    if (argc >= 4) {
      std::ofstream os(argv[3]);
      process(os);
    } else {
      process(std::cerr);
    }
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestCopositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
