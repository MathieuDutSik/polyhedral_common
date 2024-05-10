// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "StrictPositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc == 1) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCompletePositivity [eMat]\n";
      std::cerr << "CP_TestCompletePositivity [eMat] [OutFormat]\n";
      std::cerr << "CP_TestCompletePositivity [eMat] [OutFormat] [OutFile]\n";
      std::cerr << "CP_TestCompletePositivity [eMat] [OutFormat] [OutFile] "
                   "[InitialBasis]\n";
      std::cerr << "CP_TestCompletePositivity [eMat] [InitialBasis]\n";
      std::cerr << "\n";
      std::cerr << "eMat: the symmetric matrix which we want to test\n";
      std::cerr << "OutFormat: classic or GAP. Default value is classic\n";
      std::cerr << "OutFile: if present, output goes to OutFile, otherwise to "
                   "std::cerr\n";
      std::cerr << "InitialBasis: IdentityMat if absent, otherwise the one in "
                   "the named file\n";
      std::cerr << "\n";
      std::cerr << "If completely positive, we return an expression of it "
                   "using integral vector\n";
      std::cerr << "If not completely positive, we return a copositive matrix "
                   "having non-negative scalar product with it\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    //
    std::string OutFormat = "classic";
    if (argc >= 3)
      OutFormat = argv[2];
    //
    auto process = [&](std::ostream &os) -> void {
      MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
      if (argc >= 5)
        InitialBasis = ReadMatrixFile<Tint>(argv[4]);
      TestStrictPositivity<T, Tint> StrictPos =
          TestingAttemptStrictPositivity<T, Tint>(eSymmMat, InitialBasis,
                                                  std::cerr);
      WriteStrictPositivityResult(os, OutFormat, StrictPos);
    };
    if (argc >= 4) {
      std::string file = argv[3];
      std::ofstream os(file);
      process(os);
    } else {
      process(std::cerr);
    }
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestCompletePositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
