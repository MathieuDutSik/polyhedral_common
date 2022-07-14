// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "StrictPositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCompletePositivity [eMat]\n";
      std::cerr << "or\n";
      std::cerr << "CP_TestCompletePositivity [eMat] [InitialBasis]\n";
      std::cerr << "\n";
      std::cerr << "eMat: the symmetric matrix which we want to test\n";
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
    MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
    if (argc == 3)
      InitialBasis = ReadMatrixFile<Tint>(argv[2]);
    //
    TestStrictPositivity<T, Tint> StrictPos =
        TestingAttemptStrictPositivity<T, Tint>(eSymmMat, InitialBasis);
    WriteStrictPositivityResult(std::cerr, StrictPos);
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
