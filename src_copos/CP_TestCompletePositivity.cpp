// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "StrictPositivity.h"
// clang-format on

template<typename T, typename Tint>
void compute(std::string const& FileI, std::string const& OutFormat, std::ostream & os) {
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  TestStrictPositivity<T, Tint> StrictPos =
    TestingAttemptStrictPositivity<T, Tint>(eSymmMat, InitialBasis,
                                            std::cerr);
  WriteStrictPositivityResult(os, OutFormat, StrictPos);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCompletePositivity [arith] [eMat]\n";
      std::cerr << "or\n";
      std::cerr << "CP_TestCompletePositivity [arith] [eMat] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: The chosen arithmetic\n";
      std::cerr << "eMat: the symmetric matrix which we want to test\n";
      std::cerr << "OutFormat: classic or GAP. Default value is classic\n";
      std::cerr << "OutFile: if present, output goes to OutFile, otherwise to "
                   "std::cerr\n";
      std::cerr << "\n";
      std::cerr << "If completely positive, we return an expression of it "
                   "using integral vector\n";
      std::cerr << "If not completely positive, we return a copositive matrix "
                   "having non-negative scalar product with it\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute<T,Tint>(FileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    //
    if (OutFile == "stderr") {
      f(std::cerr);
    } else {
      if (OutFile == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(OutFile);
        f(os);
      }
    }
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestCompletePositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
