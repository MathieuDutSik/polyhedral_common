// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Copositivity.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void compute(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::cerr << "Reading input\n";
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  std::cerr << "eSymmMat=\n";
  WriteMatrix(std::cerr, eSymmMat);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  //
  CopositivityTestResult<Tint> eResult =
      TestCopositivity<T, Tint>(eSymmMat, InitialBasis, std::cerr);
  //
  WriteCopositivityTestResult(os_out, OutFormat, eResult);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileI;
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      compute<T>(FileI, OutFormat, *os_ptr);
    };
    add_rational_arithmetics_by_library(dispatch, f_arith);
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCopositivity [arith] [DATASYMM]\n";
      std::cerr << "or\n";
      std::cerr
          << "CP_TestCopositivity [arith] [DATASYMM] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: " << dispatch.choices() << "\n";
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
    FileI = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f_print = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    //
    FILE_PrintStderrStdoutFile(OutFile, f_print);
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestCopositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
