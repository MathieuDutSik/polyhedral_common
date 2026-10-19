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
void compute(std::string const &FileI, std::string const &strMaxNorm,
             std::string const &OutFormat, std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  T MaxNorm = ParseScalar<T>(strMaxNorm);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  //
  CopositivityEnumResult<Tint> CopoRes =
      EnumerateCopositiveShortVector<T, Tint>(eSymmMat, InitialBasis, MaxNorm,
                                              std::cerr);
  //
  WriteCopositivityEnumResult(os_out, OutFormat, eSymmMat, CopoRes);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileI;
    std::string strMaxNorm;
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      compute<T>(FileI, strMaxNorm, OutFormat, *os_ptr);
    };
    add_rational_arithmetics_by_library(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMaxNorm [arith] [DATASYMM] [MaxNorm]\n";
      std::cerr << "or\n";
      std::cerr
          << "CP_CopositiveMaxNorm [arith] [DATASYMM] [MaxNorm] [OutFormat] "
             "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: " << dispatch.choices() << "\n";
      std::cerr << "DATASYMM: The symmetric matrix on input\n";
      std::cerr << "MaxNorm: The maximm norm considered\n";
      std::cerr << "OutFormat: classic or GAP. Default is classic\n";
      std::cerr << "OutFile: if assigned output goes to OutFile. Otherwise to "
                   "std::cerr\n";
      std::cerr << "\n";
      std::cerr << "It returns the list of integer vectors v in Z^n_{>= 0} "
                   "such that A[v] <= MaxNorm\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::string arith = argv[1];
    FileI = argv[2];
    strMaxNorm = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
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
    std::cerr << "Error in CP_CopositiveMaxNorm\n";
    exit(e.eVal);
  }
  runtime(time1);
}
