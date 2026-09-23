// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "StrictPositivity.h"
// clang-format on

template <typename T>
void compute(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  TestStrictPositivity<T, Tint> StrictPos =
      TestingAttemptStrictPositivity<T, Tint>(eSymmMat, InitialBasis,
                                              std::cerr);
  WriteStrictPositivityResult(os, OutFormat, StrictPos);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCompletePositivity [arith] [eMat]\n";
      std::cerr << "or\n";
      std::cerr
          << "CP_TestCompletePositivity [arith] [eMat] [OutFormat] [OutFile]\n";
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
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f_print = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return compute<T>(FileI, OutFormat, os_out);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return compute<T>(FileI, OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return compute<T>(FileI, OutFormat, os_out);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return compute<T>(FileI, OutFormat, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Allowed values: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    //
    FILE_PrintStderrStdoutFile(OutFile, f_print);
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_TestCompletePositivity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
