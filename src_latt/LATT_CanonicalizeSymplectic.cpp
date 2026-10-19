// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticeStabEquiCan.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeCanonicalSymplectic(std::string const &FileI,
                                std::string const &OutFormat,
                                std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
  MyMatrix<Tint> B = ComputeCanonicalFormSymplectic<T, Tint>(eMat, std::cerr);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> eMat_red = B_T * eMat * B_T.transpose();
  if (OutFormat == "CPP") {
    WriteMatrix(os, eMat_red);
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, B);
    os << ", eG:=";
    WriteMatrixGAP(os, eMat_red);
    os << ");\n";
    return;
  }
  std::cerr << "LATT_CanonicalizeSymplectic: No matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_CanonicalizeSymplectic [arith] [GramI] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_CanonicalizeSymplectic [arith] [GramI]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
      std::cerr << "OutFormat values:\n";
      std::cerr << "  CPP : only the reduced matrix is in output\n";
      std::cerr << "  GAP : the basis and the reduced matrix in GAP format\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeCanonicalSymplectic<T>(FileI, OutFormat, os_out);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_CanonicalizeSymplectic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CanonicalizeSymplectic\n";
    exit(e.eVal);
  }
  runtime(time);
}
