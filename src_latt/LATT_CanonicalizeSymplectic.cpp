// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LatticeStabEquiCan.h"
// clang-format on

template <typename T, typename Tint>
void ComputeCanonicalSymplectic(std::string const &FileI,
                                std::string const &OutFormat,
                                std::ostream &os) {
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
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeCanonicalSymplectic<T, Tint>(FileI, OutFormat, os_out);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeCanonicalSymplectic<T, Tint>(FileI, OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeCanonicalSymplectic<T, Tint>(FileI, OutFormat, os_out);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return ComputeCanonicalSymplectic<T, Tint>(FileI, OutFormat, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_CanonicalizeSymplectic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CanonicalizeSymplectic\n";
    exit(e.eVal);
  }
  runtime(time);
}
