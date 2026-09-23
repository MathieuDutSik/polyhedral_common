// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "ClassicLLL.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileI);
  //
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(GramMat, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(GramMat:=";
    WriteMatrixGAP(os, recLLL.GramMatRed);
    os << ", Pmat:=";
    WriteMatrixGAP(os, recLLL.Pmat);
    os << ");\n";
    return;
  }
  if (OutFormat == "CPP_G") {
    WriteMatrix(os, recLLL.GramMatRed);
    return;
  }
  if (OutFormat == "CPP_P") {
    WriteMatrix(os, recLLL.Pmat);
    return;
  }
  std::cerr << "Failed to find a matching OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "LATT_lll arithmetic [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_lll arithmetic [FileI]\n";
      std::cerr << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "arithmetic : The arithmetic, e.g. gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "arithmetic : The arithmetic, e.g. gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      std::cerr << "FileI      : The Gram matrix on input\n";
      std::cerr << "OutFormat  : Possible values:\n";
      std::cerr << "             GAP   : the reduced Gram matrix and the "
                   "transformation matrix (GAP readable)\n";
      std::cerr << "             CPP_G : just the reduced Gram matrix\n";
      std::cerr << "             CPP_P : just the transformation matrix\n";
      std::cerr << "             Default: GAP\n";
      std::cerr << "FileO      : The output file\n";
      std::cerr << "             stdout: write to std::cout\n";
      std::cerr << "             stderr: write to std::cerr\n";
      std::cerr << "             Other filename get written to files\n";
      std::cerr << "             Default: stderr\n";
      throw TerminalException{1};
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
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileI, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(FileI, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(FileI, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return process<T, Tint>(FileI, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find a matching type for arith\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Possibilities: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_lll\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_lll\n";
    exit(e.eVal);
  }
  runtime(time);
}
