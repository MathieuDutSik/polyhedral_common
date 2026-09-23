// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryCommon.h"
#include "Positivity.h"
// clang-format on

bool ParseBoolean(std::string const &strI) {
  if (strI == "T") {
    return true;
  }
  if (strI == "F") {
    return false;
  }
  if (strI == "true") {
    return true;
  }
  if (strI == "false") {
    return false;
  }
  std::cerr << "ParseBoolean error: possible input, T, F, true, false\n";
  throw TerminalException{1};
}

template <typename T,
          typename Tint = typename underlying_z_ring<T>::ring_type>
void process(std::string const &FileI, std::string const &strCritNorm,
             std::string const &strStrictIneq, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  T CritNorm = ParseScalar<T>(strCritNorm);
  bool StrictIneq = ParseBoolean(strStrictIneq);
  //
  MyVector<Tint> V =
      GetIntegralVector_allmeth<T, Tint>(M, CritNorm, StrictIneq, std::cerr);
  if (OutFormat == "GAP") {
    os_out << "return " << StringVectorGAP(V) << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "LATT_FindPositiveVector [arith] [FileI] [CritNorm] "
                   "[StrictIneq] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr
          << "LATT_FindPositiveVector [arith] [FileI] [CritNorm] [StrictIneq]\n";
      std::cerr << "\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Possibilities for arith: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Possibilities for arith: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string strCritNorm = argv[3];
    std::string strStrictIneq = argv[4];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileO = argv[6];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return process<T>(FileI, strCritNorm, strStrictIneq, OutFormat,
                                os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileI, strCritNorm, strStrictIneq, OutFormat,
                                os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, strCritNorm, strStrictIneq, OutFormat,
                                os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(FileI, strCritNorm, strStrictIneq, OutFormat,
                                os);
      }
#endif
      std::cerr << "Failed to find matching type for arith\n";
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
    std::cerr << "Normal termination of LATT_FindIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
