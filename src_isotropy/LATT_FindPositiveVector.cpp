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
#include "Arithmetic_dispatch_types.h"
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

template <typename T>
void process(std::string const &FileI, std::string const &strCritNorm,
             std::string const &strStrictIneq, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
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
    std::string FileI;
    std::string strCritNorm;
    std::string strStrictIneq;
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(FileI, strCritNorm, strStrictIneq, OutFormat, *os_ptr);
    });
    if (argc != 5 && argc != 7) {
      std::cerr << "LATT_FindPositiveVector [arith] [FileI] [CritNorm] "
                   "[StrictIneq] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr
          << "LATT_FindPositiveVector [arith] [FileI] [CritNorm] [StrictIneq]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    FileI = argv[2];
    strCritNorm = argv[3];
    strStrictIneq = argv[4];
    if (argc == 7) {
      OutFormat = argv[5];
      FileO = argv[6];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_FindIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
