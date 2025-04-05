// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheorySafeInt.h"
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

template <typename T, typename Tint>
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
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "LATT_FindPositiveVector arith [FileI] [CritNorm] "
                   "[StrictIneq] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr
          << "LATT_FindPositiveVector arith [FileI] [CritNorm] [StrictIneq]\n";
      std::cerr << "\n";
      std::cerr << "Possibilities for arith: gmp\n";
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
        using Tint = mpz_class;
        return process<T, Tint>(FileI, strCritNorm, strStrictIneq, OutFormat,
                                os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, f);
    std::cerr << "Normal termination of LATT_FindIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
