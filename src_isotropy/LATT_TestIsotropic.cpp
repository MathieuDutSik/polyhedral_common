// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheorySafeInt.h"
#include "Isotropic.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Q = ReadMatrixFile<T>(FileI);
  //
  bool test = is_isotropic(Q, std::cerr);
  if (OutFormat == "GAP") {
    if (test) {
      os_out << "return rec(has_isotropic:=true);\n";
    } else {
      os_out << "return rec(has_isotropic:=false);\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr
          << "LATT_TestIsotropic arithmetic [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_TestIsotropic arithmetic [FileI]\n";
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
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(FileI, OutFormat, os);
      }

      /*
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileI, OutFormat, os);
      }
      */

      /*
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, OutFormat, os);
      }
      */

      /*
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileI, OutFormat, os);
      }
      */
      std::cerr << "Failed to find matching type for arith. Possibility: rational\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, f);
    std::cerr << "Normal termination of LATT_TestIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_TestIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
