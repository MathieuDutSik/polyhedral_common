// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
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
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "LATT_TestIsotropic [arith] [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_TestIsotropic [arith] [FileI]\n";
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
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(FileI, OutFormat, os);
      }
#endif

      /*
      */

      /*
#ifdef ENABLE_BOOST_TYPES
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, OutFormat, os);
      }
#endif
      */

      /*
#ifdef ENABLE_BOOST_TYPES
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileI, OutFormat, os);
      }
#endif
      */
      std::cerr << "Failed to find matching type for arith.\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Possibilities: rational, flint\n";
#else
      std::cerr << "Possibilities: rational (build with "
                << "ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_TestIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_TestIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
