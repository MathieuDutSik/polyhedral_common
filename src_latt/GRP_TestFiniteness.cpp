// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "FiniteMatrixGroupTest.h"
// clang-format on

template <typename T,
          typename Tint = typename underlying_z_ring<T>::ring_type>
void process(std::string const& FileListMat, std::string const& OutFormat, std::string const& OutFile) {
  std::vector<MyMatrix<Tint>> ListMat = ReadListMatrixFile<Tint>(FileListMat);

  bool is_finite = test_finiteness_group<T,Tint>(ListMat, std::cerr);
  auto f_print=[&](std::ostream& os) -> void {
    if (OutFormat == "Raw") {
      os << "is_finite=" << is_finite << "\n";
      return;
    }
    if (OutFormat == "GAP") {
      os << "return rec(is_finite:=" << GAP_logical(is_finite) << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_TestFiniteness [arith] [ListMat] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "GRP_TestFiniteness [arith] [ListMat]\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileListMat = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f=[&]() -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return process<T>(FileListMat, OutFormat, OutFile);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileListMat, OutFormat, OutFile);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileListMat, OutFormat, OutFile);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(FileListMat, OutFormat, OutFile);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
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
    f();
    std::cerr << "Normal termination of GRP_TestFiniteness\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_TestFiniteness\n";
    exit(e.eVal);
  }
  runtime(time);
}
