// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "edgewalk.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &MatFile, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> LorMat = ReadMatrixFile<T>(MatFile);
  //
  ResultEdgewalk<T, Tint> re =
      StandardEdgewalkAnalysis<T, Tint, Tgroup>(LorMat, std::cerr);
  bool ComputeAllSimpleRoots = true;
  PrintResultEdgewalk(LorMat, re, os_out, OutFormat, ComputeAllSimpleRoots,
                      std::cerr);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile]\n";
      std::cerr << "     or\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile] [OutFormat] "
                   "[OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(MatFile, OutFormat, os_out);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(MatFile, OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(MatFile, OutFormat, os_out);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return process<T, Tint>(MatFile, OutFormat, os_out);
      }
#endif
      std::cerr << "Failed to find matching entry for arith=" << arith << "\n";
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
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of LORENTZ_FundDomain_AllcockEdgewalk\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_FundDomain_AllcockEdgewalk\n";
    exit(e.eVal);
  }
  runtime(time);
}
