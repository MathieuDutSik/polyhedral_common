// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &QFile, std::string const &PlaneFile,
             std::string const &choice, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = Tint;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Qmat = ReadMatrixFile<T>(QFile);
  MyMatrix<Tint> Plane = ReadMatrixFile<Tint>(PlaneFile);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  auto f_get = [&]() -> size_t {
    if (choice == "plane") {
      return comb.INDEF_FORM_Invariant_IsotropicKplane(Qmat, Plane);
    }
    if (choice == "flag") {
      return comb.INDEF_FORM_Invariant_IsotropicKflag(Qmat, Plane);
    }
    std::cerr << "No correct choice. choice=" << choice << "\n";
    throw TerminalException{1};
  };
  size_t e_inv = f_get();
  if (OutFormat == "CPP") {
    os_out << e_inv << "\n";
    return;
  }
  if (OutFormat == "PYTHON") {
    os_out << e_inv << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return " << e_inv << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_AutomorphismGroup [arith] [QFile] [PlaneFile] "
                   "[choice]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_AutomorphismGroup [arith] [QFile] [PlaneFile] "
                   "[choice] [OutFormat] "
                   "[OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string QFile = argv[2];
    std::string PlaneFile = argv[3];
    std::string choice = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(QFile, PlaneFile, choice, OutFormat, os);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(QFile, PlaneFile, choice, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(QFile, PlaneFile, choice, OutFormat, os);
      }
      if (arith == "safe") {
        using T = Rational<SafeInt64>;
        using Tint = SafeInt64;
        return process<T, Tint>(QFile, PlaneFile, choice, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of INDEF_FORM_InvariantIsotropicPlane\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_InvariantIsotropicPlane, runtime=" << time
              << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
