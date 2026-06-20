// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "rational.h"
#include "Group.h"
#include "Permutation.h"
#include "ApproxAutoEquiv.h"
// clang-format on

template <typename T, typename Tint>
void compute_approx_equivalence_kernel(std::string const &eFile1,
                                       std::string const &eFile2,
                                       std::string const &strTol,
                                       std::string const &OutFormat,
                                       std::string const &OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ifstream is1(eFile1);
  MyMatrix<T> eG1 = ReadMatrix<T>(is1);
  std::ifstream is2(eFile2);
  MyMatrix<T> eG2 = ReadMatrix<T>(is2);
  T tol = ParseScalar<T>(strTol);
  std::cerr << "LATT_ApproxEquivalence: dim1=" << eG1.rows()
            << " dim2=" << eG2.rows() << " tol=" << tol << "\n";
  std::optional<MyMatrix<Tint>> opt =
      ApproximateEquivalence<T, Tint, Tgroup>(eG1, eG2, tol, std::cerr);
  auto prt = [&](std::ostream &os) -> void {
    if (OutFormat == "GAP") {
      if (opt) {
        os << "return ";
        WriteMatrixGAP(os, *opt);
        os << ";\n";
      } else {
        os << "return false;\n";
      }
      return;
    }
    if (OutFormat == "Oscar") {
      if (opt) {
        WriteMatrix(os, *opt);
      } else {
        os << "0 0\n";
      }
      return;
    }
    std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, prt);
}

void compute_approx_equivalence(std::string const &arithmetic,
                                std::string const &eFile1,
                                std::string const &eFile2,
                                std::string const &strTol,
                                std::string const &OutFormat,
                                std::string const &OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return compute_approx_equivalence_kernel<T, Tint>(eFile1, eFile2, strTol,
                                                      OutFormat, OutFile);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return compute_approx_equivalence_kernel<T, Tint>(eFile1, eFile2, strTol,
                                                      OutFormat, OutFile);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return compute_approx_equivalence_kernel<T, Tint>(eFile1, eFile2, strTol,
                                                      OutFormat, OutFile);
  }
  if (arithmetic == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return compute_approx_equivalence_kernel<T, Tint>(eFile1, eFile2, strTol,
                                                      OutFormat, OutFile);
  }
  std::cerr << "Failed to find a matching arithmetic=" << arithmetic << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ApproxEquivalence [arithmetic] [GramMat1] [GramMat2] "
                   "[tol] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ApproxEquivalence [arithmetic] [GramMat1] [GramMat2] "
                   "[tol]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile1 = argv[2];
    std::string eFile2 = argv[3];
    std::string strTol = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    compute_approx_equivalence(arithmetic, eFile1, eFile2, strTol, OutFormat,
                               OutFile);
    std::cerr << "Normal termination of LATT_ApproxEquivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ApproxEquivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}
