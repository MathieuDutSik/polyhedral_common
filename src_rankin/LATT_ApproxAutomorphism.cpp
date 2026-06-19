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
void compute_approx_automorphism_kernel(std::string const &eFile,
                                        std::string const &strTol,
                                        std::string const &OutFormat,
                                        std::string const &OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ifstream is(eFile);
  MyMatrix<T> eG = ReadMatrix<T>(is);
  T tol = ParseScalar<T>(strTol);
  std::cerr << "LATT_ApproxAutomorphism: dim=" << eG.rows() << " tol=" << tol
            << "\n";
  std::optional<std::vector<MyMatrix<Tint>>> opt =
      ApproximateAutomorphismGroup<T, Tint, Tgroup>(eG, tol, std::cerr);
  auto prt = [&](std::ostream &os) -> void {
    if (OutFormat == "GAP") {
      os << "return ";
      if (opt) {
        WriteListMatrixGAP(os, *opt);
      } else {
        os << "fail";
      }
      os << ";\n";
      return;
    }
    if (OutFormat == "Oscar") {
      if (opt) {
        std::vector<MyMatrix<Tint>> ListGen = *opt;
        os << ListGen.size() << "\n";
        for (auto &eMat : ListGen) {
          WriteMatrix(os, eMat);
        }
      } else {
        os << "-1\n";
      }
      return;
    }
    std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, prt);
}

void compute_approx_automorphism(std::string const &arithmetic,
                                 std::string const &eFile,
                                 std::string const &strTol,
                                 std::string const &OutFormat,
                                 std::string const &OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return compute_approx_automorphism_kernel<T, Tint>(eFile, strTol, OutFormat,
                                                       OutFile);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return compute_approx_automorphism_kernel<T, Tint>(eFile, strTol, OutFormat,
                                                       OutFile);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return compute_approx_automorphism_kernel<T, Tint>(eFile, strTol, OutFormat,
                                                       OutFile);
  }
  if (arithmetic == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return compute_approx_automorphism_kernel<T, Tint>(eFile, strTol, OutFormat,
                                                       OutFile);
  }
  std::cerr << "Failed to find a matching arithmetic=" << arithmetic << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ApproxAutomorphism [arithmetic] [GramMat] [tol] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ApproxAutomorphism [arithmetic] [GramMat] [tol]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile = argv[2];
    std::string strTol = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    compute_approx_automorphism(arithmetic, eFile, strTol, OutFormat, OutFile);
    std::cerr << "Normal termination of LATT_ApproxAutomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ApproxAutomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
