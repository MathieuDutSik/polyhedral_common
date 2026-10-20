// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "rational.h"
#include "Group.h"
#include "Permutation.h"
#include "ApproxAutoEquiv.h"
#include "Arithmetic_dispatch_types.h"
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
  FILE_PrintStderrStdoutFile(OutFile, prt);
}


int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string eFile1;
    std::string eFile2;
    std::string strTol;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T, typename Tint>() -> void {
      compute_approx_equivalence_kernel<T, Tint>(eFile1, eFile2, strTol,
                                                 OutFormat, OutFile);
    };
    add_rational_integer_arithmetics_by_library(dispatch, f_arith);
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ApproxEquivalence [arithmetic] [GramMat1] [GramMat2] "
                   "[tol] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ApproxEquivalence [arithmetic] [GramMat1] [GramMat2] "
                   "[tol]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arithmetic -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arithmetic = argv[1];
    eFile1 = argv[2];
    eFile2 = argv[3];
    strTol = argv[4];
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    dispatch.run(arithmetic);
    std::cerr << "Normal termination of LATT_ApproxEquivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ApproxEquivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}
