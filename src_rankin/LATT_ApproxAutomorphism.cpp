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
  FILE_PrintStderrStdoutFile(OutFile, prt);
}


int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string eFile;
    std::string strTol;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T, typename Tint>() -> void {
      compute_approx_automorphism_kernel<T, Tint>(eFile, strTol, OutFormat,
                                                  OutFile);
    };
    add_rational_integer_arithmetics_by_library(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ApproxAutomorphism [arithmetic] [GramMat] [tol] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ApproxAutomorphism [arithmetic] [GramMat] [tol]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arithmetic -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arithmetic = argv[1];
    eFile = argv[2];
    strTol = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    dispatch.run(arithmetic);
    std::cerr << "Normal termination of LATT_ApproxAutomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ApproxAutomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
