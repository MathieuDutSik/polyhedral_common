// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "rational.h"
#include "Enumeration_k_space.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint>
void compute_rankin_k_min_kernel(int const &k, std::string const &eFile,
                                 std::string const &strTol) {
  std::ifstream is(eFile);
  MyMatrix<T> A = ReadMatrix<T>(is);
  std::cerr << "A=\n";
  WriteMatrix(std::cerr, A);
  T tol = ParseScalar<T>(strTol);
  std::cerr << "k=" << k << " tol=" << tol << "\n";
  ResultKRankinMin<T, Tint> result =
      Rankin_k_minimum<T, Tint>(A, k, tol, std::cerr);
  std::cerr << "min=" << result.min << "\n";
  std::cerr << "|l_spaces|=" << result.l_spaces.size() << "\n";
}


int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  try {
    std::string strK;
    std::string eFile;
    std::string strTol;
    int k = 0;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T, typename Tint>() -> void {
      compute_rankin_k_min_kernel<T, Tint>(k, eFile, strTol);
    };
    add_rational_integer_arithmetics_by_library(dispatch, f_arith);
    dispatch.add("double", "floating point double and int64_t",
                 [&]() -> void {
                   f_arith.template operator()<double, int64_t>();
                 });
    if (argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "RANKIN_Compute_k_min [arithmetic] [k] [inputMat] [tol]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arithmetic -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arithmetic = argv[1];
    strK = argv[2];
    eFile = argv[3];
    strTol = argv[4];
    k = ParseScalar<int>(strK);
    dispatch.run(arithmetic);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
