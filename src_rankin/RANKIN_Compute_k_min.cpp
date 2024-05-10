// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "rational.h"
#include "Enumeration_k_space.h"
// clang-format on

template <typename T>
void compute_rankin_k_min_kernel(int const &k, std::string const &eFile,
                                 std::string const &strTol) {
  using Tint = int64_t;
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

void compute_k_min(std::string const &arithmetic, int const &k,
                   std::string const &eFile, std::string const &strTol) {
  if (arithmetic == "rational") {
    using T = mpq_class;
    return compute_rankin_k_min_kernel<T>(k, eFile, strTol);
  }
  if (arithmetic == "double") {
    using T = double;
    return compute_rankin_k_min_kernel<T>(k, eFile, strTol);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "Rankin_k_min [arithmetic] [k] [inputMat] [tol]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string strK = argv[2];
    std::string eFile = argv[3];
    std::string strTol = argv[4];
    int k = ParseScalar<int>(strK);
    compute_k_min(arithmetic, k, eFile, strTol);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
