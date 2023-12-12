// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T>
void compute_rankin_k_min_kernel(std::string const &eFile) {
  using Tint = int64_t;
  std::ifstream is(eFile);
  MyMatrix<T> A = ReadMatrix<T>(is);
  std::cerr << "A=\n";
  WriteMatrix(std::cerr, A);
  int k;
  T tol;
  is >> k;
  is >> tol;
  std::cerr << "k=" << k << " tol=" << tol << "\n";
  ResultKRankinMin<T, Tint> result = Rankin_k_minimum(A, k, tol);
  std::cerr << "min=" << result.min << "\n";
  std::cerr << "|l_spaces|=" << result.l_spaces.size() << "\n";
}

void compute_k_level(std::string const &arithmetic,
                     std::string const &eFile) {
  if (arithmetic == "rational") {
    using T = mpq_class;
    return compute_rankin_k_min_kernel<T>(eFile);
  }
  if (arithmetic == "double") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return compute_rankin_k_min_kernel<T>(eFile);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "Rankin_k_min [arithmetic] [inputMat]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile = argv[2];
    compute_determinant(arithmetic, eFile);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
