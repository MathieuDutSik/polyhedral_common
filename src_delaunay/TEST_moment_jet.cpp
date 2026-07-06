// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Compare the jet-based moment derivative with the 2n+1-sample interpolation.
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationDeformation.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "TEST_moment_jet [Gram.file]\n";
    return -1;
  }
  MyMatrix<T> Q = ReadMatrixFile<T>(argv[1]);
  int n = Q.rows();
  // A rank-one direction B = v v^T with v = e_0 - e_1.
  MyVector<T> v = ZeroVector<T>(n);
  v(0) = 1;
  v(1) = -1;
  MyMatrix<T> B = v * v.transpose();

  std::cerr << "=== computing interpolation moment derivative ===\n";
  MyMatrix<T> DM_interp =
      compute_moment_derivative<T, Tint, Tgroup>(Q, B, std::cerr);
  std::cerr << "DM_interp =\n" << DM_interp << "\n";
  std::cerr << "=== computing jet moment derivative ===\n";
  MyMatrix<T> DM_jet =
      compute_moment_derivative_jet<T, Tint, Tgroup>(Q, B, std::cerr);
  std::cerr << "DM_jet    =\n" << DM_jet << "\n";
  bool match = (DM_interp == DM_jet);
  std::cerr << (match ? "*** MATCH ***" : "*** MISMATCH ***") << "\n";
  return match ? 0 : 1;
}
