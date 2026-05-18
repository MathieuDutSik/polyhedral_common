// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_copos: copositivity check on a small matrix.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Copositivity.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // The identity matrix has positive diagonal and non-negative off-diagonal,
  // so the cheap sufficient test in TestCopositivityByPositivityCoeff
  // returns true.
  MyMatrix<cpp_rational> M = IdentityMat<cpp_rational>(3);
  bool is_copo = TestCopositivityByPositivityCoeff(M);
  if (!is_copo) {
    std::cerr << "Test_wasm_copos: I_3 should be copositive\n";
    return 1;
  }
  // And a clearly non-copositive matrix: negative diagonal entry rules it out.
  MyMatrix<cpp_rational> Mneg = IdentityMat<cpp_rational>(3);
  Mneg(0, 0) = -1;
  if (TestCopositivityByPositivityCoeff(Mneg)) {
    std::cerr << "Test_wasm_copos: matrix with -1 on diag flagged copositive\n";
    return 1;
  }
  std::cerr << "Test_wasm_copos: OK\n";
  return 0;
}
