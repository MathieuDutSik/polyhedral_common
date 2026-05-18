// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_copos: copositivity check on a small matrix.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Copositivity.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // The identity matrix is copositive (x^T I x = ||x||^2 >= 0 always).
  MyMatrix<cpp_rational> M = IdentityMat<cpp_rational>(3);
  RequestCopositivity<cpp_rational> request{0, true};
  CopositivityInfoReduction<cpp_rational> info =
      EnumerateCopositiveShortVector<cpp_rational, cpp_rational>(
          M, request, std::cerr);
  std::cerr << "Test_wasm_copos: |reductions|=" << info.nbReduction << "\n";
  return 0;
}
