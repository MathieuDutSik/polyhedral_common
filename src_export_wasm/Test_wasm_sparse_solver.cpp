// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_sparse_solver.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "MAT_SparseMatrix.h"
#include "GampMatlab.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // Compilation-only smoke: build a trivial sparse identity-style matrix.
  MySparseMatrix<cpp_rational> A(2, 2);
  A.coeffRef(0, 0) = 1;
  A.coeffRef(1, 1) = 1;
  std::cerr << "Test_wasm_sparse_solver: nnz=" << A.nonZeros() << "\n";
  return 0;
}
