// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_polygen.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "generalized_polytopes.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // Compilation-only smoke test: instantiate one template from the header.
  // 4-vertex polytope (a tetrahedron) in homogeneous dim 4.
  MyMatrix<cpp_rational> EXT(4, 4);
  EXT << 1, 0, 0, 0,
         1, 1, 0, 0,
         1, 0, 1, 0,
         1, 0, 0, 1;
  std::vector<MyMatrix<cpp_rational>> ListEXT{EXT};
  PolyGen<cpp_rational> gp = ListMat_to_PolyGen(ListEXT);
  std::cerr << "Test_wasm_polygen: |gp|=" << gp.ListVect.rows() << "\n";
  return 0;
}
