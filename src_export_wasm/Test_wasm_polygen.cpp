// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_polygen.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "generalized_polytopes.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // 3-simplex with one inequality per facet (FAC) and one row per vertex (EXT),
  // homogeneous dim 4. The four facets are x_i >= 0 (i=1..3) and 1-x_1-x_2-x_3 >= 0.
  MyMatrix<cpp_rational> EXT(4, 4);
  EXT << 1, 0, 0, 0,
         1, 1, 0, 0,
         1, 0, 1, 0,
         1, 0, 0, 1;
  MyMatrix<cpp_rational> FAC(4, 4);
  FAC << 0,  1,  0,  0,    // x_1 >= 0
         0,  0,  1,  0,    // x_2 >= 0
         0,  0,  0,  1,    // x_3 >= 0
         1, -1, -1, -1;    // 1 - x_1 - x_2 - x_3 >= 0
  SinglePolytope<cpp_rational> sp = get_single_polytope(FAC, EXT);
  if (sp.facets.size() != 4) {
    std::cerr << "Test_wasm_polygen: |facets|=" << sp.facets.size()
              << " expected 4\n";
    return 1;
  }
  std::cerr << "Test_wasm_polygen: OK\n";
  return 0;
}
