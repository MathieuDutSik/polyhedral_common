// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_poly. Pulls in the four public CDD entry points
// (dual description, linear programming, redundancy, skeletons) under
// WASM_PLATFORM so the external-program branches are excluded.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "POLY_cddlib.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_LinearProgramming.h"
#include "POLY_RedundancyElimination.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // 3-simplex V-representation, homogeneous (one + 3 vertex columns).
  MyMatrix<cpp_rational> EXT(4, 4);
  EXT << 1, 0, 0, 0,
         1, 1, 0, 0,
         1, 0, 1, 0,
         1, 0, 0, 1;
  // cdd dual description: expect 4 facets for a 3-simplex.
  MyMatrix<cpp_rational> FAC = cdd::DualDescription(EXT, std::cerr);
  if (FAC.rows() != 4) {
    std::cerr << "Test_wasm_poly: |FAC|=" << FAC.rows() << " expected 4\n";
    return 1;
  }
  std::cerr << "Test_wasm_poly: dual_desc OK\n";
  return 0;
}
