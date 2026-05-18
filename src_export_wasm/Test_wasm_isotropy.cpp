// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_isotropy.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Isotropic.h"
#include <iostream>

using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // The hyperbolic plane Q(x,y) = 2xy. Its Gram matrix [[0,1],[1,0]] is
  // non-degenerate and obviously isotropic (e_1 is a zero of Q).
  MyMatrix<cpp_rational> Q(2, 2);
  Q << 0, 1,
       1, 0;
  std::optional<MyVector<cpp_rational>> opt = FindIsotropic(Q, std::cerr);
  if (!opt) {
    std::cerr << "Test_wasm_isotropy: FindIsotropic returned none on [[0,1],[1,0]]\n";
    return 1;
  }
  std::cerr << "Test_wasm_isotropy: OK\n";
  return 0;
}
