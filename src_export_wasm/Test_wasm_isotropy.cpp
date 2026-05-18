// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_isotropy.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Isotropic.h"
#include "Indefinite_LLL.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // 2x2 hyperbolic form [[0,1],[1,0]] is obviously isotropic (e.g. e_1).
  MyMatrix<cpp_int> Q(2, 2);
  Q << 0, 1,
       1, 0;
  std::optional<MyVector<cpp_int>> opt =
      INDEFINITE_GetShortIsotropic<cpp_int, cpp_rational>(Q, std::cerr);
  if (!opt) {
    std::cerr << "Test_wasm_isotropy: no isotropic vector found\n";
    return 1;
  }
  std::cerr << "Test_wasm_isotropy: OK\n";
  return 0;
}
