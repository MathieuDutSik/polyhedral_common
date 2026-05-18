// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_latt: shortest-vector enumeration on Z^n.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Shvec_exact.h"
#include "InvariantVectorFamily.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // Standard lattice Z^3, shortest non-zero vectors have norm 1
  // and there are 6 of them (+/- e_i).
  MyMatrix<cpp_rational> Gram = IdentityMat<cpp_rational>(3);
  T_shvec_info<cpp_rational, cpp_int> info =
      T_computeShvec_minimum<cpp_rational, cpp_int>(Gram, std::cerr);
  if (info.short_vectors.size() != 6) {
    std::cerr << "Test_wasm_latt: |short|=" << info.short_vectors.size()
              << " expected 6\n";
    return 1;
  }
  std::cerr << "Test_wasm_latt: OK\n";
  return 0;
}
