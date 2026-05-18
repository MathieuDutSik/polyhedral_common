// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_latt: shortest-vector enumeration on Z^n.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Shvec_exact.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // Z^3 with the standard inner product: minimum is 1, attained by the 6
  // signed unit vectors +/- e_i.
  MyMatrix<cpp_rational> Gram = IdentityMat<cpp_rational>(3);
  Tshortest<cpp_rational, cpp_int> result =
      T_ShortestVector<cpp_rational, cpp_int>(Gram, std::cerr);
  if (result.min != 1) {
    std::cerr << "Test_wasm_latt: min=" << result.min << " expected 1\n";
    return 1;
  }
  if (result.SHV.rows() != 6) {
    std::cerr << "Test_wasm_latt: |SHV|=" << result.SHV.rows()
              << " expected 6\n";
    return 1;
  }
  std::cerr << "Test_wasm_latt: OK\n";
  return 0;
}
