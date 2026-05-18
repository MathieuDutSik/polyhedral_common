// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_short.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "SHORT_Realizability.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // {+/- e_1, +/- e_2} is the shortest-vector configuration of Z^2 with
  // the standard inner product, so realizability should hold.
  MyMatrix<cpp_int> SHV(4, 2);
  SHV << 1, 0,
        -1, 0,
         0, 1,
         0, -1;
  std::cerr << "Test_wasm_short: SHV rows=" << SHV.rows() << "\n";
  return 0;
}
