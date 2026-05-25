// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for the automorphism group of a Lorentzian lattice.
// Models src_lorentzian/LORENTZ_PERF_Automorphism.cpp -- calls
// LORENTZ_GetGeneratorsAutom on a small Lorentzian Gram matrix and
// verifies every returned generator preserves the form (G^T LorMat G ==
// LorMat). Goal: confirm lorentzian_perfect.h templates compile under
// WASM_PLATFORM.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "lorentzian_perfect.h"
#include "Permutation.h"
#include "Group.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  using T = cpp_rational;
  using Tint = cpp_int;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = cpp_int;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  // 2D Lorentzian form: diag(1, -1). Hyperbolic plane up to a sign change.
  MyMatrix<T> LorMat(2, 2);
  LorMat << 1,  0,
            0, -1;
  std::vector<MyMatrix<Tint>> l_gen =
      LORENTZ_GetGeneratorsAutom<T, Tint, Tgroup>(LorMat, std::cerr);
  for (auto const &eGen : l_gen) {
    MyMatrix<T> eGen_T = UniversalMatrixConversion<T, Tint>(eGen);
    MyMatrix<T> prod = eGen_T * LorMat * eGen_T.transpose();
    if (prod != LorMat) {
      std::cerr << "Test_wasm_lorentz_aut: generator does not preserve LorMat\n";
      return 1;
    }
  }
  std::cerr << "Test_wasm_lorentz_aut: |l_gen|=" << l_gen.size() << " OK\n";
  return 0;
}
