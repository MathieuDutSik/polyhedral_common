// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for the perfect-form complex enumeration.
// Models src_perfect/PERF_SerialPerfectComputation.cpp -- builds a small
// canonical T-space, configures the perfect-complex options, and runs the
// serial enumeration (EnumerateAndStore_Serial via DataPerfectTspaceFunc).
// This exercises the recursive-adjacency code path that previously was not
// WASM-clean.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "perfect_complex.h"
#include "Tspace_Generation.h"
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

  // Canonical sym n x n T-space, n = 2: the perfect form is unique (A_2 root
  // lattice up to scale), so the enumeration terminates quickly.
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(2);
  PerfectComplexOptions pco;
  pco.compute_boundary = false;
  pco.compute_contracting_homotopy = false;
  pco.only_well_rounded = false;
  FullComplexEnumeration<T, Tint, Tgroup> fce =
      get_full_complex_enumeration_kernel<T, Tint, Tgroup>(LinSpa, pco,
                                                           std::cerr);
  std::cerr << "Test_wasm_perfect: |fce.l_topdims|=" << fce.l_topdims.size() << " OK\n";
  return 0;
}
