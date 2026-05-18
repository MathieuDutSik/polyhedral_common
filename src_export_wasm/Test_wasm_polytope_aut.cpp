// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for polytope automorphism (src_group / PolytopeEquiStab).
// Mirrors src_group/GRP_LinPolytope_Automorphism.cpp.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // 3-simplex (4 vertices in homogeneous dim 4). Its linear automorphism
  // group acts as the full symmetric group S_4 (size 24).
  MyMatrix<cpp_rational> EXT(4, 4);
  EXT << 1, 0, 0, 0,
         1, 1, 0, 0,
         1, 0, 1, 0,
         1, 0, 0, 1;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, cpp_int>;
  Tgroup GRP = LinPolytope_Automorphism<cpp_rational, Tgroup>(EXT, std::cerr);
  cpp_int expected(24);
  if (GRP.size() != expected) {
    std::cerr << "Test_wasm_polytope_aut: |Aut(3-simplex)|=" << GRP.size()
              << " expected 24\n";
    return 1;
  }
  std::cerr << "Test_wasm_polytope_aut: OK\n";
  return 0;
}
