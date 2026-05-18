// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for src_group: build a tiny permutation group via permutalib.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;

int main() {
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, cpp_int>;
  // The transposition (1 2) acting on 3 points: it generates Z/2.
  std::vector<Tidx> v = {1, 0, 2};
  Telt p(std::move(v));
  Tgroup G({p}, Tidx(3));
  if (G.size() != cpp_int(2)) {
    std::cerr << "Test_wasm_group: |G|=" << G.size() << " expected 2\n";
    return 1;
  }
  std::cerr << "Test_wasm_group: OK\n";
  return 0;
}
