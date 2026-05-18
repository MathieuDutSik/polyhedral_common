// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for Gram-matrix arithmetic automorphism (src_latt).
// Mirrors src_latt/LATT_Automorphism.cpp.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include <iostream>
#include <vector>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  // Aut(I_2, Z) = signed permutations of two axes = D_4 of order 8.
  // ArithmeticAutomorphismGroupMultiple returns a set of integer-matrix
  // generators; verify we got some, and that each is unimodular and
  // preserves I_2 (G^T I G = I, i.e. G^T G = I).
  std::vector<MyMatrix<cpp_rational>> ListMat{IdentityMat<cpp_rational>(2)};
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, cpp_int>;
  std::vector<MyMatrix<cpp_int>> ListGen =
      ArithmeticAutomorphismGroupMultiple<cpp_rational, cpp_int, Tgroup>(
          ListMat, std::cerr);
  if (ListGen.empty()) {
    std::cerr << "Test_wasm_gram_aut: no generators returned\n";
    return 1;
  }
  for (auto const &G : ListGen) {
    cpp_int det = DeterminantMat(G);
    if (det != 1 && det != -1) {
      std::cerr << "Test_wasm_gram_aut: generator with det=" << det
                << " (not unimodular)\n";
      return 1;
    }
    MyMatrix<cpp_int> Gt = G.transpose();
    MyMatrix<cpp_int> prod = Gt * G;
    MyMatrix<cpp_int> Id = IdentityMat<cpp_int>(2);
    if (!TestEqualityMatrix(prod, Id)) {
      std::cerr << "Test_wasm_gram_aut: generator does not preserve I_2\n";
      return 1;
    }
  }
  std::cerr << "Test_wasm_gram_aut: |ListGen|=" << ListGen.size() << " OK\n";
  return 0;
}
