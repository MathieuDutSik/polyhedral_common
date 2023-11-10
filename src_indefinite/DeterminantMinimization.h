// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
#define SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_

#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include "WeightMatrix.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

template<typename T>
struct ResultDetMin {
  MyMatrix<T> P;
  MyMatrix<T> Mred;
};


/*
  We apply a numbr of ideas from the preprint
  "Quadratic equations in dimensions 4, 5 and more" (P1)
  Denis Simon
  ---
  The following notions are used
  * v_p(Q) to be the multiplicity of p as a prime factor of the determinant of Q.
  * d = dim Ker_{F_p}(Q).
  We have the basic result d <= v.
  ---
  Lemma 4 is clear.

  The reduction of the matrix Q prior to can be done via NullspaceMatMod function.
  So, for Lemma 5, we need to apply that function another time in order to reduce tilde(Q).
  Lemma 5 is clear.

  
 */


template<typename T>
ResultDetMin<T> DeterminantMinimization(MyMatrix<T> const& M) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  using Tring = typename underlying_ring<T>::ring_type;
  if (!IsIntegralMatrix(M)) {
    std::cerr << "The matrix M should be integral\n";
    throw TerminalException{1};
  }
  int n_row = M.rows();
  T det = DeterminantMat(M);
  T det_abs = T_abs(det);
  Tring det_ai = UniversalScalarConversion<Tring,T>(det_abs);
  std::map<Tring, size_t> map = FactorsIntMap(det_ai);
  MyMatrix<T> Mwork = M;
  // Apply Lemma 4
  std::vector<T> cases_lemma4;
  for (auto & kv : map) {
    int q = kv.second / n_row;
    for (int u=0; u<q; u++) {
      cases_lemma4.push_back(kv.first);
    }
  }
  for (auto & p : cases_lemma4) {
    map[p] -= n_row;
    if (map[p] == 0) {
      map.erase(p);
    }
    Mwork = Mwork / p;
  }
  
  while(true) {
  }

}

// clang-format off
#endif  //  SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
// clang-format on
