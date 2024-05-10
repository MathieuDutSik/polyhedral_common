// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_GRAVERBASIS_H_
#define SRC_LATT_GRAVERBASIS_H_

// clang-format off
#include "COMB_Combinatorics_buildset.h"
#include "COMB_Combinatorics_elem.h"
#include <vector>
// clang-format on

// A Graver basis is a set of vectors that can be used to
// optimize solution to an optimization problem:
// * It does not have to be optimal in any way.
// * It should be fast to compute.
// * It should be good in general, that is be "good" in some
// way and make things faster.
//
// This applies to many different context and the
// applicability is broad.

template <typename T, typename Tint>
MyMatrix<Tint> GetGraverBasis(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
  MyMatrix<Tint> M = ZeroMatrix<Tint>(2 * n, n);
  for (int i = 0; i < n; i++) {
    M(i, i) = 1;
    M(i + n, i) = -1;
  }
  return M;
}

template <typename Tint>
std::vector<MyVector<Tint>> GetGraverKbasis(int const &n, int const &k) {
  std::vector<MyVector<Tint>> ListVect;
  for (int u = 1; u <= k; u++) {
    using Tidx = SetCppIterator::Tidx;
    SetCppIterator cont(n, u);
    MyMatrix<int> LVal = BuildSet(k, 2);
    int npow = LVal.rows();
    for (auto &vect : cont) {
      MyVector<Tint> V = ZeroVector<Tint>(n);
      for (int ipow = 0; ipow < npow; ipow++) {
        for (int pos = 0; pos < k; pos++) {
          Tidx idx = vect[pos];
          Tint sign = 2 * LVal(ipow, pos) - 1;
          V(idx) = sign;
        }
        ListVect.push_back(V);
      }
    }
  }
  return ListVect;
}

// clang-format off
#endif  // SRC_LATT_GRAVERBASIS_H_
// clang-format on
