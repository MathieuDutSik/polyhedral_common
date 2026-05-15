// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_CANONICAL_H_
#define SRC_LATT_TSPACE_CANONICAL_H_

// clang-format off
#include "MAT_Matrix.h"
#include <vector>
// clang-format on

template<typename T>
std::vector<MyMatrix<T>> TSPACE_canonical_get_list_matrices(int const& n) {
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<T> eMat(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      ZeroAssignation(eMat);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListMat.push_back(eMat);
    }
  }
  return ListMat;
}

template<typename T>
MyVector<T> TSPACE_canonical_get_expression(MyMatrix<T> const& M) {
  int n = M.rows();
  int dim = n * (n + 1) / 2;
  MyVector<T> V(dim);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T val = M(i,j);
      V(pos) = val;
      pos += 1;
    }
  }
  return V;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_CANONICAL_H_
// clang-format on
