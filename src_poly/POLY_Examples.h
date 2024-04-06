// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_EXAMPLES_H_
#define SRC_POLY_POLY_EXAMPLES_H_

// clang-format off
#include "MAT_Matrix.h"
// clang-format on

template <typename T> MyMatrix<T> CyclicPolytope(int n, int k) {
  int i, j, b;
  MyMatrix<T> TheEXT(n, k + 1);
  for (i = 1; i <= n; i++) {
    T a = 1;
    b = i + 1;
    for (j = 0; j <= k; j++) {
      TheEXT(i - 1, j) = a;
      a = a * b;
    }
  }
  return TheEXT;
}

// clang-format off
#endif  // SRC_POLY_POLY_EXAMPLES_H_
// clang-format on
