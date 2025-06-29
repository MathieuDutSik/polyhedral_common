// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ISOTROPY_NORMS_H_
#define SRC_ISOTROPY_NORMS_H_

#include "MAT_Matrix.h"

// Compute an upper bound on the determinant of maximal minor
// The Hadamard bound is
// We have det(A)^2 <= Pi_{i=1}^n (sum_j x_{ij}^2)
template <typename T> T sqr_estimate_maximal_determinant(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::vector<T> ListSqr;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T sum = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      sum += M(iRow, iCol) * M(iRow, iCol);
    ListSqr.push_back(sum);
  }
  std::sort(ListSqr.begin(), ListSqr.end());
  T eProd = 1;
  for (int iRow = nbRow - nbCol; iRow < nbRow; iRow++) {
    eProd *= ListSqr[iRow];
  }
  return eProd;
}

template <typename T> T sqr_estimate_facet_coefficients(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  T max_coeff = 0;
  MyVector<T> SqrCoeffFacet(nbCol);
  for (int jCol = 0; jCol < nbCol; jCol++) {
    MyMatrix<T> Mret(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != jCol) {
          Mret(iRow, pos) = M(iRow, iCol);
          pos++;
        }
      }
    }
    T est = sqr_estimate_maximal_determinant(Mret);
    if (est > max_coeff) {
      max_coeff = est;
    }
  }
  return max_coeff;
}

// clang-format off
#endif  // SRC_ISOTROPY_NORMS_H_
// clang-format on
