// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_JET_H_
#define SRC_DELAUNAY_JET_H_

// clang-format off
#include "MAT_Matrix.h"
#include <cmath>
#include <utility>
#include <vector>
// clang-format on

// A jet (truncated Taylor expansion) c0 + c1 t + ... + c_order t^order, stored
// by its coefficients. The arithmetic below truncates every result back to the
// jet order, so a chain of operations stays at a fixed order. This replaces the
// order-specific hand-written formulas (the old std::array<T,3> code).

template <typename T> struct Jet {
  int order;
  std::vector<T> coeffs; // size order + 1
};

// The truncation to a jet of the given order of the polynomial with coefficient
// vector C (c0 + c1 t + ...). Missing coefficients are zero.
template <typename T> Jet<T> jet_from_poly(MyVector<T> const &C, int order) {
  std::vector<T> coeffs(order + 1, T(0));
  for (int i = 0; i <= order && i < C.size(); i++) {
    coeffs[i] = C(i);
  }
  return {order, std::move(coeffs)};
}

// Product of two jets of the same order, truncated to that order.
template <typename T> Jet<T> operator*(Jet<T> const &a, Jet<T> const &b) {
  int order = a.order;
  std::vector<T> coeffs(order + 1, T(0));
  for (int i = 0; i <= order; i++) {
    for (int j = 0; i + j <= order; j++) {
      coeffs[i + j] += a.coeffs[i] * b.coeffs[j];
    }
  }
  return {order, std::move(coeffs)};
}

// Scalar multiple of a jet.
template <typename T> Jet<T> jet_scalar_mult(T const &scal, Jet<T> const &a) {
  std::vector<T> coeffs(a.order + 1);
  for (int i = 0; i <= a.order; i++) {
    coeffs[i] = scal * a.coeffs[i];
  }
  return {a.order, std::move(coeffs)};
}

// Inverse 1/f of a jet f with non-zero constant term, from the recurrence
// (f * f^{-1} = 1): b_0 = 1/c_0 and b_k = -b_0 sum_{j=1}^{k} c_j b_{k-j}.
template <typename T> Jet<T> jet_inverse(Jet<T> const &f) {
  int order = f.order;
  std::vector<T> b(order + 1, T(0));
  b[0] = T(1) / f.coeffs[0];
  for (int k = 1; k <= order; k++) {
    T s(0);
    for (int j = 1; j <= k; j++) {
      s += f.coeffs[j] * b[k - j];
    }
    b[k] = -b[0] * s;
  }
  return {order, std::move(b)};
}

// f^p of a jet f with non-zero constant term and a scalar exponent p, from the
// J.C.P. Miller recurrence: b_0 = c_0^p and
// b_k = (1/(k c_0)) sum_{j=1}^{k} (p j - (k - j)) c_j b_{k-j}.
// The leading coefficient uses std::pow(c_0, p), so this is meant for a
// floating-point T (it is used for the irrational power det(t)^(-1/n)).
template <typename T> Jet<T> jet_pow(Jet<T> const &f, T const &p) {
  int order = f.order;
  std::vector<T> b(order + 1, T(0));
  T const &c0 = f.coeffs[0];
  b[0] = std::pow(c0, p);
  for (int k = 1; k <= order; k++) {
    T s(0);
    for (int j = 1; j <= k; j++) {
      s += (p * T(j) - T(k - j)) * f.coeffs[j] * b[k - j];
    }
    b[k] = s / (T(k) * c0);
  }
  return {order, std::move(b)};
}

// clang-format off
#endif  // SRC_DELAUNAY_JET_H_
// clang-format on
