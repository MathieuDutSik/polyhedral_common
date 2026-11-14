// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_CLASSICALGROUPORDERS_H_
#define SRC_GENERA_CLASSICALGROUPORDERS_H_

// clang-format off
#include "GeneraArithmetic.h"
#include <vector>
// clang-format on

/*
  The orders of the finite classical groups (symplectic, orthogonal,
  unitary, general linear) over a residue field of q elements, in the
  normalizations of the local factors of Hecke.jl.
 */

template <typename T>
inline T HeckeGroupOrderSp(int m, T const &q) {
  T res(1);
  for (int j = 2; j <= m; j += 2) {
    T pw(1);
    for (int i = 0; i < j; i++) {
      pw /= q;
    }
    res *= (1 - pw);
  }
  return res;
}

template <typename T>
inline T HeckeGroupOrderO(int m, int eps, T const &q) {
  auto qinv_pow = [&](int j) -> T {
    T pw(1);
    for (int i = 0; i < j; i++) {
      pw /= q;
    }
    return pw;
  };
  T res(2);
  if (m % 2 == 0) {
    int k = m / 2;
    res *= (1 - T(eps) * qinv_pow(k));
    for (int j = 2; j <= m - 2; j += 2) {
      res *= (1 - qinv_pow(j));
    }
  } else {
    for (int j = 2; j <= m; j += 2) {
      res *= (1 - qinv_pow(j));
    }
  }
  return res;
}

template <typename T>
inline T HermitianGroupOrderNormalized(int type, int m, T const &qK) {
  // GL (split) or U (inert): prod_{j=1}^{m} (1 - (+-qK)^{-j})
  T res(1);
  T base = qK;
  if (type != 0) {
    base = -base;
  }
  T pw(1);
  for (int j = 1; j <= m; j++) {
    pw /= base;
    res *= (1 - pw);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_CLASSICALGROUPORDERS_H_
// clang-format on
