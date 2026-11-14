// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_RATIONALFUNCTIONS_H_
#define SRC_GENERA_RATIONALFUNCTIONS_H_

// clang-format off
#include "NumberTheory.h"
#include <iostream>
#include <vector>
// clang-format on

// Univariate polynomials and rational functions over T

template <typename T> struct UnivariatePolynomial {
  std::vector<T> coeff;
  void normalize() {
    while (!coeff.empty() && coeff.back() == 0) {
      coeff.pop_back();
    }
  }
  int degree() const { return static_cast<int>(coeff.size()) - 1; }
  bool is_zero() const { return coeff.empty(); }
};

template <typename T> UnivariatePolynomial<T> PolyConstant(T const &c) {
  UnivariatePolynomial<T> P{{c}};
  P.normalize();
  return P;
}

// c * X^d
template <typename T> UnivariatePolynomial<T> PolyMonomial(T const &c, int d) {
  UnivariatePolynomial<T> P;
  if (c == 0) {
    return P;
  }
  P.coeff.assign(d + 1, T(0));
  P.coeff[d] = c;
  return P;
}

template <typename T>
UnivariatePolynomial<T> operator+(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  size_t len = std::max(A.coeff.size(), B.coeff.size());
  C.coeff.assign(len, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    C.coeff[i] += A.coeff[i];
  }
  for (size_t i = 0; i < B.coeff.size(); i++) {
    C.coeff[i] += B.coeff[i];
  }
  C.normalize();
  return C;
}

template <typename T>
UnivariatePolynomial<T> operator-(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  size_t len = std::max(A.coeff.size(), B.coeff.size());
  C.coeff.assign(len, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    C.coeff[i] += A.coeff[i];
  }
  for (size_t i = 0; i < B.coeff.size(); i++) {
    C.coeff[i] -= B.coeff[i];
  }
  C.normalize();
  return C;
}

template <typename T>
UnivariatePolynomial<T> operator*(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  if (A.is_zero() || B.is_zero()) {
    return C;
  }
  C.coeff.assign(A.coeff.size() + B.coeff.size() - 1, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    for (size_t j = 0; j < B.coeff.size(); j++) {
      C.coeff[i + j] += A.coeff[i] * B.coeff[j];
    }
  }
  C.normalize();
  return C;
}

// P(a X)
template <typename T>
UnivariatePolynomial<T> PolyScaleVariable(UnivariatePolynomial<T> const &P,
                                          T const &a) {
  UnivariatePolynomial<T> Q = P;
  T pw(1);
  for (auto &c : Q.coeff) {
    c *= pw;
    pw *= a;
  }
  Q.normalize();
  return Q;
}

template <typename T>
T PolyEvaluate(UnivariatePolynomial<T> const &P, T const &x) {
  T val(0);
  for (size_t i = P.coeff.size(); i > 0; i--) {
    val = val * x + P.coeff[i - 1];
  }
  return val;
}

// A / B for a B that divides A exactly; the remainder is checked to vanish.
template <typename T>
UnivariatePolynomial<T> PolyExactDivision(UnivariatePolynomial<T> const &A,
                                          UnivariatePolynomial<T> const &B) {
  if (B.is_zero()) {
    std::cerr << "SSER: polynomial division by zero\n";
    throw TerminalException{1};
  }
  UnivariatePolynomial<T> R = A;
  UnivariatePolynomial<T> Q;
  int dB = B.degree();
  if (R.degree() >= dB) {
    Q.coeff.assign(R.degree() - dB + 1, T(0));
  }
  T lead = B.coeff[dB];
  while (!R.is_zero() && R.degree() >= dB) {
    int shift = R.degree() - dB;
    T c = R.coeff[R.degree()] / lead;
    Q.coeff[shift] = c;
    for (int i = 0; i <= dB; i++) {
      R.coeff[shift + i] -= c * B.coeff[i];
    }
    R.normalize();
  }
  if (!R.is_zero()) {
    std::cerr << "SSER: a division that should be exact left a remainder\n";
    throw TerminalException{1};
  }
  Q.normalize();
  return Q;
}

template <typename T> struct RationalFunction {
  UnivariatePolynomial<T> num;
  UnivariatePolynomial<T> den;
};

template <typename T>
RationalFunction<T> RatFromPoly(UnivariatePolynomial<T> const &P) {
  return {P, PolyConstant<T>(T(1))};
}

template <typename T>
RationalFunction<T> operator*(RationalFunction<T> const &A,
                              RationalFunction<T> const &B) {
  return {A.num * B.num, A.den * B.den};
}

template <typename T>
RationalFunction<T> operator+(RationalFunction<T> const &A,
                              RationalFunction<T> const &B) {
  return {A.num * B.den + B.num * A.den, A.den * B.den};
}

template <typename T>
RationalFunction<T> RatScaleVariable(RationalFunction<T> const &A,
                                     T const &a) {
  return {PolyScaleVariable(A.num, a), PolyScaleVariable(A.den, a)};
}

template <typename T>
UnivariatePolynomial<T> RatToPolynomial(RationalFunction<T> const &A) {
  return PolyExactDivision(A.num, A.den);
}
// 1 - c X^d
template <typename T> UnivariatePolynomial<T> OneMinus(T const &c, int d) {
  return PolyConstant<T>(T(1)) - PolyMonomial<T>(c, d);
}

// clang-format off
#endif  // SRC_GENERA_RATIONALFUNCTIONS_H_
// clang-format on
