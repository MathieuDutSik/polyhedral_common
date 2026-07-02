// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_JET_NUMBER_H_
#define SRC_DELAUNAY_JET_NUMBER_H_

// clang-format off
#include "ExceptionsFunc.h"
#include <array>
#include <ostream>
#include <utility>
// clang-format on

// A jet number: a truncated Taylor expansion c0 + c1 t + ... + cN t^N in an
// infinitesimal parameter t, kept at a fixed compile-time order N. Unlike the
// passive Jet<T> holder in jet.h, this is a full numeric TYPE meant to be used
// as the scalar of generic templated code (matrices, determinants, CVP, the
// quantization integral, ...): a matrix Q + t H over jet<T, N> propagates the
// order-N expansion of every derived quantity in a single exact computation, so
// the derivatives fall out of the coefficients (jet_deriv below) instead of
// being interpolated from many concrete-t samples.
//
// Every operation truncates back to order N, so a chain of operations stays at
// order N (truncation is a ring homomorphism onto T[t]/(t^{N+1})).
//
// ORDER (the total order making t an infinitesimal positive): the sign of a jet
// as t -> 0^+ is the sign of its first non-zero coefficient, so a is >= 0 iff
// its first non-zero coefficient is >= 0 (an all-zero jet is 0 >= 0). All the
// comparison operators derive from this leading-coefficient rule applied to the
// difference a - b. This is exactly the ordering needed for the geometric
// predicates (positivity, CVP inequalities, pivoting) on the segment t -> 0^+.

#ifdef SANITY_CHECK
#define SANITY_CHECK_JET_NUMBER
#endif

// The extensive checks are DELIBERATELY not enabled by SANITY_CHECK: they are
// far more expensive than the ordinary invariant checks and must be requested
// on their own, e.g. make CHOICE_COMPILATION=-DSANITY_CHECK_EXTENSIVE_JET_NUMBER
// (or via the umbrella -DSANITY_CHECK_EXTENSIVE). Each such block announces
// itself with a print at the beginning of the extra computation.
#ifdef SANITY_CHECK_EXTENSIVE
#define SANITY_CHECK_EXTENSIVE_JET_NUMBER
#endif

template <typename T, int N> class jet {
public:
  std::array<T, N + 1> c; // c[k] = coefficient of t^k, 0 <= k <= N

  // Zero.
  jet() { c.fill(T(0)); }
  // Constant jets (non-explicit so generic code's T(0), T(1), 0, ... promote).
  jet(int v) {
    c.fill(T(0));
    c[0] = T(v);
  }
  jet(T const &v) {
    c.fill(T(0));
    c[0] = v;
  }
  // The infinitesimal t itself (c1 = 1, all others 0).
  static jet var() {
    jet r;
    r.c[1] = T(1);
    return r;
  }
  // The jet whose coefficient array is exactly `coeffs`.
  static jet from_coeffs(std::array<T, N + 1> coeffs) {
    jet r;
    r.c = std::move(coeffs);
    return r;
  }

  T const &operator[](int k) const { return c[k]; }
  T &operator[](int k) { return c[k]; }

  // ---- additive group ----
  jet operator-() const {
    jet r;
    for (int k = 0; k <= N; k++)
      r.c[k] = -c[k];
    return r;
  }
  jet &operator+=(jet const &o) {
    for (int k = 0; k <= N; k++)
      c[k] += o.c[k];
    return *this;
  }
  jet &operator-=(jet const &o) {
    for (int k = 0; k <= N; k++)
      c[k] -= o.c[k];
    return *this;
  }
  jet &operator*=(jet const &o) {
    *this = (*this) * o;
    return *this;
  }
  jet &operator/=(jet const &o) {
    *this = (*this) / o;
    return *this;
  }

  // Index of the first non-zero coefficient, or N+1 if the jet is zero.
  int leading_index() const {
    for (int k = 0; k <= N; k++)
      if (c[k] != T(0))
        return k;
    return N + 1;
  }
  // Sign as t -> 0^+ : sign of the first non-zero coefficient (0 if all zero).
  int sign() const {
    int k = leading_index();
    if (k > N)
      return 0;
    if (c[k] > T(0))
      return 1;
    return -1;
  }

  // Evaluate the truncated series at a concrete value t (Horner). Used by the
  // extensive checks and by the higher-level cross-validation code.
  T eval(T const &t) const {
    T s = c[N];
    for (int k = N - 1; k >= 0; k--)
      s = s * t + c[k];
    return s;
  }
};

// ---- ring operations ----
template <typename T, int N>
jet<T, N> operator+(jet<T, N> const &a, jet<T, N> const &b) {
  jet<T, N> r;
  for (int k = 0; k <= N; k++)
    r.c[k] = a.c[k] + b.c[k];
  return r;
}

template <typename T, int N>
jet<T, N> operator-(jet<T, N> const &a, jet<T, N> const &b) {
  jet<T, N> r;
  for (int k = 0; k <= N; k++)
    r.c[k] = a.c[k] - b.c[k];
  return r;
}

template <typename T, int N>
jet<T, N> operator*(jet<T, N> const &a, jet<T, N> const &b) {
  jet<T, N> r; // zero
  for (int i = 0; i <= N; i++)
    for (int j = 0; i + j <= N; j++)
      r.c[i + j] += a.c[i] * b.c[j];
  return r;
}

// Inverse 1/f of a jet with non-zero constant term, from f * f^{-1} = 1:
// b_0 = 1/c_0 and b_k = -b_0 sum_{j=1}^{k} c_j b_{k-j}.
template <typename T, int N> jet<T, N> inverse(jet<T, N> const &f) {
#ifdef SANITY_CHECK_JET_NUMBER
  if (f.c[0] == T(0)) {
    std::cerr << "jet_number: inverse of a jet with zero constant term "
                 "(would be a Laurent series, unsupported)\n";
    throw TerminalException{1};
  }
#endif
  jet<T, N> r;
  r.c[0] = T(1) / f.c[0];
  for (int k = 1; k <= N; k++) {
    T s(0);
    for (int j = 1; j <= k; j++)
      s += f.c[j] * r.c[k - j];
    r.c[k] = -r.c[0] * s;
  }
#ifdef SANITY_CHECK_EXTENSIVE_JET_NUMBER
  std::cerr << "SANITY_CHECK_EXTENSIVE_JET_NUMBER: verifying inverse via "
               "f * f^{-1} == 1\n";
  jet<T, N> prod = f * r;
  if (prod.c[0] != T(1)) {
    std::cerr << "jet_number: inverse EXTENSIVE check failed (c0 != 1)\n";
    throw TerminalException{1};
  }
  for (int k = 1; k <= N; k++)
    if (prod.c[k] != T(0)) {
      std::cerr << "jet_number: inverse EXTENSIVE check failed (c" << k
                << " != 0)\n";
      throw TerminalException{1};
    }
#endif
  return r;
}

template <typename T, int N>
jet<T, N> operator/(jet<T, N> const &a, jet<T, N> const &b) {
  jet<T, N> q = a * inverse(b);
#ifdef SANITY_CHECK_EXTENSIVE_JET_NUMBER
  std::cerr << "SANITY_CHECK_EXTENSIVE_JET_NUMBER: verifying division via "
               "(a / b) * b == a\n";
  jet<T, N> chk = q * b;
  for (int k = 0; k <= N; k++)
    if (chk.c[k] != a.c[k]) {
      std::cerr << "jet_number: division EXTENSIVE check failed at c" << k
                << "\n";
      throw TerminalException{1};
    }
#endif
  return q;
}

// ---- total order by the leading coefficient of the difference ----
template <typename T, int N>
bool operator==(jet<T, N> const &a, jet<T, N> const &b) {
  for (int k = 0; k <= N; k++)
    if (a.c[k] != b.c[k])
      return false;
  return true;
}
template <typename T, int N>
bool operator!=(jet<T, N> const &a, jet<T, N> const &b) {
  return !(a == b);
}
template <typename T, int N>
bool operator<(jet<T, N> const &a, jet<T, N> const &b) {
  return (a - b).sign() < 0;
}
template <typename T, int N>
bool operator>(jet<T, N> const &a, jet<T, N> const &b) {
  return (a - b).sign() > 0;
}
template <typename T, int N>
bool operator<=(jet<T, N> const &a, jet<T, N> const &b) {
  return (a - b).sign() <= 0;
}
template <typename T, int N>
bool operator>=(jet<T, N> const &a, jet<T, N> const &b) {
  return (a - b).sign() >= 0;
}

// The k-th derivative at t = 0: k! times the coefficient of t^k.
template <typename T, int N> T jet_deriv(jet<T, N> const &j, int k) {
  T fact(1);
  for (int i = 2; i <= k; i++)
    fact *= T(i);
  return fact * j.c[k];
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, jet<T, N> const &j) {
  os << "[";
  for (int k = 0; k <= N; k++) {
    if (k > 0)
      os << ", ";
    os << j.c[k];
  }
  os << "]";
  return os;
}

// clang-format off
#endif  // SRC_DELAUNAY_JET_NUMBER_H_
// clang-format on
