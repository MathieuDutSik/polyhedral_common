// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_JET_NUMBER_H_
#define SRC_DELAUNAY_JET_NUMBER_H_

// clang-format off
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
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

  // ---- additive group (member unary and compound assignment) ----
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
  // higher-level cross-validation code.
  T eval(T const &t) const {
    T s = c[N];
    for (int k = N - 1; k >= 0; k--)
      s = s * t + c[k];
    return s;
  }

  // ---- ring operations and total order, as hidden friends so that a mixed
  // expression like `x != 0` or `alpha * x` (with an int / T literal) promotes
  // the literal to a constant jet. Plain function templates would not apply the
  // implicit conversion during argument deduction. ----
  friend jet operator+(jet a, jet const &b) {
    a += b;
    return a;
  }
  friend jet operator-(jet a, jet const &b) {
    a -= b;
    return a;
  }
  friend jet operator*(jet const &a, jet const &b) {
    jet r; // zero
    for (int i = 0; i <= N; i++)
      for (int j = 0; i + j <= N; j++)
        r.c[i + j] += a.c[i] * b.c[j];
    return r;
  }

  // Inverse 1/f of a jet with non-zero constant term, from f * f^{-1} = 1:
  // b_0 = 1/c_0 and b_k = -b_0 sum_{j=1}^{k} c_j b_{k-j}.
  friend jet inverse(jet const &f) {
#ifdef SANITY_CHECK_JET_NUMBER
    if (f.c[0] == T(0)) {
      std::cerr << "jet_number: inverse of a jet with zero constant term "
                   "(would be a Laurent series, unsupported)\n";
      throw TerminalException{1};
    }
#endif
    jet r;
    r.c[0] = T(1) / f.c[0];
    for (int k = 1; k <= N; k++) {
      T s(0);
      for (int j = 1; j <= k; j++)
        s += f.c[j] * r.c[k - j];
      r.c[k] = -r.c[0] * s;
    }
#ifdef SANITY_CHECK_JET_NUMBER
    // Cheap invariant: f * f^{-1} == 1 (one order-N convolution).
    jet prod = f * r;
    if (prod.c[0] != T(1)) {
      std::cerr << "jet_number: inverse check failed (c0 != 1)\n";
      throw TerminalException{1};
    }
    for (int k = 1; k <= N; k++)
      if (prod.c[k] != T(0)) {
        std::cerr << "jet_number: inverse check failed (c" << k << " != 0)\n";
        throw TerminalException{1};
      }
#endif
    return r;
  }

  friend jet operator/(jet const &a, jet const &b) {
    jet q = a * inverse(b);
#ifdef SANITY_CHECK_JET_NUMBER
    // Cheap invariant: (a / b) * b == a.
    jet chk = q * b;
    for (int k = 0; k <= N; k++)
      if (chk.c[k] != a.c[k]) {
        std::cerr << "jet_number: division check failed at c" << k << "\n";
        throw TerminalException{1};
      }
#endif
    return q;
  }

  // Total order by the leading coefficient of the difference.
  friend bool operator==(jet const &a, jet const &b) {
    for (int k = 0; k <= N; k++)
      if (a.c[k] != b.c[k])
        return false;
    return true;
  }
  friend bool operator!=(jet const &a, jet const &b) { return !(a == b); }
  friend bool operator<(jet const &a, jet const &b) {
    return (a - b).sign() < 0;
  }
  friend bool operator>(jet const &a, jet const &b) {
    return (a - b).sign() > 0;
  }
  friend bool operator<=(jet const &a, jet const &b) {
    return (a - b).sign() <= 0;
  }
  friend bool operator>=(jet const &a, jet const &b) {
    return (a - b).sign() >= 0;
  }
};

// The constant term (value at t = 0) of a scalar. Generic scalars are their own
// constant term; a jet returns c0. This is the bridge that lets the
// combinatorial / canonical-form subroutines of a scalar-templated computation
// (which need a concrete field element) recover the t = 0 data from a jet
// Gram matrix, while the numeric parts keep the full expansion.
template <typename T> T const &constant_term(T const &x) { return x; }
template <typename T, int N> T const &constant_term(jet<T, N> const &j) {
  return j.c[0];
}

// The k-th derivative at t = 0: k! times the coefficient of t^k.
template <typename T, int N> T jet_deriv(jet<T, N> const &j, int k) {
  T fact(1);
  for (int i = 2; i <= k; i++)
    fact *= T(i);
  return fact * j.c[k];
}

// ---- number-type traits so jet<T, N> flows through the generic matrix code
// (MyMatrix, DeterminantMat, Inverse, SolutionMat, ...). The truncated series
// ring is a local ring, not a field, but the hand-written Gaussian-elimination
// kernels behave as if over a field as long as the pivots have non-zero
// constant term (guaranteed for Q + t H with Q non-singular, at leading order).
template <typename T, int N> struct is_ring_field<jet<T, N>> {
  static const bool value = is_ring_field<T>::value;
};
template <typename T, int N> struct is_totally_ordered<jet<T, N>> {
  static const bool value = true;
};
template <typename T, int N> struct is_exact_arithmetic<jet<T, N>> {
  static const bool value = is_exact_arithmetic<T>::value;
};
template <typename T, int N> struct overlying_field<jet<T, N>> {
  typedef jet<typename overlying_field<T>::field_type, N> field_type;
};
template <typename T, int N> struct underlying_ring<jet<T, N>> {
  typedef jet<typename underlying_ring<T>::ring_type, N> ring_type;
};

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
