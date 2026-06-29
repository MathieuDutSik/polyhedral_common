// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_POLYNOMIAL_H_
#define SRC_DELAUNAY_POLYNOMIAL_H_

// clang-format off
#include "MAT_Matrix.h"
#include "ExceptionsFunc.h"
#include <ostream>
#include <utility>
#include <vector>
// clang-format on

// Polynomials, represented by their coefficient vector (c0 + c1 t + ...), and
// rational functions P(t)/D(t), together with exact interpolation and
// reconstruction-from-samples utilities. Self-contained (no dependency on the
// quantization code that uses it).

#ifdef DEBUG
#define DEBUG_POLYNOMIAL
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_CORRECT_POLYNOMIAL
#endif

// Evaluate a polynomial (coefficients c0 + c1 t + ...) at t.
template <typename T> T eval_poly(MyVector<T> const &C, T const &t) {
  T s(0);
  T tp(1);
  for (int i = 0; i < C.size(); i++) {
    s += C(i) * tp;
    tp *= t;
  }
  return s;
}

// Product of two polynomials given by their coefficient vectors (c0 + c1 t +
// ...). The result has degree deg(A) + deg(B).
template <typename T>
MyVector<T> poly_product(MyVector<T> const &A, MyVector<T> const &B) {
  int da = A.size();
  int db = B.size();
  MyVector<T> C = ZeroVector<T>(da + db - 1);
  for (int a = 0; a < da; a++) {
    for (int b = 0; b < db; b++) {
      C(a + b) += A(a) * B(b);
    }
  }
  return C;
}

// Interpolate the polynomial of degree <= deg through the points (t_k, v_k).
template <typename T>
MyVector<T> poly_from_values(std::vector<T> const &ts, std::vector<T> const &vs,
                             int deg) {
  int N = deg + 1;
  MyMatrix<T> Vand(N, N);
  MyVector<T> rhs(N);
  for (int i = 0; i < N; i++) {
    T tp(1);
    for (int j = 0; j < N; j++) {
      Vand(i, j) = tp;
      tp *= ts[i];
    }
    rhs(i) = vs[i];
  }
  MyMatrix<T> VandInv = Inverse(Vand);
  return VandInv * rhs;
}

template <typename T> struct RationalFunc {
  MyVector<T> P;
  MyVector<T> D;
  int degree;
};

// Reconstruct a rational function N(t)/den(t) where the denominator den (with
// den(0) = 1) is known exactly and the numerator N has degree <= max_degree. The
// pool must hold exactly max_degree + 1 points: N is the unique interpolant of
// degree max_degree through them (N(t_k) = SecMoment(t_k) * den(t_k)). Its true
// degree may be lower, so trailing zero coefficients are trimmed and the trimmed
// degree is reported.
//
// Under SANITY_CHECK_CORRECT_POLYNOMIAL the reconstruction is checked against the
// sampler at extra points subdividing each pool interval [t_i, t_{i+1}] into
// tenths (9 interior points per interval). A mismatch exposes a wrong
// denominator, a numerator degree above max_degree, or the segment crossing an
// iso-Delaunay wall, and throws. Without the check the interpolant is trusted.
template <typename T, typename Fsampler>
RationalFunc<T>
reconstruct_rational_known_denominator(std::vector<T> const &tpool,
                                       Fsampler sampler, MyVector<T> const &den,
                                       int max_degree,
                                       [[maybe_unused]] std::ostream &os) {
  int N = max_degree + 1;
  std::vector<T> ts(N), ns(N);
  for (int i = 0; i < N; i++) {
    ts[i] = tpool[i];
    ns[i] = sampler(tpool[i]) * eval_poly(den, tpool[i]);
  }
  MyVector<T> P = poly_from_values(ts, ns, max_degree);
  int degree = max_degree;
  while (degree > 0 && P(degree) == T(0)) {
    degree--;
  }
  MyVector<T> Ptrim(degree + 1);
  for (int i = 0; i <= degree; i++) {
    Ptrim(i) = P(i);
  }
#ifdef SANITY_CHECK_CORRECT_POLYNOMIAL
  for (int i = 0; i + 1 < N; i++) {
    for (int j = 1; j <= 9; j++) {
      T tt = (ts[i] * T(j) + ts[i + 1] * T(10 - j)) / T(10);
      T expected = sampler(tt) * eval_poly(den, tt);
      if (eval_poly(Ptrim, tt) != expected) {
        std::cerr << "POLY: reconstruct_rational_known_denominator "
                     "SANITY_CHECK_CORRECT_POLYNOMIAL FAILED at t="
                  << tt << " (degree " << degree << ")\n";
        throw TerminalException{1};
      }
    }
  }
#endif
#ifdef DEBUG_POLYNOMIAL
  os << "POLY: numerator (known denominator) degree=" << degree << "\n";
#endif
  return RationalFunc<T>{Ptrim, den, degree};
}

// clang-format off
#endif  // SRC_DELAUNAY_POLYNOMIAL_H_
// clang-format on
