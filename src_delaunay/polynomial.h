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
#define SANITY_CHECK_POLYNOMIAL
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
// den(0) = 1) is known exactly and only the numerator polynomial N (of degree
// <= max_degree) has to be interpolated. The numerator degree is increased until
// extra held-out samples validate the fit. Throws if no numerator of degree
// <= max_degree fits, i.e. either the supplied denominator is wrong or the
// numerator degree exceeds max_degree.
template <typename T, typename Fsampler>
RationalFunc<T>
reconstruct_rational_known_denominator(std::vector<T> const &tpool,
                                       Fsampler sampler, MyVector<T> const &den,
                                       int max_degree,
                                       [[maybe_unused]] std::ostream &os) {
  std::vector<std::pair<T, T>> cache; // (t, N = value * den)
  auto get = [&](int i) -> std::pair<T, T> {
    while (static_cast<int>(cache.size()) <= i) {
      T tt = tpool[cache.size()];
      T ss = sampler(tt);
      T nn = ss * eval_poly(den, tt);
      cache.push_back({tt, nn});
    }
    return cache[i];
  };
  int nval = 4;
  for (int d = 0; d <= max_degree; d++) {
    int N = d + 1;
    if (N + nval > static_cast<int>(tpool.size())) {
      break;
    }
    std::vector<T> ts, ns;
    for (int i = 0; i < N; i++) {
      std::pair<T, T> pr = get(i);
      ts.push_back(pr.first);
      ns.push_back(pr.second);
    }
    MyVector<T> P = poly_from_values(ts, ns, d);
    bool ok = true;
    for (int i = N; i < N + nval; i++) {
      std::pair<T, T> pr = get(i);
      if (eval_poly(P, pr.first) != pr.second) {
        ok = false;
        break;
      }
    }
    if (ok) {
#ifdef SANITY_CHECK_POLYNOMIAL
      // Further validation: the interpolated numerator must reproduce every
      // remaining sample of the pool exactly, not only the nval held out for
      // acceptance. Each further value is an extra (expensive) evaluation.
      for (int i = N + nval; i < static_cast<int>(tpool.size()); i++) {
        std::pair<T, T> pr = get(i);
        if (eval_poly(P, pr.first) != pr.second) {
          std::cerr << "POLY: known-denominator interpolation SANITY_CHECK "
                       "FAILED at t="
                    << pr.first << " (degree " << d << ")\n";
          throw TerminalException{1};
        }
      }
#endif
#ifdef DEBUG_POLYNOMIAL
      os << "POLY: numerator (known denominator) converged at degree=" << d
         << "\n";
#endif
      return RationalFunc<T>{P, den, d};
    }
  }
  std::cerr << "POLY: known-denominator reconstruction failed: no numerator of "
               "degree <= "
            << max_degree << " fits (wrong denominator or degree too high)\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_DELAUNAY_POLYNOMIAL_H_
// clang-format on
