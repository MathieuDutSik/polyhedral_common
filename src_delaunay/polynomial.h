// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_POLYNOMIAL_H_
#define SRC_DELAUNAY_POLYNOMIAL_H_

// clang-format off
#include "MAT_Matrix.h"
#include "ExceptionsFunc.h"
#include <optional>
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

// Fit a rational function S(t) = P(t)/D(t) with deg P <= p, deg D <= q and
// D(0) = 1, through the given (t, s) points. Returns (P, D) or nothing if the
// linear system is singular.
template <typename T>
std::optional<std::pair<MyVector<T>, MyVector<T>>>
fit_rational(std::vector<std::pair<T, T>> const &pts, int p, int q) {
  int N = p + 1 + q;
  if (static_cast<int>(pts.size()) < N) {
    return {};
  }
  MyMatrix<T> Eq(N, N);
  MyVector<T> rhs(N);
  for (int i = 0; i < N; i++) {
    T ti = pts[i].first;
    T si = pts[i].second;
    T tp(1);
    for (int j = 0; j <= p; j++) {
      Eq(i, j) = tp;
      tp *= ti;
    }
    T tk = ti;
    for (int k = 1; k <= q; k++) {
      Eq(i, p + k) = -si * tk;
      tk *= ti;
    }
    rhs(i) = si;
  }
  if (RankMat(Eq) < N) {
    return {};
  }
  MyMatrix<T> EqInv = Inverse(Eq);
  MyVector<T> x = EqInv * rhs;
  MyVector<T> P(p + 1);
  MyVector<T> D(q + 1);
  for (int j = 0; j <= p; j++) {
    P(j) = x(j);
  }
  D(0) = T(1);
  for (int k = 1; k <= q; k++) {
    D(k) = x(p + k);
  }
  return std::make_pair(P, D);
}

template <typename T> struct RationalFunc {
  MyVector<T> P;
  MyVector<T> D;
  int degree;
};

// Reconstruct a rational function N(t)/den(t) where the denominator den (with
// den(0) = 1) is known exactly (e.g. it equals det(Q + t H)/det(Q)) and only
// the numerator polynomial N has to be interpolated. The numerator degree is
// increased until extra held-out samples validate the fit. Returns nothing if
// no polynomial numerator of degree <= max_degree fits (i.e. the denominator is
// not the expected one), so the caller can fall back to a full rational fit.
template <typename T, typename Fsampler>
std::optional<RationalFunc<T>>
reconstruct_rational_known_denominator(std::vector<T> const &tpool,
                                       Fsampler sampler, MyVector<T> const &den,
                                       int max_degree, std::ostream &os) {
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
  return {};
}

// Reconstruct a rational function P(t)/D(t) from samples produced lazily by
// sampler, at the t-values in tpool. The degree is increased until extra
// held-out samples validate the fit.
template <typename T, typename Fsampler>
RationalFunc<T> reconstruct_rational(std::vector<T> const &tpool,
                                     Fsampler sampler, int max_degree,
                                     std::ostream &os) {
  std::vector<std::pair<T, T>> cache;
  auto get = [&](int i) -> std::pair<T, T> {
    while (static_cast<int>(cache.size()) <= i) {
      T tt = tpool[cache.size()];
      T ss = sampler(tt);
      cache.push_back({tt, ss});
    }
    return cache[i];
  };
  int nval = 4;
  for (int d = 1; d <= max_degree; d++) {
    int N = 2 * d + 1;
    if (N + nval > static_cast<int>(tpool.size())) {
      break;
    }
    std::vector<std::pair<T, T>> fitpts;
    for (int i = 0; i < N; i++) {
      fitpts.push_back(get(i));
    }
    std::optional<std::pair<MyVector<T>, MyVector<T>>> opt =
        fit_rational(fitpts, d, d);
    if (!opt) {
      continue;
    }
    MyVector<T> const &P = opt->first;
    MyVector<T> const &D = opt->second;
    bool ok = true;
    for (int i = N; i < N + nval; i++) {
      std::pair<T, T> pr = get(i);
      if (eval_poly(P, pr.first) != pr.second * eval_poly(D, pr.first)) {
        ok = false;
        break;
      }
    }
    if (ok) {
#ifdef SANITY_CHECK_POLYNOMIAL
      // Further validation: the interpolated rational function must reproduce
      // every remaining sample of the pool exactly, not only the nval held out
      // for acceptance. Each further value is an extra (expensive) evaluation.
      for (int i = N + nval; i < static_cast<int>(tpool.size()); i++) {
        std::pair<T, T> pr = get(i);
        if (eval_poly(P, pr.first) != pr.second * eval_poly(D, pr.first)) {
          std::cerr << "POLY: rational interpolation SANITY_CHECK FAILED at t="
                    << pr.first << " (degree " << d << ")\n";
          throw TerminalException{1};
        }
      }
#endif
#ifdef DEBUG_POLYNOMIAL
      os << "POLY: reconstruct_rational converged at degree=" << d << "\n";
#endif
      return RationalFunc<T>{P, D, d};
    }
  }
  std::cerr << "POLY: reconstruct_rational failed up to degree=" << max_degree
            << "\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_DELAUNAY_POLYNOMIAL_H_
// clang-format on
