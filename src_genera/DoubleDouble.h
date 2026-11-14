// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_DOUBLEDOUBLE_H_
#define SRC_GENERA_DOUBLEDOUBLE_H_

// clang-format off
#include "MprealConversion.h"
#include <cmath>
#include <cstdint>
// clang-format on

/*
  Double-double numbers hi + lo, |lo| <= ulp(hi) / 2, about 106 bits, for
  the terms of the sum that need more than double precision but much less
  than the working precision (Dekker, Knuth; Hida, Li and Bailey).
 */
struct QuadDD {
  double hi;
  double lo;
};

inline QuadDD QuadDDTwoSum(double a, double b) {
  double s = a + b;
  double bb = s - a;
  double e = (a - (s - bb)) + (b - bb);
  return {s, e};
}

inline QuadDD QuadDDQuickTwoSum(double a, double b) {
  double s = a + b;
  return {s, b - (s - a)};
}

inline QuadDD QuadDDAdd(QuadDD a, QuadDD b) {
  QuadDD s = QuadDDTwoSum(a.hi, b.hi);
  QuadDD t = QuadDDTwoSum(a.lo, b.lo);
  s.lo += t.hi;
  s = QuadDDQuickTwoSum(s.hi, s.lo);
  s.lo += t.lo;
  return QuadDDQuickTwoSum(s.hi, s.lo);
}

inline QuadDD QuadDDMul(QuadDD a, QuadDD b) {
  double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  return QuadDDQuickTwoSum(p, e);
}

inline QuadDD QuadDDFromMpreal(mpfr::mpreal const &x) {
  double hi = x.toDouble();
  double lo = (x - mpfr::mpreal(hi)).toDouble();
  return QuadDDQuickTwoSum(hi, lo);
}

inline mpfr::mpreal QuadDDToMpreal(QuadDD const &x) {
  return mpfr::mpreal(x.hi) + mpfr::mpreal(x.lo);
}

// exp(x) for |x| < 700: x = k log 2 + r, exp(r / 2^10) by Taylor, then 10
// squarings kept as expm1 (2 e + e^2), and the scaling by 2^k
inline QuadDD QuadDDExp(QuadDD x) {
  const QuadDD ln2{0.6931471805599453, 2.3190468138462996e-17};
  double k = std::nearbyint(x.hi / ln2.hi);
  QuadDD r = QuadDDAdd(x, QuadDDMul(ln2, QuadDD{-k, 0.0}));
  r.hi = std::ldexp(r.hi, -10);
  r.lo = std::ldexp(r.lo, -10);
  // expm1(r) = r + r^2/2 + ... + r^9/9!
  QuadDD e = r;
  QuadDD term = r;
  for (int i = 2; i <= 9; i++) {
    term = QuadDDMul(term, r);
    term.hi /= i;
    term.lo /= i;
    term = QuadDDQuickTwoSum(term.hi, term.lo);
    e = QuadDDAdd(e, term);
  }
  for (int i = 0; i < 10; i++) {
    // expm1(2y) = 2 expm1(y) + expm1(y)^2
    QuadDD two_e{2 * e.hi, 2 * e.lo};
    e = QuadDDAdd(two_e, QuadDDMul(e, e));
  }
  QuadDD res = QuadDDAdd(QuadDD{1.0, 0.0}, e);
  int ki = static_cast<int>(k);
  return {std::ldexp(res.hi, ki), std::ldexp(res.lo, ki)};
}
// log q to double-double precision, q < 2^53: y = log q in double, then
// log q = y + log(q e^{-y}) with q e^{-y} = 1 + delta, |delta| < 2^{-50}
inline QuadDD QuadDDLogInteger(int64_t q) {
  double y = std::log(static_cast<double>(q));
  QuadDD e = QuadDDExp(QuadDD{-y, 0.0});
  QuadDD t = QuadDDMul(e, QuadDD{static_cast<double>(q), 0.0});
  QuadDD delta = QuadDDAdd(t, QuadDD{-1.0, 0.0});
  return QuadDDAdd(QuadDD{y, 0.0}, delta);
}

// clang-format off
#endif  // SRC_GENERA_DOUBLEDOUBLE_H_
// clang-format on
