// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SCALARPOLYNOMIALMODQ_H_
#define SRC_GENERA_SCALARPOLYNOMIALMODQ_H_

// clang-format off
#include "ArithmeticFunctions.h"
#include "ResidueArithmetic.h"
#include "ResidueArithmetic.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_POLYQ
#endif

/*
  Polynomials with coefficients modulo a prime q (PolyQ) and with integer
  coefficients (PolyZ): the ring operations, division, gcd and inverses
  modulo a polynomial, powers modulo a polynomial, and the factorization by
  distinct and equal degrees. The elements of the residue fields
  F_q[x] / (mu) are another type (FqElt, ScalarFiniteField.h).
 */

/*
  A polynomial modulo the prime q: its coefficients in [0, q) by increasing
  degree, trimmed (the leading coefficient is nonzero, the zero polynomial
  has none), so that equality is that of the coefficients. The
  constructors reduce and trim; the polynomials of an operation have the
  same q (checked under SANITY_CHECK_POLYQ). The default one, of q = 0, is
  only a placeholder to be assigned.
 */
struct PolyQ {
  T_small_prime q = 0;

  PolyQ() = default;
  // the zero polynomial
  explicit PolyQ(T_small_prime q_in) : q(q_in) {}
  PolyQ(T_small_prime q_in, std::vector<T_residue> coeffs)
      : q(q_in), c(std::move(coeffs)) {
    for (auto &x : c) {
      x %= q;
    }
    trim();
  }
  static PolyQ X(T_small_prime q) { return PolyQ(q, {0, 1}); }

  // -1 for the zero polynomial
  int degree() const { return static_cast<int>(c.size()) - 1; }
  bool is_zero() const { return c.empty(); }
  size_t size() const { return c.size(); }
  T_residue operator[](size_t i) const { return c[i]; }
  // the coefficient of x^i, 0 beyond the degree
  T_residue coeff(size_t i) const { return (i < c.size()) ? c[i] : 0; }
  T_residue lead() const { return c.back(); }
  std::vector<T_residue> const &coeffs() const { return c; }
  auto begin() const { return c.begin(); }
  auto end() const { return c.end(); }
  bool operator==(PolyQ const &b) const { return q == b.q && c == b.c; }

private:
  std::vector<T_residue> c;

  void trim() {
    while (!c.empty() && c.back() == 0) {
      c.pop_back();
    }
  }
};

inline void PolyCheckSameQ([[maybe_unused]] PolyQ const &a,
                           [[maybe_unused]] PolyQ const &b) {
#ifdef SANITY_CHECK_POLYQ
  if (a.q == 0 || a.q != b.q) {
    std::cerr << "POLYQ: polynomials modulo " << a.q << " and " << b.q << "\n";
    throw TerminalException{1};
  }
#endif
}

// a polynomial with integer coefficients modulo q
inline PolyQ PolyReduceSigned(std::vector<int64_t> const &a, T_small_prime q) {
  std::vector<T_residue> r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResModSigned(a[i], q);
  }
  return PolyQ(q, std::move(r));
}

inline PolyQ operator+(PolyQ const &a, PolyQ const &b) {
  PolyCheckSameQ(a, b);
  T_small_prime q = a.q;
  std::vector<T_residue> r(std::max(a.size(), b.size()));
  for (size_t i = 0; i < r.size(); i++) {
    r[i] = ResAdd(a.coeff(i), b.coeff(i), q);
  }
  return PolyQ(q, std::move(r));
}

inline PolyQ operator-(PolyQ const &a) {
  std::vector<T_residue> r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResNeg(a[i], a.q);
  }
  return PolyQ(a.q, std::move(r));
}

inline PolyQ operator-(PolyQ const &a, PolyQ const &b) {
  PolyCheckSameQ(a, b);
  T_small_prime q = a.q;
  std::vector<T_residue> r(std::max(a.size(), b.size()));
  for (size_t i = 0; i < r.size(); i++) {
    r[i] = ResSub(a.coeff(i), b.coeff(i), q);
  }
  return PolyQ(q, std::move(r));
}

inline PolyQ operator*(PolyQ const &a, PolyQ const &b) {
  PolyCheckSameQ(a, b);
  T_small_prime q = a.q;
  if (a.is_zero() || b.is_zero()) {
    return PolyQ(q);
  }
  std::vector<T_residue> r(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      r[i + j] = ResAdd(r[i + j], ResMul(a[i], b[j], q), q);
    }
  }
  return PolyQ(q, std::move(r));
}

// s a
inline PolyQ PolyScaleQ(PolyQ const &a, T_residue s) {
  std::vector<T_residue> r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResMul(a[i], s, a.q);
  }
  return PolyQ(a.q, std::move(r));
}

// quotient and remainder, b nonzero
inline std::pair<PolyQ, PolyQ> PolyDivModQ(PolyQ const &a, PolyQ const &b) {
  PolyCheckSameQ(a, b);
  T_small_prime q = a.q;
#ifdef SANITY_CHECK_POLYQ
  if (b.is_zero()) {
    std::cerr << "POLYQ: division by zero\n";
    throw TerminalException{1};
  }
#endif
  if (a.size() < b.size()) {
    return {PolyQ(q), a};
  }
  std::vector<T_residue> r = a.coeffs();
  T_residue inv = ResInv(b.lead(), q);
  std::vector<T_residue> quo(a.size() - b.size() + 1, 0);
  for (size_t k = quo.size(); k-- > 0;) {
    T_residue c = ResMul(r[k + b.size() - 1], inv, q);
    quo[k] = c;
    for (size_t j = 0; j < b.size(); j++) {
      r[k + j] = ResSub(r[k + j], ResMul(c, b[j], q), q);
    }
  }
  return {PolyQ(q, std::move(quo)), PolyQ(q, std::move(r))};
}

inline PolyQ operator/(PolyQ const &a, PolyQ const &b) {
  return PolyDivModQ(a, b).first;
}

inline PolyQ operator%(PolyQ const &a, PolyQ const &b) {
  return PolyDivModQ(a, b).second;
}

inline PolyQ PolyMonicQ(PolyQ const &a) {
  if (a.is_zero()) {
    return a;
  }
  return PolyScaleQ(a, ResInv(a.lead(), a.q));
}

inline PolyQ PolyGcdQ(PolyQ a, PolyQ b) {
  while (!b.is_zero()) {
    PolyQ r = a % b;
    a = b;
    b = r;
  }
  return PolyMonicQ(a);
}

// the inverse of a modulo b, a and b coprime
inline PolyQ PolyInverseModQ(PolyQ const &a, PolyQ const &b) {
  T_small_prime q = b.q;
  PolyQ r0 = b, r1 = a % b;
  PolyQ s0(q), s1(q, {1});
  while (!r1.is_zero()) {
    std::pair<PolyQ, PolyQ> qr = PolyDivModQ(r0, r1);
    PolyQ s2 = s0 - qr.first * s1;
    r0 = r1;
    r1 = qr.second;
    s0 = s1;
    s1 = s2;
  }
  // r0 is a nonzero constant
#ifdef SANITY_CHECK_POLYQ
  if (r0.size() != 1) {
    std::cerr << "POLYQ: the polynomials are not coprime\n";
    throw TerminalException{1};
  }
#endif
  return PolyScaleQ(s0, ResInv(r0[0], q)) % b;
}

inline PolyQ PolyPowModQ(PolyQ base, uint64_t e, PolyQ const &m) {
  PolyQ res = PolyQ(m.q, {1}) % m;
  base = base % m;
  while (e > 0) {
    if (e & 1) {
      res = (res * base) % m;
    }
    base = (base * base) % m;
    e >>= 1;
  }
  return res;
}

/*
  The irreducible factors (monic) of the squarefree g modulo q, all of
  degree f, by the equal degree splitting of Cantor and Zassenhaus (q odd)
  or by the trace map (q = 2). The random polynomials are replaced by a
  deterministic sequence, which splits as well.
 */
inline std::vector<PolyQ> EqualDegreeFactorization(PolyQ const &g, int f) {
  T_small_prime q = g.q;
  int deg = g.degree();
  if (deg == f) {
    return {PolyMonicQ(g)};
  }
  PolyQ one(q, {1});
  for (uint64_t c = 1;; c++) {
    // r = x^{1 + c mod deg} + c, a varied family of polynomials
    std::vector<T_residue> rc(2 + static_cast<size_t>(c % deg), 0);
    rc.back() = 1;
    rc[0] = c % q;
    PolyQ r(q, std::move(rc));
    PolyQ h(q);
    if (q == 2) {
      PolyQ t = r % g;
      h = t;
      for (int i = 1; i < f; i++) {
        t = (t * t) % g;
        h = h + t;
      }
    } else {
      // r^{(q^f - 1)/2} = prod_{i < f} (r^{(q-1)/2})^{q^i}, which keeps the
      // exponents within 64 bits
      PolyQ t = PolyPowModQ(r, (q - 1) / 2, g);
      h = one;
      for (int i = 0; i < f; i++) {
        h = (h * t) % g;
        t = PolyPowModQ(t, q, g);
      }
      h = h - one;
    }
    PolyQ dgc = PolyGcdQ(g, h);
    int dd = dgc.degree();
    if (dd > 0 && dd < deg) {
      std::vector<PolyQ> res = EqualDegreeFactorization(dgc, f);
      std::vector<PolyQ> res2 = EqualDegreeFactorization(g / dgc, f);
      res.insert(res.end(), res2.begin(), res2.end());
      return res;
    }
    if (c > 100000) {
      std::cerr << "POLYQ: the equal degree splitting failed\n";
      throw TerminalException{1};
    }
  }
}

// x^{deg g} g(1/x), made monic
inline PolyQ ReciprocalPolyQ(PolyQ const &g) {
  std::vector<T_residue> r(g.coeffs().rbegin(), g.coeffs().rend());
  return PolyMonicQ(PolyQ(g.q, std::move(r)));
}

/*
  Polynomials over Z (mpz), reduced modulo Phi_d (monic) and with the
  coefficients modulo Q = q^N.
 */
template <typename Tint> using PolyZ = std::vector<Tint>;

// a polynomial modulo q lifted to its residues in [0, q)
template <typename Tint> PolyZ<Tint> PolyZOfPolyQ(PolyQ const &a) {
  PolyZ<Tint> r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = UniversalScalarConversion<Tint, int64_t>(ResToSigned(a[i]));
  }
  return r;
}

template <typename Tint>
inline PolyZ<Tint> PolyZReduce(PolyZ<Tint> a, std::vector<int64_t> const &Phi,
                               Tint const &Q) {
  int deg = static_cast<int>(Phi.size()) - 1;
  for (int k = static_cast<int>(a.size()) - 1; k >= deg; k--) {
    Tint c = a[k];
    if (c != 0) {
      for (int j = 0; j <= deg; j++) {
        a[k - deg + j] -= c * UniversalScalarConversion<Tint, int64_t>(Phi[j]);
      }
    }
  }
  a.resize(deg, Tint(0));
  for (auto &x : a) {
    x = ResInt(x, Q);
    if (x < 0) {
      x += Q;
    }
  }
  return a;
}

template <typename Tint>
inline PolyZ<Tint> PolyZMul(PolyZ<Tint> const &a, PolyZ<Tint> const &b) {
  PolyZ<Tint> c(a.size() + b.size() - 1, Tint(0));
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

// The distinct irreducible factors of f modulo q with their multiplicities.
inline std::vector<std::pair<PolyQ, int>>
FactorModQWithMultiplicity(PolyQ const &f_in) {
  T_small_prime q = f_in.q;
  PolyQ f = PolyMonicQ(f_in);
  std::vector<std::pair<PolyQ, int>> res;
  PolyQ rem = f;
  int deg = f.degree();
  PolyQ xpol = PolyQ::X(q);
  for (int k = 1; k <= deg && rem.degree() > 0; k++) {
    // x^{q^k} - x modulo rem: the product of the irreducible factors of
    // degree dividing k, the smaller ones being already removed
    PolyQ y = xpol;
    for (int i = 0; i < k; i++) {
      y = PolyPowModQ(y, q, rem);
    }
    PolyQ g = PolyGcdQ(y - xpol, rem);
    if (g.degree() <= 0) {
      continue;
    }
    for (auto const &fac : EqualDegreeFactorization(g, k)) {
      int mult = 0;
      while (true) {
        std::pair<PolyQ, PolyQ> qr = PolyDivModQ(rem, fac);
        if (!qr.second.is_zero()) {
          break;
        }
        rem = qr.first;
        mult++;
      }
      res.push_back({fac, mult});
    }
  }
  return res;
}

// Rabin's test of irreducibility over F_p.
inline bool IsIrreducibleModP(PolyQ const &f) {
  T_small_prime p = f.q;
  int d = f.degree();
  if (d <= 0) {
    return false;
  }
  if (d == 1) {
    return true;
  }
  PolyQ x = PolyQ::X(p);
  auto frob_power = [&](int e) -> PolyQ {
    PolyQ y = x;
    for (int i = 0; i < e; i++) {
      y = PolyPowModQ(y, p, f);
    }
    return y;
  };
  // f divides x^{p^d} - x
  if (!((frob_power(d) - x) % f).is_zero()) {
    return false;
  }
  // gcd(x^{p^{d/r}} - x, f) = 1 for the primes r dividing d
  for (int r = 2; r <= d; r++) {
    if (d % r != 0) {
      continue;
    }
    bool is_prime = true;
    for (int s = 2; s * s <= r; s++) {
      if (r % s == 0) {
        is_prime = false;
      }
    }
    if (!is_prime) {
      continue;
    }
    if (PolyGcdQ(frob_power(d / r) - x, f).degree() > 0) {
      return false;
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_GENERA_SCALARPOLYNOMIALMODQ_H_
// clang-format on
