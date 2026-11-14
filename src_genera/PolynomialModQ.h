// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_POLYNOMIALMODQ_H_
#define SRC_GENERA_POLYNOMIALMODQ_H_

// clang-format off
#include "GeneraArithmetic.h"
#include "ResidueArithmetic.h"
#include <cstdint>
#include <random>
#include <utility>
#include <vector>
// clang-format on

/*
  Polynomials with coefficients modulo a prime q (PolyQ, the coefficients
  in increasing degree) and with integer coefficients (PolyZ): products,
  division, gcd and inverses modulo q, powers modulo a polynomial, and the
  factorization modulo q by distinct and equal degrees.
 */

using PolyQ = std::vector<T_residue>;

inline void PolyTrim(PolyQ &a) {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
}

inline PolyQ PolyReduce(PolyQ a, T_small_prime q) {
  for (auto &x : a) {
    x %= q;
  }
  PolyTrim(a);
  return a;
}

// a polynomial with integer coefficients modulo q
inline PolyQ PolyReduceSigned(std::vector<int64_t> const &a, T_small_prime q) {
  PolyQ r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResModSigned(a[i], q);
  }
  PolyTrim(r);
  return r;
}

inline PolyQ PolyMulQ(PolyQ const &a, PolyQ const &b, T_small_prime q) {
  if (a.empty() || b.empty()) {
    return {};
  }
  PolyQ c(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] = ResAdd(c[i + j], ResMul(a[i], b[j], q), q);
    }
  }
  PolyTrim(c);
  return c;
}

inline PolyQ PolySubQ(PolyQ a, PolyQ const &b, T_small_prime q) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (size_t i = 0; i < b.size(); i++) {
    a[i] = ResSub(a[i], b[i], q);
  }
  PolyTrim(a);
  return a;
}

// quotient and remainder modulo q, b nonzero
inline std::pair<PolyQ, PolyQ> PolyDivModQ(PolyQ a, PolyQ const &b, T_small_prime q) {
  PolyTrim(a);
  if (a.size() < b.size()) {
    return {{}, a};
  }
  T_residue inv = ResInv(b.back(), q);
  PolyQ quo(a.size() - b.size() + 1, 0);
  for (size_t k = a.size() - b.size() + 1; k-- > 0;) {
    T_residue c = ResMul(a[k + b.size() - 1], inv, q);
    quo[k] = c;
    for (size_t j = 0; j < b.size(); j++) {
      a[k + j] = ResSub(a[k + j], ResMul(c, b[j], q), q);
    }
  }
  PolyTrim(a);
  PolyTrim(quo);
  return {quo, a};
}

inline PolyQ PolyMonicQ(PolyQ a, T_small_prime q) {
  if (a.empty()) {
    return a;
  }
  T_residue inv = ResInv(a.back(), q);
  for (auto &x : a) {
    x = ResMul(x, inv, q);
  }
  return a;
}

inline PolyQ PolyGcdQ(PolyQ a, PolyQ b, T_small_prime q) {
  PolyTrim(a);
  PolyTrim(b);
  while (!b.empty()) {
    PolyQ r = PolyDivModQ(a, b, q).second;
    a = b;
    b = r;
  }
  return PolyMonicQ(a, q);
}

// u with u a = gcd modulo b (a and b coprime): the inverse of a modulo b
inline PolyQ PolyInverseModQ(PolyQ const &a, PolyQ const &b, T_small_prime q) {
  PolyQ r0 = b, r1 = PolyDivModQ(a, b, q).second;
  PolyQ s0 = {}, s1 = {1};
  while (!r1.empty()) {
    std::pair<PolyQ, PolyQ> qr = PolyDivModQ(r0, r1, q);
    PolyQ s2 = PolySubQ(s0, PolyMulQ(qr.first, s1, q), q);
    r0 = r1;
    r1 = qr.second;
    s0 = s1;
    s1 = s2;
  }
  // r0 is a nonzero constant
  if (r0.size() != 1) {
    std::cerr << "CYCLO: the polynomials are not coprime\n";
    throw TerminalException{1};
  }
  T_residue inv = ResInv(r0[0], q);
  for (auto &x : s0) {
    x = ResMul(x, inv, q);
  }
  return PolyDivModQ(s0, b, q).second;
}

inline PolyQ PolyPowModQ(PolyQ base, uint64_t e, PolyQ const &m, T_small_prime q) {
  PolyQ res = {1};
  base = PolyDivModQ(base, m, q).second;
  while (e > 0) {
    if (e & 1) {
      res = PolyDivModQ(PolyMulQ(res, base, q), m, q).second;
    }
    base = PolyDivModQ(PolyMulQ(base, base, q), m, q).second;
    e >>= 1;
  }
  return res;
}

/*
  The quadratic character of u != 0 in F_q[x] / g, q odd: that of its norm
  to F_q, the product of its conjugates u^{q^i}, i < deg g (the norm is
  onto and takes the squares to the squares), with no exponent beyond q.
 */
inline int ResidueQuadraticCharacter(PolyQ const &u, PolyQ const &g, T_small_prime q) {
  int f = static_cast<int>(g.size()) - 1;
  PolyQ prod = {1};
  PolyQ conj = PolyDivModQ(u, g, q).second;
  for (int i = 0; i < f; i++) {
    prod = PolyDivModQ(PolyMulQ(prod, conj, q), g, q).second;
    conj = PolyPowModQ(conj, q, g, q);
  }
  PolyTrim(prod);
  if (prod.size() != 1) {
    std::cerr << "CYCLO: the norm is not in F_q\n";
    throw TerminalException{1};
  }
  return KroneckerResidue(prod[0], q);
}

/*
  The irreducible factors (monic) of the squarefree g modulo q, all of
  degree f, by the equal degree splitting of Cantor and Zassenhaus (q odd)
  or by the trace map (q = 2). The random polynomials are replaced by a
  deterministic sequence, which splits as well.
 */
inline std::vector<PolyQ> EqualDegreeFactorization(PolyQ const &g, int f,
                                                   T_small_prime q) {
  int deg = static_cast<int>(g.size()) - 1;
  if (deg == f) {
    return {PolyMonicQ(g, q)};
  }
  for (uint64_t c = 1;; c++) {
    // r = x^{1 + c mod deg} + c, a varied family of polynomials
    PolyQ r(2 + static_cast<size_t>(c % deg), 0);
    r.back() = 1;
    r[0] = c % q;
    r = PolyReduce(r, q);
    PolyQ h;
    if (q == 2) {
      PolyQ t = PolyDivModQ(r, g, q).second;
      h = t;
      for (int i = 1; i < f; i++) {
        t = PolyDivModQ(PolyMulQ(t, t, q), g, q).second;
        h = PolySubQ(h, PolySubQ({}, t, q), q);
      }
    } else {
      // r^{(q^f - 1)/2} = prod_{i < f} (r^{(q-1)/2})^{q^i}, which keeps the
      // exponents within 64 bits
      PolyQ t = PolyPowModQ(r, (q - 1) / 2, g, q);
      h = {1};
      for (int i = 0; i < f; i++) {
        h = PolyDivModQ(PolyMulQ(h, t, q), g, q).second;
        t = PolyPowModQ(t, q, g, q);
      }
      h = PolySubQ(h, {1}, q);
    }
    PolyQ dgc = PolyGcdQ(g, h, q);
    int dd = static_cast<int>(dgc.size()) - 1;
    if (dd > 0 && dd < deg) {
      std::vector<PolyQ> res = EqualDegreeFactorization(dgc, f, q);
      std::vector<PolyQ> res2 =
          EqualDegreeFactorization(PolyDivModQ(g, dgc, q).first, f, q);
      res.insert(res.end(), res2.begin(), res2.end());
      return res;
    }
    if (c > 100000) {
      std::cerr << "CYCLO: the equal degree splitting failed\n";
      throw TerminalException{1};
    }
  }
}

inline PolyQ ReciprocalPolyQ(PolyQ const &g, T_small_prime q) {
  PolyQ r(g.rbegin(), g.rend());
  return PolyMonicQ(PolyReduce(r, q), q);
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
inline std::vector<std::pair<PolyQ, int>> FactorModQWithMultiplicity(PolyQ f,
                                                                      T_small_prime q) {
  f = PolyMonicQ(PolyReduce(f, q), q);
  std::vector<std::pair<PolyQ, int>> res;
  PolyQ rem = f;
  int deg = static_cast<int>(f.size()) - 1;
  PolyQ xpol{0, 1};
  for (int k = 1; k <= deg && rem.size() > 1; k++) {
    // x^{q^k} - x modulo rem: the product of the irreducible factors of
    // degree dividing k, the smaller ones being already removed
    PolyQ y = xpol;
    for (int i = 0; i < k; i++) {
      y = PolyPowModQ(y, q, rem, q);
    }
    PolyQ g = PolyGcdQ(PolySubQ(y, xpol, q), rem, q);
    PolyTrim(g);
    if (g.size() <= 1) {
      continue;
    }
    for (auto const &fac : EqualDegreeFactorization(PolyMonicQ(g, q), k, q)) {
      int mult = 0;
      while (true) {
        std::pair<PolyQ, PolyQ> qr = PolyDivModQ(rem, fac, q);
        PolyQ r = qr.second;
        PolyTrim(r);
        if (!r.empty()) {
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
inline bool IsIrreducibleModP(PolyQ const &f, T_small_prime p) {
  int d = static_cast<int>(f.size()) - 1;
  if (d <= 0) {
    return false;
  }
  if (d == 1) {
    return true;
  }
  PolyQ x{0, 1};
  auto frob_power = [&](int e) -> PolyQ {
    PolyQ y = x;
    for (int i = 0; i < e; i++) {
      y = PolyPowModQ(y, p, f, p);
    }
    return y;
  };
  // f divides x^{p^d} - x
  PolyQ y = PolySubQ(frob_power(d), x, p);
  PolyTrim(y);
  if (!y.empty() && !(y.size() == 1 && y[0] == 0)) {
    std::pair<PolyQ, PolyQ> qr = PolyDivModQ(y, f, p);
    PolyQ rem = qr.second;
    PolyTrim(rem);
    if (!rem.empty() && !(rem.size() == 1 && rem[0] == 0)) {
      return false;
    }
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
    PolyQ z = PolySubQ(frob_power(d / r), x, p);
    PolyQ g = PolyGcdQ(z, f, p);
    PolyTrim(g);
    if (g.size() > 1) {
      return false;
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_GENERA_POLYNOMIALMODQ_H_
// clang-format on
