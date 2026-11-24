// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SCALARFINITEFIELD_H_
#define SRC_GENERA_SCALARFINITEFIELD_H_

// clang-format off
#include "ScalarPolynomialModQ.h"
#include "ResidueArithmetic.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
// clang-format on

/*
  The finite fields F_q = F_p[x] / (mu), mu monic irreducible of degree f,
  q = p^f: the scalar type of the residue fields of the number fields and
  of the eigenvalues of the isometries modulo p. The constructions over
  F_q (linear algebra, projective spaces, forms) are in the headers that
  use them.

  An element (FqElt) is a polynomial of degree < f in x, its coefficients
  in [0, p) by increasing degree, trimmed, so that 0 has no coefficient
  and equality is that of the coefficients. It does not carry the field:
  it is made and operated on through the FiniteField F (FqFromPoly,
  FqFromCoeffs, FqAdd, ...), which reduces modulo p and mu.

  The codes: an element a_0 + a_1 x + ... is coded by the integer a_0 + a_1
  p + ... in [0, q), when q <= kMaxResidueModulus (q = 0 otherwise, and
  the codes are refused). The operations on the codes do not allocate for
  f <= 2; for a small q the tables FqTables give them by lookup.
 */

struct FiniteField {
  T_small_prime p;
  int f;
  PolyQ mu;
  // p^f, or 0 beyond kMaxResidueModulus
  uint64_t q;
};

struct FqElt {
  FqElt() = default;

  bool is_zero() const { return c.empty(); }
  size_t size() const { return c.size(); }
  T_residue operator[](size_t i) const { return c[i]; }
  // the coefficient of x^i, 0 beyond the degree
  T_residue coeff(size_t i) const { return (i < c.size()) ? c[i] : 0; }
  auto begin() const { return c.begin(); }
  auto end() const { return c.end(); }
  bool operator==(FqElt const &b) const { return c == b.c; }

private:
  std::vector<T_residue> c;

  friend FqElt FqFromPoly(FiniteField const &F, PolyQ const &a);
};

inline FiniteField MakeFiniteField(T_small_prime p, PolyQ const &mu_in) {
  if (mu_in.q != p) {
    std::cerr << "FQ: the modulus is not a polynomial modulo " << p << "\n";
    throw TerminalException{1};
  }
  if (mu_in.degree() < 1) {
    std::cerr << "FQ: the modulus has degree < 1\n";
    throw TerminalException{1};
  }
  PolyQ mu = PolyMonicQ(mu_in);
  int f = mu.degree();
  uint64_t q = 1;
  for (int i = 0; i < f && q != 0; i++) {
    q = (q > kMaxResidueModulus / p) ? 0 : q * p;
  }
  return {p, f, mu, q};
}

// F_p as F_p[x] / (x)
inline FiniteField MakePrimeField(T_small_prime p) {
  return MakeFiniteField(p, PolyQ::X(p));
}

// the class of a modulo mu
inline FqElt FqFromPoly(FiniteField const &F, PolyQ const &a) {
  if (a.q != F.p) {
    std::cerr << "FQ: a polynomial modulo " << a.q << " in F_" << F.p
              << "^" << F.f << "\n";
    throw TerminalException{1};
  }
  PolyQ r = a;
  if (r.degree() >= F.f) {
    r = r % F.mu;
  }
  FqElt res;
  res.c = r.coeffs();
  return res;
}

// the class of sum_i a_i x^i, the a_i taken modulo p
template <typename Tcoef>
inline FqElt FqFromCoeffs(FiniteField const &F, std::vector<Tcoef> const &a) {
  std::vector<T_residue> r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResModSigned(static_cast<int64_t>(a[i]), F.p);
  }
  return FqFromPoly(F, PolyQ(F.p, std::move(r)));
}

// the polynomial of degree < f of a
inline PolyQ FqToPoly(FiniteField const &F, FqElt const &a) {
  return PolyQ(F.p, std::vector<T_residue>(a.begin(), a.end()));
}

inline FqElt FqOne(FiniteField const &F) {
  return FqFromPoly(F, PolyQ(F.p, {1}));
}

// x^k
inline FqElt FqXPower(FiniteField const &F, int k) {
  std::vector<T_residue> r(k + 1, 0);
  r[k] = 1;
  return FqFromPoly(F, PolyQ(F.p, std::move(r)));
}

inline FqElt FqAdd(FiniteField const &F, FqElt const &a, FqElt const &b) {
  return FqFromPoly(F, FqToPoly(F, a) + FqToPoly(F, b));
}

inline FqElt FqNeg(FiniteField const &F, FqElt const &a) {
  return FqFromPoly(F, -FqToPoly(F, a));
}

inline FqElt FqSub(FiniteField const &F, FqElt const &a, FqElt const &b) {
  return FqFromPoly(F, FqToPoly(F, a) - FqToPoly(F, b));
}

inline FqElt FqMul(FiniteField const &F, FqElt const &a, FqElt const &b) {
  return FqFromPoly(F, FqToPoly(F, a) * FqToPoly(F, b));
}

inline FqElt FqInv(FiniteField const &F, FqElt const &a) {
  if (a.is_zero()) {
    std::cerr << "FQ: the inverse of 0\n";
    throw TerminalException{1};
  }
  FqElt res = FqFromPoly(F, PolyInverseModQ(FqToPoly(F, a), F.mu));
#ifdef SANITY_CHECK_FINITE_FIELD
  if (FqMul(F, a, res) != FqOne(F)) {
    std::cerr << "FQ: the inverse is wrong\n";
    throw TerminalException{1};
  }
#endif
  return res;
}

inline FqElt FqPow(FiniteField const &F, FqElt const &a, uint64_t e) {
  return FqFromPoly(F, PolyPowModQ(FqToPoly(F, a), e, F.mu));
}

// a^(p^k)
inline FqElt FqFrobenius(FiniteField const &F, FqElt a, int k) {
  for (int i = 0; i < k; i++) {
    a = FqPow(F, a, F.p);
  }
  return a;
}

// The trace to F_p, the sum of the f conjugates.
inline T_residue FqTrace(FiniteField const &F, FqElt const &a) {
  FqElt t;
  FqElt x = a;
  for (int i = 0; i < F.f; i++) {
    t = FqAdd(F, t, x);
    x = FqFrobenius(F, x, 1);
  }
  if (t.size() > 1) {
    std::cerr << "FQ: the trace is not in the prime field\n";
    throw TerminalException{1};
  }
  return t.coeff(0);
}

/*
  The quadratic character of u != 0, p odd: that of its norm to F_p, the
  product of its conjugates u^{p^i}, i < f (the norm is onto and takes the
  squares to the squares), with no exponent beyond p.
 */
inline int FqQuadraticCharacter(FiniteField const &F, FqElt const &u) {
  if (u.is_zero()) {
    std::cerr << "FQ: the quadratic character of 0\n";
    throw TerminalException{1};
  }
  FqElt prod = FqOne(F);
  FqElt conj = u;
  for (int i = 0; i < F.f; i++) {
    prod = FqMul(F, prod, conj);
    conj = FqFrobenius(F, conj, 1);
  }
  if (prod.size() != 1) {
    std::cerr << "FQ: the norm is not in F_p\n";
    throw TerminalException{1};
  }
  return ResKronecker(prod[0], F.p);
}

inline void FqCheckCodes(FiniteField const &F) {
  if (F.q == 0) {
    std::cerr << "FQ: the field F_" << F.p << "^" << F.f
              << " is too large for the codes of its elements\n";
    throw TerminalException{1};
  }
}

inline uint64_t FqEncode(FiniteField const &F, FqElt const &a) {
  FqCheckCodes(F);
  uint64_t c = 0;
  for (size_t i = a.size(); i-- > 0;) {
    c = c * F.p + a[i];
  }
  return c;
}

inline FqElt FqDecode(FiniteField const &F, uint64_t c) {
  FqCheckCodes(F);
  std::vector<T_residue> a(F.f);
  for (int i = 0; i < F.f; i++) {
    a[i] = c % F.p;
    c /= F.p;
  }
  return FqFromPoly(F, PolyQ(F.p, std::move(a)));
}

// The operations on the codes, without allocation for f <= 2 (mu = x^2 +
// c_1 x + c_0, x^2 = -c_1 x - c_0).
inline uint64_t FqCodeAddFast(FiniteField const &F, uint64_t x, uint64_t y) {
  if (F.f == 1) {
    return ResAdd(x, y, F.p);
  }
  if (F.f == 2) {
    T_residue a = ResAdd(x % F.p, y % F.p, F.p);
    T_residue b = ResAdd(x / F.p, y / F.p, F.p);
    return a + F.p * b;
  }
  return FqEncode(F, FqAdd(F, FqDecode(F, x), FqDecode(F, y)));
}

inline uint64_t FqCodeMulFast(FiniteField const &F, uint64_t x, uint64_t y) {
  if (F.f == 1) {
    return ResMul(x, y, F.p);
  }
  if (F.f == 2) {
    T_residue p = F.p;
    T_residue a1 = x % p, b1 = x / p, a2 = y % p, b2 = y / p;
    T_residue bb = ResMul(b1, b2, p);
    T_residue a = ResSub(ResMul(a1, a2, p), ResMul(F.mu.coeff(0), bb, p), p);
    T_residue b = ResSub(ResAdd(ResMul(a1, b2, p), ResMul(a2, b1, p), p),
                         ResMul(F.mu.coeff(1), bb, p), p);
    return a + p * b;
  }
  return FqEncode(F, FqMul(F, FqDecode(F, x), FqDecode(F, y)));
}

inline uint64_t FqCodeInvFast(FiniteField const &F, uint64_t x) {
  if (x == 0) {
    std::cerr << "FQ: the inverse of 0\n";
    throw TerminalException{1};
  }
  if (F.f == 1) {
    return ResInv(x, F.p);
  }
  if (F.f == 2) {
    // (a + b x)(a + b x') = a^2 - c_1 a b + c_0 b^2, x' = -c_1 - x
    T_residue p = F.p;
    T_residue a = x % p, b = x / p;
    T_residue norm =
        ResAdd(ResSub(ResMul(a, a, p), ResMul(F.mu.coeff(1), ResMul(a, b, p), p), p),
               ResMul(F.mu.coeff(0), ResMul(b, b, p), p), p);
    T_residue ni = ResInv(norm, p);
    T_residue ca = ResMul(ResSub(a, ResMul(F.mu.coeff(1), b, p), p), ni, p);
    T_residue cb = ResMul(ResNeg(b, p), ni, p);
    return ca + p * cb;
  }
  return FqEncode(F, FqInv(F, FqDecode(F, x)));
}

/*
  The operations on the codes, checked against the operations on the
  elements under SANITY_CHECK_FINITE_FIELD (the formulas of f <= 2).
 */
inline uint64_t FqCodeAdd(FiniteField const &F, uint64_t x, uint64_t y) {
  uint64_t z = FqCodeAddFast(F, x, y);
#ifdef SANITY_CHECK_FINITE_FIELD
  if (z != FqEncode(F, FqAdd(F, FqDecode(F, x), FqDecode(F, y)))) {
    std::cerr << "FQ: the sum of the codes is wrong\n";
    throw TerminalException{1};
  }
#endif
  return z;
}

inline uint64_t FqCodeMul(FiniteField const &F, uint64_t x, uint64_t y) {
  uint64_t z = FqCodeMulFast(F, x, y);
#ifdef SANITY_CHECK_FINITE_FIELD
  if (z != FqEncode(F, FqMul(F, FqDecode(F, x), FqDecode(F, y)))) {
    std::cerr << "FQ: the product of the codes is wrong\n";
    throw TerminalException{1};
  }
#endif
  return z;
}

inline uint64_t FqCodeInv(FiniteField const &F, uint64_t x) {
  uint64_t z = FqCodeInvFast(F, x);
#ifdef SANITY_CHECK_FINITE_FIELD
  if (z != FqEncode(F, FqInv(F, FqDecode(F, x)))) {
    std::cerr << "FQ: the inverse of the code is wrong\n";
    throw TerminalException{1};
  }
#endif
  return z;
}

/*
  F_q for a small q as tables on the codes: the products of the orbit
  enumerations are read there, without the polynomial arithmetic and the
  allocations of every product (the orbits of P^11(F_4), 5.6 10^6 points,
  under the centralizer of (A2^12, w) took about ten minutes per class).
 */
struct FqTables {
  uint64_t q;
  std::vector<int> add;
  std::vector<int> mul;
  std::vector<int> inv;
};

inline size_t const fq_tables_max_q = 4096;

inline FqTables BuildFqTables(FiniteField const &F) {
  FqCheckCodes(F);
  uint64_t q = F.q;
  FqTables tb{q, std::vector<int>(q * q), std::vector<int>(q * q),
              std::vector<int>(q, 0)};
  for (uint64_t a = 0; a < q; a++) {
    for (uint64_t b = 0; b < q; b++) {
      tb.add[a * q + b] = static_cast<int>(FqCodeAdd(F, a, b));
      int c = static_cast<int>(FqCodeMul(F, a, b));
      tb.mul[a * q + b] = c;
      if (c == 1) {
        tb.inv[a] = static_cast<int>(b);
      }
    }
  }
  return tb;
}

// clang-format off
#endif  // SRC_GENERA_SCALARFINITEFIELD_H_
// clang-format on
