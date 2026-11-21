// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_FINITEFIELD_H_
#define SRC_GENERA_FINITEFIELD_H_

// clang-format off
#include "PolynomialModQ.h"
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

  An element is a polynomial of degree < f in x (a PolyQ, coefficients in
  [0, p) by increasing degree), kept trimmed, so that 0 is the empty
  vector and equality is that of the vectors. The operations accept any
  PolyQ (for instance a vector of length f with zeros at the end) and
  return the trimmed one.

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

using FqElt = PolyQ;

inline FiniteField MakeFiniteField(T_small_prime p, PolyQ const &mu_in) {
  PolyQ mu = PolyReduce(mu_in, p);
  PolyTrim(mu);
  if (mu.size() < 2) {
    std::cerr << "FQ: the modulus has degree < 1\n";
    throw TerminalException{1};
  }
  mu = PolyMonicQ(mu, p);
  int f = static_cast<int>(mu.size()) - 1;
  uint64_t q = 1;
  for (int i = 0; i < f && q != 0; i++) {
    q = (q > kMaxResidueModulus / p) ? 0 : q * p;
  }
  return {p, f, mu, q};
}

// F_p as F_p[x] / (x)
inline FiniteField MakePrimeField(T_small_prime p) {
  return MakeFiniteField(p, PolyQ{0, 1});
}

inline FqElt FqNormalize(FiniteField const &F, FqElt a) {
  a = PolyReduce(a, F.p);
  PolyTrim(a);
  if (static_cast<int>(a.size()) > F.f) {
    a = PolyDivModQ(a, F.mu, F.p).second;
    PolyTrim(a);
  }
  return a;
}

inline bool FqIsZero(FqElt const &a) {
  return std::all_of(a.begin(), a.end(), [](T_residue x) { return x == 0; });
}

inline FqElt FqOne() { return FqElt{1}; }

inline FqElt FqAdd(FiniteField const &F, FqElt const &a, FqElt const &b) {
  FqElt r(std::max(a.size(), b.size()), 0);
  for (size_t i = 0; i < r.size(); i++) {
    T_residue x = (i < a.size()) ? a[i] % F.p : 0;
    T_residue y = (i < b.size()) ? b[i] % F.p : 0;
    r[i] = ResAdd(x, y, F.p);
  }
  return FqNormalize(F, r);
}

inline FqElt FqNeg(FiniteField const &F, FqElt const &a) {
  FqElt r(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    r[i] = ResNeg(a[i] % F.p, F.p);
  }
  return FqNormalize(F, r);
}

inline FqElt FqSub(FiniteField const &F, FqElt const &a, FqElt const &b) {
  return FqAdd(F, a, FqNeg(F, b));
}

inline FqElt FqMul(FiniteField const &F, FqElt const &a, FqElt const &b) {
  return FqNormalize(F, PolyMulQ(FqNormalize(F, a), FqNormalize(F, b), F.p));
}

inline FqElt FqInv(FiniteField const &F, FqElt const &a) {
  FqElt a2 = FqNormalize(F, a);
  if (a2.empty()) {
    std::cerr << "FQ: the inverse of 0\n";
    throw TerminalException{1};
  }
  return FqNormalize(F, PolyInverseModQ(a2, F.mu, F.p));
}

inline FqElt FqPow(FiniteField const &F, FqElt const &a, uint64_t e) {
  return FqNormalize(F, PolyPowModQ(FqNormalize(F, a), e, F.mu, F.p));
}

// a^(p^k)
inline FqElt FqFrobenius(FiniteField const &F, FqElt a, int k) {
  a = FqNormalize(F, a);
  for (int i = 0; i < k; i++) {
    a = FqPow(F, a, F.p);
  }
  return a;
}

// The trace to F_p, the sum of the f conjugates.
inline T_residue FqTrace(FiniteField const &F, FqElt const &a) {
  FqElt t;
  FqElt x = FqNormalize(F, a);
  for (int i = 0; i < F.f; i++) {
    t = FqAdd(F, t, x);
    x = FqFrobenius(F, x, 1);
  }
  if (t.size() > 1) {
    std::cerr << "FQ: the trace is not in the prime field\n";
    throw TerminalException{1};
  }
  return t.empty() ? 0 : t[0];
}

inline void FqCheckCodes(FiniteField const &F) {
  if (F.q == 0) {
    std::cerr << "FQ: the field F_" << F.p << "^" << F.f
              << " is too large for the codes of its elements\n";
    throw TerminalException{1};
  }
}

inline uint64_t FqEncode(FiniteField const &F, FqElt const &a_in) {
  FqCheckCodes(F);
  FqElt a = FqNormalize(F, a_in);
  uint64_t c = 0;
  for (size_t i = a.size(); i-- > 0;) {
    c = c * F.p + a[i];
  }
  return c;
}

inline FqElt FqDecode(FiniteField const &F, uint64_t c) {
  FqCheckCodes(F);
  FqElt a(F.f);
  for (int i = 0; i < F.f; i++) {
    a[i] = c % F.p;
    c /= F.p;
  }
  PolyTrim(a);
  return a;
}

// The operations on the codes, without allocation for f <= 2 (mu = x^2 +
// c_1 x + c_0, x^2 = -c_1 x - c_0).
inline uint64_t FqCodeAdd(FiniteField const &F, uint64_t x, uint64_t y) {
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

inline uint64_t FqCodeMul(FiniteField const &F, uint64_t x, uint64_t y) {
  if (F.f == 1) {
    return ResMul(x, y, F.p);
  }
  if (F.f == 2) {
    T_residue p = F.p;
    T_residue a1 = x % p, b1 = x / p, a2 = y % p, b2 = y / p;
    T_residue bb = ResMul(b1, b2, p);
    T_residue a = ResSub(ResMul(a1, a2, p), ResMul(F.mu[0], bb, p), p);
    T_residue b = ResSub(ResAdd(ResMul(a1, b2, p), ResMul(a2, b1, p), p),
                         ResMul(F.mu[1], bb, p), p);
    return a + p * b;
  }
  return FqEncode(F, FqMul(F, FqDecode(F, x), FqDecode(F, y)));
}

inline uint64_t FqCodeInv(FiniteField const &F, uint64_t x) {
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
        ResAdd(ResSub(ResMul(a, a, p), ResMul(F.mu[1], ResMul(a, b, p), p), p),
               ResMul(F.mu[0], ResMul(b, b, p), p), p);
    T_residue ni = ResInv(norm, p);
    T_residue ca = ResMul(ResSub(a, ResMul(F.mu[1], b, p), p), ni, p);
    T_residue cb = ResMul(ResNeg(b, p), ni, p);
    return ca + p * cb;
  }
  return FqEncode(F, FqInv(F, FqDecode(F, x)));
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
#endif  // SRC_GENERA_FINITEFIELD_H_
// clang-format on
