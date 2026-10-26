// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALDENSITY_H_
#define SRC_GENERA_LOCALDENSITY_H_

// clang-format off
#include "SiegelSpecialValues.h"
#include <cstdint>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_LOCAL_DENSITY
#endif

/*
  Representation masses from local densities.

  For a positive definite form S of rank m and a genus G of rank n, with
  n - m >= 1, Siegel's theorem reads

      M_S / mass(G) = c * alpha_oo(S, T) * prod_p alpha_p(S, T),

      M_S = sum_{L in G} r(S, L) / |Aut L|,   r(S, L) = #{X : X^t L X = S},

  with T any lattice of G. The conventions are those of the Sage prototype
  in Paper_strongly_perfect_lattices/scripts (density_harness.py and
  friends), where every factor was validated against exact representation
  masses of completely enumerated genera: S and T are EVEN bilinear Gram
  matrices, the local density at p counts solutions of X^t T X = S modulo
  p^k off the diagonal and 2 p^k on it, and with that normalisation the
  archimedean constant is c = 2 for n - m >= 2 and c = 1 for n - m = 1
  (Siegel's epsilon; checked on the certified masses of ranks 12 and 13).

  This file holds the part of the product that does not depend on the
  p-adic structure beyond the determinants: the archimedean factor times
  the product over ALL the good primes p (those not dividing 2 det S det T).
  At a good prime both forms are unimodular and alpha_p is the number of
  embeddings over F_p normalised, a rational function of p that factors
  into terms (1 - chi_D(p) p^{-s})^e; the product over all p of such terms
  is 1 / L(s, chi_D)^e, and the transcendental parts of the L-values cancel
  against alpha_oo.
 */

/*
  One factor (1 - chi_D(p) p^{-s})^e of the good-prime density; D = 1
  stands for the trivial character.
 */
struct EulerTerm {
  int s;
  int64_t D;
  int e;
};

/*
  The Euler factors of alpha_p(S, T) at a good prime, for S of rank m and
  determinant detS and T of rank n and determinant detT. By Witt, O(V_T)
  acts transitively on the embeddings of the space of S, with stabiliser
  O(W), W the orthogonal complement, so alpha_p = |O(V_T)| / |O(W)| / p^d;
  the classical orders of the orthogonal groups over F_p factor as below.
  (Derived in the Sage prototype and checked there against the group-order
  formula at many primes, for both parities of n.)
 */
inline std::vector<EulerTerm> GoodPrimeEulerTerms(int n, int m, int64_t detS,
                                                  int64_t detT) {
  std::vector<EulerTerm> terms;
  if (n % 2 == 0) {
    int r = n / 2;
    terms.push_back({r, FundamentalDiscriminant(SignedPower(r, detT)), 1});
    if (m % 2 == 0) {
      int rp = (n - m) / 2;
      int64_t DW = FundamentalDiscriminant(SignedPower(rp, detS * detT));
      for (int i = rp; i < r; i++) {
        terms.push_back({2 * i, 1, 1});
      }
      terms.push_back({rp, DW, -1});
    } else {
      int rp = (n - m - 1) / 2;
      for (int i = rp + 1; i < r; i++) {
        terms.push_back({2 * i, 1, 1});
      }
    }
    return terms;
  }
  int r = (n - 1) / 2;
  if (m % 2 == 1) {
    int rp = (n - m) / 2;
    int64_t DW = FundamentalDiscriminant(SignedPower(rp, detS * detT));
    for (int i = rp; i <= r; i++) {
      terms.push_back({2 * i, 1, 1});
    }
    terms.push_back({rp, DW, -1});
  } else {
    int rp = (n - m - 1) / 2;
    for (int i = rp + 1; i <= r; i++) {
      terms.push_back({2 * i, 1, 1});
    }
  }
  return terms;
}

// prod (1 - chi_D(p) p^{-s})^e at one prime, exactly.
template <typename T>
T EulerValueAtPrime(std::vector<EulerTerm> const &terms, int64_t p) {
  T val(1);
  T p_T = ScalarFromInt64<T>(p);
  for (auto const &term : terms) {
    int chi = QuadraticCharacter(term.D, p);
    T ps(1);
    for (int i = 0; i < term.s; i++) {
      ps *= p_T;
    }
    T factor = T(1) - T(chi) / ps;
    if (term.e >= 0) {
      for (int i = 0; i < term.e; i++) {
        val *= factor;
      }
    } else {
      for (int i = 0; i < -term.e; i++) {
        val /= factor;
      }
    }
  }
  return val;
}

// The primes dividing 2 det S det T, in increasing order.
inline std::vector<int64_t> BadPrimes(int64_t detS, int64_t detT) {
  int64_t N = 2 * detS * detT;
  std::map<int64_t, size_t> fact = FactorsIntMap<int64_t>(N);
  std::vector<int64_t> list;
  for (auto &kv : fact) {
    list.push_back(kv.first);
  }
  return list;
}

/*
  The archimedean factor, without the constant c:

      alpha_oo(S, T) = pi^{m(2n-m+1)/4} / prod_{j=n-m+1}^{n} Gamma(j/2)
                       * det(S)^{(n-m-1)/2} * det(T)^{-m/2}.
 */
template <typename T>
PiRadical<T> ArchimedeanFactor(int n, int m, int64_t detS, int64_t detT) {
  PiRadical<T> val = {T(1), T(m * (2 * n - m + 1)) / T(4), T(1)};
  for (int j = n - m + 1; j <= n; j++) {
    val = val * PiRadicalInverse(GammaHalfInteger<T>(j));
  }
  val = val * RationalHalfPower<T>(ScalarFromInt64<T>(detS), n - m - 1);
  val = val * RationalHalfPower<T>(ScalarFromInt64<T>(detT), -m);
  return val;
}

/*
  alpha_oo(S, T) times the product of alpha_p(S, T) over the good primes,
  without the constant c: the infinite product of the Euler terms over all
  primes, as L-values, divided by their values at the bad primes.
 */
template <typename T>
T GenericDensityPart(int n, int m, int64_t detS, int64_t detT) {
  if (n - m < 1) {
    std::cerr << "LDEN: GenericDensityPart needs n - m >= 1, got n=" << n
              << " m=" << m << "\n";
    throw TerminalException{1};
  }
  std::vector<EulerTerm> terms = GoodPrimeEulerTerms(n, m, detS, detT);
  PiRadical<T> val = ArchimedeanFactor<T>(n, m, detS, detT);
  for (auto const &term : terms) {
    PiRadical<T> L = LQuadraticValue<T>(term.s, term.D);
    val = val * PiRadicalPower(L, -term.e);
  }
  T rat = PiRadicalToRational(val);
  for (int64_t p : BadPrimes(detS, detT)) {
    rat /= EulerValueAtPrime<T>(terms, p);
  }
  return rat;
}

// clang-format off
#endif  // SRC_GENERA_LOCALDENSITY_H_
// clang-format on
