// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ABELIANQUADRATICMASS_H_
#define SRC_GENERA_ABELIANQUADRATICMASS_H_

// clang-format off
#include "AbelianHermitianMass.h"
#include <optional>
#include <set>
#include <vector>
// clang-format on

/*
  Masses of quadratic lattices over a totally real abelian field K, the
  subfield of Q(zeta_n) fixed by a subgroup H containing -1 (Hecke.jl,
  QuadForm/Quad/MassQuad.jl, after Siegel and Kirschmer). The standard mass
  of a lattice of rank m, r = floor(m / 2), is

      2^{-[K:Q] r} prod_{i odd, i <= m - 2} zeta_K(-i)            (m odd)
      2^{-[K:Q] r} prod_{i odd, i <= m - 3} zeta_K(-i) L_r        (m even)

  with L_r = zeta_K(1 - r) when the discriminant (-1)^{m(m-1)/2} det is a
  square in K, and L(1 - r, chi_{E/K}) for E = K(sqrt disc) otherwise. As
  K is abelian, zeta_K(1 - k) is the product of the L(1 - k, chi) over the
  Dirichlet characters trivial on H; when E is abelian too, given as the
  subfield of Q(zeta_n) fixed by H_E (of index 2 in H), L(1 - r, chi_{E/K})
  is the product over the characters trivial on H_E but not on H.
 */

inline mpq_class DedekindZetaAbelian(AbelianCMField const &F, int k) {
  // zeta_K(1 - k): all the characters of K, even as K is real
  return ProductDirichletLValues(F.n, k, 0, F.H);
}

// The quadratic field data: a totally real field as (n, H) with -1 in H.
inline AbelianCMField MakeAbelianRealField(int64_t n,
                                           std::vector<int64_t> const &gens) {
  std::vector<int64_t> g = gens;
  if (n > 2) {
    g.push_back(n - 1);
  }
  // MakeAbelianCMField refuses -1: build the subgroup directly
  std::set<int64_t> Hs{1 % std::max<int64_t>(n, 1)};
  std::vector<int64_t> queue{1 % std::max<int64_t>(n, 1)};
  while (!queue.empty()) {
    int64_t a = queue.back();
    queue.pop_back();
    for (int64_t x : g) {
      int64_t b = ((a * (x % n)) % n + n) % n;
      if (Hs.insert(b).second) {
        queue.push_back(b);
      }
    }
  }
  return {n, std::vector<int64_t>(Hs.begin(), Hs.end())};
}

/*
  The standard mass; for m even with a non-square discriminant, FE gives
  E = K(sqrt disc) as the subfield of Q(zeta_n) fixed by H_E (same n).
  Nothing when E is needed and not given.
 */
inline std::optional<mpq_class>
StandardQuadraticMassAbelian(AbelianCMField const &FK, int m, bool disc_square,
                             std::optional<AbelianCMField> const &FE) {
  if (m == 0) {
    return mpq_class(1);
  }
  if (m == 1) {
    return mpq_class(1, 2);
  }
  int64_t deg = AbelianFieldDegree(FK);
  int r = m / 2;
  mpq_class res(1);
  for (int64_t i = 0; i < deg * r; i++) {
    res /= 2;
  }
  int top = (m % 2 == 1) ? m - 2 : m - 3;
  for (int i = 1; i <= top; i += 2) {
    res *= DedekindZetaAbelian(FK, i + 1);
  }
  if (m % 2 == 0) {
    if (disc_square) {
      res *= DedekindZetaAbelian(FK, r);
    } else {
      if (!FE || FE->n != FK.n) {
        return {};
      }
      // the characters of E not of K, whatever their parity
      mpq_class prod(1);
      for (int parity = 0; parity <= 1; parity++) {
        prod *= ProductDirichletLValues(FK.n, r, parity, FE->H, FK.H);
      }
      res *= prod;
    }
  }
  return (res < 0) ? mpq_class(-res) : res;
}

// clang-format off
#endif  // SRC_GENERA_ABELIANQUADRATICMASS_H_
// clang-format on
