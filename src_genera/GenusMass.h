// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSMASS_H_
#define SRC_GENERA_GENUSMASS_H_

// clang-format off
#include "LocalDensityOverlattice.h"
#include <array>
#include <cstdint>
#include <vector>
// clang-format on

/*
  The mass sum_{L in G} 1 / |Aut L| of the genus G of a positive definite
  even lattice, by the formula of Conway and Sloane (Low-dimensional
  lattices IV: the mass formula, Proc. R. Soc. Lond. A 419 (1988) 259-286):

      mass(G) = std(G) * prod_{p | 2 det} m_p(G) / std_p(G),

  std(G) the standard mass (gamma, zeta and L values), m_p the p-mass read
  off the p-adic symbol (diagonal factors M(species), cross terms, and at
  p = 2 the type factor), std_p the standard p-mass. The code follows the
  sage.quadratic_forms.genera.genus implementation (Genus_Symbol_p_adic_ring
  .mass, ._standard_mass, ._species_list and GenusSymbol_global_ring
  ._standard_mass); the Gram matrix is the bilinear form, as Sage's
  Genus(M).
 */

// The diagonal factor M_p(species) of Conway-Sloane, Table 2.
template <typename T> T DiagonalFactorMp(int species, int64_t p) {
  if (species == 0) {
    return T(1);
  }
  int n = (species > 0) ? species : -species;
  int s = (n + 1) / 2;
  T p_T = ScalarFromInt64<T>(p);
  T mp(2);
  T pk(1);
  for (int k = 1; k < s; k++) {
    pk *= p_T * p_T;
    mp *= T(1) - T(1) / pk;
  }
  if (n % 2 == 0) {
    T ps(1);
    for (int i = 0; i < s; i++) {
      ps *= p_T;
    }
    T sign = (species > 0) ? T(1) : T(-1);
    mp *= T(1) - sign / ps;
  }
  return T(1) / mp;
}

// Sage's local symbol at an odd p: (scale exponent, rank, Legendre symbol
// of the determinant of the component of the Gram matrix itself).
template <typename T, typename Tint>
std::vector<std::array<int, 3>> OddSymbolOfGram(MyMatrix<T> const &Gram,
                                                int64_t p) {
  std::vector<std::array<int, 3>> sym = OddLocalSymbol<T, Tint>(Gram, p);
  int leg2 = LegendreUnit(2, p);
  for (auto &c : sym) {
    if (c[1] % 2 == 1) {
      c[2] *= leg2;
    }
  }
  return sym;
}

inline std::vector<int> SpeciesListOdd(std::vector<std::array<int, 3>> const &sym,
                                       int64_t p) {
  int leg_m1 = (p % 4 == 1) ? 1 : -1;
  std::vector<int> res;
  for (auto const &c : sym) {
    int n = c[1];
    int expected = (n % 4 == 2) ? leg_m1 : 1;
    if (n % 2 == 0 && c[2] != expected) {
      res.push_back(-n);
    } else {
      res.push_back(n);
    }
  }
  return res;
}

// The species at p = 2 (Conway-Sloane, Table 1): a component is bound when
// a neighbouring scale carries an odd component, and the octane value is
// the oddity, plus 4 when the determinant is 3 or 5 mod 8.
inline std::vector<int>
SpeciesListTwo(std::vector<std::array<int, 5>> const &sym) {
  std::vector<std::array<int, 5>> dense;
  dense.push_back({-2, 0, 1, 0, 0});
  dense.push_back({-1, 0, 1, 0, 0});
  size_t pos = 0;
  int s_max = sym.back()[0];
  for (int k = 0; k <= s_max; k++) {
    if (pos < sym.size() && sym[pos][0] == k) {
      dense.push_back(sym[pos]);
      pos++;
    } else {
      dense.push_back({k, 0, 1, 0, 0});
    }
  }
  dense.push_back({s_max + 1, 0, 1, 0, 0});
  dense.push_back({s_max + 2, 0, 1, 0, 0});
  std::vector<int> res;
  for (size_t k = 1; k + 1 < dense.size(); k++) {
    bool free = dense[k - 1][3] != 1 && dense[k + 1][3] != 1;
    int n = dense[k][1];
    int o = dense[k][4];
    int d = ModPositive(dense[k][2], 8);
    if (d == 3 || d == 5) {
      o = (o + 4) % 8;
    }
    int t = (dense[k][3] == 0 || n % 2 == 1) ? n / 2 : n / 2 - 1;
    if (free && (o == 0 || o == 1 || o == 7)) {
      res.push_back(2 * t);
    } else if (free && (o == 3 || o == 4 || o == 5)) {
      res.push_back(-2 * t);
    } else {
      res.push_back(2 * t + 1);
    }
  }
  return res;
}

// The p-mass m_p, as q * sqrt(r): the cross terms may be half-integral
// powers of p.
template <typename T>
PiRadical<T> PMass(std::vector<int> const &species,
                   std::vector<std::array<int, 2>> const &scale_rank,
                   int64_t p) {
  T val(1);
  for (int s : species) {
    val *= DiagonalFactorMp<T>(s, p);
  }
  int ct2 = 0;
  for (size_t j = 0; j < scale_rank.size(); j++) {
    for (size_t i = 0; i < j; i++) {
      ct2 += (scale_rank[j][0] - scale_rank[i][0]) * scale_rank[i][1] *
             scale_rank[j][1];
    }
  }
  return PiRadicalRational<T>(val) *
         RationalHalfPower<T>(ScalarFromInt64<T>(p), ct2);
}

// The standard p-mass std_p of a lattice of rank n and determinant det.
template <typename T> T StandardPMass(int n, int64_t det, int64_t p) {
  int s = (n + 1) / 2;
  T p_T = ScalarFromInt64<T>(p);
  T std_inv(2);
  T pk(1);
  for (int k = 1; k < s; k++) {
    pk *= p_T * p_T;
    std_inv *= T(1) - T(1) / pk;
  }
  if (n % 2 == 0) {
    int64_t D = SignedPower(s, det);
    // the Kronecker symbol (4 D / p): 0 at p = 2 and at p | D
    int eps = (p == 2) ? 0 : KroneckerSymbol(D, p);
    T ps(1);
    for (int i = 0; i < s; i++) {
      ps *= p_T;
    }
    std_inv *= T(1) - ScalarFromInt64<T>(eps) / ps;
  }
  return T(1) / std_inv;
}

template <typename T, typename Tint> T GenusMass(MyMatrix<T> const &Gram) {
  int n = Gram.rows();
  if (n == 1) {
    return T(1) / T(2);
  }
  int64_t det = DeterminantInt64<T, Tint>(Gram);
  int s = (n + 1) / 2;
  std::vector<int64_t> primes = BadPrimes(1, det);
  // the standard mass
  PiRadical<T> std_mass = {T(2), T(-n * (n + 1)) / T(4), T(1)};
  for (int j = 1; j <= n; j++) {
    std_mass = std_mass * GammaHalfInteger<T>(j);
  }
  for (int k = 1; k < s; k++) {
    std_mass = std_mass * ZetaEvenValue<T>(2 * k);
  }
  if (n % 2 == 0) {
    int64_t D = SignedPower(s, det);
    int64_t d = FundamentalDiscriminant(D);
    std_mass = std_mass * LQuadraticValue<T>(s, d);
    // the Euler factors at the primes of the local symbols
    for (int64_t p : primes) {
      T ps(1);
      for (int i = 0; i < s; i++) {
        ps *= ScalarFromInt64<T>(p);
      }
      T fac = T(1) - ScalarFromInt64<T>(KroneckerSymbol(d, p)) / ps;
      std_mass = std_mass * PiRadicalRational<T>(fac);
    }
  }
  PiRadical<T> mass = std_mass;
  for (int64_t p : primes) {
    PiRadical<T> m_p = PiRadicalRational<T>(T(1));
    if (p == 2) {
      std::vector<std::array<int, 5>> sym = TwoAdicSymbolRaw<T, Tint>(Gram);
      std::vector<std::array<int, 2>> scale_rank;
      int n_II = 0, n_I_I = 0;
      for (size_t k = 0; k < sym.size(); k++) {
        scale_rank.push_back({sym[k][0], sym[k][1]});
        if (sym[k][3] == 0) {
          n_II += sym[k][1];
        }
        if (k + 1 < sym.size() && sym[k][3] == 1 && sym[k + 1][3] == 1 &&
            sym[k][0] + 1 == sym[k + 1][0]) {
          n_I_I++;
        }
      }
      m_p = PMass<T>(SpeciesListTwo(sym), scale_rank, p);
      T type_factor(1);
      int e = n_I_I - n_II;
      for (int i = 0; i < (e >= 0 ? e : -e); i++) {
        type_factor *= T(2);
      }
      if (e < 0) {
        type_factor = T(1) / type_factor;
      }
      m_p = m_p * PiRadicalRational<T>(type_factor);
    } else {
      std::vector<std::array<int, 3>> sym = OddSymbolOfGram<T, Tint>(Gram, p);
      std::vector<std::array<int, 2>> scale_rank;
      for (auto const &c : sym) {
        scale_rank.push_back({c[0], c[1]});
      }
      m_p = PMass<T>(SpeciesListOdd(sym, p), scale_rank, p);
    }
    mass = mass * m_p *
           PiRadicalRational<T>(T(1) / StandardPMass<T>(n, det, p));
  }
  return PiRadicalToRational(mass);
}

// clang-format off
#endif  // SRC_GENERA_GENUSMASS_H_
// clang-format on
