// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENUS_H_
#define SRC_GENERA_HERMITIANGENUS_H_

// clang-format off
#include "NumberTheory.h"
#include "GeneraArithmetic.h"
#include "CyclotomicArithmetic.h"
#include "PolynomialModQ.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
// clang-format on

/*
  The genera of positive definite hermitian lattices over E = Q(zeta_d)
  relative to K_0 = Q(zeta_d)^+, d >= 3, as in hermitian_genera of Hecke.jl
  (Kirschmer's local symbols), for the d that are not powers of 2: there
  E / K_0 is not ramified at a dyadic prime.

  The primes. For a prime l, write d = l^a m with m prime to l (d = 2 mod 4
  is replaced by d / 2). In E the primes above l have ramification index
  phi(l^a), residue degree f = ord(l mod m) and number phi(m) / f. E / K_0
  is ramified above l exactly when d = l^a (a prime power); otherwise the
  primes of K_0 above l are inert in E when -1 is a power of l mod m, split
  when not. So the primes of K_0 above l all behave alike, and they are
  phi(m) / f in number when inert or ramified, phi(m) / (2 f) when split.

  The local symbols at a prime p of K_0, of rank n, scale valuations in
  [min, max] (in the prime P of E above p):
  - p split or inert: the Jordan blocks (s, r, det), the det forced: 1, or
    -1 at an inert prime when s r is odd (the norm is not surjective);
  - p ramified, not dyadic: s r even, the det free in {1, -1} for s even,
    hyperbolic_det^(r / 2) for s odd, hyperbolic_det = 1 when -1 is a local
    norm, that is (-1)^((l - 1) / 2) for d = l^a, l odd (K_0 is totally
    ramified at l, with residue field F_l).
  The local det is the product of the dets of the blocks. A family of local
  symbols is the genus of a global positive definite lattice when the
  number of non split primes of det -1 is even (Hecke's
  _check_global_genus with all signatures 0).

  The primes taken are those of K_0 above 2, the ramified ones and those
  where the scale bounds are not trivial; at the others the symbol is
  unimodular.
 */

struct HermPrimeData {
  int64_t l;
  // 0 split, 1 inert, 2 ramified
  int kind;
  bool dyadic;
  // the number of primes of K_0 above l, and their ramification index
  // over l
  int64_t n_prime;
  int64_t e0;
};

inline HermPrimeData HermPrimesAbove(int64_t d_in, int64_t l) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int64_t la = 1;
  int64_t m = d;
  while (m % l == 0) {
    m /= l;
    la *= l;
  }
  HermPrimeData res;
  res.l = l;
  res.dyadic = (l == 2);
  int64_t e_E = EulerPhi(la);
  int64_t f = MultiplicativeOrder(l, m);
  int64_t g_E = EulerPhi(m) / f;
  if (la > 1 && m == 1) {
    res.kind = 2;
    res.n_prime = g_E;
    res.e0 = e_E / 2;
  } else if (IsPowerModulo(-1, l, m)) {
    res.kind = 1;
    res.n_prime = g_E;
    res.e0 = e_E;
  } else {
    res.kind = 0;
    res.n_prime = g_E / 2;
    res.e0 = e_E;
  }
  return res;
}

// A local symbol: the Jordan blocks (scale, rank, det).
struct HermLocalSymbol {
  std::vector<std::array<int, 3>> blocks;
  int det() const {
    int res = 1;
    for (auto const &b : blocks) {
      res *= b[2];
    }
    return res;
  }
  bool operator<(HermLocalSymbol const &o) const { return blocks < o.blocks; }
  bool operator==(HermLocalSymbol const &o) const {
    return blocks == o.blocks;
  }
};

// The compositions of n into k non-negative parts.
inline void HermCompositions(int n, int k, std::vector<int> &cur,
                             std::vector<std::vector<int>> &out) {
  if (k == 1) {
    cur.push_back(n);
    out.push_back(cur);
    cur.pop_back();
    return;
  }
  for (int i = 0; i <= n; i++) {
    cur.push_back(i);
    HermCompositions(n - i, k - 1, cur, out);
    cur.pop_back();
  }
}

/*
  The local symbols of rank n at a prime of the given kind (not dyadic
  ramified), scales in [min_scale, max_scale]; hyperbolic_det for the
  ramified case.
 */
inline std::vector<HermLocalSymbol>
HermitianLocalGenera(HermPrimeData const &pd, int n, int min_scale,
                     int max_scale, int hyperbolic_det) {
  if (pd.kind == 2 && pd.dyadic) {
    std::cerr << "HERM_GENUS: the dyadic ramified case is not done\n";
    throw TerminalException{1};
  }
  std::vector<HermLocalSymbol> res;
  if (max_scale < min_scale) {
    return res;
  }
  std::vector<std::vector<int>> comps;
  std::vector<int> cur;
  HermCompositions(n, max_scale - min_scale + 1, cur, comps);
  for (auto const &c : comps) {
    std::vector<std::array<int, 2>> sr;
    for (size_t i = 0; i < c.size(); i++) {
      if (c[i] != 0) {
        sr.push_back({min_scale + static_cast<int>(i), c[i]});
      }
    }
    if (pd.kind != 2) {
      HermLocalSymbol g;
      for (auto const &[s, r] : sr) {
        int det = (pd.kind == 1 && (s * r) % 2 != 0) ? -1 : 1;
        g.blocks.push_back({s, r, det});
      }
      res.push_back(g);
      continue;
    }
    bool ok = true;
    for (auto const &[s, r] : sr) {
      ok = ok && ((s * r) % 2 == 0);
    }
    if (!ok) {
      continue;
    }
    // the dets of the blocks: free for an even scale, forced for an odd one
    std::vector<std::vector<int>> choices;
    for (auto const &[s, r] : sr) {
      if (s % 2 == 0) {
        choices.push_back({1, -1});
      } else {
        int h = 1;
        for (int k = 0; k < r / 2; k++) {
          h *= hyperbolic_det;
        }
        choices.push_back({h});
      }
    }
    std::vector<size_t> idx(sr.size(), 0);
    while (true) {
      HermLocalSymbol g;
      for (size_t k = 0; k < sr.size(); k++) {
        g.blocks.push_back({sr[k][0], sr[k][1], choices[k][idx[k]]});
      }
      res.push_back(g);
      size_t k = 0;
      while (k < idx.size()) {
        idx[k]++;
        if (idx[k] < choices[k].size()) {
          break;
        }
        idx[k] = 0;
        k++;
      }
      if (k == idx.size()) {
        break;
      }
    }
  }
  return res;
}

/*
  A prime of K_0 to treat, with the scale bounds (in the prime of E above
  it) and its local symbols.
 */
struct HermPrimeSymbols {
  HermPrimeData pd;
  int min_scale;
  int max_scale;
  std::vector<HermLocalSymbol> symbols;
};

/*
  The primes of K_0 to treat for hermitian_genera(E, n, all signatures 0;
  min_scale = O_E, max_scale = N O_E): above 2, the ramified one and those
  dividing N. The scale of N O_E at a prime p of K_0 above l: the norm N^2
  has p-valuation 2 e0 v_l(N), halved (Hecke's div) when p is not ramified.
  Each of the n_prime primes above l is listed.
 */
inline std::vector<HermPrimeSymbols> HermitianPrimeSymbols(int64_t d, int n,
                                                           int64_t N) {
  int64_t dn = NormalizeCyclotomicIndex(d);
  std::set<int64_t> primes{2};
  int64_t l_ram = PrimeOfPrimePower(dn);
  if (l_ram != 0) {
    primes.insert(l_ram);
  }
  {
    if (N > 1) {
      for (auto const &[p, mult] : FactorsIntMap(N)) {
        primes.insert(p);
      }
    }
  }
  std::vector<HermPrimeSymbols> res;
  for (int64_t l : primes) {
    HermPrimeData pd = HermPrimesAbove(d, l);
    int64_t vN = 0;
    for (int64_t x = N; x % l == 0; x /= l) {
      vN++;
    }
    int max_scale = static_cast<int>(2 * pd.e0 * vN);
    if (pd.kind != 2) {
      max_scale /= 2;
    }
    int hyperbolic_det = 1;
    if (pd.kind == 2 && !pd.dyadic) {
      hyperbolic_det = (((l - 1) / 2) % 2 == 0) ? 1 : -1;
    }
    std::vector<HermLocalSymbol> syms =
        HermitianLocalGenera(pd, n, 0, max_scale, hyperbolic_det);
    for (int64_t k = 0; k < pd.n_prime; k++) {
      res.push_back({pd, 0, max_scale, syms});
    }
  }
  return res;
}

/*
  The global genera: one local symbol per prime with an even number of non
  split primes of det -1. Each genus as the indices of its local symbols.
 */
inline std::vector<std::vector<size_t>>
HermitianGlobalGenera(std::vector<HermPrimeSymbols> const &ps) {
  std::vector<std::vector<size_t>> res;
  for (auto const &p : ps) {
    if (p.symbols.empty()) {
      return res;
    }
  }
  std::vector<size_t> idx(ps.size(), 0);
  while (true) {
    int n_minus = 0;
    for (size_t i = 0; i < ps.size(); i++) {
      if (ps[i].pd.kind != 0 && ps[i].symbols[idx[i]].det() == -1) {
        n_minus++;
      }
    }
    if (n_minus % 2 == 0) {
      res.push_back(idx);
    }
    size_t k = 0;
    while (k < idx.size()) {
      idx[k]++;
      if (idx[k] < ps[k].symbols.size()) {
        break;
      }
      idx[k] = 0;
      k++;
    }
    if (k == idx.size()) {
      break;
    }
  }
  return res;
}

// The factors of Phi_m mod l (l prime to m), all of degree ord(l mod m).
inline std::vector<PolyQ> HermCyclotomicFactorsModL(int64_t m, int64_t l) {
  T_small_prime q = static_cast<T_small_prime>(l);
  return EqualDegreeFactorization(
      PolyReduceSigned(CyclotomicPolynomial(static_cast<int>(m)), q),
      MultiplicativeOrder(l, m), q);
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENUS_H_
// clang-format on
