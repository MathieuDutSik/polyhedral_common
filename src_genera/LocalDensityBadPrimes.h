// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALDENSITYBADPRIMES_H_
#define SRC_GENERA_LOCALDENSITYBADPRIMES_H_

// clang-format off
#include "SiegelSeries.h"
#include <cstdint>
#include <map>
#include <optional>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_LOCAL_DENSITY_BAD_PRIMES
#endif

/*
  Local densities at a prime where T is not unimodular.

  RankOneDensity computes alpha_p(a, T) = lim_k p^{-k(n-1)} N_k(a),
  N_k(a) = #{x in (Z/p^k)^n : Q(x) = a mod p^k},  Q(x) = x^t T x / 2,
  exactly, by counting. The counting is cheap because N_k depends on a only
  through its class under multiplication by unit squares (x -> u x is a
  bijection), and there are few such classes modulo p^k: 0, and p^v u with
  u up to squares -- its Legendre symbol for odd p, u mod 8 (or mod 4, 2
  when p^{k-v} is small) for p = 2. The distribution of Q is the
  convolution of the distributions of the Jordan blocks of T, and a
  convolution of class functions only needs the number of ways a fixed
  element of a class is the sum of an element of a second class and one of
  a third, counted once per (p, k). The level k is raised until two
  consecutive levels give the same density, so the stabilisation is checked
  rather than assumed.

  UnitChainDensity computes alpha_p(S, T) at an odd prime by splitting off
  the unit diagonal entries of S one at a time: a vector of T whose norm is
  a p-adic unit splits off, and by Witt cancellation its orthogonal
  complement is determined -- the unimodular Jordan component loses one
  dimension and its determinant is divided by the norm. Hence
      alpha_p(s_1 + ... + s_m, T) = alpha_p(s_1, T) alpha_p(s_2 + ..., T'),
  exact as long as every entry but the last is a unit. This is the method
  of the Sage prototype (density_oddbad.py), which reproduced 54 of 54
  extracted values and every reference value it accepted; it refuses S
  with two or more non-unit entries at p rather than guess.
 */

//
// Class functions on Z / p^k under unit squares
//

struct ResidueClasses {
  int64_t p;
  int k;
  int64_t modulus;
  // class index of every residue, and the size of each class
  std::vector<int> class_of;
  std::vector<int64_t> size;
  // a representative of each class
  std::vector<int64_t> rep;
  // structure constants: n_sum[A][B][C] = #{(c1, c2) in A x B :
  // c1 + c2 = rep[C]}
  std::vector<std::vector<std::vector<int64_t>>> n_sum;
  int n_class() const { return static_cast<int>(size.size()); }
};

inline int64_t PowerInt64(int64_t a, int e) {
  int64_t r = 1;
  for (int i = 0; i < e; i++) {
    r *= a;
  }
  return r;
}

// The class key of a residue c modulo p^k: (v, unit key), v = k for c = 0.
inline std::pair<int, int64_t> ResidueClassKey(int64_t c, int64_t p, int k,
                                                int64_t modulus) {
  c = ModPositive(c, modulus);
  if (c == 0) {
    return {k, 0};
  }
  int v = 0;
  while (c % p == 0) {
    c /= p;
    v++;
  }
  if (p != 2) {
    return {v, LegendreUnit(c, p)};
  }
  int64_t rem = PowerInt64(2, k - v);
  int64_t mod_u = std::min<int64_t>(8, rem);
  return {v, c % mod_u};
}

inline ResidueClasses BuildResidueClasses(int64_t p, int k) {
  ResidueClasses rc;
  rc.p = p;
  rc.k = k;
  rc.modulus = PowerInt64(p, k);
  std::map<std::pair<int, int64_t>, int> index;
  rc.class_of.resize(rc.modulus);
  for (int64_t c = 0; c < rc.modulus; c++) {
    auto key = ResidueClassKey(c, p, k, rc.modulus);
    auto iter = index.find(key);
    int idx;
    if (iter == index.end()) {
      idx = static_cast<int>(rc.size.size());
      index[key] = idx;
      rc.size.push_back(0);
      rc.rep.push_back(c);
    } else {
      idx = iter->second;
    }
    rc.class_of[c] = idx;
    rc.size[idx]++;
  }
  int nc = rc.n_class();
  rc.n_sum.assign(nc, std::vector<std::vector<int64_t>>(
                          nc, std::vector<int64_t>(nc, 0)));
  for (int C = 0; C < nc; C++) {
    int64_t c0 = rc.rep[C];
    for (int64_t c1 = 0; c1 < rc.modulus; c1++) {
      int64_t c2 = ModPositive(c0 - c1, rc.modulus);
      rc.n_sum[rc.class_of[c1]][rc.class_of[c2]][C]++;
    }
  }
  return rc;
}

/*
  The caches of the bad-prime densities. In a sweep over the root lattices
  S for a fixed genus, the same residue classes, block distributions and
  rank one densities come back again and again: T is fixed, and the
  targets take few p-adic classes.
 */
template <typename T> struct BadPrimeDensityCache {
  std::map<std::pair<int64_t, int>, ResidueClasses> classes;
  std::map<std::tuple<int64_t, int, int, int, int64_t>, std::vector<T>>
      block_dist;
  std::map<std::tuple<std::tuple<int64_t,
                                 std::vector<std::tuple<int, int, int64_t>>>,
                      int, int64_t>,
           T>
      rank_one;
};

inline ResidueClasses const &
GetResidueClasses(std::map<std::pair<int64_t, int>, ResidueClasses> &classes,
                  int64_t p, int k) {
  auto key = std::make_pair(p, k);
  auto iter = classes.find(key);
  if (iter != classes.end()) {
    return iter->second;
  }
  return classes.emplace(key, BuildResidueClasses(p, k)).first->second;
}

/*
  A distribution of values modulo p^k as a class function: count[A] is the
  number of vectors whose value is a GIVEN element of class A (the same for
  every element of A).
 */
template <typename T>
std::vector<T> ConvolveClassFunctions(ResidueClasses const &rc,
                                      std::vector<T> const &f,
                                      std::vector<T> const &g) {
  int nc = rc.n_class();
  std::vector<T> h(nc, T(0));
  for (int A = 0; A < nc; A++) {
    if (f[A] == 0) {
      continue;
    }
    for (int B = 0; B < nc; B++) {
      if (g[B] == 0) {
        continue;
      }
      T fg = f[A] * g[B];
      for (int C = 0; C < nc; C++) {
        int64_t cnt = rc.n_sum[A][B][C];
        if (cnt != 0) {
          h[C] += fg * ScalarFromInt64<T>(cnt);
        }
      }
    }
  }
  return h;
}

// The class function of the values of one Jordan block modulo p^k.
template <typename T>
std::vector<T> BlockValueDistribution(ResidueClasses const &rc,
                                      SiegelJordanBlock const &b) {
  int64_t M = rc.modulus;
  std::vector<int64_t> hist(M, 0);
  int64_t pe = 1;
  for (int i = 0; i < b.e && pe < M; i++) {
    pe *= rc.p;
  }
  if (b.e >= rc.k) {
    pe = 0;
  }
  if (b.kind == SiegelBlockKind::Unit) {
    for (int64_t x = 0; x < M; x++) {
      hist[MulMod(MulMod(pe, b.unit, M), MulMod(x, x, M), M)]++;
    }
  } else {
    // 2^e x y, or 2^e (x^2 + x y + y^2)
    for (int64_t x = 0; x < M; x++) {
      for (int64_t y = 0; y < M; y++) {
        int64_t q = MulMod(x, y, M);
        if (b.kind == SiegelBlockKind::Y) {
          q = ModPositive(q + MulMod(x, x, M) + MulMod(y, y, M), M);
        }
        hist[MulMod(pe, q, M)]++;
      }
    }
  }
  int nc = rc.n_class();
  std::vector<T> f(nc, T(0));
  std::vector<bool> seen(nc, false);
  for (int64_t c = 0; c < M; c++) {
    int A = rc.class_of[c];
    if (!seen[A]) {
      f[A] = ScalarFromInt64<T>(hist[c]);
      seen[A] = true;
    }
#ifdef SANITY_CHECK_LOCAL_DENSITY_BAD_PRIMES
    else if (f[A] != ScalarFromInt64<T>(hist[c])) {
      std::cerr << "LDBP: the distribution of a block is not a class "
                   "function\n";
      throw TerminalException{1};
    }
#endif
  }
  return f;
}

// N_k(a) for the blocks q, a given by its p-adic class.
template <typename T>
T CountRepresentationsModPk(ResidueClasses const &rc,
                            SiegelJordanBlocks const &q, PadicClass const &a,
                            BadPrimeDensityCache<T> &cache) {
  int nc = rc.n_class();
  std::vector<T> total(nc, T(0));
  // the distribution of the empty form: only the value 0
  total[rc.class_of[0]] = T(1);
  for (auto const &b : q.blocks) {
    auto key = std::make_tuple(rc.p, rc.k, b.e, static_cast<int>(b.kind),
                               b.unit);
    auto iter = cache.block_dist.find(key);
    if (iter == cache.block_dist.end()) {
      iter = cache.block_dist.emplace(key, BlockValueDistribution<T>(rc, b))
                 .first;
    }
    total = ConvolveClassFunctions(rc, total, iter->second);
  }
  int64_t a_res;
  if (a.val >= rc.k) {
    a_res = 0;
  } else {
    a_res = MulMod(PowerInt64(rc.p, a.val), a.unit, rc.modulus);
    // a.unit is only known modulo UnitModulus(p), which is all the class
    // of a depends on
  }
  return total[rc.class_of[a_res]];
}

/*
  alpha_p(a, T) for the Jordan blocks q of T/2 and a target a (a p-adic
  class), raising k until two consecutive levels agree.
 */
template <typename T>
T RankOneDensity(SiegelJordanBlocks const &q, PadicClass const &a,
                 BadPrimeDensityCache<T> &cache) {
  auto cache_key = std::make_tuple(SiegelBlocksKey(q), a.val, a.unit);
  auto cache_iter = cache.rank_one.find(cache_key);
  if (cache_iter != cache.rank_one.end()) {
    return cache_iter->second;
  }
  int64_t p = q.p;
  int n = q.dim();
  int e_max = 0;
  for (auto const &b : q.blocks) {
    e_max = std::max(e_max, b.e);
  }
  int k = a.val + e_max + ((p == 2) ? 4 : 2);
  T p_T = ScalarFromInt64<T>(p);
  auto density_at = [&](int kk) -> T {
    ResidueClasses const &rc = GetResidueClasses(cache.classes, p, kk);
    return CountRepresentationsModPk<T>(rc, q, a, cache) /
           PowerInt(p_T, kk * (n - 1));
  };
  T prev = density_at(k);
  while (true) {
    k++;
    T next = density_at(k);
    if (next == prev) {
      cache.rank_one[cache_key] = next;
      return next;
    }
    prev = next;
    if (k > a.val + e_max + 40) {
      std::cerr << "LDBP: the rank one density did not stabilise\n";
      throw TerminalException{1};
    }
  }
}

//
// The unit chain at an odd prime
//

// The diagonal Jordan entries of the half-integral S / 2 at an odd prime,
// as p-adic classes in the order produced by the decomposition.
template <typename T, typename Tint>
std::vector<PadicClass> OddDiagonalClasses(MyMatrix<T> const &S, int64_t p) {
  std::vector<T> diag = JordanDiagonalOdd<T, Tint>(S / T(2), p);
  std::vector<PadicClass> res;
  for (auto const &d : diag) {
    res.push_back(PadicFromRational<T, Tint>(d, p));
  }
  return res;
}

inline SiegelJordanBlocks BlocksOfDiagonal(std::vector<PadicClass> const &D,
                                           int64_t p) {
  SiegelJordanBlocks q;
  q.p = p;
  for (auto const &c : D) {
    q.blocks.push_back({c.val, SiegelBlockKind::Unit, c.unit});
  }
  return q;
}

/*
  alpha_p(S, T), p odd, by the unit chain; nothing when S has two or more
  non-unit diagonal entries at p, or when the unimodular part of T is too
  small for a unit to split off.
 */
template <typename T, typename Tint>
std::optional<T> UnitChainDensity(MyMatrix<T> const &S,
                                  std::vector<PadicClass> const &DT, int64_t p,
                                  BadPrimeDensityCache<T> &cache) {
  if (p == 2) {
    std::cerr << "LDBP: the unit chain is for odd primes\n";
    throw TerminalException{1};
  }
  std::vector<PadicClass> D = DT;
  std::vector<PadicClass> DS = OddDiagonalClasses<T, Tint>(S, p);
  std::vector<PadicClass> chain;
  std::vector<PadicClass> nonunits;
  for (auto const &c : DS) {
    if (c.val == 0) {
      chain.push_back(c);
    } else {
      nonunits.push_back(c);
    }
  }
  if (nonunits.size() > 1) {
    return {};
  }
  chain.insert(chain.end(), nonunits.begin(), nonunits.end());
  T alpha(1);
  for (size_t i = 0; i < chain.size(); i++) {
    PadicClass const &s = chain[i];
    alpha *= RankOneDensity<T>(BlocksOfDiagonal(D, p), s, cache);
    if (i + 1 == chain.size()) {
      break;
    }
    // split off a unit-norm vector: the unimodular part loses one entry and
    // its determinant is divided by s
    std::vector<PadicClass> uni, rest;
    for (auto const &c : D) {
      if (c.val == 0) {
        uni.push_back(c);
      } else {
        rest.push_back(c);
      }
    }
    if (uni.size() < 2) {
      return {};
    }
    int64_t prod = 1;
    for (auto const &c : uni) {
      prod = MulMod(prod, c.unit, p);
    }
    std::vector<PadicClass> keep(uni.begin(), uni.end() - 2);
    int64_t prod_keep = 1;
    for (auto const &c : keep) {
      prod_keep = MulMod(prod_keep, c.unit, p);
    }
    int64_t last = MulMod(prod, InverseModInt64(MulMod(s.unit, prod_keep, p), p), p);
    keep.push_back({0, last});
    D = keep;
    D.insert(D.end(), rest.begin(), rest.end());
  }
  return alpha;
}


//
// The representation mass M_S / mass(G)
//

/*
  M_S / mass(G) for the genus G of the even lattice Tgram, of even rank n
  and odd determinant, and a positive definite even S of rank m with
  n - m >= 3, from the local densities alone:

      2 * (generic part) * prod_{p | 2 det S det T} alpha_p(S, T),

  alpha_p from the Siegel series where T is unimodular at p (p = 2 always,
  T_2 being even unimodular), and from the unit chain at the odd primes
  dividing det T. Nothing when the unit chain refuses S.
 */
/*
  What the densities need of the genus, computed once: its rank and
  determinant, and the diagonal Jordan classes of T / 2 at every odd prime
  dividing the determinant.
 */
struct GenusDensityData {
  int n;
  int64_t detT;
  std::map<int64_t, std::vector<PadicClass>> odd_diagonal;
};

template <typename T, typename Tint>
int64_t DeterminantInt64(MyMatrix<T> const &M) {
  Tint d = UniversalScalarConversion<Tint, T>(DeterminantMat(M));
  return UniversalScalarConversion<int64_t, Tint>(d);
}

template <typename T, typename Tint>
GenusDensityData ComputeGenusDensityData(MyMatrix<T> const &Tgram) {
  GenusDensityData data;
  data.n = Tgram.rows();
  data.detT = DeterminantInt64<T, Tint>(Tgram);
  for (int64_t p : BadPrimes(1, data.detT)) {
    if (p != 2) {
      data.odd_diagonal[p] = OddDiagonalClasses<T, Tint>(Tgram, p);
    }
  }
  return data;
}

template <typename T> struct LocalDensityCache {
  SiegelSeriesCache<T> siegel;
  BadPrimeDensityCache<T> bad;
};

template <typename T, typename Tint>
std::optional<T> MassRatioFromDensities(MyMatrix<T> const &S,
                                        GenusDensityData const &data,
                                        LocalDensityCache<T> &cache) {
  int n = data.n;
  int m = S.rows();
  if (n % 2 != 0) {
    std::cerr << "LDBP: MassRatioFromDensities needs an even rank\n";
    throw TerminalException{1};
  }
  if (n - m < 3) {
    std::cerr << "LDBP: MassRatioFromDensities needs n - m >= 3\n";
    throw TerminalException{1};
  }
  int64_t detS = DeterminantInt64<T, Tint>(S);
  int64_t detT = data.detT;
  if (detT % 2 == 0) {
    std::cerr << "LDBP: MassRatioFromDensities needs an odd determinant\n";
    throw TerminalException{1};
  }
  int k = n / 2;
  T val = T(2) * GenericDensityPart<T>(n, m, detS, detT);
  for (int64_t p : BadPrimes(detS, detT)) {
    if (detT % p == 0) {
      std::optional<T> opt = UnitChainDensity<T, Tint>(
          S, data.odd_diagonal.at(p), p, cache.bad);
      if (!opt) {
        return {};
      }
      val *= *opt;
    } else {
      int xi_T = KroneckerSymbol(SignedPower(k, detT), p);
      val *= AlphaUnimodularTarget<T, Tint>(S, p, k, xi_T, cache.siegel);
    }
  }
  return val;
}

// clang-format off
#endif  // SRC_GENERA_LOCALDENSITYBADPRIMES_H_
// clang-format on
