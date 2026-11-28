// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENERA_H_
#define SRC_GENERA_HERMITIANGENERA_H_

// clang-format off
#include "HermitianGenusSeed.h"
#include <array>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
# define DEBUG_HERMITIAN_GENERA
#endif

/*
  All the genera of Hermitian lattices over E / K_0 of a given rank m,
  determinant and signatures, as Hecke.jl's hermitian_genera does.

  The local genera at a prime p of K_0 (P a prime of E above it, the scales
  and the determinant in P-valuations): Jordan blocks (s, r) with scales in
  [min_scale, max_scale], of total rank m and sum of s r the valuation of
  the determinant, and
  * p split: the determinant class 1 (every unit is a norm);
  * p inert: the class (-1)^{s r}, the norm being onto the units;
  * p ramified, not dyadic: r even for an odd s, the class of a block of odd
    scale being that of the hyperbolic planes, h^{r/2} with h = 1 if -1 is a
    norm and -1 otherwise; both classes for an even s;
  * p ramified and dyadic, e = v_P(D) the different: r even for an odd s; a
    block of odd rank has an even scale s and the norm valuation s/2 and
    both classes; a block of even rank has the norm valuations k with
    s <= 2k < s + e and both classes, and when s + e is even also k = (s +
    e)/2 with the class h^{r/2}.
  Equivalent local symbols (HermLocalSymbolEquivalent) are kept once.

  A global genus is a choice of one local genus at each of the relevant
  primes (those over 2, the ramified ones and those of the determinant and
  of the scale bounds), subject to the product formula: the number of
  non-split primes where the determinant is not a norm and of the real
  places of K_0 with an odd number of negative entries is even. In the
  output, as in GENUS_HermitianEnumerate's input, the primes over 2 and the
  ramified ones are always written, the others only when the local genus is
  not unimodular.
 */

// The bounds at a prime p of K_0: the P-valuation of the determinant and
// of the scales.
struct HermGeneraBounds {
  int det_val;
  int min_scale;
  int max_scale;
};

/*
  The rank sequences: (s, r) blocks, scales in [smin, smax], total rank m,
  sum s r = dv.
 */
inline void HermGeneraScaleRanks(int m, int smin, int smax, int dv,
                                 std::vector<std::vector<std::array<int, 2>>> &out) {
  std::vector<std::array<int, 2>> cur;
  // rem_rank, rem_det, the next scale
  std::function<void(int, int, int)> rec = [&](int rr, int rd, int s) {
    if (rr == 0) {
      if (rd == 0) {
        out.push_back(cur);
      }
      return;
    }
    if (s > smax) {
      return;
    }
    for (int r = rr; r >= 0; r--) {
      int dd = s * r;
      // the remaining blocks have scales above s: their contribution is at
      // least (s + 1) (rr - r) and at most smax (rr - r)
      int rest = rd - dd;
      if (rr - r > 0) {
        // the scales above s are in [s + 1, smax]
        if (s + 1 > smax || rest < (s + 1) * (rr - r) || rest > smax * (rr - r)) {
          continue;
        }
      } else if (rest != 0) {
        continue;
      }
      if (r > 0) {
        cur.push_back({s, r});
      }
      rec(rr - r, rest, s + 1);
      if (r > 0) {
        cur.pop_back();
      }
    }
  };
  rec(m, dv, smin);
}

/*
  The local genera at the prime p of K_0 represented by l_prime[i_P]
  (HermK0PrimesAbove).
 */
template <typename T, typename Tint>
std::vector<HermLocalSymbol>
HermLocalGenera(CMOrder<T, Tint> const &ord,
                std::vector<HermIdealPrime<Tint>> const &l_prime, size_t i_P,
                int m, HermGeneraBounds const &bd) {
  HermIdealPrime<Tint> const &P = l_prime[i_P];
  bool ramified = P.type == HermIdealType::Ramified;
  bool inert = P.type == HermIdealType::Inert;
  bool dyadic = ramified && P.ell == 2;
  std::vector<std::vector<std::array<int, 2>>> l_sr;
  HermGeneraScaleRanks(m, bd.min_scale, bd.max_scale, bd.det_val, l_sr);
  std::vector<HermLocalSymbol> res;
  if (!ramified) {
    for (auto const &g : l_sr) {
      HermLocalSymbol sym;
      for (auto const &b : g) {
        int d = (inert && (b[0] * b[1]) % 2 != 0) ? -1 : 1;
        sym.blocks.push_back({b[0], b[1], d});
      }
      res.push_back(sym);
    }
    return res;
  }
  // h = 1 if -1 is a norm at P, -1 otherwise
  MyVector<T> minus_one = ZeroVector<T>(ord.deg);
  minus_one(0) = T(-1);
  int h = HermLocalNormClass(ord, l_prime, i_P, minus_one);
  int e = dyadic ? HermDifferentVal(ord, l_prime, i_P) : 0;
  auto add_unique = [&](HermLocalSymbol const &sym) -> void {
    for (auto const &o : res) {
      if (HermLocalSymbolEquivalent(o, sym, true, dyadic, e)) {
        return;
      }
    }
    res.push_back(sym);
  };
  for (auto const &g : l_sr) {
    bool ok = true;
    for (auto const &b : g) {
      if ((b[0] * b[1]) % 2 != 0) {
        ok = false;
      }
    }
    if (!ok) {
      continue;
    }
    // per block the possible (class, norm valuation)
    std::vector<std::vector<std::array<int, 2>>> l_choice;
    for (auto const &b : g) {
      std::vector<std::array<int, 2>> ch;
      int s = b[0], r = b[1];
      if (!dyadic) {
        if (s % 2 == 0) {
          ch = {{1, 0}, {-1, 0}};
        } else {
          // the hyperbolic planes: h^{r / 2}
          ch = {{(r / 2) % 2 == 0 ? 1 : h, 0}};
        }
      } else if (r % 2 == 1) {
        ch = {{1, s / 2}, {-1, s / 2}};
      } else {
        // k from ceil(s / 2) while 2 k < s + e
        int k = (s >= 0) ? (s + 1) / 2 : -((-s) / 2);
        while (2 * k < s + e) {
          ch.push_back({1, k});
          ch.push_back({-1, k});
          k++;
        }
        if ((s + e) % 2 == 0) {
          int hd = ((r / 2) % 2 == 0 || h == 1) ? 1 : -1;
          ch.push_back({hd, k});
        }
      }
      l_choice.push_back(ch);
    }
    // the cartesian product of the choices
    std::vector<size_t> idx(g.size(), 0);
    while (true) {
      HermLocalSymbol sym;
      for (size_t i = 0; i < g.size(); i++) {
        auto const &c = l_choice[i][idx[i]];
        sym.blocks.push_back({g[i][0], g[i][1], c[0]});
        if (dyadic) {
          sym.norms.push_back(c[1]);
        }
      }
      add_unique(sym);
      size_t pos = 0;
      while (pos < g.size()) {
        idx[pos]++;
        if (idx[pos] < l_choice[pos].size()) {
          break;
        }
        idx[pos] = 0;
        pos++;
      }
      if (pos == g.size()) {
        break;
      }
    }
  }
  return res;
}

// A relevant prime of K_0 for the listing.
template <typename Tint> struct HermGeneraPrime {
  int64_t ell;
  std::vector<HermIdealPrime<Tint>> l_prime;
  HermK0Prime<Tint> p;
  // always written (over 2 or ramified)
  bool always;
  std::vector<HermLocalSymbol> l_local;
};

/*
  The genera. primes: the relevant primes of K_0 with their bounds (a prime
  not given has the bounds (0, 0, 0): unimodular); n_odd_signature: the
  number of real places with an odd number of negative entries.
 */
template <typename T, typename Tint>
std::vector<std::string>
HermitianGenera(CMOrder<T, Tint> const &ord, int m,
                std::vector<std::pair<int64_t, std::map<size_t, HermGeneraBounds>>> const &l_bound,
                int n_odd_signature, [[maybe_unused]] std::ostream &os) {
  std::vector<HermGeneraPrime<Tint>> l_p;
  for (auto const &[ell, bmap] : l_bound) {
    std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
    for (auto const &p : HermK0PrimesAbove(ord, l_prime)) {
      HermGeneraBounds bd{0, 0, 0};
      auto iter = bmap.find(p.i_P);
      if (iter != bmap.end()) {
        bd = iter->second;
      }
      HermGeneraPrime<Tint> gp;
      gp.ell = ell;
      gp.l_prime = l_prime;
      gp.p = p;
      gp.always = (ell == 2) || p.type == HermIdealType::Ramified;
      gp.l_local = HermLocalGenera(ord, l_prime, p.i_P, m, bd);
#ifdef DEBUG_HERMITIAN_GENERA
      os << "HERMGENERA: at " << p.name << ", " << gp.l_local.size()
         << " local genera\n";
#endif
      if (gp.l_local.empty()) {
        return {};
      }
      l_p.push_back(gp);
    }
  }
  HermLocalSymbol unimod;
  unimod.blocks.push_back({0, m, 1});
  std::vector<std::string> res;
  std::vector<size_t> idx(l_p.size(), 0);
  while (true) {
    int n_non_norm = n_odd_signature;
    for (size_t i = 0; i < l_p.size(); i++) {
      HermLocalSymbol const &sym = l_p[i].l_local[idx[i]];
      if (l_p[i].p.type != HermIdealType::Split && sym.det() == -1) {
        n_non_norm++;
      }
    }
    if (n_non_norm % 2 == 0) {
      HermGlobalSymbol G;
      G.m = m;
      for (size_t i = 0; i < l_p.size(); i++) {
        HermLocalSymbol const &sym = l_p[i].l_local[idx[i]];
        if (l_p[i].always || !(sym == unimod)) {
          G.entries.push_back({l_p[i].ell, l_p[i].p.name, sym});
        }
      }
      res.push_back(HermGlobalSymbolString(G));
    }
    size_t pos = 0;
    while (pos < l_p.size()) {
      idx[pos]++;
      if (idx[pos] < l_p[pos].l_local.size()) {
        break;
      }
      idx[pos] = 0;
      pos++;
    }
    if (pos == l_p.size()) {
      break;
    }
  }
  return res;
}

/*
  The comparison of genus symbols written differently (another program's
  output): a symbol as its local symbols by prime of K_0, and the
  equivalence of two genera, prime by prime (HermLocalSymbolEquivalent), a
  prime absent from a symbol being unimodular.
 */
template <typename T, typename Tint>
std::map<HermK0Ref, HermLocalSymbol>
HermGenusLocalMap(HermPrimeCache<T, Tint> &cache, std::string const &str,
                  int &m) {
  HermSymbolInput<T> input = ParseHermSymbol<T>(str);
  m = input.m;
  std::map<HermK0Ref, HermLocalSymbol> res;
  for (size_t i = 0; i < input.l_ell.size(); i++) {
    res[HermResolvePrime(cache, input.l_ell[i], input.l_elt[i])] =
        input.l_sym[i];
  }
  return res;
}

template <typename T, typename Tint>
bool HermGenusEquivalent(HermPrimeCache<T, Tint> &cache, int m,
                         std::map<HermK0Ref, HermLocalSymbol> const &a,
                         std::map<HermK0Ref, HermLocalSymbol> const &b) {
  HermLocalSymbol unimod;
  unimod.blocks.push_back({0, m, 1});
  std::set<HermK0Ref> keys;
  for (auto const &kv : a) {
    keys.insert(kv.first);
  }
  for (auto const &kv : b) {
    keys.insert(kv.first);
  }
  for (auto const &ref : keys) {
    std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ref.ell);
    bool ramified = l_prime[ref.i_P].type == HermIdealType::Ramified;
    bool dyadic = ramified && ref.ell == 2;
    int e = dyadic ? HermDifferentVal(*cache.ord, l_prime, ref.i_P) : 0;
    auto ia = a.find(ref);
    auto ib = b.find(ref);
    HermLocalSymbol const &sa = (ia == a.end()) ? unimod : ia->second;
    HermLocalSymbol const &sb = (ib == b.end()) ? unimod : ib->second;
    if (!HermLocalSymbolEquivalent(sa, sb, ramified, dyadic, e)) {
      return false;
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENERA_H_
// clang-format on
