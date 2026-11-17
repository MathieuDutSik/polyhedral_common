// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENUS_H_
#define SRC_GENERA_HERMITIANGENUS_H_

// clang-format off
#include "NumberTheory.h"
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

inline int64_t HermEulerPhi(int64_t n) {
  int64_t res = n;
  int64_t m = n;
  for (int64_t p = 2; p * p <= m; p++) {
    if (m % p == 0) {
      while (m % p == 0) {
        m /= p;
      }
      res -= res / p;
    }
  }
  if (m > 1) {
    res -= res / m;
  }
  return res;
}

// The multiplicative order of x modulo m (gcd(x, m) = 1), 1 for m <= 2.
inline int64_t HermOrderMod(int64_t x, int64_t m) {
  if (m <= 2) {
    return 1;
  }
  int64_t y = x % m;
  int64_t k = 1;
  while (y != 1) {
    y = (y * (x % m)) % m;
    k++;
  }
  return k;
}

inline bool HermIsPowerMod(int64_t target, int64_t x, int64_t m) {
  if (m <= 2) {
    return true;
  }
  int64_t t = ((target % m) + m) % m;
  int64_t y = 1;
  int64_t ord = HermOrderMod(x, m);
  for (int64_t k = 0; k < ord; k++) {
    if (y == t) {
      return true;
    }
    y = (y * (x % m)) % m;
  }
  return false;
}

inline int64_t HermNormalizedConductor(int64_t d) {
  if (d % 4 == 2) {
    return d / 2;
  }
  return d;
}

inline bool HermIsPrimePower(int64_t d, int64_t &l) {
  for (int64_t p = 2; p <= d; p++) {
    if (d % p == 0) {
      int64_t m = d;
      while (m % p == 0) {
        m /= p;
      }
      l = p;
      return m == 1;
    }
  }
  return false;
}

inline HermPrimeData HermPrimesAbove(int64_t d_in, int64_t l) {
  int64_t d = HermNormalizedConductor(d_in);
  int64_t la = 1;
  int64_t m = d;
  while (m % l == 0) {
    m /= l;
    la *= l;
  }
  HermPrimeData res;
  res.l = l;
  res.dyadic = (l == 2);
  int64_t e_E = HermEulerPhi(la);
  int64_t f = HermOrderMod(l, m);
  int64_t g_E = HermEulerPhi(m) / f;
  if (la > 1 && m == 1) {
    res.kind = 2;
    res.n_prime = g_E;
    res.e0 = e_E / 2;
  } else if (HermIsPowerMod(-1, l, m)) {
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
  int64_t dn = HermNormalizedConductor(d);
  std::set<int64_t> primes{2};
  int64_t l_ram = 0;
  if (HermIsPrimePower(dn, l_ram)) {
    primes.insert(l_ram);
  }
  {
    int64_t x = N;
    for (int64_t p = 2; p * p <= x; p++) {
      while (x % p == 0) {
        primes.insert(p);
        x /= p;
      }
    }
    if (x > 1) {
      primes.insert(x);
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

/*
  Polynomials over F_l (l prime, < 2^31), the coefficients in [0, l), of
  increasing degree, without leading zeros (the zero polynomial is empty).
 */
using HermPolyFl = std::vector<int64_t>;

inline void HermFlTrim(HermPolyFl &a) {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
}

inline int64_t HermFlInv(int64_t x, int64_t l) {
  // x^(l - 2) mod l
  int64_t res = 1, b = x % l;
  int64_t e = l - 2;
  while (e > 0) {
    if (e & 1) {
      res = static_cast<int64_t>((static_cast<__int128>(res) * b) % l);
    }
    b = static_cast<int64_t>((static_cast<__int128>(b) * b) % l);
    e >>= 1;
  }
  return res;
}

inline HermPolyFl HermFlMul(HermPolyFl const &a, HermPolyFl const &b,
                            int64_t l) {
  if (a.empty() || b.empty()) {
    return {};
  }
  HermPolyFl c(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] = (c[i + j] + a[i] * b[j]) % l;
    }
  }
  HermFlTrim(c);
  return c;
}

// The remainder of a modulo b (b non zero).
inline HermPolyFl HermFlRem(HermPolyFl a, HermPolyFl const &b, int64_t l) {
  int64_t inv = HermFlInv(b.back(), l);
  while (a.size() >= b.size()) {
    int64_t c = (a.back() * inv) % l;
    size_t shift = a.size() - b.size();
    for (size_t i = 0; i < b.size(); i++) {
      a[shift + i] = ((a[shift + i] - c * b[i]) % l + l) % l;
    }
    HermFlTrim(a);
  }
  return a;
}

inline HermPolyFl HermFlMonic(HermPolyFl a, int64_t l) {
  if (a.empty()) {
    return a;
  }
  int64_t inv = HermFlInv(a.back(), l);
  for (auto &x : a) {
    x = (x * inv) % l;
  }
  return a;
}

inline HermPolyFl HermFlGcd(HermPolyFl a, HermPolyFl b, int64_t l) {
  while (!b.empty()) {
    HermPolyFl r = HermFlRem(a, b, l);
    a = b;
    b = r;
  }
  return HermFlMonic(a, l);
}

// The quotient of a by b (exact or not, the remainder dropped).
inline HermPolyFl HermFlQuo(HermPolyFl a, HermPolyFl const &b, int64_t l) {
  int64_t inv = HermFlInv(b.back(), l);
  HermPolyFl q;
  if (a.size() >= b.size()) {
    q.assign(a.size() - b.size() + 1, 0);
  }
  while (a.size() >= b.size()) {
    int64_t c = (a.back() * inv) % l;
    size_t shift = a.size() - b.size();
    q[shift] = c;
    for (size_t i = 0; i < b.size(); i++) {
      a[shift + i] = ((a[shift + i] - c * b[i]) % l + l) % l;
    }
    HermFlTrim(a);
  }
  HermFlTrim(q);
  return q;
}

inline HermPolyFl HermFlPowMod(HermPolyFl base, unsigned __int128 e,
                               HermPolyFl const &f, int64_t l) {
  HermPolyFl res{1};
  base = HermFlRem(base, f, l);
  while (e > 0) {
    if (e & 1) {
      res = HermFlRem(HermFlMul(res, base, l), f, l);
    }
    base = HermFlRem(HermFlMul(base, base, l), f, l);
    e >>= 1;
  }
  return res;
}

/*
  The monic irreducible factors of the squarefree monic f over F_l, all of
  degree deg (Cantor-Zassenhaus; for l = 2 by the trace map), with a fixed
  sequence of trial polynomials.
 */
inline std::vector<HermPolyFl> HermFlEqualDegreeFactors(HermPolyFl const &f,
                                                        int deg, int64_t l) {
  int n = static_cast<int>(f.size()) - 1;
  if (n == deg) {
    return {HermFlMonic(f, l)};
  }
  uint64_t seed = 12345;
  for (int iter = 0; iter < 100000; iter++) {
    HermPolyFl a(n, 0);
    for (int i = 0; i < n; i++) {
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      a[i] = static_cast<int64_t>((seed >> 33) % static_cast<uint64_t>(l));
    }
    HermFlTrim(a);
    if (a.empty()) {
      continue;
    }
    HermPolyFl b;
    if (l == 2) {
      // a + a^2 + ... + a^(2^(deg - 1))
      HermPolyFl t = HermFlRem(a, f, l), pw = t;
      for (int k = 1; k < deg; k++) {
        pw = HermFlRem(HermFlMul(pw, pw, l), f, l);
        HermPolyFl s(std::max(t.size(), pw.size()), 0);
        for (size_t i = 0; i < s.size(); i++) {
          s[i] = ((i < t.size() ? t[i] : 0) + (i < pw.size() ? pw[i] : 0)) % 2;
        }
        HermFlTrim(s);
        t = s;
      }
      b = t;
    } else {
      unsigned __int128 q = 1;
      for (int k = 0; k < deg; k++) {
        q *= static_cast<unsigned __int128>(l);
      }
      b = HermFlPowMod(a, (q - 1) / 2, f, l);
      if (b.empty()) {
        b = HermPolyFl{static_cast<int64_t>(l - 1)};
      } else {
        b[0] = (b[0] + l - 1) % l;
        HermFlTrim(b);
      }
    }
    HermPolyFl g = HermFlGcd(f, b, l);
    int dg = static_cast<int>(g.size()) - 1;
    if (dg > 0 && dg < n) {
      std::vector<HermPolyFl> res = HermFlEqualDegreeFactors(g, deg, l);
      for (auto const &h :
           HermFlEqualDegreeFactors(HermFlQuo(f, g, l), deg, l)) {
        res.push_back(h);
      }
      std::sort(res.begin(), res.end());
      return res;
    }
  }
  std::cerr << "HERM_GENUS: the equal degree factorization failed\n";
  throw TerminalException{1};
}

// Phi_m over Z, of increasing degree.
inline std::vector<int64_t> HermCyclotomic(int64_t m) {
  std::vector<int64_t> num(static_cast<size_t>(m) + 1, 0);
  num[0] = -1;
  num[static_cast<size_t>(m)] = 1;
  for (int64_t e = 1; e < m; e++) {
    if (m % e != 0) {
      continue;
    }
    std::vector<int64_t> den = HermCyclotomic(e);
    size_t dn = num.size() - 1, dd = den.size() - 1;
    std::vector<int64_t> quo(dn - dd + 1, 0);
    for (size_t k = dn - dd + 1; k-- > 0;) {
      quo[k] = num[k + dd];
      for (size_t i = 0; i <= dd; i++) {
        num[k + i] -= quo[k] * den[i];
      }
    }
    num = quo;
  }
  return num;
}

// The factors of Phi_m mod l (l prime to m), all of degree ord(l mod m).
inline std::vector<HermPolyFl> HermCyclotomicFactorsModL(int64_t m,
                                                         int64_t l) {
  std::vector<int64_t> phi = HermCyclotomic(m);
  HermPolyFl f;
  for (auto c : phi) {
    f.push_back(((c % l) + l) % l);
  }
  HermFlTrim(f);
  return HermFlEqualDegreeFactors(f, static_cast<int>(HermOrderMod(l, m)), l);
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENUS_H_
// clang-format on
