// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SYMMETRICORBITCOUNTING_H_
#define SRC_GENERA_SYMMETRICORBITCOUNTING_H_

// clang-format off
#include "SiegelSeries.h"
#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <tuple>
#include <utility>
#include <vector>
// clang-format on

inline int64_t PowerInt64(int64_t a, int e) {
  int64_t r = 1;
  for (int i = 0; i < e; i++) {
    r *= a;
  }
  return r;
}

//
// Counting for several non-unit columns: orbit class functions on Sym_j
//

/*
  N_k(Sigma) = #{X in M_{n x j}(Z/p^k) : X^t C X = Sigma}, C = diag(c_1, ...,
  c_n), for p odd. X -> X g (g in GL_j(Z/p^k)) is a bijection between the
  solutions for Sigma and for g^t Sigma g, so the distribution of X^t C X is
  a class function on the GL_j-orbits of Sym_j(Z/p^k). An orbit is given by
  its Jordan invariants: for every level v < k the number of diagonal
  entries p^v u of a diagonalisation and the Legendre symbol of the product
  of their units; the entries of level >= k are 0. Adding a row x of
  coefficient c adds c x x^t, so
      new[C] = sum_A H_c[C][A] old[A],
      H_c[C][A] = #{x in (Z/p^k)^j : Sigma_C - c x x^t in orbit A},
  with Sigma_C a fixed element of the orbit C. H_c is counted once per
  (p, k, j, class of c) by running over x, which costs p^{kj} times the
  number of orbits; beyond a budget the counting declines.
 */

// The Jordan invariants of a symmetric j x j matrix modulo p^k, p odd, as a
// key: for v = 0, ..., k-1 the pair (count, Legendre symbol of the product),
// then the number of entries of level >= k.
inline std::vector<int> SymmetricOrbitKey(std::vector<int64_t> a, int j,
                                          int64_t p, int k, int64_t mod) {
  auto val = [&](int64_t x) -> int {
    x = ModPositive(x, mod);
    if (x == 0) {
      return k;
    }
    int v = 0;
    while (x % p == 0) {
      x /= p;
      v++;
    }
    return v;
  };
  std::vector<int> count(k + 1, 0);
  std::vector<int64_t> prod(k, 1);
  std::vector<bool> active(j, true);
  int n_active = j;
  while (n_active > 0) {
    // the entry of minimal valuation, a diagonal one if possible
    int best_i = -1, best_j = -1, best_v = k + 1;
    bool best_diag = false;
    for (int r = 0; r < j; r++) {
      if (!active[r]) {
        continue;
      }
      for (int c = r; c < j; c++) {
        if (!active[c]) {
          continue;
        }
        int v = val(a[r * j + c]);
        bool diag = (r == c);
        if (v < best_v || (v == best_v && diag && !best_diag)) {
          best_v = v;
          best_i = r;
          best_j = c;
          best_diag = diag;
        }
      }
    }
    if (best_v >= k) {
      count[k] += n_active;
      break;
    }
    int i0 = best_i;
    if (!best_diag) {
      // make the diagonal entry (i0, i0) of minimal valuation: add row and
      // column best_j to row and column i0
      int j0 = best_j;
      for (int c = 0; c < j; c++) {
        a[i0 * j + c] = ModPositive(a[i0 * j + c] + a[j0 * j + c], mod);
      }
      for (int r = 0; r < j; r++) {
        a[r * j + i0] = ModPositive(a[r * j + i0] + a[r * j + j0], mod);
      }
    }
    int64_t piv = ModPositive(a[i0 * j + i0], mod);
    int v = val(piv);
    int64_t pv = PowerInt64(p, v);
    int64_t u = piv / pv;
    int64_t u_inv = InverseModInt64(u, mod);
    // the Schur complement: a_rs <- a_rs - (a_ri / a_ii) a_is
    for (int r = 0; r < j; r++) {
      if (!active[r] || r == i0) {
        continue;
      }
      int64_t ari = ModPositive(a[r * j + i0], mod);
      int64_t f = MulMod(ari / pv, u_inv, mod);
      for (int c = 0; c < j; c++) {
        if (!active[c] || c == i0) {
          continue;
        }
        a[r * j + c] = ModPositive(a[r * j + c] - MulMod(f, a[i0 * j + c], mod), mod);
      }
    }
    count[v]++;
    prod[v] = MulMod(prod[v], u, p);
    active[i0] = false;
    n_active--;
  }
  std::vector<int> key;
  for (int v = 0; v < k; v++) {
    key.push_back(count[v]);
    key.push_back(count[v] > 0 ? LegendreUnit(prod[v], p) : 0);
  }
  key.push_back(count[k]);
  return key;
}

struct SymmetricOrbits {
  int64_t p;
  int k;
  int j;
  int64_t mod;
  std::map<std::vector<int>, int> index;
  std::vector<std::vector<int64_t>> rep;
  // H[class of c][C][A], keyed by the class of c
  std::map<std::pair<int, int64_t>, std::vector<std::vector<int64_t>>> H;
  int n_orbit() const { return static_cast<int>(rep.size()); }
};

// All the orbits, with a diagonal representative of each.
inline SymmetricOrbits BuildSymmetricOrbits(int64_t p, int k, int j) {
  SymmetricOrbits so;
  so.p = p;
  so.k = k;
  so.j = j;
  so.mod = PowerInt64(p, k);
  int64_t nonres = 2;
  while (LegendreUnit(nonres, p) != -1) {
    nonres++;
  }
  // the levels of the j diagonal entries, non-decreasing, level k for 0
  std::vector<int> levels(j, 0);
  std::function<void(int, int)> rec = [&](int pos, int min_level) -> void {
    if (pos == j) {
      // the levels present, and a choice of Legendre symbol for each
      std::vector<int> present;
      for (int v = 0; v < k; v++) {
        if (std::count(levels.begin(), levels.end(), v) > 0) {
          present.push_back(v);
        }
      }
      int n_present = static_cast<int>(present.size());
      for (int mask = 0; mask < (1 << n_present); mask++) {
        std::vector<int64_t> M(j * j, 0);
        for (int ip = 0; ip < n_present; ip++) {
          int v = present[ip];
          bool minus = ((mask >> ip) & 1) != 0;
          int last = -1;
          for (int e = 0; e < j; e++) {
            if (levels[e] == v) {
              M[e * j + e] = PowerInt64(p, v);
              last = e;
            }
          }
          if (minus) {
            M[last * j + last] = MulMod(PowerInt64(p, v), nonres, so.mod);
          }
        }
        std::vector<int> key = SymmetricOrbitKey(M, j, p, k, so.mod);
        if (so.index.count(key) == 0) {
          so.index[key] = so.n_orbit();
          so.rep.push_back(M);
        }
      }
      return;
    }
    for (int v = min_level; v <= k; v++) {
      levels[pos] = v;
      rec(pos + 1, v);
    }
  };
  rec(0, 0);
  return so;
}

// H_c for the class (e, u) of a coefficient c = p^e u.
inline std::vector<std::vector<int64_t>> const &
SymmetricOrbitTransitions(SymmetricOrbits &so, int e, int64_t u) {
  auto key = std::make_pair(e, u);
  auto iter = so.H.find(key);
  if (iter != so.H.end()) {
    return iter->second;
  }
  int j = so.j;
  int64_t mod = so.mod;
  int64_t c = (e >= so.k) ? 0 : MulMod(PowerInt64(so.p, e), u, mod);
  int no = so.n_orbit();
  std::vector<std::vector<int64_t>> H(no, std::vector<int64_t>(no, 0));
  std::vector<int64_t> x(j, 0);
  int64_t n_vec = 1;
  for (int i = 0; i < j; i++) {
    n_vec *= mod;
  }
  for (int C = 0; C < no; C++) {
    std::vector<int64_t> const &SC = so.rep[C];
    for (int64_t idx = 0; idx < n_vec; idx++) {
      int64_t w = idx;
      for (int i = 0; i < j; i++) {
        x[i] = w % mod;
        w /= mod;
      }
      std::vector<int64_t> M = SC;
      if (c != 0) {
        for (int r = 0; r < j; r++) {
          int64_t cx = MulMod(c, x[r], mod);
          for (int s = 0; s < j; s++) {
            M[r * j + s] = ModPositive(M[r * j + s] - MulMod(cx, x[s], mod), mod);
          }
        }
      }
      auto it = so.index.find(SymmetricOrbitKey(M, j, so.p, so.k, mod));
      if (it == so.index.end()) {
        std::cerr << "LDBP: a symmetric matrix outside the listed orbits\n";
        throw TerminalException{1};
      }
      H[C][it->second]++;
    }
  }
  return so.H.emplace(key, std::move(H)).first->second;
}


// The orbit tables, one per (p, k, j), filled lazily and shared by all the
// counts of a computation.
struct SymmetricOrbitsCache {
  std::map<std::tuple<int64_t, int, int>, SymmetricOrbits> map;
  SymmetricOrbits &get(int64_t p, int k, int j) {
    auto key = std::make_tuple(p, k, j);
    auto iter = map.find(key);
    if (iter != map.end()) {
      return iter->second;
    }
    return map.emplace(key, BuildSymmetricOrbits(p, k, j)).first->second;
  }
};

// clang-format off
#endif  // SRC_GENERA_SYMMETRICORBITCOUNTING_H_
// clang-format on
