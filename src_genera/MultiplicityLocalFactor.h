// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
#define SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_

// clang-format off
#include "GroupNeighborModP.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_MULTIPLICITY_LOCAL_FACTOR
#endif

/*
  The local factors of the mass of a component V_rho = W (x) U of a finite
  group G with End_G(V_rho) = Q (IsometryMass.h, FiniteGroupComponents).

  X = L cap (l (x) U) is the multiplicity lattice, N = Z[G] X subset
  L_rho. The centraliser C_rho = O(U) commutes with G, so the stabiliser
  of L_rho in C_rho(Q_p) fixes X_p = L_p cap (l (x) U) and N_p, and

      [O(X_p) : Stab(L_p)] = |O(X_p) . (L_p / N_p)|,

  the orbit of L_p / N_p in the finite p-group N_p^# / N_p. At p not
  dividing |G| the order Z_p[G] is maximal, L_p = N_p and the factor is 1.

  O(X_p) acts on N_p^# / N_p through X / p^c X, p^c the exponent of
  N_p^# / N_p: if g' = g mod p^c then (g' - g) N^# subset p^{-c} (g' - g) N
  subset N. The image of O(X_p) in GL(X / p^c X) is found from the Newton
  step g -> g - (1/2) A^{-1} g^{-T} E, E = g^T A g - A, which leaves an
  error (1/4) E g^{-1} A^{-1} g^{-T} E: if v_p(E) = k > e + 2t (p^e A^{-1}
  integral, t = v_p(2)), there is an isometry h = g mod p^{k-e-t}. So the
  image modulo p^c is the set of the solutions g of g^T A g = A modulo
  p^c that lift to solutions modulo p^K, K = max(c + e + t, e + 2t + 1).
  The solutions are lifted level by level, g -> g + p^j Delta, Delta
  solving the linear system Delta^T A g + g^T A Delta = -E / p^j modulo p.

  On V_rho the element g of O(X_p) acts as 1 (x) g; in the basis
  {gamma_a x_i} of V_rho (gamma_a in G with gamma_a x_1 a basis of the
  W-side, x_i a basis of X) it is block diagonal.
 */

struct MultiplicityFactorBounds {
  size_t max_solution = 2000000;
  size_t max_module = 1000000;
  // the nodes of the backtrack of the solutions modulo p
  size_t max_node = 5000000;
};

/*
  The solutions modulo p^{j+1} lifting g (a solution modulo p^j), j >= 1:
  g + p^j Delta for the Delta modulo p of the affine solution space.
 */
inline std::vector<ModMat> LiftSolutions(ModMat const &g, ModMat const &A,
                                         int64_t p, int j, int64_t pj) {
  size_t m = A.size();
  int64_t pj1 = pj * p;
  (void)j;
  // E = g^T A g - A, exactly over the integers of the representatives
  std::vector<std::vector<__int128>> Ag(m, std::vector<__int128>(m, 0));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      __int128 s = 0;
      for (size_t c = 0; c < m; c++) {
        s += static_cast<__int128>(A[a][c]) * g[c][b];
      }
      Ag[a][b] = s;
    }
  }
  std::vector<std::vector<int64_t>> Ered(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      __int128 s = 0;
      for (size_t c = 0; c < m; c++) {
        s += static_cast<__int128>(g[c][a]) * Ag[c][b];
      }
      s -= A[a][b];
      __int128 r = s % pj1;
      if (r < 0) {
        r += pj1;
      }
      if (r % pj != 0) {
        std::cerr << "MULT_FACTOR: g is not a solution modulo p^j\n";
        throw TerminalException{1};
      }
      Ered[a][b] = static_cast<int64_t>(r / pj);
    }
  }
  // the linear system over the m^2 entries of Delta, equations a <= b:
  // sum_c Delta(c, a) (A g)(c, b) + (A g)(c, a) Delta(c, b) = -E(a,b)/p^j
  size_t n_var = m * m;
  ModMat coef;
  std::vector<int64_t> rhs;
  for (size_t a = 0; a < m; a++) {
    for (size_t b = a; b < m; b++) {
      std::vector<int64_t> row(n_var, 0);
      for (size_t c = 0; c < m; c++) {
        int64_t x = ModP(static_cast<int64_t>(Ag[c][b] % p), p);
        row[c * m + a] = ModP(row[c * m + a] + x, p);
        int64_t y = ModP(static_cast<int64_t>(Ag[c][a] % p), p);
        row[c * m + b] = ModP(row[c * m + b] + y, p);
      }
      coef.push_back(row);
      rhs.push_back(ModP(-Ered[a][b], p));
    }
  }
  // the affine solutions: x with coef x = rhs
  size_t n_eq = coef.size();
  ModMat aug(n_eq, std::vector<int64_t>(n_var + 1));
  for (size_t e = 0; e < n_eq; e++) {
    for (size_t v = 0; v < n_var; v++) {
      aug[e][v] = coef[e][v];
    }
    aug[e][n_var] = rhs[e];
  }
  ModRowSpace rs = ModRowReduce(aug, p);
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    if (rs.pivots[i] == n_var) {
      return {};
    }
  }
  std::vector<int64_t> part(n_var, 0);
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    part[rs.pivots[i]] = rs.basis[i][n_var];
  }
  std::vector<bool> is_pivot(n_var, false);
  for (size_t piv : rs.pivots) {
    is_pivot[piv] = true;
  }
  ModMat ker;
  for (size_t f = 0; f < n_var; f++) {
    if (is_pivot[f]) {
      continue;
    }
    std::vector<int64_t> v(n_var, 0);
    v[f] = 1;
    for (size_t i = 0; i < rs.pivots.size(); i++) {
      v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
    }
    ker.push_back(v);
  }
  std::vector<ModMat> res;
  size_t n_sol = 1;
  for (size_t i = 0; i < ker.size(); i++) {
    n_sol *= static_cast<size_t>(p);
  }
  for (size_t s = 0; s < n_sol; s++) {
    std::vector<int64_t> x = part;
    size_t t = s;
    for (auto const &kv : ker) {
      int64_t c = static_cast<int64_t>(t % static_cast<size_t>(p));
      t /= static_cast<size_t>(p);
      for (size_t v = 0; v < n_var; v++) {
        x[v] = (x[v] + c * kv[v]) % p;
      }
    }
    ModMat g2 = g;
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        g2[a][b] = ModP(g[a][b] + pj * x[a * m + b], pj1);
      }
    }
    res.push_back(g2);
  }
  return res;
}

// The solutions modulo p of g^T A g = A with g invertible, column by column.
inline std::optional<std::vector<ModMat>>
SolutionsModP(ModMat const &A, int64_t p, size_t max_solution,
              size_t max_node) {
  size_t m = A.size();
  std::vector<std::vector<int64_t>> vects;
  size_t n_vect = 1;
  for (size_t i = 0; i < m; i++) {
    n_vect *= static_cast<size_t>(p);
  }
  for (size_t u = 0; u < n_vect; u++) {
    std::vector<int64_t> v(m);
    size_t t = u;
    for (size_t i = 0; i < m; i++) {
      v[i] = static_cast<int64_t>(t % static_cast<size_t>(p));
      t /= static_cast<size_t>(p);
    }
    vects.push_back(v);
  }
  auto bil = [&](std::vector<int64_t> const &x,
                 std::vector<int64_t> const &y) -> int64_t {
    int64_t s = 0;
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        s = (s + x[a] * ModP(A[a][b], p) % p * y[b]) % p;
      }
    }
    return s;
  };
  std::vector<ModMat> res;
  std::vector<std::vector<int64_t>> cols(m);
  bool overflow = false;
  size_t n_node = 0;
  std::function<void(size_t)> rec = [&](size_t j) -> void {
    if (overflow) {
      return;
    }

    if (j == m) {
      ModMat g(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          g[a][b] = cols[b][a];
        }
      }
      if (ModRank(g, p) == m) {
        res.push_back(g);
        if (res.size() > max_solution) {
          overflow = true;
        }
      }
      return;
    }
    for (auto const &v : vects) {
      // the budget counts the candidate columns tried
      if (++n_node > max_node) {
        overflow = true;
        return;
      }
      bool ok = true;
      for (size_t i = 0; i <= j && ok; i++) {
        std::vector<int64_t> const &w = (i == j) ? v : cols[i];
        if (bil(w, v) != ModP(A[i][j], p)) {
          ok = false;
        }
      }
      if (ok) {
        cols[j] = v;
        rec(j + 1);
      }
    }
  };
  rec(0);
  if (overflow) {
    return {};
  }
  return res;
}

/*
  The image of O(X_p) in GL_m(Z / p^c), X with Gram matrix A (integral);
  nothing beyond the bounds.
 */
inline std::optional<std::vector<ModMat>>
LocalOrthogonalImage(ModMat const &A, int64_t p, int c, int e,
                     MultiplicityFactorBounds const &bounds) {
  int t = (p == 2) ? 1 : 0;
  int K = std::max(c + e + t, e + 2 * t + 1);
  std::optional<std::vector<ModMat>> opt =
      SolutionsModP(A, p, bounds.max_solution, bounds.max_node);
  if (!opt) {
    return {};
  }
  std::vector<ModMat> level = *opt;
  int64_t pj = p;
  for (int j = 1; j < c; j++) {
    std::vector<ModMat> next;
    for (auto const &g : level) {
      std::vector<ModMat> lifts = LiftSolutions(g, A, p, j, pj);
      next.insert(next.end(), lifts.begin(), lifts.end());
      if (next.size() > bounds.max_solution) {
        return {};
      }
    }
    level = std::move(next);
    pj *= p;
  }
  // keep those that lift to level K
  std::function<bool(ModMat const &, int, int64_t)> liftable =
      [&](ModMat const &g, int j, int64_t pj_loc) -> bool {
    if (j >= K) {
      return true;
    }
    std::vector<ModMat> lifts = LiftSolutions(g, A, p, j, pj_loc);
    for (auto const &g2 : lifts) {
      if (liftable(g2, j + 1, pj_loc * p)) {
        return true;
      }
    }
    return false;
  };
  std::vector<ModMat> res;
  for (auto const &g : level) {
    if (liftable(g, c, pj)) {
      res.push_back(g);
    }
  }
  return res;
}

/*
  The image modulo p^c of the stabiliser {g in Lambda (x) Z_p : g A g^T =
  A} (rows convention), Lambda the Z-span of the m x m integer matrices of
  basis: the centraliser of a group acting on a lattice of Gram matrix A,
  when Lambda is its commutant. The Newton step g -> g - (1/2) E g^{-T}
  A^{-1}, E = g A g^T - A, stays in the commutant when g does, since the
  group preserves A and E; so Lemma (multfactor) holds within Lambda, and
  the image is the set of the solutions modulo p^c in Lambda that lift to
  solutions modulo p^K, lifted level by level by Delta in Lambda / p
  Lambda. Nothing beyond the bounds.
 */
inline std::optional<std::vector<ModMat>>
LocalStabilizerImage(ModMat const &A, std::vector<ModMat> const &basis,
                     int64_t p, int c, int e,
                     MultiplicityFactorBounds const &bounds) {
  // A and basis are exact integer matrices: a solution is kept by its
  // coordinates in the basis, so that it stays in Lambda (the entrywise
  // reduction of an element of Lambda need not be in Lambda, and the
  // lifts g + p^j Delta, Delta in Lambda, would miss the isometries)
  size_t m = A.size();
  size_t dim = basis.size();
  int t = (p == 2) ? 1 : 0;
  int K = std::max(c + e + t, e + 2 * t + 1);
  using Coord = std::vector<int64_t>;
  auto matrix_of = [&](Coord const &x) -> std::vector<std::vector<__int128>> {
    std::vector<std::vector<__int128>> g(m, std::vector<__int128>(m, 0));
    for (size_t i = 0; i < dim; i++) {
      if (x[i] == 0) {
        continue;
      }
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          g[a][b] += static_cast<__int128>(x[i]) * basis[i][a][b];
        }
      }
    }
    return g;
  };
  // E = g A g^T - A, exactly
  auto error = [&](std::vector<std::vector<__int128>> const &g)
      -> std::vector<std::vector<__int128>> {
    std::vector<std::vector<__int128>> gA(m, std::vector<__int128>(m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        __int128 sum = 0;
        for (size_t u = 0; u < m; u++) {
          sum += g[a][u] * A[u][b];
        }
        gA[a][b] = sum;
      }
    }
    std::vector<std::vector<__int128>> E(m, std::vector<__int128>(m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        __int128 sum = 0;
        for (size_t u = 0; u < m; u++) {
          sum += gA[a][u] * g[b][u];
        }
        E[a][b] = sum - A[a][b];
      }
    }
    return E;
  };
  auto mod_p = [&](__int128 v) -> int64_t {
    __int128 r = v % p;
    if (r < 0) {
      r += p;
    }
    return static_cast<int64_t>(r);
  };
  // level 1: all x in {0, ..., p-1}^dim
  size_t n_x = 1;
  for (size_t i = 0; i < dim; i++) {
    n_x *= static_cast<size_t>(p);
    if (n_x > bounds.max_solution * 10) {
      return {};
    }
  }
  std::vector<Coord> level;
  for (size_t u = 0; u < n_x; u++) {
    Coord x(dim);
    size_t tt = u;
    for (size_t i = 0; i < dim; i++) {
      x[i] = static_cast<int64_t>(tt % static_cast<size_t>(p));
      tt /= static_cast<size_t>(p);
    }
    std::vector<std::vector<__int128>> g = matrix_of(x);
    std::vector<std::vector<__int128>> E = error(g);
    bool ok = true;
    for (size_t a = 0; a < m && ok; a++) {
      for (size_t b = 0; b < m && ok; b++) {
        if (E[a][b] % p != 0) {
          ok = false;
        }
      }
    }
    if (!ok) {
      continue;
    }
    ModMat gp(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        gp[a][b] = mod_p(g[a][b]);
      }
    }
    if (ModRank(gp, p) == m) {
      level.push_back(x);
      if (level.size() > bounds.max_solution) {
        return {};
      }
    }
  }
  // the lifts x + p^j y, y modulo p solving the linear system
  auto lifts = [&](Coord const &x, int64_t pj) -> std::vector<Coord> {
    int64_t pj1 = pj * p;
    std::vector<std::vector<__int128>> g = matrix_of(x);
    std::vector<std::vector<__int128>> E = error(g);
    std::vector<std::vector<int64_t>> gA(m, std::vector<int64_t>(m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        __int128 sum = 0;
        for (size_t u = 0; u < m; u++) {
          sum += g[a][u] * A[u][b];
        }
        gA[a][b] = mod_p(sum);
      }
    }
    ModMat aug;
    for (size_t a = 0; a < m; a++) {
      for (size_t b = a; b < m; b++) {
        std::vector<int64_t> row(dim + 1, 0);
        for (size_t i = 0; i < dim; i++) {
          // (Delta A g^T + g A Delta^T)(a, b), A symmetric
          int64_t v = 0;
          for (size_t u = 0; u < m; u++) {
            v += mod_p(basis[i][a][u]) * gA[b][u] % p;
            v += gA[a][u] * mod_p(basis[i][b][u]) % p;
          }
          row[i] = ModP(v, p);
        }
        __int128 r = E[a][b] % pj1;
        if (r < 0) {
          r += pj1;
        }
        if (r % pj != 0) {
          std::cerr << "MULT_FACTOR: not a solution modulo p^j\n";
          throw TerminalException{1};
        }
        row[dim] = ModP(-static_cast<int64_t>(r / pj), p);
        aug.push_back(row);
      }
    }
    ModRowSpace rs = ModRowReduce(aug, p);
    for (size_t piv : rs.pivots) {
      if (piv == dim) {
        return {};
      }
    }
    Coord part(dim, 0);
    for (size_t i = 0; i < rs.pivots.size(); i++) {
      part[rs.pivots[i]] = rs.basis[i][dim];
    }
    std::vector<bool> is_pivot(dim, false);
    for (size_t piv : rs.pivots) {
      is_pivot[piv] = true;
    }
    std::vector<Coord> ker;
    for (size_t f = 0; f < dim; f++) {
      if (is_pivot[f]) {
        continue;
      }
      Coord v(dim, 0);
      v[f] = 1;
      for (size_t i = 0; i < rs.pivots.size(); i++) {
        v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
      }
      ker.push_back(v);
    }
    size_t n_sol = 1;
    for (size_t i = 0; i < ker.size(); i++) {
      n_sol *= static_cast<size_t>(p);
    }
    std::vector<Coord> res;
    for (size_t u = 0; u < n_sol; u++) {
      Coord y = part;
      size_t tt = u;
      for (auto const &kv : ker) {
        int64_t cc = static_cast<int64_t>(tt % static_cast<size_t>(p));
        tt /= static_cast<size_t>(p);
        for (size_t i = 0; i < dim; i++) {
          y[i] = (y[i] + cc * kv[i]) % p;
        }
      }
      Coord x2 = x;
      for (size_t i = 0; i < dim; i++) {
        x2[i] = x[i] + pj * y[i];
      }
      res.push_back(x2);
    }
    return res;
  };
  int64_t pj = p;
  for (int j = 1; j < c; j++) {
    std::vector<Coord> next;
    for (auto const &x : level) {
      std::vector<Coord> l = lifts(x, pj);
      next.insert(next.end(), l.begin(), l.end());
      if (next.size() > bounds.max_solution) {
        return {};
      }
    }
    level = std::move(next);
    pj *= p;
  }
  std::function<bool(Coord const &, int, int64_t)> liftable =
      [&](Coord const &x, int j, int64_t pj_loc) -> bool {
    if (j >= K) {
      return true;
    }
    for (auto const &x2 : lifts(x, pj_loc)) {
      if (liftable(x2, j + 1, pj_loc * p)) {
        return true;
      }
    }
    return false;
  };
  std::vector<ModMat> res;
  for (auto const &x : level) {
    if (liftable(x, std::max(c, 1), pj)) {
      std::vector<std::vector<__int128>> g = matrix_of(x);
      ModMat gm(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          gm[a][b] = static_cast<int64_t>(g[a][b]);
        }
      }
      res.push_back(gm);
    }
  }
  return res;
}

template <typename T> int ValuationOfRational(T const &x, int64_t p) {
  // x nonzero
  T y = x;
  T p_T = UniversalScalarConversion<T, int64_t>(p);
  int v = 0;
  while (IsInteger(y / p_T)) {
    y /= p_T;
    v++;
  }
  return v;
}

// clang-format off
#endif  // SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
// clang-format on
