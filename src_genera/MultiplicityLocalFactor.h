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
SolutionsModP(ModMat const &A, int64_t p, size_t max_solution) {
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
      SolutionsModP(A, p, bounds.max_solution);
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

/*
  The product over the primes of [O(X_p) : Stab(L_p)].
  A: the Gram matrix of L; Xb: a basis of X (rows, coordinates of L); E: a
  basis of L_rho = L cap V_rho; elts: the elements of G (on the rows).
  Nothing when the bounds are exceeded.
 */
template <typename T, typename Tint>
std::optional<T> MultiplicityLocalFactor(MyMatrix<T> const &A,
                                         MyMatrix<T> const &Xb,
                                         MyMatrix<T> const &E,
                                         std::vector<MyMatrix<T>> const &elts,
                                         MultiplicityFactorBounds const &bounds,
                                         [[maybe_unused]] std::ostream &os) {
  int m = Xb.rows();
  int k = E.rows();
  int n = A.cols();
  if (k % m != 0) {
    return {};
  }
  int w = k / m;
  // the gamma_a: gamma x_1 independent
  std::vector<size_t> ListGamma;
  {
    MyMatrix<T> rows(0, n);
    for (size_t ig = 0; ig < elts.size() && static_cast<int>(ListGamma.size()) < w; ig++) {
      MyVector<T> v = elts[ig].transpose() * GetMatrixRow(Xb, 0);
      MyMatrix<T> test(rows.rows() + 1, n);
      for (int i = 0; i < rows.rows(); i++) {
        for (int j = 0; j < n; j++) {
          test(i, j) = rows(i, j);
        }
      }
      for (int j = 0; j < n; j++) {
        test(rows.rows(), j) = v(j);
      }
      if (RankMat(test) == test.rows()) {
        rows = test;
        ListGamma.push_back(ig);
      }
    }
  }
  if (static_cast<int>(ListGamma.size()) != w) {
    return {};
  }
  // the basis B_V = {gamma_a x_i}, index a * m + i
  MyMatrix<T> BV(k, n);
  for (int a = 0; a < w; a++) {
    MyMatrix<T> Y = Xb * elts[ListGamma[a]];
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        BV(a * m + i, j) = Y(i, j);
      }
    }
  }
  if (RankMat(BV) != k) {
    return {};
  }
  // coordinates in B_V of vectors of V_rho
  auto coords = [&](MyMatrix<T> const &M) -> MyMatrix<T> {
    // M = C BV, C = M BV^T (BV BV^T)^{-1}
    MyMatrix<T> BB = BV * BV.transpose();
    return M * BV.transpose() * Inverse(BB);
  };
  // N = Z[G] X
  MyMatrix<T> Ngen(static_cast<int>(elts.size()) * m, n);
  for (size_t ig = 0; ig < elts.size(); ig++) {
    MyMatrix<T> Y = Xb * elts[ig];
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        Ngen(static_cast<int>(ig) * m + i, j) = Y(i, j);
      }
    }
  }
  MyMatrix<Tint> Ngen_i = UniversalMatrixConversion<Tint, T>(Ngen);
  MyMatrix<T> Nb = UniversalMatrixConversion<T, Tint>(GetZbasis(Ngen_i));
  if (Nb.rows() != k) {
    return {};
  }
  MyMatrix<T> GramN = Nb * A * Nb.transpose();
  MyMatrix<T> Nsharp = Inverse(GramN) * Nb;
  MyMatrix<T> N_c = coords(Nb);
  MyMatrix<T> Ns_c = coords(Nsharp);
  MyMatrix<T> L_c = coords(E);
  // N in the coordinates of N^#: R = Nb Nsharp^{-1}, integral
  MyMatrix<T> Rm = N_c * Inverse(Ns_c);
  MyMatrix<Tint> R_i = UniversalMatrixConversion<Tint, T>(Rm);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(R_i);
  MyMatrix<Tint> Dm = snf.ROW * R_i * snf.COL;
  MyMatrix<T> COL_T = UniversalMatrixConversion<T, Tint>(snf.COL);
  // f = COL^{-1} Nsharp: the coordinates of x in N^# are x_F COL
  MyMatrix<T> F_c = Inverse(COL_T) * Ns_c;
  MyMatrix<T> F_inv = Inverse(F_c);
  // the index [L_rho : N]
  T index = T_abs(DeterminantMat(MyMatrix<T>(N_c * Inverse(L_c))));
  T factor(1);
  if (index == T(1)) {
    return factor;
  }
  std::vector<Tint> primes =
      FactorsInt(UniversalScalarConversion<Tint, T>(index));
  std::sort(primes.begin(), primes.end());
  primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
  MyMatrix<T> GramX = Xb * A * Xb.transpose();
  for (auto const &p_i : primes) {
    int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
    // the p-part of N^# / N: coordinates i with p | d_i, orders p^{v_i}
    std::vector<int> ListIdx;
    std::vector<int64_t> ListOrd;
    int c = 0;
    for (int i = 0; i < k; i++) {
      T d = T_abs(UniversalScalarConversion<T, Tint>(Dm(i, i)));
      if (d == T(0)) {
        return {};
      }
      int v = ValuationOfRational(d, p);
      if (v > 0) {
        ListIdx.push_back(i);
        int64_t o = 1;
        for (int u = 0; u < v; u++) {
          o *= p;
        }
        ListOrd.push_back(o);
        c = std::max(c, v);
      }
    }
    size_t r = ListIdx.size();
    // a vector of V_rho (B_V coordinates, p-integral in N^#) to its class
    auto residue = [&](MyVector<T> const &x) -> std::vector<int64_t> {
      MyVector<T> y = F_inv.transpose() * x;
      std::vector<int64_t> res(r);
      for (size_t s = 0; s < r; s++) {
        T val = y(ListIdx[s]);
        // p-integral: numerator times the inverse of the denominator
        Tint num = val.get_num();
        Tint den = val.get_den();
        Tint o_T = UniversalScalarConversion<Tint, int64_t>(ListOrd[s]);
        Tint nm = ((num % o_T) + o_T) % o_T;
        Tint dn = ((den % o_T) + o_T) % o_T;
        int64_t nm64 = UniversalScalarConversion<int64_t, Tint>(nm);
        int64_t dn64 = UniversalScalarConversion<int64_t, Tint>(dn);
        res[s] = MulModQ(nm64, InvModQ(dn64, ListOrd[s]), ListOrd[s]);
      }
      return res;
    };
    size_t q_size = 1;
    for (auto o : ListOrd) {
      q_size *= static_cast<size_t>(o);
      if (q_size > bounds.max_module) {
        return {};
      }
    }
    auto encode = [&](std::vector<int64_t> const &x) -> size_t {
      size_t idx = 0;
      for (size_t s = r; s-- > 0;) {
        idx = idx * static_cast<size_t>(ListOrd[s]) + static_cast<size_t>(x[s]);
      }
      return idx;
    };
    // the subgroup generated by vectors, as the sorted set of its elements
    auto span = [&](std::vector<std::vector<int64_t>> const &gens)
        -> std::vector<size_t> {
      std::set<size_t> seen;
      std::vector<std::vector<int64_t>> queue{std::vector<int64_t>(r, 0)};
      seen.insert(0);
      while (!queue.empty()) {
        std::vector<int64_t> x = queue.back();
        queue.pop_back();
        for (auto const &g : gens) {
          std::vector<int64_t> y(r);
          for (size_t s = 0; s < r; s++) {
            y[s] = (x[s] + g[s]) % ListOrd[s];
          }
          if (seen.insert(encode(y)).second) {
            queue.push_back(y);
          }
        }
      }
      return std::vector<size_t>(seen.begin(), seen.end());
    };
    // L / N
    std::vector<std::vector<int64_t>> Lgens;
    for (int i = 0; i < k; i++) {
      Lgens.push_back(residue(GetMatrixRow(L_c, i)));
    }
    // the image of O(X_p) modulo p^c
    int e = 0;
    {
      MyMatrix<T> Ainv = Inverse(GramX);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          if (Ainv(i, j) != T(0)) {
            T den = GetDenominator(Ainv(i, j));
            e = std::max(e, ValuationOfRational(den, p));
          }
        }
      }
    }
    ModMat A_mod(m, std::vector<int64_t>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        A_mod[i][j] = UniversalScalarConversion<int64_t, T>(GramX(i, j));
      }
    }
    std::optional<std::vector<ModMat>> opt =
        LocalOrthogonalImage(A_mod, p, c, e, bounds);
    if (!opt) {
      return {};
    }
    // the action of an element on the generators f_s of the p-part: the
    // rows f_s (B_V coordinates) moved by 1 (x) g
    auto act_gens = [&](ModMat const &g) -> std::vector<std::vector<int64_t>> {
      std::vector<std::vector<int64_t>> res;
      for (size_t s = 0; s < r; s++) {
        MyVector<T> f = GetMatrixRow(F_c, ListIdx[s]);
        MyVector<T> img = ZeroVector<T>(k);
        for (int a = 0; a < w; a++) {
          for (int i = 0; i < m; i++) {
            T x = f(a * m + i);
            if (x == T(0)) {
              continue;
            }
            for (int j = 0; j < m; j++) {
              // rows convention: x_i -> sum_j g(j, i) x_j, g acting on the
              // columns of coordinates of X
              img(a * m + j) += x * UniversalScalarConversion<T, int64_t>(g[j][i]);
            }
          }
        }
        res.push_back(residue(img));
      }
      return res;
    };
    std::vector<std::vector<std::vector<int64_t>>> ListAct;
    for (auto const &g : *opt) {
      ListAct.push_back(act_gens(g));
    }
    auto apply = [&](std::vector<std::vector<int64_t>> const &act,
                     std::vector<int64_t> const &x) -> std::vector<int64_t> {
      std::vector<int64_t> y(r, 0);
      for (size_t s = 0; s < r; s++) {
        if (x[s] == 0) {
          continue;
        }
        for (size_t u = 0; u < r; u++) {
          y[u] = (y[u] + MulModQ(x[s], act[s][u], ListOrd[u])) % ListOrd[u];
        }
      }
      return y;
    };
    auto decode = [&](size_t idx) -> std::vector<int64_t> {
      std::vector<int64_t> x(r);
      for (size_t s = 0; s < r; s++) {
        x[s] = static_cast<int64_t>(idx % static_cast<size_t>(ListOrd[s]));
        idx /= static_cast<size_t>(ListOrd[s]);
      }
      return x;
    };
    std::vector<size_t> L0 = span(Lgens);
    std::set<std::vector<size_t>> orbit{L0};
    std::vector<std::vector<size_t>> queue{L0};
    while (!queue.empty()) {
      std::vector<size_t> H = queue.back();
      queue.pop_back();
      for (auto const &act : ListAct) {
        std::vector<size_t> img;
        for (size_t idx : H) {
          img.push_back(encode(apply(act, decode(idx))));
        }
        std::sort(img.begin(), img.end());
        if (orbit.insert(img).second) {
          queue.push_back(img);
        }
      }
    }
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    os << "MULT_FACTOR: p=" << p << " c=" << c << " e=" << e
       << " |image of O(X_p)|=" << opt->size() << " |N^#/N|_p=" << q_size
       << " |L/N|_p=" << L0.size() << " orbit=" << orbit.size() << "\n";
#endif
    factor *= UniversalScalarConversion<T, int64_t>(
        static_cast<int64_t>(orbit.size()));
  }
  return factor;
}

// clang-format off
#endif  // SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
// clang-format on
