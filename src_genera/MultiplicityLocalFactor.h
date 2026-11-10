// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
#define SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_

// clang-format off
#include "GroupNeighborModP.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <tuple>
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
  // the order of the p-part of N^# / N, which only enters an encoding of
  // its elements (the orbits are of the subgroup L / N, spanned)
  size_t max_module = size_t(1) << 62;
  // the nodes of the backtrack of the solutions modulo p
  size_t max_node = 5000000;
};

/*
  The solutions modulo p^{j+1} lifting g (a solution modulo p^j), j >= 1:
  g + p^j Delta for the Delta modulo p of an affine space, given by a
  particular solution and a basis of the kernel (nothing when it is
  empty).
 */
struct LiftAffineSpace {
  std::vector<int64_t> part;
  ModMat ker;
};

// the linear system of the lifts of g from p^j to p^{j+1}: coef Delta =
// rhs modulo p, the equations indexed by a <= b
struct LiftSystem {
  ModMat coef;
  std::vector<int64_t> rhs;
};

inline LiftSystem MakeLiftSystem(ModMat const &g, ModMat const &A, int64_t p,
                                 int64_t pj) {
  size_t m = A.size();
  int64_t pj1 = pj * p;
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
  return LiftSystem{coef, rhs};
}

inline std::optional<LiftAffineSpace> LiftAffine(ModMat const &g,
                                                 ModMat const &A, int64_t p,
                                                 int64_t pj) {
  size_t m = A.size();
  size_t n_var = m * m;
  LiftSystem sys = MakeLiftSystem(g, A, p, pj);
  ModMat const &coef = sys.coef;
  std::vector<int64_t> const &rhs = sys.rhs;
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
  return LiftAffineSpace{part, ker};
}

// the element g + p^j Delta, Delta = part + sum_i coef_i ker_i
inline ModMat LiftAffineElement(ModMat const &g, LiftAffineSpace const &sp,
                                std::vector<int64_t> const &coef, int64_t p,
                                int64_t pj) {
  size_t m = g.size();
  int64_t pj1 = pj * p;
  std::vector<int64_t> x = sp.part;
  for (size_t i = 0; i < sp.ker.size(); i++) {
    if (coef[i] != 0) {
      for (size_t v = 0; v < x.size(); v++) {
        x[v] = (x[v] + coef[i] * sp.ker[i][v]) % p;
      }
    }
  }
  ModMat g2 = g;
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      g2[a][b] = ModP(g[a][b] + pj * x[a * m + b], pj1);
    }
  }
  return g2;
}

inline std::vector<ModMat> LiftSolutions(ModMat const &g, ModMat const &A,
                                         int64_t p, int j, int64_t pj) {
  (void)j;
  std::optional<LiftAffineSpace> sp = LiftAffine(g, A, p, pj);
  if (!sp) {
    return {};
  }
  std::vector<ModMat> res;
  size_t n_sol = 1;
  for (size_t i = 0; i < sp->ker.size(); i++) {
    n_sol *= static_cast<size_t>(p);
  }
  std::vector<int64_t> coef(sp->ker.size());
  for (size_t s = 0; s < n_sol; s++) {
    size_t t = s;
    for (size_t i = 0; i < coef.size(); i++) {
      coef[i] = static_cast<int64_t>(t % static_cast<size_t>(p));
      t /= static_cast<size_t>(p);
    }
    res.push_back(LiftAffineElement(g, *sp, coef, p, pj));
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
  The searches below keep the residues modulo p^K in machine integers and
  form g A g^T exactly in __int128: this holds when m^2 max|A| p^{2K} stays
  below 2^120. Beyond that they return nothing (the factor is then not
  computed), rather than overflow.
 */
inline bool ResiduesFitInt128(ModMat const &A, int64_t p, int K) {
  double bits = 0;
  int64_t mx = 1;
  for (auto const &row : A) {
    for (int64_t x : row) {
      mx = std::max(mx, (x >= 0) ? x : -x);
    }
  }
  bits += std::log2(static_cast<double>(mx));
  bits += 2.0 * K * std::log2(static_cast<double>(p));
  bits += 2.0 * std::log2(static_cast<double>(A.size()) + 1);
  return bits < 120.0;
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
  if (!ResiduesFitInt128(A, p, K)) {
    return {};
  }
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

// the precision W of the Newton steps: p^W < 2^61, so that the products
// of residues fit in 128 bits
inline int MultiplicityNewtonPrecision(int64_t p) {
  int W = 1;
  int64_t pw = p;
  while (pw <= (int64_t(1) << 61) / p) {
    pw *= p;
    W++;
  }
  return W;
}

/*
  Generators of the image G of O(X_p) in GL_m(Z / p^c) (columns: g^T A g =
  A), complete by construction and without listing G, which grows like
  p^{m(m-1)/2}:

  - G_1 = G mod p by its stabiliser chain on the basis e_0, ..., e_{m-1}.
    From the last level up, the orbit of e_i under the elements found so
    far is closed, and for each admissible image y of e_i (B(e_k, y) =
    A_{k i} for k < i, B(y, y) = A_{i i}) outside it, an element fixing
    e_0, ..., e_{i-1} with e_i -> y is searched: the further columns one
    at a time on their affine spaces, the matrix invertible, and a lift to
    p^K (Lemma: then it comes from an isometry). Once the orbit holds every
    y with such an element, the elements found generate the stabiliser H
    of e_0, ..., e_{i-1} in G_1, as they contain the stabiliser of e_i and
    move e_i over its whole H-orbit.
  - The layers G_j / G_{j+1}, 1 <= j < c, of the elements congruent to 1
    modulo p^j: F_p-spaces of the Delta with 1 + p^j Delta lifting to p^K,
    given by a basis of the solutions when its elements all lift (the
    liftable ones form a subspace), else by enumerating the space.

  An element of G congruent to a solution g modulo p^j is found from a
  solution modulo p^k, k = max(j + e + t, e + 2t + 1), extending g (an
  affine system at each level, a single one for p odd and e <= 1) and the
  Newton steps of the Lemma, Be = p^e A^{-1} being given modulo p^W, W as
  large as 64-bit products allow; g mod p is in G_1 exactly when such a
  solution exists.

  The lifts modulo p^c of these elements generate G (by the filtration
  G_1 = G mod p, G_j / G_{j+1}), and |G| is the product of the orbit
  lengths and of the orders of the layers (order, when asked). Nothing
  beyond the bounds.
 */
inline std::optional<std::vector<ModMat>> LocalOrthogonalImageGeneratorsCore(
    ModMat const &A, ModMat const &Be, std::vector<ModMat> const &Om,
    std::vector<ModMat> const &Lam, int64_t p, int c, int e, size_t rad_start,
    MultiplicityFactorBounds const &bounds, double *order) {
  size_t m = A.size();
  int t = (p == 2) ? 1 : 0;
  size_t n_pts = 1;
  for (size_t i = 0; i < m; i++) {
    n_pts *= static_cast<size_t>(p);
    if (n_pts > (size_t(1) << 24)) {
      return {};
    }
  }
  int64_t pc = 1;
  for (int i = 0; i < c; i++) {
    pc *= p;
  }
  ModMat Ap(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      Ap[a][b] = ModP(A[a][b], p);
    }
  }
  auto bil = [&](std::vector<int64_t> const &x,
                 std::vector<int64_t> const &y) -> int64_t {
    int64_t s = 0;
    for (size_t a = 0; a < m; a++) {
      if (x[a] == 0) {
        continue;
      }
      for (size_t b = 0; b < m; b++) {
        s = (s + x[a] * Ap[a][b] % p * y[b]) % p;
      }
    }
    return s;
  };
  size_t n_node = 0;
  bool over = false;
  int W = MultiplicityNewtonPrecision(p);
  if (W < c + 4 * e + 2 * t + 4) {
    return {};
  }
  // the inverse of a modulo M, gcd(a, M) = 1
  auto inv_mod = [](int64_t a, int64_t M) -> int64_t {
    int64_t r0 = M, r1 = ((a % M) + M) % M, s0 = 0, s1 = 1;
    while (r1 != 0) {
      int64_t q = r0 / r1;
      int64_t r2 = r0 - q * r1;
      r0 = r1;
      r1 = r2;
      int64_t s2 = static_cast<int64_t>(
          (static_cast<__int128>(s0) - static_cast<__int128>(q) * s1) % M);
      s0 = s1;
      s1 = s2;
    }
    return ((s0 % M) + M) % M;
  };
  auto pow_p = [&](int k) -> int64_t {
    int64_t r = 1;
    for (int i = 0; i < k; i++) {
      r *= p;
    }
    return r;
  };
  auto mulmod = [](int64_t x, int64_t y, int64_t M) -> int64_t {
    __int128 r = static_cast<__int128>(x) * y % M;
    return static_cast<int64_t>((r < 0) ? r + M : r);
  };
  auto matmul = [&](ModMat const &X, ModMat const &Y, int64_t M) -> ModMat {
    ModMat Z(m, std::vector<int64_t>(m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t k = 0; k < m; k++) {
        if (X[a][k] == 0) {
          continue;
        }
        for (size_t b = 0; b < m; b++) {
          Z[a][b] = (Z[a][b] + mulmod(X[a][k], Y[k][b], M)) % M;
        }
      }
    }
    return Z;
  };
  auto transpose = [&](ModMat const &X) -> ModMat {
    ModMat Z(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Z[a][b] = X[b][a];
      }
    }
    return Z;
  };
  auto reduce = [&](ModMat X, int64_t M) -> ModMat {
    for (auto &row : X) {
      for (auto &x : row) {
        x = ModP(x, M);
      }
    }
    return X;
  };
  // the inverse modulo p^k of g, invertible modulo p: modulo p by Gauss,
  // then X <- X (2 - g X)
  auto inverse_mod = [&](ModMat const &g, int k) -> std::optional<ModMat> {
    ModMat aug(m, std::vector<int64_t>(2 * m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        aug[a][b] = ModP(g[a][b], p);
      }
      aug[a][m + a] = 1;
    }
    for (size_t col = 0; col < m; col++) {
      size_t piv = col;
      while (piv < m && aug[piv][col] == 0) {
        piv++;
      }
      if (piv == m) {
        return {};
      }
      std::swap(aug[piv], aug[col]);
      int64_t iv = inv_mod(aug[col][col], p);
      for (auto &x : aug[col]) {
        x = x * iv % p;
      }
      for (size_t a = 0; a < m; a++) {
        if (a != col && aug[a][col] != 0) {
          int64_t f0 = aug[a][col];
          for (size_t b = 0; b < 2 * m; b++) {
            aug[a][b] = ModP(aug[a][b] - f0 * aug[col][b], p);
          }
        }
      }
    }
    ModMat X(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        X[a][b] = aug[a][m + b];
      }
    }
    int prec = 1;
    while (prec < k) {
      prec = std::min(2 * prec, k);
      int64_t M = pow_p(prec);
      ModMat gX = matmul(reduce(g, M), X, M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          gX[a][b] = ModP(((a == b) ? 2 : 0) - gX[a][b], M);
        }
      }
      X = matmul(X, gX, M);
    }
    return X;
  };
  // the valuation of g^T A g - A, up to k
  auto error_val = [&](ModMat const &g, int k) -> int {
    int64_t M = pow_p(k);
    ModMat E = matmul(transpose(g), matmul(reduce(A, M), g, M), M);
    int v = k;
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        int64_t x = ModP(E[a][b] - A[a][b], M);
        int vx = 0;
        while (vx < k && x % p == 0) {
          x /= p;
          vx++;
        }
        v = std::min(v, vx);
      }
    }
    return v;
  };
  // Newton (Lemma): from a solution modulo p^k, k > e + 2t, to an
  // isometry known modulo p^c, the steps g -> g - (1/2) A^{-1} g^{-T} E,
  // A^{-1} = Be / p^e, the precision of the representatives decreasing by
  // e + t at each step
  auto newton = [&](ModMat g) -> std::optional<ModMat> {
    int prec = W;
    for (int it = 0; it < 64; it++) {
      int v = error_val(g, prec);
      if (v >= c + e + t) {
        return reduce(g, pc);
      }
      if (prec - e - t <= c + e + t) {
        return {};
      }
      int64_t M = pow_p(prec);
      std::optional<ModMat> gi = inverse_mod(g, prec);
      if (!gi) {
        return {};
      }
      ModMat E = matmul(transpose(g), matmul(reduce(A, M), g, M), M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          E[a][b] = ModP(E[a][b] - A[a][b], M);
        }
      }
      ModMat F = matmul(reduce(Be, M), matmul(transpose(*gi), E, M), M);
      // F / p^e (/ 2 when p = 2) holds modulo p^{prec - e - t}
      int64_t div = pow_p(e + t);
      int nprec = prec - e - t;
      int64_t NM = pow_p(nprec);
      int64_t inv2 = (p == 2) ? 1 : inv_mod(2, NM);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          int64_t x = F[a][b];
          if (x % div != 0) {
            return {};
          }
          x = (x / div) % NM;
          x = mulmod(x, inv2, NM);
          g[a][b] = ModP(g[a][b] - x, NM);
        }
      }
      prec = nprec;
    }
    return {};
  };
  // With a commutant Lam (a centre K), the lifts g + p^j Delta keep Delta
  // in Lam mod p: the system of MakeLiftSystem restricted to the
  // coordinates of Lam, its solutions written back as matrices.
  size_t d_lam = Lam.size();
  ModMat Mlam;  // m^2 x d_lam, the vec of the Lam_lambda modulo p
  if (d_lam > 0) {
    Mlam.assign(m * m, std::vector<int64_t>(d_lam, 0));
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t a2 = 0; a2 < m; a2++) {
        for (size_t b2 = 0; b2 < m; b2++) {
          Mlam[a2 * m + b2][l] = ModP(Lam[l][a2][b2], p);
        }
      }
    }
  }
  auto restrict_coef = [&](ModMat const &coef) -> ModMat {
    if (d_lam == 0) {
      return coef;
    }
    ModMat res(coef.size(), std::vector<int64_t>(d_lam, 0));
    for (size_t r2 = 0; r2 < coef.size(); r2++) {
      for (size_t v = 0; v < m * m; v++) {
        if (coef[r2][v] == 0) {
          continue;
        }
        for (size_t l = 0; l < d_lam; l++) {
          res[r2][l] = (res[r2][l] + coef[r2][v] * Mlam[v][l]) % p;
        }
      }
    }
    return res;
  };
  auto to_matrix_vec = [&](std::vector<int64_t> const &dv) -> std::vector<int64_t> {
    std::vector<int64_t> x(m * m, 0);
    for (size_t l = 0; l < d_lam; l++) {
      if (dv[l] == 0) {
        continue;
      }
      for (size_t v = 0; v < m * m; v++) {
        x[v] = (x[v] + dv[l] * Mlam[v][l]) % p;
      }
    }
    return x;
  };
  auto lift_aff = [&](ModMat const &g,
                      int64_t pjl) -> std::optional<LiftAffineSpace> {
    if (d_lam == 0) {
      return LiftAffine(g, A, p, pjl);
    }
    LiftSystem sys = MakeLiftSystem(g, A, p, pjl);
    ModMat cf = restrict_coef(sys.coef);
    ModMat aug;
    for (size_t r2 = 0; r2 < cf.size(); r2++) {
      std::vector<int64_t> row = cf[r2];
      row.push_back(sys.rhs[r2]);
      aug.push_back(row);
    }
    ModRowSpace rs = ModRowReduce(aug, p);
    for (size_t piv : rs.pivots) {
      if (piv == d_lam) {
        return {};
      }
    }
    std::vector<int64_t> part(d_lam, 0);
    for (size_t i2 = 0; i2 < rs.pivots.size(); i2++) {
      part[rs.pivots[i2]] = rs.basis[i2][d_lam];
    }
    std::vector<bool> is_piv(d_lam, false);
    for (size_t pv : rs.pivots) {
      is_piv[pv] = true;
    }
    LiftAffineSpace res;
    res.part = to_matrix_vec(part);
    for (size_t f = 0; f < d_lam; f++) {
      if (is_piv[f]) {
        continue;
      }
      std::vector<int64_t> v(d_lam, 0);
      v[f] = 1;
      for (size_t i2 = 0; i2 < rs.pivots.size(); i2++) {
        v[rs.pivots[i2]] = ModP(-rs.basis[i2][f], p);
      }
      res.ker.push_back(to_matrix_vec(v));
    }
    return res;
  };
  // the element of Lam congruent modulo p^k to g, which is in Lam mod p^k:
  // the coordinates from d_lam rows of Mlam independent modulo p, solved
  // modulo p^k, then the combination of the Lam_lambda modulo p^W
  std::vector<size_t> lam_rows;
  if (d_lam > 0) {
    ModMat test;
    for (size_t v = 0; v < m * m && lam_rows.size() < d_lam; v++) {
      test.push_back(Mlam[v]);
      if (ModRank(test, p) == test.size()) {
        lam_rows.push_back(v);
      } else {
        test.pop_back();
      }
    }
    if (lam_rows.size() != d_lam) {
      return {};
    }
  }
  int64_t pWfull = pow_p(W);
  auto project = [&](ModMat const &g, int k) -> ModMat {
    if (d_lam == 0) {
      return g;
    }
    int64_t M = pow_p(k);
    // S: the d x d matrix of the chosen rows, its inverse modulo p^k
    ModMat S(d_lam, std::vector<int64_t>(d_lam));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      size_t v = lam_rows[r2];
      for (size_t l = 0; l < d_lam; l++) {
        S[r2][l] = ModP(Lam[l][v / m][v % m], M);
      }
    }
    // the inverse modulo p, then Newton
    ModMat aug(d_lam, std::vector<int64_t>(2 * d_lam, 0));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      for (size_t l = 0; l < d_lam; l++) {
        aug[r2][l] = ModP(S[r2][l], p);
      }
      aug[r2][d_lam + r2] = 1;
    }
    ModRowSpace rr = ModRowReduce(aug, p);
    ModMat Sinv(d_lam, std::vector<int64_t>(d_lam));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      for (size_t l = 0; l < d_lam; l++) {
        Sinv[r2][l] = rr.basis[r2][d_lam + l];
      }
    }
    int prec2 = 1;
    while (prec2 < k) {
      prec2 = std::min(2 * prec2, k);
      int64_t M2 = pow_p(prec2);
      ModMat SX(d_lam, std::vector<int64_t>(d_lam, 0));
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t k2 = 0; k2 < d_lam; k2++) {
          for (size_t b2 = 0; b2 < d_lam; b2++) {
            SX[a2][b2] = (SX[a2][b2] + mulmod(ModP(S[a2][k2], M2), Sinv[k2][b2], M2)) % M2;
          }
        }
      }
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t b2 = 0; b2 < d_lam; b2++) {
          SX[a2][b2] = ModP(((a2 == b2) ? 2 : 0) - SX[a2][b2], M2);
        }
      }
      ModMat NS(d_lam, std::vector<int64_t>(d_lam, 0));
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t k2 = 0; k2 < d_lam; k2++) {
          for (size_t b2 = 0; b2 < d_lam; b2++) {
            NS[a2][b2] = (NS[a2][b2] + mulmod(Sinv[a2][k2], SX[k2][b2], M2)) % M2;
          }
        }
      }
      Sinv = NS;
    }
    std::vector<int64_t> coefs(d_lam, 0);
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t r2 = 0; r2 < d_lam; r2++) {
        size_t v = lam_rows[r2];
        coefs[l] = (coefs[l] + mulmod(Sinv[l][r2], ModP(g[v / m][v % m], M), M)) % M;
      }
    }
    ModMat h(m, std::vector<int64_t>(m, 0));
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t a2 = 0; a2 < m; a2++) {
        for (size_t b2 = 0; b2 < m; b2++) {
          h[a2][b2] = (h[a2][b2] + mulmod(coefs[l], ModP(Lam[l][a2][b2], pWfull), pWfull)) % pWfull;
        }
      }
    }
    return h;
  };
  // a solution modulo p^target extending the solution g modulo p^j, the
  // first found by a depth-first search over the affine spaces
  std::function<std::optional<ModMat>(ModMat const &, int, int64_t, int)>
      extend = [&](ModMat const &g, int j, int64_t pjl,
                   int target) -> std::optional<ModMat> {
    if (j >= target) {
      return g;
    }
    if (over || ++n_node > bounds.max_node) {
      over = true;
      return {};
    }
    // with a commutant, g in Lam modulo p^{j+1} (it is modulo p^j), so that
    // the lifts g + p^j Delta, Delta in Lam, are in Lam modulo p^{j+1}
    ModMat gl = (d_lam > 0) ? reduce(project(g, j), pjl * p) : g;
    std::optional<LiftAffineSpace> sp = lift_aff(gl, pjl);
    if (!sp) {
      return {};
    }
    size_t dim = sp->ker.size();
    std::vector<int64_t> coef(dim, 0);
    while (!over) {
      std::optional<ModMat> r = extend(LiftAffineElement(gl, *sp, coef, p, pjl),
                                       j + 1, pjl * p, target);
      if (r) {
        return r;
      }
      size_t i = 0;
      while (i < dim && coef[i] == p - 1) {
        coef[i] = 0;
        i++;
      }
      if (i == dim) {
        return {};
      }
      coef[i]++;
    }
    return {};
  };
  // the element of G (modulo p^c) congruent to the solution g modulo p^j,
  // when there is one: a solution modulo p^k, k = max(j + e + t, e + 2t +
  // 1), extending g, then Newton
  auto lift_path = [&](ModMat const &g, int j,
                       int64_t pjl) -> std::optional<ModMat> {
    int target = std::max(j + e + t, e + 2 * t + 1);
    std::optional<ModMat> sol = extend(g, j, pjl, target);
    if (!sol) {
      return {};
    }
    std::optional<ModMat> h = newton(project(*sol, target));
    if (!h) {
      std::cerr << "MULT_FACTOR: Newton did not converge\n";
      throw TerminalException{1};
    }
    return h;
  };
  // the affine space {x : (fixed_k^T A_p) x = rhs_k}
  auto affine = [&](std::vector<std::vector<int64_t>> const &fixed,
                    std::vector<int64_t> const &rhs)
      -> std::optional<std::pair<std::vector<int64_t>, ModMat>> {
    if (fixed.empty()) {
      ModMat ker;
      for (size_t v = 0; v < m; v++) {
        std::vector<int64_t> e1(m, 0);
        e1[v] = 1;
        ker.push_back(e1);
      }
      return std::make_pair(std::vector<int64_t>(m, 0), ker);
    }
    ModMat aug;
    for (size_t k = 0; k < fixed.size(); k++) {
      std::vector<int64_t> row(m + 1, 0);
      for (size_t b = 0; b < m; b++) {
        int64_t sacc = 0;
        for (size_t a = 0; a < m; a++) {
          sacc = (sacc + fixed[k][a] * Ap[a][b]) % p;
        }
        row[b] = sacc;
      }
      row[m] = ModP(rhs[k], p);
      aug.push_back(row);
    }
    ModRowSpace rs = ModRowReduce(aug, p);
    for (size_t piv : rs.pivots) {
      if (piv == m) {
        return {};
      }
    }
    std::vector<int64_t> part(m, 0);
    for (size_t i = 0; i < rs.pivots.size(); i++) {
      part[rs.pivots[i]] = rs.basis[i][m];
    }
    std::vector<bool> is_piv(m, false);
    for (size_t pv : rs.pivots) {
      is_piv[pv] = true;
    }
    ModMat ker;
    for (size_t f = 0; f < m; f++) {
      if (is_piv[f]) {
        continue;
      }
      std::vector<int64_t> v(m, 0);
      v[f] = 1;
      for (size_t i = 0; i < rs.pivots.size(); i++) {
        v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
      }
      ker.push_back(v);
    }
    return std::make_pair(part, ker);
  };
  // f(x) for the x of the affine space with B(x, x) = q, until f is true
  auto for_each_vector =
      [&](std::pair<std::vector<int64_t>, ModMat> const &sp, int64_t q,
          std::function<bool(std::vector<int64_t> const &)> const &f) -> void {
    size_t dim = sp.second.size();
    std::vector<int64_t> coef(dim, 0);
    while (!over) {
      std::vector<int64_t> x = sp.first;
      for (size_t i = 0; i < dim; i++) {
        if (coef[i] != 0) {
          for (size_t v = 0; v < m; v++) {
            x[v] = (x[v] + coef[i] * sp.second[i][v]) % p;
          }
        }
      }
      if (bil(x, x) == q && f(x)) {
        return;
      }
      size_t i = 0;
      while (i < dim && coef[i] == p - 1) {
        coef[i] = 0;
        i++;
      }
      if (i == dim) {
        return;
      }
      coef[i]++;
    }
  };
  // an element of G_1 with the given first columns, lifted: the further
  // columns searched one at a time
  std::vector<std::vector<int64_t>> cols(m);
  // e = 1, basis vectors rad_start, ... spanning the radical R of A mod p:
  // the image of such a vector lies in R, and the form A / p on R is kept
  // modulo p (representatives in [0, p), which are in X cap p X^#)
  int64_t p2 = p * p;
  auto radical_ok = [&](std::vector<int64_t> const &x, size_t i) -> bool {
    if (i < rad_start) {
      return true;
    }
    for (size_t a = 0; a < rad_start; a++) {
      if (x[a] != 0) {
        return false;
      }
    }
    auto pair_val = [&](std::vector<int64_t> const &u,
                        std::vector<int64_t> const &w) -> int64_t {
      __int128 sacc = 0;
      for (size_t a = 0; a < m; a++) {
        if (u[a] == 0) {
          continue;
        }
        for (size_t b = 0; b < m; b++) {
          sacc += static_cast<__int128>(u[a]) * A[a][b] * w[b];
        }
      }
      return ModP(static_cast<int64_t>(sacc % p2), p2);
    };
    if (pair_val(x, x) != ModP(A[i][i], p2)) {
      return false;
    }
    for (size_t k = rad_start; k < i; k++) {
      if (pair_val(x, cols[k]) != ModP(A[k][i], p2)) {
        return false;
      }
    }
    return true;
  };
  // With an action Om of O_K (a centre K): the image of e_j is forced when
  // e_j is in the span modulo p of the w e_k, k < j, w in Om; otherwise the
  // conditions B(w y_k, y) = B(w e_k, e_j) for all k < j and w, and B(w y,
  // y) = B(w e_j, e_j) for all w.
  size_t n_om = Om.size();
  std::vector<std::vector<std::vector<int64_t>>> Omc(n_om);
  for (size_t w = 0; w < n_om; w++) {
    Omc[w].assign(m, std::vector<int64_t>(m));
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        Omc[w][a2][b2] = ModP(Om[w][a2][b2], p);
      }
    }
  }
  auto om_apply = [&](size_t w, std::vector<int64_t> const &x) -> std::vector<int64_t> {
    std::vector<int64_t> y(m, 0);
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        y[a2] = (y[a2] + Omc[w][a2][b2] * x[b2]) % p;
      }
    }
    return y;
  };
  // forced[j]: the coefficients (k, w, coefficient) of e_j
  std::vector<bool> forced(m, false);
  std::vector<std::vector<std::tuple<size_t, size_t, int64_t>>> forced_combo(m);
  if (n_om > 0) {
    std::vector<std::pair<size_t, size_t>> labels;
    ModMat span;
    for (size_t j = 0; j < m; j++) {
      std::vector<int64_t> ej(m, 0);
      ej[j] = 1;
      if (!span.empty()) {
        // solve sum c_r span_r = e_j
        ModMat aug(m, std::vector<int64_t>(span.size() + 1, 0));
        for (size_t r2 = 0; r2 < span.size(); r2++) {
          for (size_t a2 = 0; a2 < m; a2++) {
            aug[a2][r2] = span[r2][a2];
          }
        }
        for (size_t a2 = 0; a2 < m; a2++) {
          aug[a2][span.size()] = ej[a2];
        }
        ModRowSpace rs = ModRowReduce(aug, p);
        bool solvable = true;
        for (size_t pv : rs.pivots) {
          if (pv == span.size()) {
            solvable = false;
          }
        }
        if (solvable) {
          forced[j] = true;
          for (size_t i2 = 0; i2 < rs.pivots.size(); i2++) {
            int64_t cf = rs.basis[i2][span.size()];
            if (cf != 0) {
              auto lb = labels[rs.pivots[i2]];
              forced_combo[j].push_back({lb.first, lb.second, cf});
            }
          }
        }
      }
      // e_j itself (label n_om), and the w e_j
      span.push_back(ej);
      labels.push_back({j, n_om});
      for (size_t w = 0; w < n_om; w++) {
        span.push_back(om_apply(w, ej));
        labels.push_back({j, w});
      }
    }
  }
  // B(w e_k, e_j) modulo p
  auto om_entry = [&](size_t w, size_t k, size_t j) -> int64_t {
    int64_t sacc = 0;
    for (size_t a2 = 0; a2 < m; a2++) {
      sacc = (sacc + Omc[w][a2][k] * Ap[a2][j]) % p;
    }
    return sacc;
  };
  // the conditions on the image y of e_j beyond B(y, y) and B(y_k, y)
  auto om_ok = [&](std::vector<int64_t> const &y, size_t j) -> bool {
    for (size_t w = 0; w < n_om; w++) {
      if (bil(om_apply(w, y), y) != om_entry(w, j, j)) {
        return false;
      }
    }
    return true;
  };
  // the full check of a candidate modulo p: isometry, invertible, and
  // commuting with Om
  auto full_ok = [&](ModMat const &g) -> bool {
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        std::vector<int64_t> ca(m), cb(m);
        for (size_t u = 0; u < m; u++) {
          ca[u] = g[u][a2];
          cb[u] = g[u][b2];
        }
        if (bil(ca, cb) != Ap[a2][b2]) {
          return false;
        }
      }
    }
    for (size_t w = 0; w < n_om; w++) {
      for (size_t a2 = 0; a2 < m; a2++) {
        for (size_t b2 = 0; b2 < m; b2++) {
          int64_t x1 = 0, x2 = 0;
          for (size_t u = 0; u < m; u++) {
            x1 = (x1 + g[a2][u] * Omc[w][u][b2]) % p;
            x2 = (x2 + Omc[w][a2][u] * g[u][b2]) % p;
          }
          if (x1 != x2) {
            return false;
          }
        }
      }
    }
    return ModRank(g, p) == m;
  };
  // the affine conditions on the image of e_j given the columns k < j
  auto conditions = [&](size_t j, std::vector<std::vector<int64_t>> &fixed,
                        std::vector<int64_t> &rhs) -> void {
    fixed.clear();
    rhs.clear();
    for (size_t k = 0; k < j; k++) {
      fixed.push_back(cols[k]);
      rhs.push_back(Ap[k][j]);
      for (size_t w = 0; w < n_om; w++) {
        fixed.push_back(om_apply(w, cols[k]));
        rhs.push_back(om_entry(w, k, j));
      }
    }
  };
  std::function<std::optional<ModMat>(size_t)> complete =
      [&](size_t j) -> std::optional<ModMat> {
    if (j == m) {
      ModMat g(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          g[a][b] = cols[b][a];
        }
      }
      if (!full_ok(g)) {
        return {};
      }
      return lift_path(g, 1, p);
    }
    if (forced[j]) {
      std::vector<int64_t> y(m, 0);
      for (auto const &[k, w, cf] : forced_combo[j]) {
        std::vector<int64_t> v = (w == n_om) ? cols[k] : om_apply(w, cols[k]);
        for (size_t a2 = 0; a2 < m; a2++) {
          y[a2] = (y[a2] + cf * v[a2]) % p;
        }
      }
      if (!radical_ok(y, j)) {
        return {};
      }
      cols[j] = y;
      return complete(j + 1);
    }
    std::vector<std::vector<int64_t>> fixed;
    std::vector<int64_t> rhs;
    conditions(j, fixed, rhs);
    std::optional<std::pair<std::vector<int64_t>, ModMat>> sp =
        affine(fixed, rhs);
    if (!sp) {
      return {};
    }
    std::optional<ModMat> res;
    for_each_vector(*sp, Ap[j][j], [&](std::vector<int64_t> const &x) -> bool {
      if (++n_node > bounds.max_node) {
        over = true;
        return true;
      }
      if (!radical_ok(x, j) || !om_ok(x, j)) {
        return false;
      }
      cols[j] = x;
      res = complete(j + 1);
      return res.has_value();
    });
    return res;
  };
  auto encode = [&](std::vector<int64_t> const &x) -> size_t {
    size_t idx = 0;
    for (size_t a = m; a-- > 0;) {
      idx = idx * static_cast<size_t>(p) + static_cast<size_t>(x[a]);
    }
    return idx;
  };
  auto apply = [&](ModMat const &g,
                   std::vector<int64_t> const &x) -> std::vector<int64_t> {
    std::vector<int64_t> y(m, 0);
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        y[a] = (y[a] + ModP(g[a][b], p) * x[b]) % p;
      }
    }
    return y;
  };
  std::vector<ModMat> gens1;
  std::vector<ModMat> gens;
  double log_order = 0;
  for (size_t i = m; i-- > 0;) {
    std::vector<int64_t> ei(m, 0);
    ei[i] = 1;
    // the orbit of e_i under gens1
    std::vector<char> in_orbit(n_pts, 0);
    auto close = [&]() -> size_t {
      std::vector<std::vector<int64_t>> queue;
      size_t cnt = 0;
      for (size_t idx = 0; idx < n_pts; idx++) {
        if (in_orbit[idx]) {
          cnt++;
        }
      }
      if (cnt == 0) {
        in_orbit[encode(ei)] = 1;
        cnt = 1;
      }
      // from every point of the orbit
      for (size_t idx = 0; idx < n_pts; idx++) {
        if (in_orbit[idx]) {
          std::vector<int64_t> x(m);
          size_t tt = idx;
          for (size_t a = 0; a < m; a++) {
            x[a] = static_cast<int64_t>(tt % static_cast<size_t>(p));
            tt /= static_cast<size_t>(p);
          }
          queue.push_back(x);
        }
      }
      while (!queue.empty()) {
        std::vector<int64_t> x = queue.back();
        queue.pop_back();
        for (auto const &g : gens1) {
          std::vector<int64_t> y = apply(g, x);
          size_t iy = encode(y);
          if (!in_orbit[iy]) {
            in_orbit[iy] = 1;
            cnt++;
            queue.push_back(y);
          }
        }
      }
      return cnt;
    };
    size_t orbit_len = close();
    if (forced[i]) {
      // fixing e_0, ..., e_{i-1}, the K-linear maps fix e_i
      continue;
    }
    std::vector<std::vector<int64_t>> idcols(i);
    for (size_t k = 0; k < i; k++) {
      std::vector<int64_t> ek(m, 0);
      ek[k] = 1;
      idcols[k] = ek;
      cols[k] = ek;
    }
    std::vector<std::vector<int64_t>> fixed;
    std::vector<int64_t> rhs;
    conditions(i, fixed, rhs);
    std::optional<std::pair<std::vector<int64_t>, ModMat>> sp =
        affine(fixed, rhs);
    if (!sp) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: no admissible image of e_i, n_node=" << n_node << "\n";
#endif
      return {};
    }
    size_t n_cand = 0, n_fail = 0;
    size_t node0 = n_node;
    for_each_vector(*sp, Ap[i][i], [&](std::vector<int64_t> const &y) -> bool {
      if (!radical_ok(y, i) || !om_ok(y, i)) {
        return false;
      }
      n_cand++;
      if (in_orbit[encode(y)]) {
        return false;
      }
      for (size_t k = 0; k < i; k++) {
        cols[k] = idcols[k];
      }
      cols[i] = y;
      std::optional<ModMat> g = complete(i + 1);
      if (over) {
        return true;
      }
      if (g) {
        ModMat g1 = *g;
        for (auto &row : g1) {
          for (auto &x : row) {
            x = ModP(x, p);
          }
        }
        gens1.push_back(g1);
        gens.push_back(*g);
        orbit_len = close();
      } else {
        n_fail++;
      }
      return false;
    });
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: level " << i << " candidates " << n_cand
              << " failed " << n_fail << " orbit " << orbit_len << " nodes "
              << n_node - node0 << "\n";
#endif
    if (over) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: budget at level 1, n_node=" << n_node << "\n";
#endif
      return {};
    }
    log_order += std::log(static_cast<double>(orbit_len));
  }
  // the layers
  ModMat Id(m, std::vector<int64_t>(m, 0));
  for (size_t a = 0; a < m; a++) {
    Id[a][a] = 1;
  }
  int64_t pj = p;
  for (int j = 1; j < c; j++) {
    std::optional<LiftAffineSpace> sp = lift_aff(Id, pj);
    if (!sp) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: no layer, n_node=" << n_node << "\n";
#endif
      return {};
    }
    std::vector<ModMat> lifts;
    bool all = true;
    for (auto const &kv : sp->ker) {
      std::vector<int64_t> coef(sp->ker.size(), 0);
      for (size_t i = 0; i < sp->ker.size(); i++) {
        coef[i] = (&kv == &sp->ker[i]) ? 1 : 0;
      }
      ModMat g = LiftAffineElement(Id, *sp, coef, p, pj);
      std::optional<ModMat> lf = lift_path(g, j + 1, pj * p);
      if (over) {
        return {};
      }
      if (!lf) {
        all = false;
        break;
      }
      lifts.push_back(*lf);
    }
    // p odd and e = 1: one level decides the lift, by a system whose
    // matrix is constant on the layer; its solvability is f(Delta) = L rhs
    // = 0, L a basis of the left kernel, f = l + q with l linear and q
    // quadratic. The liftable Delta form a subspace V (an image of a
    // group); f(x Delta) = x l + x^2 q for all x in F_p gives V in W0 = ker
    // l, and V = W0 when f vanishes on a basis of W0 and on the sums of its
    // pairs (q then vanishes on W0).
    bool done = false;
    if (!all && p != 2 && e == 1) {
      size_t dim = sp->ker.size();
      int64_t pj1 = pj * p;
      auto layer_elt = [&](std::vector<int64_t> const &coefv) -> ModMat {
        return LiftAffineElement(Id, *sp, coefv, p, pj);
      };
      LiftSystem sys0 = MakeLiftSystem(layer_elt(std::vector<int64_t>(dim, 0)),
                                       A, p, pj1);
      ModMat coef0 = restrict_coef(sys0.coef);
      size_t n_eq = coef0.size();
      size_t n_var = coef0.empty() ? 0 : coef0[0].size();
      // the left kernel: the y with y^T coef = 0
      ModMat coefT(n_var, std::vector<int64_t>(n_eq));
      for (size_t a = 0; a < n_eq; a++) {
        for (size_t b = 0; b < n_var; b++) {
          coefT[b][a] = coef0[a][b];
        }
      }
      auto nullspace = [&](ModMat const &M, size_t ncol) -> ModMat {
        ModRowSpace rs = ModRowReduce(M, p);
        std::vector<bool> is_piv(ncol, false);
        for (size_t pv : rs.pivots) {
          is_piv[pv] = true;
        }
        ModMat ker;
        for (size_t f = 0; f < ncol; f++) {
          if (is_piv[f]) {
            continue;
          }
          std::vector<int64_t> v(ncol, 0);
          v[f] = 1;
          for (size_t i = 0; i < rs.pivots.size(); i++) {
            v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
          }
          ker.push_back(v);
        }
        return ker;
      };
      ModMat L = nullspace(coefT, n_eq);
      auto fval = [&](std::vector<int64_t> const &coefv) -> std::vector<int64_t> {
        LiftSystem sys = MakeLiftSystem(layer_elt(coefv), A, p, pj1);
        std::vector<int64_t> r(L.size(), 0);
        for (size_t k = 0; k < L.size(); k++) {
          for (size_t a = 0; a < n_eq; a++) {
            r[k] = (r[k] + L[k][a] * sys.rhs[a]) % p;
          }
        }
        return r;
      };
      int64_t inv2 = (p + 1) / 2;
      // l on the basis of the layer: l(b) = f(b) - q(b), q(b) = (f(2b) -
      // 2 f(b)) / 2
      ModMat lmat(L.size(), std::vector<int64_t>(dim, 0));
      for (size_t i = 0; i < dim; i++) {
        std::vector<int64_t> c1(dim, 0), c2(dim, 0);
        c1[i] = 1;
        c2[i] = 2 % p;
        std::vector<int64_t> f1 = fval(c1), f2 = fval(c2);
        for (size_t k = 0; k < L.size(); k++) {
          int64_t q = ModP(f2[k] - 2 * f1[k], p) * inv2 % p;
          lmat[k][i] = ModP(f1[k] - q, p);
        }
      }
      ModMat W0 = nullspace(lmat, dim);
      bool q_zero = true;
      for (size_t a = 0; a < W0.size() && q_zero; a++) {
        for (size_t b = a; b < W0.size() && q_zero; b++) {
          std::vector<int64_t> cv = W0[a];
          if (b != a) {
            for (size_t i = 0; i < dim; i++) {
              cv[i] = (cv[i] + W0[b][i]) % p;
            }
          }
          for (int64_t x : fval(cv)) {
            if (x != 0) {
              q_zero = false;
            }
          }
        }
      }
      if (q_zero) {
        for (auto const &w : W0) {
          std::optional<ModMat> lf = lift_path(layer_elt(w), j + 1, pj1);
          if (over) {
            return {};
          }
          if (!lf) {
            std::cerr << "MULT_FACTOR: an element of the layer does not lift\n";
            throw TerminalException{1};
          }
          gens.push_back(*lf);
        }
        log_order += static_cast<double>(W0.size()) *
                     std::log(static_cast<double>(p));
        done = true;
      }
    }
    if (done) {
      // found above
    } else if (all) {
      gens.insert(gens.end(), lifts.begin(), lifts.end());
      log_order += static_cast<double>(sp->ker.size()) *
                   std::log(static_cast<double>(p));
    } else {
      // the liftable subspace, by enumerating the space
      size_t dim = sp->ker.size();
      size_t n_el = 1;
      for (size_t i = 0; i < dim; i++) {
        n_el *= static_cast<size_t>(p);
        if (n_el > 100000) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: layer too large, n_node=" << n_node << "\n";
#endif
          return {};
        }
      }
      ModMat liftable_coefs;
      std::vector<ModMat> liftable_elts;
      std::vector<int64_t> coef(dim, 0);
      for (size_t s2 = 0; s2 < n_el; s2++) {
        size_t tt = s2;
        for (size_t i = 0; i < dim; i++) {
          coef[i] = static_cast<int64_t>(tt % static_cast<size_t>(p));
          tt /= static_cast<size_t>(p);
        }
        ModMat g = LiftAffineElement(Id, *sp, coef, p, pj);
        std::optional<ModMat> lf = lift_path(g, j + 1, pj * p);
        if (over) {
          return {};
        }
        if (lf) {
          ModMat test = liftable_coefs;
          test.push_back(coef);
          if (ModRank(test, p) == test.size()) {
            liftable_coefs = test;
            liftable_elts.push_back(*lf);
          }
        }
      }
      gens.insert(gens.end(), liftable_elts.begin(), liftable_elts.end());
      log_order += static_cast<double>(liftable_coefs.size()) *
                   std::log(static_cast<double>(p));
    }
    pj *= p;
  }
  if (order) {
    *order = std::exp(log_order);
  }
  return gens;
}

/*
  LocalOrthogonalImageGeneratorsCore in a basis adapted to the radical R of
  A mod p when e = 1: the complement first, then lifts in [0, p) of a basis
  of R, P of determinant prime to p; the search runs on A' = P^T A P with
  the pruning of the radical, and the generators g' come back as P g'
  P^{-1}. For e != 1 the standard basis, without that pruning. Om (the
  action of a Z-basis of O_K, a centre K) and Lam (a Z-basis of the
  integral commutant of O_K) are empty when End_G = Q; all the matrices
  act on the columns of coordinates.
 */
inline std::optional<std::vector<ModMat>>
LocalOrthogonalImageGenerators(ModMat const &A, ModMat const &Be,
                               std::vector<ModMat> const &Om,
                               std::vector<ModMat> const &Lam, int64_t p,
                               int c, int e,
                               MultiplicityFactorBounds const &bounds,
                               double *order) {
  size_t m = A.size();
  if (e != 1) {
    return LocalOrthogonalImageGeneratorsCore(A, Be, Om, Lam, p, c, e, m,
                                              bounds, order);
  }
  // the radical of A mod p
  ModMat Ap(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      Ap[a][b] = ModP(A[a][b], p);
    }
  }
  ModRowSpace rs = ModRowReduce(Ap, p);
  std::vector<bool> is_piv(m, false);
  for (size_t pv : rs.pivots) {
    is_piv[pv] = true;
  }
  ModMat rad;
  for (size_t f = 0; f < m; f++) {
    if (is_piv[f]) {
      continue;
    }
    std::vector<int64_t> v(m, 0);
    v[f] = 1;
    for (size_t i = 0; i < rs.pivots.size(); i++) {
      v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
    }
    rad.push_back(v);
  }
  if (rad.empty()) {
    return LocalOrthogonalImageGeneratorsCore(A, Be, Om, Lam, p, c, e, m,
                                              bounds, order);
  }
  // a complement from the standard basis
  ModMat cols_P;
  ModMat test = rad;
  for (size_t a = 0; a < m && cols_P.size() + rad.size() < m; a++) {
    std::vector<int64_t> v(m, 0);
    v[a] = 1;
    test.push_back(v);
    if (ModRank(test, p) == test.size()) {
      cols_P.push_back(v);
    } else {
      test.pop_back();
    }
  }
  size_t rad_start = cols_P.size();
  for (auto const &v : rad) {
    cols_P.push_back(v);
  }
  ModMat P(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      P[a][b] = cols_P[b][a];
    }
  }
  int W = MultiplicityNewtonPrecision(p);
  int64_t pW = 1;
  for (int i = 0; i < W; i++) {
    pW *= p;
  }
  auto mulmod = [](int64_t x, int64_t y, int64_t M) -> int64_t {
    __int128 r = static_cast<__int128>(x) * y % M;
    return static_cast<int64_t>((r < 0) ? r + M : r);
  };
  auto matmul = [&](ModMat const &X, ModMat const &Y, int64_t M) -> ModMat {
    ModMat Z(m, std::vector<int64_t>(m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t k = 0; k < m; k++) {
        if (X[a][k] == 0) {
          continue;
        }
        for (size_t b = 0; b < m; b++) {
          Z[a][b] = (Z[a][b] + mulmod(X[a][k], Y[k][b], M)) % M;
        }
      }
    }
    return Z;
  };
  auto transpose = [&](ModMat const &X) -> ModMat {
    ModMat Z(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Z[a][b] = X[b][a];
      }
    }
    return Z;
  };
  // P^{-1} modulo p^W: modulo p by Gauss, then X <- X (2 - P X)
  ModMat Pinv(m, std::vector<int64_t>(m));
  {
    ModMat aug(m, std::vector<int64_t>(2 * m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        aug[a][b] = ModP(P[a][b], p);
      }
      aug[a][m + a] = 1;
    }
    ModRowSpace rr = ModRowReduce(aug, p);
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Pinv[a][b] = rr.basis[a][m + b];
      }
    }
    int64_t M = p;
    while (M < pW) {
      M = (M > pW / M) ? pW : M * M;
      ModMat PX = matmul(P, Pinv, M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          PX[a][b] = ModP(((a == b) ? 2 : 0) - PX[a][b], M);
        }
      }
      Pinv = matmul(Pinv, PX, M);
    }
  }
  // A' = P^T A P exactly, Be' = P^{-1} Be P^{-T} modulo p^W
  ModMat A2(m, std::vector<int64_t>(m, 0));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      __int128 sacc = 0;
      for (size_t k = 0; k < m; k++) {
        for (size_t l = 0; l < m; l++) {
          sacc += static_cast<__int128>(P[k][a]) * A[k][l] * P[l][b];
        }
      }
      A2[a][b] = static_cast<int64_t>(sacc);
    }
  }
  ModMat Be2 = matmul(matmul(Pinv, Be, pW), transpose(Pinv), pW);
  // the action of O_K and the commutant in the new basis: P^{-1} X P
  auto conj = [&](ModMat const &X) -> ModMat {
    ModMat Xr = X;
    for (auto &row : Xr) {
      for (auto &x : row) {
        x = ModP(x, pW);
      }
    }
    return matmul(matmul(Pinv, Xr, pW), P, pW);
  };
  std::vector<ModMat> Om2, Lam2;
  for (auto const &X : Om) {
    Om2.push_back(conj(X));
  }
  for (auto const &X : Lam) {
    Lam2.push_back(conj(X));
  }
  std::optional<std::vector<ModMat>> gens =
      LocalOrthogonalImageGeneratorsCore(A2, Be2, Om2, Lam2, p, c, e,
                                         rad_start, bounds, order);
  if (!gens) {
    return {};
  }
  int64_t pc = 1;
  for (int i = 0; i < c; i++) {
    pc *= p;
  }
  std::vector<ModMat> res;
  for (auto const &g : *gens) {
    res.push_back(matmul(matmul(P, g, pc), Pinv, pc));
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
  if (!ResiduesFitInt128(A, p, K)) {
    return {};
  }
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
