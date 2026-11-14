// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
#define SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_

// clang-format off
#include "GeneraArithmetic.h"
#include "GroupNeighborModP.h"
#include "ResidueArithmetic.h"
#include "MAT_MatrixModPrime.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_set>
#include <limits>
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
  L_rho. The centralizer C_rho = O(U) commutes with G, so the stabilizer
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
  ModMatrix ker;
};

// the linear system of the lifts of g from p^j to p^{j+1}: coef Delta =
// rhs modulo p, the equations indexed by a <= b
struct LiftSystem {
  ModMatrix coef;
  std::vector<int64_t> rhs;
};

inline LiftSystem MakeLiftSystem(ModMatrix const &g, ModMatrix const &A, int64_t p,
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
  ModMatrix coef;
  std::vector<int64_t> rhs;
  for (size_t a = 0; a < m; a++) {
    for (size_t b = a; b < m; b++) {
      std::vector<int64_t> row(n_var, 0);
      for (size_t c = 0; c < m; c++) {
        int64_t x = ResMod(static_cast<int64_t>(Ag[c][b] % p), p);
        row[c * m + a] = ResMod(row[c * m + a] + x, p);
        int64_t y = ResMod(static_cast<int64_t>(Ag[c][a] % p), p);
        row[c * m + b] = ResMod(row[c * m + b] + y, p);
      }
      coef.push_back(row);
      rhs.push_back(ResMod(-Ered[a][b], p));
    }
  }
  return LiftSystem{coef, rhs};
}
// The right-hand side of the system of MakeLiftSystem alone: -E / p^j
// modulo p for the a <= b, E = g^T A g - A (its matrix depends on g modulo
// p only, so a search that changes g by p^j Delta needs only this)
inline std::vector<int64_t> LiftRhs(ModMatrix const &g, ModMatrix const &A, int64_t p,
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
  std::vector<int64_t> rhs;
  for (size_t a = 0; a < m; a++) {
    for (size_t b = a; b < m; b++) {
      rhs.push_back(ResMod(-Ered[a][b], p));
    }
  }
  return rhs;
}

inline std::optional<LiftAffineSpace> LiftAffine(ModMatrix const &g,
                                                 ModMatrix const &A, int64_t p,
                                                 int64_t pj) {
  size_t m = A.size();
  size_t n_var = m * m;
  LiftSystem sys = MakeLiftSystem(g, A, p, pj);
  ModMatrix const &coef = sys.coef;
  std::vector<int64_t> const &rhs = sys.rhs;
  // the affine solutions: x with coef x = rhs
  size_t n_eq = coef.size();
  ModMatrix aug(n_eq, std::vector<int64_t>(n_var + 1));
  for (size_t e = 0; e < n_eq; e++) {
    for (size_t v = 0; v < n_var; v++) {
      aug[e][v] = coef[e][v];
    }
    aug[e][n_var] = rhs[e];
  }
  RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
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
  ModMatrix ker;
  for (size_t f = 0; f < n_var; f++) {
    if (is_pivot[f]) {
      continue;
    }
    std::vector<int64_t> v(n_var, 0);
    v[f] = 1;
    for (size_t i = 0; i < rs.pivots.size(); i++) {
      v[rs.pivots[i]] = ResMod(-rs.basis[i][f], p);
    }
    ker.push_back(v);
  }
  return LiftAffineSpace{part, ker};
}

// the element g + p^j Delta, Delta = part + sum_i coef_i ker_i
inline ModMatrix LiftAffineElement(ModMatrix const &g, LiftAffineSpace const &sp,
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
  ModMatrix g2 = g;
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      g2[a][b] = ResMod(g[a][b] + pj * x[a * m + b], pj1);
    }
  }
  return g2;
}

inline std::vector<ModMatrix> LiftSolutions(ModMatrix const &g, ModMatrix const &A,
                                         int64_t p, int j, int64_t pj) {
  (void)j;
  std::optional<LiftAffineSpace> sp = LiftAffine(g, A, p, pj);
  if (!sp) {
    return {};
  }
  std::vector<ModMatrix> res;
  size_t n_sol = MyPow<size_t>(static_cast<size_t>(p), sp->ker.size());
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
inline std::optional<std::vector<ModMatrix>>
SolutionsModP(ModMatrix const &A, int64_t p, size_t max_solution,
              size_t max_node) {
  size_t m = A.size();
  std::vector<std::vector<int64_t>> vects;
  size_t n_vect = MyPow<size_t>(static_cast<size_t>(p), m);
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
        s = (s + x[a] * ResMod(A[a][b], p) % p * y[b]) % p;
      }
    }
    return s;
  };
  std::vector<ModMatrix> res;
  std::vector<std::vector<int64_t>> cols(m);
  bool overflow = false;
  size_t n_node = 0;
  std::function<void(size_t)> rec = [&](size_t j) -> void {
    if (overflow) {
      return;
    }

    if (j == m) {
      ModMatrix g(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          g[a][b] = cols[b][a];
        }
      }
      if (RankModPrime(g, p) == m) {
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
        if (bil(w, v) != ResMod(A[i][j], p)) {
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
  The searches of the image of O(X_p) modulo p^c (X with p^e A^{-1}
  integral) only see A and the solutions modulo p^{K+1}, K as below; the
  modulus, or nothing beyond 2^62. Exact data (a Gram matrix, a basis of
  the commutant, Tint) is given to them as residues modulo this modulus,
  so that it never has to fit in 64 bits.
 */
inline std::optional<int64_t> LocalImageModulus(int64_t p, int c, int e) {
  int t = (p == 2) ? 1 : 0;
  int K = std::max(c + e + t, e + 2 * t + 1);
  int64_t M = 1;
  for (int i = 0; i <= K; i++) {
    if (M > (int64_t(1) << 62) / p) {
      return {};
    }
    M *= p;
  }
  return M;
}

/*
  The searches below keep the residues modulo p^K in machine integers and
  form g A g^T exactly in __int128: this holds when m^2 max|A| p^{2K} stays
  below 2^120. Beyond that they return nothing (the factor is then not
  computed), rather than overflow.
 */
inline bool ResiduesFitInt128(ModMatrix const &A, int64_t p, int K) {
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
inline std::optional<std::vector<ModMatrix>>
LocalOrthogonalImage(ModMatrix const &A_in, int64_t p, int c, int e,
                     MultiplicityFactorBounds const &bounds) {
  int t = (p == 2) ? 1 : 0;
  int K = std::max(c + e + t, e + 2 * t + 1);
  // A is only read modulo p^{j+1}, j < c: centered modulo p^{K+1}
  std::optional<int64_t> optM = LocalImageModulus(p, c, e);
  if (!optM) {
    return {};
  }
  ModMatrix A = A_in;
  for (auto &row : A) {
    for (auto &x : row) {
      x = CenteredModP(x, *optM);
    }
  }
  if (!ResiduesFitInt128(A, p, K)) {
    return {};
  }
  std::optional<std::vector<ModMatrix>> opt =
      SolutionsModP(A, p, bounds.max_solution, bounds.max_node);
  if (!opt) {
    return {};
  }
  std::vector<ModMatrix> level = *opt;
  int64_t pj = p;
  for (int j = 1; j < c; j++) {
    std::vector<ModMatrix> next;
    for (auto const &g : level) {
      std::vector<ModMatrix> lifts = LiftSolutions(g, A, p, j, pj);
      next.insert(next.end(), lifts.begin(), lifts.end());
      if (next.size() > bounds.max_solution) {
        return {};
      }
    }
    level = std::move(next);
    pj *= p;
  }
  // keep those that lift to level K
  std::function<bool(ModMatrix const &, int, int64_t)> liftable =
      [&](ModMatrix const &g, int j, int64_t pj_loc) -> bool {
    if (j >= K) {
      return true;
    }
    std::vector<ModMatrix> lifts = LiftSolutions(g, A, p, j, pj_loc);
    for (auto const &g2 : lifts) {
      if (liftable(g2, j + 1, pj_loc * p)) {
        return true;
      }
    }
    return false;
  };
  std::vector<ModMatrix> res;
  for (auto const &g : level) {
    if (liftable(g, c, pj)) {
      res.push_back(g);
    }
  }
  return res;
}

/*
  The Howell normal form of the subgroup of (Z / p^C)^r generated by the
  rows, p^C < 2^62: an echelon form whose pivots are powers p^v of p (in
  increasing columns), with the saturation rows p^{C-v} (pivot row) added
  as they arise, so that the rows with zeros in the first j columns span
  the part of the subgroup with zeros there, and the entries above a pivot
  p^v reduced into [0, p^v). Two generating sets give the same form exactly
  when they generate the same subgroup.
 */
inline ModMatrix HowellFormPrimePower(ModMatrix rows, int64_t p, int C) {
  int64_t N = MyPow<int64_t>(p, C);
  auto mulmod = [&](int64_t x, int64_t y) -> int64_t {
    return ResMul(x, y, N);
  };
  auto inv_mod = [&](int64_t a) -> int64_t { return ResInv(a, N); };
  auto val = [&](int64_t x) -> int {
    if (x == 0) {
      return C;
    }
    int v = 0;
    while (x % p == 0) {
      x /= p;
      v++;
    }
    return v;
  };
  size_t r = rows.empty() ? 0 : rows[0].size();
  for (auto &row : rows) {
    for (auto &x : row) {
      x = ((x % N) + N) % N;
    }
  }
  ModMatrix W = rows;
  ModMatrix R;
  std::vector<size_t> piv_col;
  std::vector<int> piv_val;
  for (size_t j = 0; j < r; j++) {
    // the row of least valuation in column j
    int vmin = C;
    size_t imin = 0;
    for (size_t i = 0; i < W.size(); i++) {
      int v = val(W[i][j]);
      if (v < vmin) {
        vmin = v;
        imin = i;
      }
    }
    if (vmin == C) {
      continue;
    }
    std::vector<int64_t> pr = W[imin];
    W.erase(W.begin() + static_cast<long>(imin));
    int64_t pv = MyPow<int64_t>(p, vmin);
    // pivot p^vmin
    int64_t u = pr[j] / pv;
    int64_t ui = inv_mod(u);
    for (auto &x : pr) {
      x = mulmod(x, ui);
    }
    for (auto &row : W) {
      if (row[j] == 0) {
        continue;
      }
      int64_t q = row[j] / pv;
      for (size_t b = 0; b < r; b++) {
        row[b] = ((row[b] - mulmod(q, pr[b])) % N + N) % N;
      }
    }
    // the saturation row p^{C - vmin} pr, zero in column j
    int64_t sat = N / pv;
    std::vector<int64_t> srow(r);
    bool nonzero = false;
    for (size_t b = 0; b < r; b++) {
      srow[b] = mulmod(sat, pr[b]);
      if (srow[b] != 0) {
        nonzero = true;
      }
    }
    if (nonzero) {
      W.push_back(srow);
    }
    R.push_back(pr);
    piv_col.push_back(j);
    piv_val.push_back(vmin);
  }
  // the entries above the pivots in [0, p^v)
  for (size_t i = 0; i < R.size(); i++) {
    size_t j = piv_col[i];
    int64_t pv = MyPow<int64_t>(p, piv_val[i]);
    for (size_t k = 0; k < i; k++) {
      int64_t q = R[k][j] / pv;
      if (q != 0) {
        for (size_t b = 0; b < r; b++) {
          R[k][b] = ((R[k][b] - mulmod(q, R[i][b])) % N + N) % N;
        }
      }
    }
  }
  return R;
}

// The order of an image found by LocalOrthogonalImageGeneratorsCore: the
// product of the orbit lengths of the stabilizer chain modulo p times
// p^layer_exp, the layers G_j / G_{j+1}
struct LocalImageOrder {
  std::vector<size_t> orbits;
  size_t layer_exp = 0;
};

inline double LocalImageOrderValue(LocalImageOrder const &o, int64_t p) {
  double v = 1;
  for (size_t x : o.orbits) {
    v *= static_cast<double>(x);
  }
  return v * std::pow(static_cast<double>(p), static_cast<double>(o.layer_exp));
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
  A Jordan basis of the symmetric A modulo p^W, by congruence: the columns
  of P, with P^T A P block diagonal, blocks p^l U with U unimodular, of
  size 1 (and 2 at p = 2 when an off-diagonal entry has the least
  valuation), and the level l of each column, nondecreasing. In that basis
  the lattice X^(i) = {x : A x = 0 mod p^i} is spanned modulo p by the
  columns of level >= i, which lie in it exactly. The precision decreases
  by the levels met, far above what the search uses.
 */
inline std::pair<ModMatrix, std::vector<int>>
JordanBasisModPW(ModMatrix const &A, int64_t p, int W) {
  size_t m = A.size();
  int64_t pW = MyPow<int64_t>(p, W);
  auto mulmod = [&](int64_t x, int64_t y) -> int64_t {
    return ResMul(x, y, pW);
  };
  auto inv_mod = [&](int64_t a) -> int64_t { return ResInv(a, pW); };
  auto val = [&](int64_t x) -> int {
    x = ResMod(x, pW);
    if (x == 0) {
      return W;
    }
    int v = 0;
    while (x % p == 0) {
      x /= p;
      v++;
    }
    return v;
  };
  ModMatrix M(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      M[a][b] = ResMod(A[a][b], pW);
    }
  }
  ModMatrix P(m, std::vector<int64_t>(m, 0));
  for (size_t a = 0; a < m; a++) {
    P[a][a] = 1;
  }
  // column b += f column a, and the same on the rows (a congruence)
  auto add = [&](size_t b, size_t a, int64_t f) -> void {
    if (f == 0) {
      return;
    }
    for (size_t u = 0; u < m; u++) {
      M[u][b] = ResMod(M[u][b] + mulmod(f, M[u][a]), pW);
    }
    for (size_t u = 0; u < m; u++) {
      M[b][u] = ResMod(M[b][u] + mulmod(f, M[a][u]), pW);
    }
    for (size_t u = 0; u < m; u++) {
      P[u][b] = ResMod(P[u][b] + mulmod(f, P[u][a]), pW);
    }
  };
  auto swap_idx = [&](size_t a, size_t b) -> void {
    if (a == b) {
      return;
    }
    std::swap(M[a], M[b]);
    for (size_t u = 0; u < m; u++) {
      std::swap(M[u][a], M[u][b]);
      std::swap(P[u][a], P[u][b]);
    }
  };
  std::vector<int> levels(m, 0);
  size_t start = 0;
  while (start < m) {
    int vmin = W + 1;
    size_t ia = start, ib = start;
    for (size_t a = start; a < m; a++) {
      for (size_t b = start; b < m; b++) {
        int v = val(M[a][b]);
        if (v < vmin) {
          vmin = v;
          ia = a;
          ib = b;
        }
      }
    }
    if (vmin >= W) {
      // the rest is zero to the precision: level W
      for (size_t a = start; a < m; a++) {
        levels[a] = W;
      }
      break;
    }
    size_t diag = m;
    for (size_t a = start; a < m && diag == m; a++) {
      if (val(M[a][a]) == vmin) {
        diag = a;
      }
    }
#ifdef DEBUG_JORDAN
    std::cerr << "JORDAN: start " << start << " vmin " << vmin << " ia " << ia
              << " ib " << ib << " diag " << diag << "\n";
#endif
    if (diag == m && p != 2) {
      // x_a + x_b: the diagonal entry M_aa + 2 M_ab + M_bb has valuation vmin
      add(ia, ib, 1);
      diag = ia;
    }
    int64_t pv = MyPow<int64_t>(p, vmin);
    if (diag != m) {
      swap_idx(diag, start);
      int64_t u0 = M[start][start] / pv;
      int64_t iu = inv_mod(u0);
      for (size_t b = start + 1; b < m; b++) {
        if (M[start][b] == 0) {
          continue;
        }
        int64_t q = M[start][b] / pv;
        add(b, start, ResMod(-mulmod(q, iu), pW));
      }
      levels[start] = vmin;
      start++;
      continue;
    }
    // p = 2, a 2 x 2 block on (ia, ib): [[x, y], [y, z]] / 2^vmin is
    // [[0, 1], [1, 0]] modulo 2
    if (ia == ib) {
      std::cerr << "JORDAN: no pivot found\n";
      throw TerminalException{1};
    }
    swap_idx(ia, start);
    if (ib == start) {
      ib = ia;
    }
    swap_idx(ib, start + 1);
    int64_t x = M[start][start] / pv, y = M[start][start + 1] / pv,
            z = M[start + 1][start + 1] / pv;
    int64_t det = ResMod(mulmod(x, z) - mulmod(y, y), pW);
    int64_t idet = inv_mod(det);
    for (size_t c = start + 2; c < m; c++) {
      int64_t r1 = M[start][c] / pv, r2 = M[start + 1][c] / pv;
      if (r1 == 0 && r2 == 0) {
        continue;
      }
      // (f1, f2) = -B^{-1} (r1, r2), B^{-1} = [[z, -y], [-y, x]] / det
      int64_t f1 = ResMod(-mulmod(idet, ResMod(mulmod(z, r1) - mulmod(y, r2), pW)), pW);
      int64_t f2 = ResMod(-mulmod(idet, ResMod(mulmod(x, r2) - mulmod(y, r1), pW)), pW);
      add(c, start, f1);
      add(c, start + 1, f2);
    }
    levels[start] = vmin;
    levels[start + 1] = vmin;
    start += 2;
  }
  // P^T A P is the matrix reduced
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      __int128 sacc = 0;
      for (size_t k = 0; k < m; k++) {
        for (size_t l = 0; l < m; l++) {
          sacc = (sacc + static_cast<__int128>(P[k][a]) * ResMod(A[k][l], pW) % pW * P[l][b]) % pW;
        }
      }
      if (ResMod(static_cast<int64_t>(sacc), pW) != M[a][b]) {
        std::cerr << "JORDAN: P^T A P differs from the reduced matrix at " << a
                  << "," << b << "\n";
        throw TerminalException{1};
      }
    }
  }
  return {P, levels};
}

/*
  The integral Lie lattice of the local centralizer, on the columns:

      u = {Y in Lambda : Y^T A + A Y = 0},

  Lambda the Z-span of comm (the integral commutant, saturated in M_m(Z)),
  or all of M_m(Z) when comm is empty. u is saturated in M_m(Z) (the
  condition is linear over Q and Lambda is saturated), so a Z-basis of u is
  a Z_p-basis of u (x) Z_p whose reduction modulo p is independent.

  Its use (LocalOrthogonalImageGenerators): the Cayley transform Y -> (1 +
  Y)(1 - Y)^{-1} maps p^j u onto the congruence subgroup Gamma(p^j) = {g in
  U(X_p) : g = 1 mod p^j} for p odd and j >= 1, and 2^{j-1} u onto
  Gamma(2^j) for j >= 2 (for g in Gamma(p^j), Y = (g - 1)(g + 1)^{-1} is
  integral, divisible by p^j, resp. 2^{j-1}, skew-adjoint and in Lambda).
  The elements of the image modulo p^c congruent to 1 modulo p^j are the
  images of Gamma(p^j), j < c, so the layer G_j / G_{j+1} is {1 + p^j Y0 :
  Y0 in u} modulo p^{j+1}, of order p^{rank u}, generated by the Cayley
  transforms of p^j Y_i, resp. 2^{j-1} Y_i, Y_i a Z-basis of u: linear
  algebra, without lifting, for every e.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>>
IntegralLieLatticeColumns(MyMatrix<T> const &A,
                          std::vector<MyMatrix<Tint>> const &comm) {
  int m = A.rows();
  std::vector<MyMatrix<Tint>> B = comm;
  if (B.empty()) {
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        MyMatrix<Tint> E = ZeroMatrix<Tint>(m, m);
        E(a, b) = Tint(1);
        B.push_back(E);
      }
    }
  }
  // A made integral: the condition is homogeneous
  Tint den(1);
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      den = LCMpair(den, GetDenominator_z(A(a, b)));
    }
  }
  MyMatrix<Tint> Ai(m, m);
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      Ai(a, b) = UniversalScalarConversion<Tint, T>(
          A(a, b) * UniversalScalarConversion<T, Tint>(den));
    }
  }
  int nb = static_cast<int>(B.size());
  MyMatrix<Tint> Eq(nb, m * m);
  for (int i = 0; i < nb; i++) {
    MyMatrix<Tint> Z = B[i].transpose() * Ai + Ai * B[i];
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        Eq(i, a * m + b) = Z(a, b);
      }
    }
  }
  // the coefficient vectors x with x Eq = 0
  MyMatrix<Tint> K = NullspaceIntTrMat(MyMatrix<Tint>(Eq.transpose()));
  std::vector<MyMatrix<Tint>> res;
  for (int r = 0; r < K.rows(); r++) {
    MyMatrix<Tint> Y = ZeroMatrix<Tint>(m, m);
    for (int i = 0; i < nb; i++) {
      if (K(r, i) != Tint(0)) {
        MyMatrix<Tint> term = B[i] * K(r, i);
        Y += term;
      }
    }
    res.push_back(Y);
  }
  return res;
}

// Integer matrices as residues modulo M, for the Lie lattice given to
// LocalOrthogonalImageGenerators.
template <typename Tint>
std::vector<ModMatrix> ResidueMatrices(std::vector<MyMatrix<Tint>> const &l,
                                    int64_t M) {
  std::vector<ModMatrix> res;
  for (auto const &Y : l) {
    ModMatrix R(Y.rows(), std::vector<int64_t>(Y.cols()));
    for (int a = 0; a < Y.rows(); a++) {
      for (int b = 0; b < Y.cols(); b++) {
        R[a][b] = ResModOf(Y(a, b), M);
      }
    }
    res.push_back(R);
  }
  return res;
}

/*
  The generators of the layers G_j / G_{j+1}, j0 <= j < c, by the Cayley
  transforms of the Lie lattice (IntegralLieLatticeColumns), j0 = 1 for p
  odd and 2 for p = 2, modulo p^c, on the columns; the exponents of the
  layers are added to layer_exp.
 */
inline std::vector<ModMatrix> CayleyLayerGenerators(std::vector<ModMatrix> const &Lie,
                                                 int64_t p, int c,
                                                 size_t &layer_exp) {
  std::vector<ModMatrix> res;
  int64_t pc = ResiduePrimePower(p, c);
  int j0 = (p == 2) ? 2 : 1;
  for (int j = j0; j < c; j++) {
    int64_t scale = ResiduePrimePower(p, (p == 2) ? j - 1 : j);
    for (auto const &Y0 : Lie) {
      size_t m = Y0.size();
      ModMatrix Y(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          Y[a][b] = ResMul(ResMod(Y0[a][b], pc), ResMod(scale, pc), pc);
        }
      }
      auto mul = [&](ModMatrix const &X, ModMatrix const &Z) -> ModMatrix {
        return ModMatrixProduct(X, Z, pc);
      };
      // (1 - Y)^{-1} = sum_{k < c} Y^k, Y being divisible by p
      ModMatrix Id(m, std::vector<int64_t>(m, 0));
      for (size_t a = 0; a < m; a++) {
        Id[a][a] = ResMod(1, pc);
      }
      ModMatrix inv = Id, pw = Id;
      for (int k = 1; k < c; k++) {
        pw = mul(pw, Y);
        for (size_t a = 0; a < m; a++) {
          for (size_t b = 0; b < m; b++) {
            inv[a][b] = ResAdd(inv[a][b], pw[a][b], pc);
          }
        }
      }
      ModMatrix onePlus = Id;
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          onePlus[a][b] = ResAdd(onePlus[a][b], Y[a][b], pc);
        }
      }
      res.push_back(mul(onePlus, inv));
    }
    layer_exp += Lie.size();
  }
  return res;
}

/*
  Generators of the image G of O(X_p) in GL_m(Z / p^c) (columns: g^T A g =
  A), complete by construction and without listing G, which grows like
  p^{m(m-1)/2}:

  - G_1 = G mod p by its stabilizer chain on the basis e_0, ..., e_{m-1}.
    From the last level up, the orbit of e_i under the elements found so
    far is closed, and for each admissible image y of e_i (B(e_k, y) =
    A_{k i} for k < i, B(y, y) = A_{i i}) outside it, an element fixing
    e_0, ..., e_{i-1} with e_i -> y is searched: the further columns one
    at a time on their affine spaces, the matrix invertible, and a lift to
    p^K (Lemma: then it comes from an isometry). Once the orbit holds every
    y with such an element, the elements found generate the stabilizer H
    of e_0, ..., e_{i-1} in G_1, as they contain the stabilizer of e_i and
    move e_i over its whole H-orbit.
  - The layers G_j / G_{j+1}, 1 <= j < c, of the elements congruent to 1
    modulo p^j: F_p-spaces of the Delta with 1 + p^j Delta lifting to p^K,
    given by a basis of the solutions when its elements all lift (the
    liftable ones form a subspace V), else from the basis vectors that
    lift and, for each that does not, a lifting sum with a list of
    vectors independent modulo V (at most p^{codim V} tries each).

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
inline std::optional<std::vector<ModMatrix>> LocalOrthogonalImageGeneratorsCore(
    ModMatrix const &A, ModMatrix const &Be, std::vector<ModMatrix> const &Om,
    std::vector<ModMatrix> const &Lam, std::vector<ModMatrix> const &Bil, int64_t p,
    int c, int e, std::vector<int> const &levels, int layer_max,
    MultiplicityFactorBounds const &bounds, LocalImageOrder *order) {
  size_t m = A.size();
  int t = (p == 2) ? 1 : 0;
  // the points of F_p^m are encoded in a size_t; the orbits are kept in a
  // table of all the points when there are at most 2^24 of them, else in a
  // hash set (the orbits of the image are then expected to be small)
  size_t n_pts = 1;
  bool small_space = true;
  for (size_t i = 0; i < m; i++) {
    if (n_pts > std::numeric_limits<size_t>::max() / 2 / static_cast<size_t>(p)) {
      return {};
    }
    n_pts *= static_cast<size_t>(p);
    if (n_pts > (size_t(1) << 24)) {
      small_space = false;
    }
  }
  int64_t pc = MyPow<int64_t>(p, c);
  ModMatrix Ap(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      Ap[a][b] = ResMod(A[a][b], p);
    }
  }
  auto bil = [&](std::vector<int64_t> const &x,
                 std::vector<int64_t> const &y) -> int64_t {
    int64_t s = 0;
    for (size_t a = 0; a < m; a++) {
      if (x[a] != 0) {
        for (size_t b = 0; b < m; b++) {
          s = (s + x[a] * Ap[a][b] % p * y[b]) % p;
        }
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
  auto inv_mod = [](int64_t a, int64_t M) -> int64_t { return ResInv(a, M); };
  auto pow_p = [&](int k) -> int64_t {
    int64_t r = MyPow<int64_t>(p, k);
    return r;
  };
  auto mulmod = [](int64_t x, int64_t y, int64_t M) -> int64_t {
    return ResMul(x, y, M);
  };
  auto matmul = [](ModMatrix const &X, ModMatrix const &Y, int64_t M) -> ModMatrix {
    return ModMatrixProduct(X, Y, M);
  };
  auto transpose = [&](ModMatrix const &X) -> ModMatrix {
    ModMatrix Z(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Z[a][b] = X[b][a];
      }
    }
    return Z;
  };
  auto reduce = [&](ModMatrix X, int64_t M) -> ModMatrix {
    for (auto &row : X) {
      for (auto &x : row) {
        x = ResMod(x, M);
      }
    }
    return X;
  };
  // the inverse modulo p^k of g, invertible modulo p: modulo p by Gauss,
  // then X <- X (2 - g X)
  auto inverse_mod = [&](ModMatrix const &g, int k) -> std::optional<ModMatrix> {
    ModMatrix aug(m, std::vector<int64_t>(2 * m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        aug[a][b] = ResMod(g[a][b], p);
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
            aug[a][b] = ResMod(aug[a][b] - f0 * aug[col][b], p);
          }
        }
      }
    }
    ModMatrix X(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        X[a][b] = aug[a][m + b];
      }
    }
    int prec = 1;
    while (prec < k) {
      prec = std::min(2 * prec, k);
      int64_t M = pow_p(prec);
      ModMatrix gX = matmul(reduce(g, M), X, M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          gX[a][b] = ResMod(((a == b) ? 2 : 0) - gX[a][b], M);
        }
      }
      X = matmul(X, gX, M);
    }
    return X;
  };
  // the valuation of g^T A g - A, up to k
  auto error_val = [&](ModMatrix const &g, int k) -> int {
    int64_t M = pow_p(k);
    ModMatrix E = matmul(transpose(g), matmul(reduce(A, M), g, M), M);
    int v = k;
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        int64_t x = ResMod(E[a][b] - A[a][b], M);
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
  auto newton = [&](ModMatrix g) -> std::optional<ModMatrix> {
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
      std::optional<ModMatrix> gi = inverse_mod(g, prec);
      if (!gi) {
        return {};
      }
      ModMatrix E = matmul(transpose(g), matmul(reduce(A, M), g, M), M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          E[a][b] = ResMod(E[a][b] - A[a][b], M);
        }
      }
      ModMatrix F = matmul(reduce(Be, M), matmul(transpose(*gi), E, M), M);
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
          g[a][b] = ResMod(g[a][b] - x, NM);
        }
      }
      prec = nprec;
    }
    return {};
  };
  // With a commutant Lam (a center K), the lifts g + p^j Delta keep Delta
  // in Lam mod p: the system of MakeLiftSystem restricted to the
  // coordinates of Lam, its solutions written back as matrices.
  size_t d_lam = Lam.size();
  ModMatrix Mlam;  // m^2 x d_lam, the vec of the Lam_lambda modulo p
  if (d_lam > 0) {
    Mlam.assign(m * m, std::vector<int64_t>(d_lam, 0));
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t a2 = 0; a2 < m; a2++) {
        for (size_t b2 = 0; b2 < m; b2++) {
          Mlam[a2 * m + b2][l] = ResMod(Lam[l][a2][b2], p);
        }
      }
    }
  }
  auto restrict_coef = [&](ModMatrix const &coef) -> ModMatrix {
    if (d_lam == 0) {
      return coef;
    }
    ModMatrix res(coef.size(), std::vector<int64_t>(d_lam, 0));
    for (size_t r2 = 0; r2 < coef.size(); r2++) {
      for (size_t v = 0; v < m * m; v++) {
        if (coef[r2][v] != 0) {
          for (size_t l = 0; l < d_lam; l++) {
            res[r2][l] = (res[r2][l] + coef[r2][v] * Mlam[v][l]) % p;
          }
        }
      }
    }
    return res;
  };
  auto to_matrix_vec = [&](std::vector<int64_t> const &dv) -> std::vector<int64_t> {
    std::vector<int64_t> x(m * m, 0);
    for (size_t l = 0; l < d_lam; l++) {
      if (dv[l] != 0) {
        for (size_t v = 0; v < m * m; v++) {
          x[v] = (x[v] + dv[l] * Mlam[v][l]) % p;
        }
      }
    }
    return x;
  };
  auto lift_aff = [&](ModMatrix const &g,
                      int64_t pjl) -> std::optional<LiftAffineSpace> {
    if (d_lam == 0) {
      return LiftAffine(g, A, p, pjl);
    }
    LiftSystem sys = MakeLiftSystem(g, A, p, pjl);
    ModMatrix cf = restrict_coef(sys.coef);
    ModMatrix aug;
    for (size_t r2 = 0; r2 < cf.size(); r2++) {
      std::vector<int64_t> row = cf[r2];
      row.push_back(sys.rhs[r2]);
      aug.push_back(row);
    }
    RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
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
        v[rs.pivots[i2]] = ResMod(-rs.basis[i2][f], p);
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
    ModMatrix test;
    for (size_t v = 0; v < m * m && lam_rows.size() < d_lam; v++) {
      test.push_back(Mlam[v]);
      if (RankModPrime(test, p) == test.size()) {
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
  // the inverse modulo p^k of the matrix S of the chosen rows of Lam,
  // computed once for each k (project is called at every node)
  std::map<int, ModMatrix> sinv_cache;
  auto lam_sinv = [&](int k) -> ModMatrix const & {
    auto iter = sinv_cache.find(k);
    if (iter != sinv_cache.end()) {
      return iter->second;
    }
    int64_t M = pow_p(k);
    // S: the d x d matrix of the chosen rows, its inverse modulo p^k
    ModMatrix S(d_lam, std::vector<int64_t>(d_lam));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      size_t v = lam_rows[r2];
      for (size_t l = 0; l < d_lam; l++) {
        S[r2][l] = ResMod(Lam[l][v / m][v % m], M);
      }
    }
    // the inverse modulo p, then Newton
    ModMatrix aug(d_lam, std::vector<int64_t>(2 * d_lam, 0));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      for (size_t l = 0; l < d_lam; l++) {
        aug[r2][l] = ResMod(S[r2][l], p);
      }
      aug[r2][d_lam + r2] = 1;
    }
    RowSpaceModP rr = ComputeRowSpaceModP(aug, p);
    ModMatrix Sinv(d_lam, std::vector<int64_t>(d_lam));
    for (size_t r2 = 0; r2 < d_lam; r2++) {
      for (size_t l = 0; l < d_lam; l++) {
        Sinv[r2][l] = rr.basis[r2][d_lam + l];
      }
    }
    int prec2 = 1;
    while (prec2 < k) {
      prec2 = std::min(2 * prec2, k);
      int64_t M2 = pow_p(prec2);
      ModMatrix SX(d_lam, std::vector<int64_t>(d_lam, 0));
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t k2 = 0; k2 < d_lam; k2++) {
          for (size_t b2 = 0; b2 < d_lam; b2++) {
            SX[a2][b2] = (SX[a2][b2] + mulmod(ResMod(S[a2][k2], M2), Sinv[k2][b2], M2)) % M2;
          }
        }
      }
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t b2 = 0; b2 < d_lam; b2++) {
          SX[a2][b2] = ResMod(((a2 == b2) ? 2 : 0) - SX[a2][b2], M2);
        }
      }
      ModMatrix NS(d_lam, std::vector<int64_t>(d_lam, 0));
      for (size_t a2 = 0; a2 < d_lam; a2++) {
        for (size_t k2 = 0; k2 < d_lam; k2++) {
          for (size_t b2 = 0; b2 < d_lam; b2++) {
            NS[a2][b2] = (NS[a2][b2] + mulmod(Sinv[a2][k2], SX[k2][b2], M2)) % M2;
          }
        }
      }
      Sinv = NS;
    }
    return sinv_cache.emplace(k, Sinv).first->second;
  };
  auto project = [&](ModMatrix const &g, int k) -> ModMatrix {
    if (d_lam == 0) {
      return g;
    }
    int64_t M = pow_p(k);
    ModMatrix const &Sinv = lam_sinv(k);
    std::vector<int64_t> coefs(d_lam, 0);
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t r2 = 0; r2 < d_lam; r2++) {
        size_t v = lam_rows[r2];
        coefs[l] = (coefs[l] + mulmod(Sinv[l][r2], ResMod(g[v / m][v % m], M), M)) % M;
      }
    }
    ModMatrix h(m, std::vector<int64_t>(m, 0));
    for (size_t l = 0; l < d_lam; l++) {
      for (size_t a2 = 0; a2 < m; a2++) {
        for (size_t b2 = 0; b2 < m; b2++) {
          h[a2][b2] = (h[a2][b2] + mulmod(coefs[l], ResMod(Lam[l][a2][b2], pWfull), pWfull)) % pWfull;
        }
      }
    }
    return h;
  };
  // a solution modulo p^target extending the solution g modulo p^j, the
  // first found by a depth-first search over the affine spaces
  std::function<std::optional<ModMatrix>(ModMatrix const &, int, int64_t, int)>
      extend = [&](ModMatrix const &g, int j, int64_t pjl,
                   int target) -> std::optional<ModMatrix> {
    if (j >= target) {
      return g;
    }
    if (over || ++n_node > bounds.max_node) {
      over = true;
      return {};
    }
    // with a commutant, g in Lam modulo p^{j+1} (it is modulo p^j), so that
    // the lifts g + p^j Delta, Delta in Lam, are in Lam modulo p^{j+1}
    ModMatrix gl = (d_lam > 0) ? reduce(project(g, j), pjl * p) : g;
    std::optional<LiftAffineSpace> sp = lift_aff(gl, pjl);
    if (!sp) {
      return {};
    }
    size_t dim = sp->ker.size();
    // The next level is solvable for g' = g + p^j Delta exactly when the
    // residual r = L rhs vanishes, L the left kernel of its (restricted)
    // matrix, which depends on g modulo p only. r is of degree <= 2 in
    // Delta, and affine for j >= 2: rhs = E(g')/p^{j+1} with g' = g + p^j
    // Delta has the quadratic part p^{j-1} Delta^T A Delta, and r does not
    // depend on the representative of g' modulo p^{j+2}. For j = 1 it is
    // checked on the pairs of basis vectors (and the doubles for p odd).
    // When r is affine the good Delta are an affine subspace, found by
    // solving, and the search runs on it only; otherwise every Delta is
    // tried.
    std::optional<LiftAffineSpace> good;
    if (j + 1 < target && dim > 0) {
      int64_t pj1 = pjl * p;
      auto next_sys = [&](std::vector<int64_t> const &cf0) -> LiftSystem {
        ModMatrix g2 = LiftAffineElement(gl, *sp, cf0, p, pjl);
        ModMatrix gp = (d_lam > 0) ? reduce(project(g2, j + 1), pj1 * p) : g2;
        return MakeLiftSystem(gp, A, p, pj1);
      };
      LiftSystem sys0 = next_sys(std::vector<int64_t>(dim, 0));
      ModMatrix cf = restrict_coef(sys0.coef);
      size_t n_eq = cf.size();
      size_t n_var = cf.empty() ? 0 : cf[0].size();
      ModMatrix cfT(n_var, std::vector<int64_t>(n_eq));
      for (size_t a = 0; a < n_eq; a++) {
        for (size_t b = 0; b < n_var; b++) {
          cfT[b][a] = cf[a][b];
        }
      }
      RowSpaceModP rsT = ComputeRowSpaceModP(cfT, p);
      std::vector<bool> is_piv(n_eq, false);
      for (size_t pv : rsT.pivots) {
        is_piv[pv] = true;
      }
      ModMatrix L;
      for (size_t f = 0; f < n_eq; f++) {
        if (is_piv[f]) {
          continue;
        }
        std::vector<int64_t> v(n_eq, 0);
        v[f] = 1;
        for (size_t i2 = 0; i2 < rsT.pivots.size(); i2++) {
          v[rsT.pivots[i2]] = ResMod(-rsT.basis[i2][f], p);
        }
        L.push_back(v);
      }
      auto residual = [&](std::vector<int64_t> const &cf0) -> std::vector<int64_t> {
        ModMatrix g2 = LiftAffineElement(gl, *sp, cf0, p, pjl);
        ModMatrix gp = (d_lam > 0) ? reduce(project(g2, j + 1), pj1 * p) : g2;
        std::vector<int64_t> rv = LiftRhs(gp, A, p, pj1);
        std::vector<int64_t> r(L.size(), 0);
        for (size_t k2 = 0; k2 < L.size(); k2++) {
          for (size_t a = 0; a < n_eq; a++) {
            r[k2] = (r[k2] + L[k2][a] * rv[a]) % p;
          }
        }
        return r;
      };
      std::vector<int64_t> r0 = residual(std::vector<int64_t>(dim, 0));
      std::vector<std::vector<int64_t>> ri(dim);
      for (size_t i = 0; i < dim; i++) {
        std::vector<int64_t> ci(dim, 0);
        ci[i] = 1;
        ri[i] = residual(ci);
      }
      bool affine_r = true;
#ifdef SANITY_CHECK_MULT
      bool check_pairs = true;
#else
      bool check_pairs = (j < 2);
#endif
      for (size_t i = 0; i < dim && affine_r && check_pairs; i++) {
        for (size_t k2 = i; k2 < dim && affine_r; k2++) {
          if (k2 == i && p == 2) {
            continue;
          }
          std::vector<int64_t> ci(dim, 0);
          ci[i] = (ci[i] + 1) % p;
          ci[k2] = (ci[k2] + 1) % p;
          std::vector<int64_t> rik = residual(ci);
          for (size_t a = 0; a < L.size(); a++) {
            if (ResMod(rik[a] - ri[i][a] - ri[k2][a] + r0[a], p) != 0) {
              affine_r = false;
            }
          }
        }
      }
      if (!affine_r && j >= 2) {
        std::cerr << "MULT_FACTOR: the residual is not affine at j=" << j
                  << "\n";
        throw TerminalException{1};
      }
      if (affine_r) {
        // sum_i c_i (r_i - r_0) = -r_0
        ModMatrix aug;
        for (size_t a = 0; a < L.size(); a++) {
          std::vector<int64_t> row(dim + 1);
          for (size_t i = 0; i < dim; i++) {
            row[i] = ResMod(ri[i][a] - r0[a], p);
          }
          row[dim] = ResMod(-r0[a], p);
          aug.push_back(row);
        }
        LiftAffineSpace gs;
        gs.part.assign(dim, 0);
        if (!aug.empty()) {
          RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
          for (size_t pv : rs.pivots) {
            if (pv == dim) {
              return {};
            }
          }
          std::vector<bool> piv(dim, false);
          for (size_t i2 = 0; i2 < rs.pivots.size(); i2++) {
            piv[rs.pivots[i2]] = true;
            gs.part[rs.pivots[i2]] = rs.basis[i2][dim];
          }
          for (size_t f = 0; f < dim; f++) {
            if (piv[f]) {
              continue;
            }
            std::vector<int64_t> v(dim, 0);
            v[f] = 1;
            for (size_t i2 = 0; i2 < rs.pivots.size(); i2++) {
              v[rs.pivots[i2]] = ResMod(-rs.basis[i2][f], p);
            }
            gs.ker.push_back(v);
          }
        } else {
          for (size_t f = 0; f < dim; f++) {
            std::vector<int64_t> v(dim, 0);
            v[f] = 1;
            gs.ker.push_back(v);
          }
        }
        good = gs;
      }
    }
    // the coefficients in the space of Delta from those in the good space
    auto to_coef = [&](std::vector<int64_t> const &u) -> std::vector<int64_t> {
      if (!good) {
        return u;
      }
      std::vector<int64_t> c2 = good->part;
      for (size_t i = 0; i < good->ker.size(); i++) {
        if (u[i] != 0) {
          for (size_t v = 0; v < dim; v++) {
            c2[v] = (c2[v] + u[i] * good->ker[i][v]) % p;
          }
        }
      }
      return c2;
    };
    size_t dim_s = good ? good->ker.size() : dim;
    std::vector<int64_t> coef(dim_s, 0);
    while (!over) {
      std::optional<ModMatrix> r = extend(
          LiftAffineElement(gl, *sp, to_coef(coef), p, pjl), j + 1, pjl * p,
          target);
      if (r) {
        return r;
      }
      size_t i = 0;
      while (i < dim_s && coef[i] == p - 1) {
        coef[i] = 0;
        i++;
      }
      if (i == dim_s) {
        return {};
      }
      coef[i]++;
    }
    return {};
  };
  // the element of G (modulo p^c) congruent to the solution g modulo p^j,
  // when there is one: a solution modulo p^k, k = max(j + e + t, e + 2t +
  // 1), extending g, then Newton
  auto lift_path = [&](ModMatrix const &g, int j,
                       int64_t pjl) -> std::optional<ModMatrix> {
    int target = std::max(j + e + t, e + 2 * t + 1);
    std::optional<ModMatrix> sol = extend(g, j, pjl, target);
    if (!sol) {
      return {};
    }
    std::optional<ModMatrix> h = newton(project(*sol, target));
    if (!h) {
      std::cerr << "MULT_FACTOR: Newton did not converge\n";
      throw TerminalException{1};
    }
    return h;
  };
  // the affine space {x : fixed_k x = rhs_k}, the fixed_k linear forms
  auto affine = [&](std::vector<std::vector<int64_t>> const &fixed,
                    std::vector<int64_t> const &rhs)
      -> std::optional<std::pair<std::vector<int64_t>, ModMatrix>> {
    if (fixed.empty()) {
      ModMatrix ker;
      for (size_t v = 0; v < m; v++) {
        std::vector<int64_t> e1(m, 0);
        e1[v] = 1;
        ker.push_back(e1);
      }
      return std::make_pair(std::vector<int64_t>(m, 0), ker);
    }
    ModMatrix aug;
    for (size_t k = 0; k < fixed.size(); k++) {
      std::vector<int64_t> row(m + 1, 0);
      for (size_t b = 0; b < m; b++) {
        row[b] = fixed[k][b];
      }
      row[m] = ResMod(rhs[k], p);
      aug.push_back(row);
    }
    RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
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
    ModMatrix ker;
    for (size_t f = 0; f < m; f++) {
      if (is_piv[f]) {
        continue;
      }
      std::vector<int64_t> v(m, 0);
      v[f] = 1;
      for (size_t i = 0; i < rs.pivots.size(); i++) {
        v[rs.pivots[i]] = ResMod(-rs.basis[i][f], p);
      }
      ker.push_back(v);
    }
    return std::make_pair(part, ker);
  };
  // f(x) for the x of the affine space with B(x, x) = q, until f is true
  auto for_each_vector =
      [&](std::pair<std::vector<int64_t>, ModMatrix> const &sp, int64_t q,
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
  // The levels (a Jordan basis, JordanBasisModPW): column j lies in
  // X^(l_j) = {x : A x = 0 mod p^{l_j}}, and so does the image y of e_j:
  // zero coordinates on the columns of lower level. The image differs from
  // its representative (coordinates in [0, p) on the columns of level >=
  // l_j) by p z, z in X^(l_j - 1) (not X^(l_j): the coordinates on the
  // lower columns are divisible by p, not zero), so the representatives
  // fix B(y, y) modulo p^{l_j + 1} (modulo 4 at level 0 when p = 2) and
  // B(y_k, y) modulo p^{l + 1} only when l_k = l_j = l; between different
  // levels nothing beyond the conditions modulo p.
  auto ppow = [&](int k) -> __int128 {
    __int128 r = MyPow<__int128>(p, k);
    return r;
  };
  auto pair_mod = [&](std::vector<int64_t> const &u,
                      std::vector<int64_t> const &w, __int128 M) -> __int128 {
    __int128 sacc = 0;
    for (size_t a = 0; a < m; a++) {
      if (u[a] != 0) {
        for (size_t b = 0; b < m; b++) {
          if (w[b] == 0) {
            continue;
          }
          sacc = (sacc + static_cast<__int128>(u[a]) * (A[a][b] % M) % M * w[b]) % M;
        }
      }
    }
    return ((sacc % M) + M) % M;
  };
  auto modM = [&](int64_t x, __int128 M) -> __int128 {
    __int128 r = static_cast<__int128>(x) % M;
    return (r < 0) ? r + M : r;
  };
  auto radical_ok = [&](std::vector<int64_t> const &x, size_t i) -> bool {
    int li = levels[i];
    if (li == 0) {
      return true;
    }
    for (size_t a = 0; a < m; a++) {
      if (levels[a] < li && x[a] != 0) {
        return false;
      }
    }
    __int128 Md = ppow(li + 1);
    if (pair_mod(x, x, Md) != modM(A[i][i], Md)) {
      return false;
    }
    for (size_t k = 0; k < i; k++) {
      if (levels[k] != li) {
        continue;
      }
      __int128 Mk = ppow(li + 1);
      if (pair_mod(x, cols[k], Mk) != modM(A[k][i], Mk)) {
        return false;
      }
    }
    return true;
  };
  // p = 2: B(y, y) modulo 4 depends only on y modulo 2 (B(y + 2z, y + 2z)
  // = B(y, y) + 4 (B(y, z) + B(z, z))), so the image y of e_j has B(y, y)
  // = A_jj modulo 4, which the bilinear conditions modulo 2 miss on an
  // even lattice
  auto diag4_ok = [&](std::vector<int64_t> const &y, size_t j) -> bool {
    if (p != 2) {
      return true;
    }
    __int128 sacc = 0;
    for (size_t a = 0; a < m; a++) {
      if (y[a] != 0) {
        for (size_t b = 0; b < m; b++) {
          if (y[b] != 0) {
            sacc += A[a][b];
          }
        }
      }
    }
    __int128 d = (sacc - A[j][j]) % 4;
    return d == 0;
  };
  // With an action Om of O_K (a center K): the image of e_j is forced when
  // e_j is in the span modulo p of the w e_k, k < j, w in Om; otherwise the
  // conditions B(w y_k, y) = B(w e_k, e_j) for all k < j and w, and B(w y,
  // y) = B(w e_j, e_j) for all w.
  size_t n_om = Om.size();
  std::vector<std::vector<std::vector<int64_t>>> Omc(n_om);
  for (size_t w = 0; w < n_om; w++) {
    Omc[w].assign(m, std::vector<int64_t>(m));
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        Omc[w][a2][b2] = ResMod(Om[w][a2][b2], p);
      }
    }
  }
  // further integral bilinear forms N preserved by the group (the
  // B(w x, y) for w in the trace dual of the order acting): the conditions
  // y_k^T N y = N_kj modulo p, sharper than those of Om where the order is
  // ramified
  std::vector<ModMatrix> Bilp;
  for (auto const &N : Bil) {
    ModMatrix Np(m, std::vector<int64_t>(m));
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        Np[a2][b2] = ResMod(N[a2][b2], p);
      }
    }
    Bilp.push_back(Np);
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
    ModMatrix span;
    for (size_t j = 0; j < m; j++) {
      std::vector<int64_t> ej(m, 0);
      ej[j] = 1;
      if (!span.empty()) {
        // solve sum c_r span_r = e_j
        ModMatrix aug(m, std::vector<int64_t>(span.size() + 1, 0));
        for (size_t r2 = 0; r2 < span.size(); r2++) {
          for (size_t a2 = 0; a2 < m; a2++) {
            aug[a2][r2] = span[r2][a2];
          }
        }
        for (size_t a2 = 0; a2 < m; a2++) {
          aug[a2][span.size()] = ej[a2];
        }
        RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
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

    for (auto const &N : Bilp) {
      int64_t sacc = 0;
      for (size_t a = 0; a < m; a++) {
        if (y[a] != 0) {
          for (size_t b = 0; b < m; b++) {
            sacc = (sacc + y[a] * N[a][b] % p * y[b]) % p;
          }
        }
      }
      if (sacc != N[j][j]) {
        return false;
      }
    }
    return true;
  };
  // the rows vec(Lam_l) modulo p, independent (lam_rows): g modulo p is in
  // Lam modulo p when adding vec(g) keeps the rank. Commuting with Om
  // modulo p is weaker when Om modulo p spans less than the order it
  // generates over Z_p does.
  ModMatrix lam_vecs;
  for (size_t l = 0; l < d_lam; l++) {
    std::vector<int64_t> v(m * m);
    for (size_t w = 0; w < m * m; w++) {
      v[w] = Mlam[w][l];
    }
    lam_vecs.push_back(v);
  }
  auto in_lam = [&](ModMatrix const &g) -> bool {
    if (d_lam == 0) {
      return true;
    }
    ModMatrix test = lam_vecs;
    std::vector<int64_t> v(m * m);
    for (size_t a2 = 0; a2 < m; a2++) {
      for (size_t b2 = 0; b2 < m; b2++) {
        v[a2 * m + b2] = ResMod(g[a2][b2], p);
      }
    }
    test.push_back(v);
    return RankModPrime(test, p) == d_lam;
  };
  // the full check of a candidate modulo p: isometry, invertible,
  // commuting with Om and in Lam
  auto full_ok = [&](ModMatrix const &g) -> bool {
    if (!in_lam(g)) {
      return false;
    }
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
    return RankModPrime(g, p) == m;
  };
  // the affine conditions on the image of e_j given the columns k < j
  auto conditions = [&](size_t j, std::vector<std::vector<int64_t>> &fixed,
                        std::vector<int64_t> &rhs) -> void {
    fixed.clear();
    rhs.clear();
    // the linear form x -> u^T N x modulo p
    auto form = [&](std::vector<int64_t> const &u,
                    ModMatrix const &N) -> std::vector<int64_t> {
      std::vector<int64_t> row(m, 0);
      for (size_t a = 0; a < m; a++) {
        if (u[a] != 0) {
          for (size_t b = 0; b < m; b++) {
            row[b] = (row[b] + u[a] * N[a][b]) % p;
          }
        }
      }
      return row;
    };
    int lj = levels[j];
    // the columns of lower level: zero coordinates
    for (size_t a = 0; a < m; a++) {
      if (levels[a] < lj) {
        std::vector<int64_t> row(m, 0);
        row[a] = 1;
        fixed.push_back(row);
        rhs.push_back(0);
      }
    }
    // the linear form y -> B(u, y) / p^M modulo p on the coordinates of
    // level >= l_j, and b / p^M, u in X^(M) or M <= l_j
    auto level_form = [&](std::vector<int64_t> const &u, int M,
                          int64_t bval) -> void {
      __int128 Mk = ppow(M + 1);
      __int128 pM = ppow(M);
      std::vector<int64_t> row(m, 0);
      for (size_t b = 0; b < m; b++) {
        if (levels[b] < lj) {
          continue;
        }
        __int128 sacc = 0;
        for (size_t a = 0; a < m; a++) {
          if (u[a] != 0) {
            sacc = (sacc + static_cast<__int128>(u[a]) * modM(A[a][b], Mk)) % Mk;
          }
        }
        if (sacc % pM != 0) {
          std::cerr << "MULT_FACTOR: a pairing is not divisible by its level, p="
                    << p << " j=" << j << " b=" << b << " M=" << M << " levels";
          for (int l : levels) {
            std::cerr << " " << l;
          }
          std::cerr << " u =";
          for (auto v : u) {
            std::cerr << " " << v;
          }
          std::cerr << " A row b mod p^(M+1):";
          for (size_t a = 0; a < m; a++) {
            std::cerr << " " << static_cast<int64_t>(modM(A[a][b], Mk));
          }
          std::cerr << "\n";
          throw TerminalException{1};
        }
        row[b] = static_cast<int64_t>((sacc / pM) % p);
      }
      __int128 bm = modM(bval, Mk);
      if (bm % pM != 0) {
        std::cerr << "MULT_FACTOR: an entry of A is not divisible by its level\n";
        throw TerminalException{1};
      }
      fixed.push_back(row);
      rhs.push_back(static_cast<int64_t>((bm / pM) % p));
    };
    for (size_t k = 0; k < j; k++) {
      if (lj > 0 && levels[k] == lj) {
        // B(y_k, y) / p^l modulo p, both of level l
        level_form(cols[k], lj, A[k][j]);
      } else {
        fixed.push_back(form(cols[k], Ap));
        rhs.push_back(Ap[k][j]);
      }
      for (size_t w = 0; w < n_om; w++) {
        fixed.push_back(form(om_apply(w, cols[k]), Ap));
        rhs.push_back(om_entry(w, k, j));
      }
      for (auto const &N : Bilp) {
        fixed.push_back(form(cols[k], N));
        rhs.push_back(N[k][j]);
        // and x^T N y_k = N_jk, N need not be symmetric
        std::vector<int64_t> row(m, 0);
        for (size_t a = 0; a < m; a++) {
          for (size_t b = 0; b < m; b++) {
            row[a] = (row[a] + N[a][b] * cols[k][b]) % p;
          }
        }
        fixed.push_back(row);
        rhs.push_back(N[j][k]);
      }
    }
  };
  // the columns 0, ..., j independent modulo p: g is invertible, and a
  // forced column (an image of the w y_k) can make them dependent, which
  // is better seen at once than after the further columns
  size_t n_full_fail = 0, n_lift_fail = 0, n_dual_fail = 0;
  size_t n_forced_rad = 0, n_forced_diag = 0, n_forced_dep = 0;
  // g^T A g = A is g Be g^T = Be, Be = p^e A^{-1}: modulo p for g modulo
  // p, and on the diagonal modulo 4 when p = 2 (as for diag4_ok). When A
  // is singular modulo p (e >= 1) this is not implied by the conditions on
  // the columns, and it rejects at once most of the solutions modulo p
  // that do not lift.
  int64_t mod_dual = (p == 2) ? 4 : p;
  auto dual_ok = [&](ModMatrix const &g) -> bool {
    if (e == 0) {
      return true;
    }
    for (size_t a = 0; a < m; a++) {
      for (size_t b = a; b < m; b++) {
        __int128 sacc = 0;
        for (size_t u = 0; u < m; u++) {
          if (g[a][u] != 0) {
            for (size_t v = 0; v < m; v++) {
              if (g[b][v] != 0) {
                sacc += static_cast<__int128>(g[a][u]) * (Be[u][v] % mod_dual) * g[b][v];
              }
            }
          }
        }
        int64_t md = (a == b) ? mod_dual : p;
        __int128 d = (sacc - Be[a][b]) % md;
        if (d != 0) {
          return false;
        }
      }
    }
    return true;
  };
  auto independent = [&](size_t j) -> bool {
    ModMatrix Mc(cols.begin(), cols.begin() + j + 1);
    return RankModPrime(Mc, p) == j + 1;
  };
  std::function<std::optional<ModMatrix>(size_t)> complete =
      [&](size_t j) -> std::optional<ModMatrix> {
    if (j == m) {
      ModMatrix g(m, std::vector<int64_t>(m));
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          g[a][b] = cols[b][a];
        }
      }
      if (!full_ok(g)) {
        n_full_fail++;
        return {};
      }
      if (!dual_ok(g)) {
        n_dual_fail++;
        return {};
      }
      std::optional<ModMatrix> lg = lift_path(g, 1, p);
      if (!lg) {
        n_lift_fail++;
      }
      return lg;
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
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
        if (n_forced_rad < 3) {
          std::cerr << "MULT_FACTOR: forced column " << j << " level "
                    << levels[j] << " y =";
          for (auto v : y) {
            std::cerr << " " << v;
          }
          std::cerr << " | cols:";
          for (size_t k = 0; k < j; k++) {
            std::cerr << " [";
            for (auto v : cols[k]) {
              std::cerr << v;
            }
            std::cerr << "]";
          }
          std::cerr << "\n";
        }
#endif
        n_forced_rad++;
        return {};
      }
      if (!diag4_ok(y, j)) {
        n_forced_diag++;
        return {};
      }
      cols[j] = y;
      if (!independent(j)) {
        n_forced_dep++;
        return {};
      }
      return complete(j + 1);
    }
    std::vector<std::vector<int64_t>> fixed;
    std::vector<int64_t> rhs;
    conditions(j, fixed, rhs);
    std::optional<std::pair<std::vector<int64_t>, ModMatrix>> sp =
        affine(fixed, rhs);
    if (!sp) {
      return {};
    }
    std::optional<ModMatrix> res;
    for_each_vector(*sp, Ap[j][j], [&](std::vector<int64_t> const &x) -> bool {
      if (++n_node > bounds.max_node) {
        over = true;
        return true;
      }
      if (!radical_ok(x, j) || !diag4_ok(x, j) || !om_ok(x, j)) {
        return false;
      }
      cols[j] = x;
      if (!independent(j)) {
        return false;
      }
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
  auto apply = [&](ModMatrix const &g,
                   std::vector<int64_t> const &x) -> std::vector<int64_t> {
    std::vector<int64_t> y(m, 0);
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        y[a] = (y[a] + ResMod(g[a][b], p) * x[b]) % p;
      }
    }
    return y;
  };
  std::vector<ModMatrix> gens1;
  std::vector<ModMatrix> gens;
  LocalImageOrder ord;
  for (size_t i = m; i-- > 0;) {
    std::vector<int64_t> ei(m, 0);
    ei[i] = 1;
    // the orbit of e_i under gens1
    std::vector<char> in_orbit_tab(small_space ? n_pts : 0, 0);
    std::unordered_set<size_t> in_orbit_set;
    auto orbit_has = [&](size_t idx) -> bool {
      return small_space ? (in_orbit_tab[idx] != 0) : (in_orbit_set.count(idx) > 0);
    };
    auto orbit_add = [&](size_t idx) -> void {
      if (small_space) {
        in_orbit_tab[idx] = 1;
      } else {
        in_orbit_set.insert(idx);
      }
    };
    // the orbit, closed incrementally: from e_i under all of gens1 at
    // first, then for a new generator g its images of the orbit points and
    // the closure of the new points under all of gens1
    std::vector<std::vector<int64_t>> orbit_pts;
    auto bfs = [&](std::vector<std::vector<int64_t>> queue) -> void {
      while (!queue.empty()) {
        std::vector<int64_t> x = queue.back();
        queue.pop_back();
        for (auto const &g : gens1) {
          std::vector<int64_t> y = apply(g, x);
          size_t iy = encode(y);
          if (!orbit_has(iy)) {
            orbit_add(iy);
            orbit_pts.push_back(y);
            queue.push_back(y);
          }
        }
      }
    };
    auto close_init = [&]() -> size_t {
      orbit_add(encode(ei));
      orbit_pts.push_back(ei);
      bfs({ei});
      return orbit_pts.size();
    };
    auto close_new = [&](ModMatrix const &g) -> size_t {
      std::vector<std::vector<int64_t>> queue;
      size_t n0 = orbit_pts.size();
      for (size_t i0 = 0; i0 < n0; i0++) {
        std::vector<int64_t> y = apply(g, orbit_pts[i0]);
        size_t iy = encode(y);
        if (!orbit_has(iy)) {
          orbit_add(iy);
          orbit_pts.push_back(y);
          queue.push_back(y);
        }
      }
      bfs(queue);
      return orbit_pts.size();
    };
    size_t orbit_len = close_init();
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
    std::optional<std::pair<std::vector<int64_t>, ModMatrix>> sp =
        affine(fixed, rhs);
    if (!sp) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: no admissible image of e_i, n_node=" << n_node << "\n";
#endif
      return {};
    }
    size_t n_cand = 0, n_fail = 0;
    size_t node0 = n_node;
    (void)node0;
    for_each_vector(*sp, Ap[i][i], [&](std::vector<int64_t> const &y) -> bool {
      if (!radical_ok(y, i) || !diag4_ok(y, i) || !om_ok(y, i)) {
        return false;
      }
      n_cand++;
      if (orbit_has(encode(y))) {
        return false;
      }
      for (size_t k = 0; k < i; k++) {
        cols[k] = idcols[k];
      }
      cols[i] = y;
      if (!independent(i)) {
        n_fail++;
        return false;
      }
      std::optional<ModMatrix> g = complete(i + 1);
      if (over) {
        return true;
      }
      if (g) {
        ModMatrix g1 = *g;
        for (auto &row : g1) {
          for (auto &x : row) {
            x = ResMod(x, p);
          }
        }
        gens1.push_back(g1);
        gens.push_back(*g);
        orbit_len = close_new(g1);
      } else {
        n_fail++;
      }
      return false;
    });
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: level " << i << " candidates " << n_cand
              << " failed " << n_fail << " orbit " << orbit_len << " nodes "
              << n_node - node0 << " full_fail " << n_full_fail
              << " dual_fail " << n_dual_fail << " lift_fail " << n_lift_fail
              << " forced rad/diag/dep " << n_forced_rad << "/" << n_forced_diag
              << "/" << n_forced_dep << "\n";
#endif
    if (over) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: budget at level 1, n_node=" << n_node << "\n";
#endif
      return {};
    }
    ord.orbits.push_back(orbit_len);
  }
  // the layers below layer_max (the others by the Cayley transforms of
  // the Lie lattice, in LocalOrthogonalImageGenerators)
  ModMatrix Id(m, std::vector<int64_t>(m, 0));
  for (size_t a = 0; a < m; a++) {
    Id[a][a] = 1;
  }
  int64_t pj = p;
  for (int j = 1; j < std::min(c, layer_max); j++) {
    std::optional<LiftAffineSpace> sp = lift_aff(Id, pj);
    if (!sp) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: no layer, n_node=" << n_node << "\n";
#endif
      return {};
    }
    // the basis vectors of the layer that lift, and those that do not
    std::vector<ModMatrix> lifts;
    std::vector<size_t> lifted, failed;
    for (size_t i0 = 0; i0 < sp->ker.size(); i0++) {
      std::vector<int64_t> coef(sp->ker.size(), 0);
      coef[i0] = 1;
      ModMatrix g = LiftAffineElement(Id, *sp, coef, p, pj);
      std::optional<ModMatrix> lf = lift_path(g, j + 1, pj * p);
      if (over) {
        return {};
      }
      if (lf) {
        lifts.push_back(*lf);
        lifted.push_back(i0);
      } else {
        failed.push_back(i0);
      }
    }
    bool all = failed.empty();
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: layer " << j << " dim " << sp->ker.size()
              << " lifted " << lifted.size() << " failed " << failed.size()
              << "\n";
#endif
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
      auto layer_elt = [&](std::vector<int64_t> const &coefv) -> ModMatrix {
        return LiftAffineElement(Id, *sp, coefv, p, pj);
      };
      LiftSystem sys0 = MakeLiftSystem(layer_elt(std::vector<int64_t>(dim, 0)),
                                       A, p, pj1);
      ModMatrix coef0 = restrict_coef(sys0.coef);
      size_t n_eq = coef0.size();
      size_t n_var = coef0.empty() ? 0 : coef0[0].size();
      // the left kernel: the y with y^T coef = 0
      ModMatrix coefT(n_var, std::vector<int64_t>(n_eq));
      for (size_t a = 0; a < n_eq; a++) {
        for (size_t b = 0; b < n_var; b++) {
          coefT[b][a] = coef0[a][b];
        }
      }
      auto nullspace = [&](ModMatrix const &M, size_t ncol) -> ModMatrix {
        RowSpaceModP rs = ComputeRowSpaceModP(M, p);
        std::vector<bool> is_piv(ncol, false);
        for (size_t pv : rs.pivots) {
          is_piv[pv] = true;
        }
        ModMatrix ker;
        for (size_t f = 0; f < ncol; f++) {
          if (is_piv[f]) {
            continue;
          }
          std::vector<int64_t> v(ncol, 0);
          v[f] = 1;
          for (size_t i = 0; i < rs.pivots.size(); i++) {
            v[rs.pivots[i]] = ResMod(-rs.basis[i][f], p);
          }
          ker.push_back(v);
        }
        return ker;
      };
      ModMatrix L = nullspace(coefT, n_eq);
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
      ModMatrix lmat(L.size(), std::vector<int64_t>(dim, 0));
      for (size_t i = 0; i < dim; i++) {
        std::vector<int64_t> c1(dim, 0), c2(dim, 0);
        c1[i] = 1;
        c2[i] = 2 % p;
        std::vector<int64_t> f1 = fval(c1), f2 = fval(c2);
        for (size_t k = 0; k < L.size(); k++) {
          int64_t q = ResMod(f2[k] - 2 * f1[k], p) * inv2 % p;
          lmat[k][i] = ResMod(f1[k] - q, p);
        }
      }
      ModMatrix W0 = nullspace(lmat, dim);
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
          std::optional<ModMatrix> lf = lift_path(layer_elt(w), j + 1, pj1);
          if (over) {
            return {};
          }
          if (!lf) {
            std::cerr << "MULT_FACTOR: an element of the layer does not lift\n";
            throw TerminalException{1};
          }
          gens.push_back(*lf);
        }
        ord.layer_exp += W0.size();
        done = true;
      }
    }
    if (done) {
      // found above
    } else if (all) {
      gens.insert(gens.end(), lifts.begin(), lifts.end());
      ord.layer_exp += sp->ker.size();
    } else {
      // the liftable subspace V contains the basis vectors that lift. The
      // failed ones f are taken in turn with a list T of failed vectors
      // independent modulo V: when some f + sum_t a_t t lifts it is a new
      // element of V, otherwise f is independent of T modulo V and joins
      // it. Then W = V + span(T) with T independent modulo V, so V is
      // spanned by the elements found, at the cost of at most |failed|
      // p^{codim V} lifts.
      size_t dim = sp->ker.size();
      ModMatrix liftable_coefs;
      std::vector<ModMatrix> liftable_elts = lifts;
      for (size_t i0 : lifted) {
        std::vector<int64_t> coef(dim, 0);
        coef[i0] = 1;
        liftable_coefs.push_back(coef);
      }
      ModMatrix reps;
      for (size_t i0 : failed) {
        size_t n_el = 1;
        for (size_t i = 0; i < reps.size(); i++) {
          n_el *= static_cast<size_t>(p);
          if (n_el > 100000) {
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
    std::cerr << "MULT_FACTOR: generators stop: layer too large, n_node=" << n_node << "\n";
#endif
            return {};
          }
        }
        bool in_v = false;
        for (size_t s2 = 0; s2 < n_el && !in_v; s2++) {
          std::vector<int64_t> coef(dim, 0);
          coef[i0] = 1;
          size_t tt = s2;
          for (size_t i = 0; i < reps.size(); i++) {
            int64_t a = static_cast<int64_t>(tt % static_cast<size_t>(p));
            tt /= static_cast<size_t>(p);
            for (size_t v = 0; v < dim; v++) {
              coef[v] = (coef[v] + a * reps[i][v]) % p;
            }
          }
          if (s2 == 0) {
            // f itself is known not to lift
            continue;
          }
          ModMatrix g = LiftAffineElement(Id, *sp, coef, p, pj);
          std::optional<ModMatrix> lf = lift_path(g, j + 1, pj * p);
          if (over) {
            return {};
          }
          if (lf) {
            liftable_coefs.push_back(coef);
            liftable_elts.push_back(*lf);
            in_v = true;
          }
        }
        if (!in_v) {
          std::vector<int64_t> coef(dim, 0);
          coef[i0] = 1;
          reps.push_back(coef);
        }
      }
      if (RankModPrime(liftable_coefs, p) != liftable_coefs.size()) {
        std::cerr << "MULT_FACTOR: the layer elements found are dependent\n";
        throw TerminalException{1};
      }
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
      std::cerr << "MULT_FACTOR: layer " << j << " codimension "
                << reps.size() << "\n";
#endif
      gens.insert(gens.end(), liftable_elts.begin(), liftable_elts.end());
      ord.layer_exp += liftable_coefs.size();
    }
    pj *= p;
  }
  if (order) {
    *order = ord;
  }
  return gens;
}

/*
  LocalOrthogonalImageGeneratorsCore in a Jordan basis when e >= 1
  (JordanBasisModPW, P of determinant prime to p): the search runs on A' =
  P^T A P with the levels of the columns (the images of the columns of
  level l stay in X^(l), and their pairings are known modulo higher powers
  of p), and the generators g' come back as P g' P^{-1}. For e = 0 the
  standard basis, all of level 0. Om (the
  action of a Z-basis of O_K, a center K) and Lam (a Z-basis of the
  integral commutant of O_K) are empty when End_G = Q; Bil are further
  integral bilinear forms preserved by the group (N -> P^T N P here),
  imposed modulo p on the columns; all the matrices act on the columns of
  coordinates.
 */
/*
  Lie: a Z-basis of the integral Lie lattice u (IntegralLieLatticeColumns)
  as residues modulo p^W, on the columns. When given, the layers G_j /
  G_{j+1} are those of its Cayley transforms (CayleyLayerGenerators), for j
  >= 1 at p odd and j >= 2 at p = 2, and the search only builds G mod p
  (and the first layer at p = 2); when not, all the layers are searched.
 */
inline std::optional<std::vector<ModMatrix>> LocalOrthogonalImageGeneratorsSearch(
    ModMatrix const &A, ModMatrix const &Be, std::vector<ModMatrix> const &Om,
    std::vector<ModMatrix> const &Lam, std::vector<ModMatrix> const &Bil, int64_t p,
    int c, int e, int layer_max, MultiplicityFactorBounds const &bounds,
    LocalImageOrder *order);

inline std::optional<std::vector<ModMatrix>>
LocalOrthogonalImageGenerators(ModMatrix const &A, ModMatrix const &Be,
                               std::vector<ModMatrix> const &Om,
                               std::vector<ModMatrix> const &Lam,
                               std::vector<ModMatrix> const &Bil,
                               std::optional<std::vector<ModMatrix>> const &Lie,
                               int64_t p, int c, int e,
                               MultiplicityFactorBounds const &bounds,
                               LocalImageOrder *order) {
  int layer_max = c;
  if (Lie) {
    layer_max = (p == 2) ? 2 : 1;
  }
  LocalImageOrder ord;
  std::optional<std::vector<ModMatrix>> gens = LocalOrthogonalImageGeneratorsSearch(
      A, Be, Om, Lam, Bil, p, c, e, layer_max, bounds, &ord);
  if (!gens) {
    return {};
  }
  if (Lie) {
    std::vector<ModMatrix> layer = CayleyLayerGenerators(*Lie, p, c, ord.layer_exp);
    gens->insert(gens->end(), layer.begin(), layer.end());
  }
  if (order) {
    *order = ord;
  }
  return gens;
}

inline std::optional<std::vector<ModMatrix>> LocalOrthogonalImageGeneratorsSearch(
    ModMatrix const &A, ModMatrix const &Be, std::vector<ModMatrix> const &Om,
    std::vector<ModMatrix> const &Lam, std::vector<ModMatrix> const &Bil, int64_t p,
    int c, int e, int layer_max, MultiplicityFactorBounds const &bounds,
    LocalImageOrder *order) {
  size_t m = A.size();
  if (e == 0) {
    return LocalOrthogonalImageGeneratorsCore(A, Be, Om, Lam, Bil, p, c, e,
                                              std::vector<int>(m, 0), layer_max,
                                              bounds, order);
  }
  // a Jordan basis and the levels of its columns
  auto [P, levels] = JordanBasisModPW(A, p, MultiplicityNewtonPrecision(p));
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
  std::cerr << "MULT_FACTOR: Jordan levels";
  for (int l : levels) {
    std::cerr << " " << l;
  }
  std::cerr << "\n";
#endif
  int W = MultiplicityNewtonPrecision(p);
  int64_t pW = MyPow<int64_t>(p, W);
  auto matmul = [](ModMatrix const &X, ModMatrix const &Y, int64_t M) -> ModMatrix {
    return ModMatrixProduct(X, Y, M);
  };
  auto transpose = [&](ModMatrix const &X) -> ModMatrix {
    ModMatrix Z(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Z[a][b] = X[b][a];
      }
    }
    return Z;
  };
  // P^{-1} modulo p^W: modulo p by Gauss, then X <- X (2 - P X)
  ModMatrix Pinv(m, std::vector<int64_t>(m));
  {
    ModMatrix aug(m, std::vector<int64_t>(2 * m, 0));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        aug[a][b] = ResMod(P[a][b], p);
      }
      aug[a][m + a] = 1;
    }
    RowSpaceModP rr = ComputeRowSpaceModP(aug, p);
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        Pinv[a][b] = rr.basis[a][m + b];
      }
    }
    int64_t M = p;
    while (M < pW) {
      M = (M > pW / M) ? pW : M * M;
      ModMatrix PX = matmul(P, Pinv, M);
      for (size_t a = 0; a < m; a++) {
        for (size_t b = 0; b < m; b++) {
          PX[a][b] = ResMod(((a == b) ? 2 : 0) - PX[a][b], M);
        }
      }
      Pinv = matmul(Pinv, PX, M);
    }
  }
  // A' = P^T A P and Be' = P^{-1} Be P^{-T} modulo p^W: the search uses A
  // modulo p^W only, and A may be given as residues modulo p^W, for which
  // the exact product does not fit in 64 bits
  ModMatrix A2(m, std::vector<int64_t>(m, 0));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      __int128 sacc = 0;
      for (size_t k = 0; k < m; k++) {
        if (P[k][a] != 0) {
          for (size_t l = 0; l < m; l++) {
            // reduced after each product: the entries of P (a Jordan basis)
            // are residues modulo p^W, and three of them overflow 128 bits
            __int128 t1 = static_cast<__int128>(P[k][a]) * ResMod(A[k][l], pW) % pW;
            sacc = (sacc + t1 * P[l][b]) % pW;
          }
        }
      }
      A2[a][b] = ResMod(static_cast<int64_t>(sacc), pW);
    }
  }
  ModMatrix Be2 = matmul(matmul(Pinv, Be, pW), transpose(Pinv), pW);
  // the levels: A'_ab divisible by p^{max(l_a, l_b)}
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      int64_t pl = 1;
      for (int i = 0; i < std::max(levels[a], levels[b]) && i < 20; i++) {
        pl *= p;
      }
      if (ResMod(A2[a][b], pW) % pl != 0) {
        std::cerr << "MULT_FACTOR: the Jordan basis does not have the levels "
                     "given, p=" << p << " a=" << a << " b=" << b << " levels";
        for (int l : levels) {
          std::cerr << " " << l;
        }
        std::cerr << "\nA mod p^3 =";
        for (size_t u = 0; u < m; u++) {
          for (size_t v = 0; v < m; v++) {
            std::cerr << " " << ResMod(A[u][v], p * p * p);
          }
          std::cerr << " |";
        }
        std::cerr << "\nA' mod p^3 =";
        for (size_t u = 0; u < m; u++) {
          for (size_t v = 0; v < m; v++) {
            std::cerr << " " << ResMod(A2[u][v], p * p * p);
          }
          std::cerr << " |";
        }
        std::cerr << "\n";
        throw TerminalException{1};
      }
    }
  }
#ifdef DEBUG_MULTIPLICITY_LOCAL_FACTOR
  std::cerr << "MULT_FACTOR: A' modulo p^3:";
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      std::cerr << " " << ResMod(A2[a][b], p * p * p);
    }
    std::cerr << " |";
  }
  std::cerr << "\n";
#endif
  // the action of O_K and the commutant in the new basis: P^{-1} X P
  auto conj = [&](ModMatrix const &X) -> ModMatrix {
    ModMatrix Xr = X;
    for (auto &row : Xr) {
      for (auto &x : row) {
        x = ResMod(x, pW);
      }
    }
    return matmul(matmul(Pinv, Xr, pW), P, pW);
  };
  std::vector<ModMatrix> Om2, Lam2;
  for (auto const &X : Om) {
    Om2.push_back(conj(X));
  }
  for (auto const &X : Lam) {
    Lam2.push_back(conj(X));
  }
  // the bilinear forms: x = P x', so N -> P^T N P
  std::vector<ModMatrix> Bil2;
  for (auto const &N : Bil) {
    ModMatrix Nr = N;
    for (auto &row : Nr) {
      for (auto &x : row) {
        x = ResMod(x, pW);
      }
    }
    ModMatrix PT(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        PT[a][b] = P[b][a];
      }
    }
    Bil2.push_back(matmul(matmul(PT, Nr, pW), P, pW));
  }
  std::optional<std::vector<ModMatrix>> gens =
      LocalOrthogonalImageGeneratorsCore(A2, Be2, Om2, Lam2, Bil2, p, c, e,
                                         levels, layer_max, bounds, order);
  if (!gens) {
    return {};
  }
  int64_t pc = MyPow<int64_t>(p, c);
  std::vector<ModMatrix> res;
  for (auto const &g : *gens) {
    res.push_back(matmul(matmul(P, g, pc), Pinv, pc));
  }
  return res;
}

/*
  The image modulo p^c of the stabilizer {g in Lambda (x) Z_p : g A g^T =
  A} (rows convention), Lambda the Z-span of the m x m integer matrices of
  basis: the centralizer of a group acting on a lattice of Gram matrix A,
  when Lambda is its commutant. The Newton step g -> g - (1/2) E g^{-T}
  A^{-1}, E = g A g^T - A, stays in the commutant when g does, since the
  group preserves A and E; so Lemma (multfactor) holds within Lambda, and
  the image is the set of the solutions modulo p^c in Lambda that lift to
  solutions modulo p^K, lifted level by level by Delta in Lambda / p
  Lambda. Nothing beyond the bounds.
 */
inline std::optional<std::vector<ModMatrix>>
LocalStabilizerImage(ModMatrix const &A_in, std::vector<ModMatrix> const &basis_in,
                     int64_t p, int c, int e,
                     MultiplicityFactorBounds const &bounds) {
  // A and the basis of Lambda, exact or modulo p^{K+1}: a solution is kept
  // by its coordinates in the basis, so that it stays in Lambda (the
  // entrywise reduction of an element of Lambda need not be in Lambda, and
  // the lifts g + p^j Delta, Delta in Lambda, would miss the isometries)
  size_t m = A_in.size();
  size_t dim = basis_in.size();
  int t = (p == 2) ? 1 : 0;
  int K = std::max(c + e + t, e + 2 * t + 1);
  /*
    E = g A g^T - A is only read modulo p^{j+1}, j < c <= K: A, the basis
    and g are taken modulo M = p^{K+1} (centered), so that their size does
    not depend on the exact data, and the guard below bounds the products
    m^2 |g|^2 |A| with |g| <= M / 2.
   */
  std::optional<int64_t> optM = LocalImageModulus(p, c, e);
  if (!optM) {
    return {};
  }
  int64_t const M = *optM;
  ModMatrix A(m, std::vector<int64_t>(m));
  for (size_t a = 0; a < m; a++) {
    for (size_t b = 0; b < m; b++) {
      A[a][b] = CenteredModP(A_in[a][b], M);
    }
  }
  std::vector<ModMatrix> basis = basis_in;
  for (auto &B : basis) {
    for (auto &row : B) {
      for (auto &x : row) {
        x = CenteredModP(x, M);
      }
    }
  }
  if (!ResiduesFitInt128(A, p, K + 1)) {
    return {};
  }
  using Coord = std::vector<int64_t>;
  auto matrix_of = [&](Coord const &x) -> std::vector<std::vector<__int128>> {
    std::vector<std::vector<__int128>> g(m, std::vector<__int128>(m, 0));
    for (size_t i = 0; i < dim; i++) {
      if (x[i] != 0) {
        for (size_t a = 0; a < m; a++) {
          for (size_t b = 0; b < m; b++) {
            g[a][b] = (g[a][b] + static_cast<__int128>(x[i]) * basis[i][a][b]) % M;
          }
        }
      }
    }
    for (auto &row : g) {
      for (auto &v : row) {
        if (v > M / 2) {
          v -= M;
        } else if (v < -(M / 2)) {
          v += M;
        }
      }
    }
    return g;
  };
  // E = g A g^T - A, exactly for these representatives (so modulo M)
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
    ModMatrix gp(m, std::vector<int64_t>(m));
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) {
        gp[a][b] = mod_p(g[a][b]);
      }
    }
    if (RankModPrime(gp, p) == m) {
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
    ModMatrix aug;
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
          row[i] = ResMod(v, p);
        }
        __int128 r = E[a][b] % pj1;
        if (r < 0) {
          r += pj1;
        }
        if (r % pj != 0) {
          std::cerr << "MULT_FACTOR: not a solution modulo p^j\n";
          throw TerminalException{1};
        }
        row[dim] = ResMod(-static_cast<int64_t>(r / pj), p);
        aug.push_back(row);
      }
    }
    RowSpaceModP rs = ComputeRowSpaceModP(aug, p);
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
        v[rs.pivots[i]] = ResMod(-rs.basis[i][f], p);
      }
      ker.push_back(v);
    }
    size_t n_sol = MyPow<size_t>(static_cast<size_t>(p), ker.size());
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
  std::vector<ModMatrix> res;
  for (auto const &x : level) {
    if (liftable(x, std::max(c, 1), pj)) {
      std::vector<std::vector<__int128>> g = matrix_of(x);
      ModMatrix gm(m, std::vector<int64_t>(m));
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

// clang-format off
#endif  // SRC_GENERA_MULTIPLICITYLOCALFACTOR_H_
// clang-format on
