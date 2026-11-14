// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GROUPNEIGHBORMODP_H_
#define SRC_GENERA_GROUPNEIGHBORMODP_H_

// clang-format off
#include "FiniteFieldProjective.h"
#include "PolynomialModQ.h"
#include "CyclotomicHermitianMass.h"
#include "MAT_MatrixModPrime.h"
#include <algorithm>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GROUP_NEIGHBOR_MOD_P
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GROUP_NEIGHBOR_MOD_P
#endif

/*
  The G-stable neighbors of a lattice L for a finite group G, at an odd
  prime p not dividing |G| det.

  L/pL is then a semisimple F_p[G]-module. A G-stable neighbor
  N = L_S + (1/p) S~ comes from a G-stable subspace S of L/pL that is
  totally isotropic, lifted to a G-stable S~ modulo p^2 that is totally
  isotropic modulo p^2, with L_S = {x in L : b(x, S~) in pZ}; then
  [L : L cap N] = [N : L cap N] = p^dim S. The irreducible S are taken: for
  an irreducible type tau with End_G(tau) = F_q, q = p^f, the
  tau-isotypic part is tau (x)_{F_q} U with U = Hom_G(tau, L/pL) of
  dimension m over F_q, and its irreducible submodules are the images of
  the nonzero phi in U, that is the points of P^{m-1}(F_q). For G cyclic
  and tau of dimension one these are the eigenvectors of the cyclic
  enumeration; for an eigenvalue outside F_p, tau has dimension f and the
  points are those of the Hermitian space over the residue field of an
  inert prime.

  The modules are computed by linear algebra modulo p:
  - an irreducible submodule is spun from a random vector and split until
    End_G(S) is a field, which an element whose minimal polynomial is
    irreducible of degree dim End_G(S) certifies (it then generates
    End_G(S); if End_G(S) is not a field no element does);
  - the isotypic part of the type of S0 is the sum of the images of
    Hom_G(S0, L/pL), and the next type is searched in a G-stable
    complement, obtained by averaging a projector over G.
  The lift: E~ the G-equivariant projector on S averaged over G modulo p^2
  and made idempotent by 3 E^2 - 2 E^3 gives a G-stable S~ free of rank
  dim S; the form b(s~_i, s~_j)/p = beta_ij modulo p is corrected by
  s~ -> s~ + p phi(s), with phi : S -> T G-equivariant, T a G-stable
  complement of S^perp, b(s', phi(s)) = -beta(s', s)/2, which stays
  G-stable and makes S~ isotropic modulo p^2.
 */

/*
  A G-submodule given by an echelon basis, with the matrices of the
  generators of G in that basis (rows: s_i g = sum_j act[g](i, j) s_j).
 */
struct GroupSubmodule {
  RowSpaceModP space;
  std::vector<ModMatrix> act;
};

inline GroupSubmodule MakeGroupSubmodule(RowSpaceModP const &space,
                                         std::vector<ModMatrix> const &gens,
                                         T_small_prime p) {
  GroupSubmodule sm{space, {}};
  int k = space.basis.rows();
  for (auto const &g : gens) {
    ModMatrix img = ModMatrixProduct(space.basis, g, p);
    // the coordinates of a vector of the row space are its pivot entries
    ModMatrix a(k, k);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        a(i, j) = img(i, static_cast<int>(space.pivots[j]));
      }
    }
    sm.act.push_back(a);
  }
  return sm;
}

// The G-submodule spanned by the orbit of the rows of V.
inline RowSpaceModP SpinModule(ModMatrix const &V,
                               std::vector<ModMatrix> const &gens,
                               T_small_prime p) {
  RowSpaceModP rs = ComputeRowSpaceModP(V, p);
  while (true) {
    ModMatrix rows = rs.basis;
    for (auto const &g : gens) {
      rows = Concatenate(rows, ModMatrixProduct(rs.basis, g, p));
    }
    RowSpaceModP rs2 = ComputeRowSpaceModP(rows, p);
    if (rs2.basis.rows() == rs.basis.rows()) {
      return rs;
    }
    rs = std::move(rs2);
  }
}

/*
  Hom_G(S, T) for G-submodules S, T: the k x l matrices X with
  act_S[g] X = X act_T[g], by the nullspace of the linear system.
 */
inline std::vector<ModMatrix> HomGroupModules(GroupSubmodule const &S,
                                              GroupSubmodule const &T,
                                              T_small_prime p) {
  int k = S.space.basis.rows();
  int l = T.space.basis.rows();
  int n_var = k * l;
  int n_eq = static_cast<int>(S.act.size()) * k * l;
  // the equations, as the columns of a matrix with n_var rows
  ModMatrix M = ZeroMatrix<T_residue>(n_var, n_eq);
  int e = 0;
  for (size_t ig = 0; ig < S.act.size(); ig++) {
    ModMatrix const &a = S.act[ig];
    ModMatrix const &b = T.act[ig];
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < l; j++) {
        // (a X - X b)(i, j) = sum_u a(i,u) X(u,j) - sum_v X(i,v) b(v,j)
        for (int u = 0; u < k; u++) {
          M(u * l + j, e) = ResAdd(M(u * l + j, e), a(i, u), p);
        }
        for (int v = 0; v < l; v++) {
          M(i * l + v, e) = ResSub(M(i * l + v, e), b(v, j), p);
        }
        e++;
      }
    }
  }
  ModMatrix ker = LeftNullspaceModP(M, p);
  std::vector<ModMatrix> res;
  for (int r = 0; r < ker.rows(); r++) {
    ModMatrix X(k, l);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < l; j++) {
        X(i, j) = ker(r, i * l + j);
      }
    }
    res.push_back(X);
  }
  return res;
}

// The minimal polynomial of a square matrix, by the first dependence of
// its powers.
inline PolyQ ModMinimalPolynomial(ModMatrix const &X, T_small_prime p) {
  int k = X.rows();
  ModMatrix flat(0, k * k);
  ModMatrix pw = ModMatrixIdentity(k);
  for (int d = 0; d <= k; d++) {
    ModVector v = FlattenModMatrix(pw);
    if (flat.rows() > 0) {
      std::optional<ModVector> sol = SolveRowsModP(flat, v, p);
      if (sol) {
        // X^d = sum_i sol_i X^i
        PolyQ mu(d + 1, 0);
        for (int i = 0; i < d; i++) {
          mu[i] = ResNeg((*sol)[i], p);
        }
        mu[d] = 1;
        return mu;
      }
    }
    AppendModRow(flat, v);
    pw = ModMatrixProduct(pw, X, p);
  }
  std::cerr << "GROUP_NEIGHBOR: no minimal polynomial found\n";
  throw TerminalException{1};
}

inline ModMatrix PolyOfMatrix(PolyQ const &f, ModMatrix const &X,
                              T_small_prime p) {
  int k = X.rows();
  ModMatrix res = ZeroMatrix<T_residue>(k, k);
  for (size_t i = f.size(); i-- > 0;) {
    res = ModMatrixProduct(res, X, p);
    for (int j = 0; j < k; j++) {
      res(j, j) = ResAdd(res(j, j), f[i] % p, p);
    }
  }
  return res;
}

/*
  An irreducible type: the submodule S0 (dimension k), a generator theta
  of the field End_G(S0) = F_q, q = p^f, with minimal polynomial mu, and
  an F_q-basis Phi_1, ..., Phi_m of Hom_G(S0, L/pL), each a k x n matrix
  (the images of the basis of S0), F_q acting by alpha . Phi = alpha(theta)
  Phi.
 */
struct IrreducibleType {
  GroupSubmodule S0;
  ModMatrix theta;
  PolyQ mu;
  int f;
  std::vector<ModMatrix> Phi;
  // the totally isotropic type: b vanishes on its isotypic part
  bool isotypic_isotropic;
};

/*
  An irreducible submodule of the G-submodule W (echelon basis), with the
  certificate of irreducibility: an element theta of End_G(S) whose
  minimal polynomial is irreducible of degree dim End_G(S).
 */
inline std::pair<GroupSubmodule, std::pair<ModMatrix, PolyQ>>
IrreducibleSubmodule(RowSpaceModP const &W, std::vector<ModMatrix> const &gens,
                     T_small_prime p, std::mt19937_64 &rng) {
  int n = W.basis.cols();
  std::uniform_int_distribution<T_residue> unif(0, p - 1);
  ModMatrix v(1, n);
  do {
    for (int j = 0; j < n; j++) {
      v(0, j) = 0;
    }
    for (int r = 0; r < W.basis.rows(); r++) {
      T_residue c = unif(rng);
      for (int j = 0; j < n; j++) {
        v(0, j) = ResAdd(v(0, j), ResMul(c, W.basis(r, j), p), p);
      }
    }
  } while (IsZeroMatrix(v));
  RowSpaceModP S = SpinModule(v, gens, p);
  for (int iter = 0; iter < 10000; iter++) {
    GroupSubmodule sm = MakeGroupSubmodule(S, gens, p);
    std::vector<ModMatrix> End = HomGroupModules(sm, sm, p);
    size_t d = End.size();
    int k = S.basis.rows();
    if (d == 1) {
      return {sm, {ModMatrixIdentity(k), PolyQ{p - 1, 1}}};
    }
    // a random element of End_G(S)
    ModMatrix x = ZeroMatrix<T_residue>(k, k);
    for (auto const &E : End) {
      T_residue c = unif(rng);
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          x(i, j) = ResAdd(x(i, j), ResMul(c, E(i, j), p), p);
        }
      }
    }
    PolyQ mu = ModMinimalPolynomial(x, p);
    bool irr = IsIrreducibleModP(mu, p);
    if (irr && mu.size() - 1 == d) {
      return {sm, {x, mu}};
    }
    if (irr) {
      // irreducible but of smaller degree: try another element
      continue;
    }
    // a proper factor h of mu: ker h(x) is a proper G-submodule
    PolyQ xpol{0, 1};
    PolyQ h;
    for (size_t e = 1; e <= mu.size() - 1; e++) {
      PolyQ y = xpol;
      for (size_t i = 0; i < e; i++) {
        y = PolyPowModQ(y, p, mu, p);
      }
      PolyQ g = PolyGcdQ(PolySubQ(y, xpol, p), mu, p);
      PolyTrim(g);
      if (g.size() > 1 && g.size() < mu.size()) {
        h = PolyMonicQ(g, p);
        break;
      }
    }
    if (h.empty()) {
      // mu squarefree with all its factors of one degree e, or mu = g^e:
      // the equal degree splitting of gcd(x^{p^e} - x, mu)
      for (size_t e = 1; e < mu.size() && h.empty(); e++) {
        PolyQ y = xpol;
        for (size_t i = 0; i < e; i++) {
          y = PolyPowModQ(y, p, mu, p);
        }
        PolyQ g = PolyGcdQ(PolySubQ(y, xpol, p), mu, p);
        PolyTrim(g);
        if (g.size() > 1) {
          std::vector<PolyQ> fac =
              EqualDegreeFactorization(PolyMonicQ(g, p), static_cast<int>(e),
                                       p);
          h = fac[0];
        }
      }
    }
    if (h.empty()) {
      continue;
    }
    ModMatrix hx = PolyOfMatrix(h, x, p);
    ModMatrix ker = LeftNullspaceModP(hx, p);
    if (ker.rows() == 0 || ker.rows() == k) {
      continue;
    }
    // back to the ambient coordinates
    S = ComputeRowSpaceModP(ModMatrixProduct(ker, S.basis, p), p);
  }
  std::cerr << "GROUP_NEIGHBOR: no irreducible submodule found\n";
  throw TerminalException{1};
}

/*
  The G-equivariant projector on the G-submodule U (echelon basis), on the
  row vectors, by averaging over the elements of G (p does not divide |G|).
  ListElt and ListEltInv: the elements and their inverses modulo p.
 */
inline ModMatrix EquivariantProjector(RowSpaceModP const &U,
                                      std::vector<ModMatrix> const &ListElt,
                                      std::vector<ModMatrix> const &ListEltInv,
                                      T_small_prime p) {
  int n = ListElt[0].rows();
  ModMatrix E0 = ZeroMatrix<T_residue>(n, n);
  for (size_t i = 0; i < U.pivots.size(); i++) {
    E0.row(static_cast<int>(U.pivots[i])) = U.basis.row(static_cast<int>(i));
  }
  ModMatrix E = ZeroMatrix<T_residue>(n, n);
  for (size_t ig = 0; ig < ListElt.size(); ig++) {
    ModMatrix X = ModMatrixProduct(ModMatrixProduct(ListEltInv[ig], E0, p),
                                   ListElt[ig], p);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E(i, j) = ResAdd(E(i, j), X(i, j), p);
      }
    }
  }
  T_residue inv = ResInv(ListElt.size() % p, p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      E(i, j) = ResMul(E(i, j), inv, p);
    }
  }
  return E;
}

// The row space of I - E.
inline RowSpaceModP ComplementOfProjector(ModMatrix const &E,
                                          T_small_prime p) {
  int n = E.rows();
  ModMatrix C = ModMatrixIdentity(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C(i, j) = ResSub(C(i, j), E(i, j), p);
    }
  }
  return ComputeRowSpaceModP(C, p);
}

/*
  The decomposition of F_p^n into isotypic parts: the irreducible types,
  with the F_q-bases of the Hom spaces. A = Gram matrix modulo p, for the
  test of total isotropy of an isotypic part.
 */
inline std::vector<IrreducibleType>
IrreducibleTypes(int n, std::vector<ModMatrix> const &gens,
                 std::vector<ModMatrix> const &ListElt,
                 std::vector<ModMatrix> const &ListEltInv, ModMatrix const &A,
                 T_small_prime p, std::mt19937_64 &rng) {
  std::vector<IrreducibleType> res;
  RowSpaceModP Full = ComputeRowSpaceModP(ModMatrixIdentity(n), p);
  GroupSubmodule FullMod = MakeGroupSubmodule(Full, gens, p);
  ModMatrix found(0, n);
  RowSpaceModP Rest = Full;
  while (Rest.basis.rows() > 0) {
    std::pair<GroupSubmodule, std::pair<ModMatrix, PolyQ>> irr =
        IrreducibleSubmodule(Rest, gens, p, rng);
    IrreducibleType tp;
    tp.S0 = irr.first;
    tp.theta = irr.second.first;
    tp.mu = irr.second.second;
    tp.f = static_cast<int>(tp.mu.size()) - 1;
    int k = tp.S0.space.basis.rows();
    // the k x n matrices of the maps (the full module has the standard
    // basis)
    std::vector<ModMatrix> HomMat = HomGroupModules(tp.S0, FullMod, p);
    // an F_q-basis
    std::vector<ModMatrix> theta_pow{ModMatrixIdentity(k)};
    for (int a = 1; a < tp.f; a++) {
      theta_pow.push_back(ModMatrixProduct(theta_pow.back(), tp.theta, p));
    }
    ModMatrix span(0, k * n);
    for (auto const &X : HomMat) {
      ModMatrix test = span;
      AppendModRow(test, FlattenModMatrix(X));
      if (RankModPrime(test, p) == static_cast<size_t>(span.rows())) {
        continue;
      }
      tp.Phi.push_back(X);
      for (auto const &tp_a : theta_pow) {
        AppendModRow(span, FlattenModMatrix(ModMatrixProduct(tp_a, X, p)));
      }
      span = ComputeRowSpaceModP(span, p).basis;
    }
    // the isotypic part
    ModMatrix iso_rows(0, n);
    for (auto const &X : tp.Phi) {
      for (auto const &tp_a : theta_pow) {
        iso_rows = Concatenate(iso_rows, ModMatrixProduct(tp_a, X, p));
      }
    }
    RowSpaceModP Iso = ComputeRowSpaceModP(iso_rows, p);
    // b on the isotypic part
    ModMatrix B = ModMatrixProduct(ModMatrixProduct(Iso.basis, A, p),
                                   ModMatrix(Iso.basis.transpose()), p);
    tp.isotypic_isotropic = IsZeroMatrix(B);
#ifdef DEBUG_GROUP_NEIGHBOR_MOD_P
    std::cerr << "GROUP_NEIGHBOR: p=" << p << " type dim=" << k
              << " f=" << tp.f << " m=" << tp.Phi.size()
              << " isotypic dim=" << Iso.basis.rows()
              << " isotropic=" << tp.isotypic_isotropic << "\n";
#endif
    if (Iso.basis.rows() != k * static_cast<int>(tp.Phi.size())) {
      std::cerr << "GROUP_NEIGHBOR: the isotypic part has dimension "
                << Iso.basis.rows() << " instead of " << k * tp.Phi.size()
                << "\n";
      throw TerminalException{1};
    }
    found = Concatenate(found, Iso.basis);
    res.push_back(tp);
    RowSpaceModP Found = ComputeRowSpaceModP(found, p);
    if (Found.basis.rows() == n) {
      break;
    }
    ModMatrix E = EquivariantProjector(Found, ListElt, ListEltInv, p);
    Rest = ComplementOfProjector(E, p);
  }
  return res;
}


inline std::vector<ModMatrix> ThetaPowers(IrreducibleType const &tp,
                                          T_small_prime p) {
  int k = tp.S0.space.basis.rows();
  std::vector<ModMatrix> res{ModMatrixIdentity(k)};
  for (int a = 1; a < tp.f; a++) {
    res.push_back(ModMatrixProduct(res.back(), tp.theta, p));
  }
  return res;
}

// sum_i alpha_i(theta) Phi_i, a k x n matrix whose rows span the submodule.
inline ModMatrix HomOfPoint(IrreducibleType const &tp,
                            std::vector<ModMatrix> const &theta_pow,
                            std::vector<ModVector> const &alpha,
                            T_small_prime p) {
  int k = tp.S0.space.basis.rows();
  int n = tp.Phi[0].cols();
  ModMatrix res = ZeroMatrix<T_residue>(k, n);
  for (size_t i = 0; i < tp.Phi.size(); i++) {
    for (int a = 0; a < tp.f; a++) {
      T_residue c = alpha[i][a];
      if (c != 0) {
        ModMatrix Y = ModMatrixProduct(theta_pow[a], tp.Phi[i], p);
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < n; v++) {
            res(u, v) = ResAdd(res(u, v), ResMul(c, Y(u, v), p), p);
          }
        }
      }
    }
  }
  return res;
}

/*
  The action of h, commuting with G (on the rows, modulo p), on
  Hom_G(S0, L/pL): Phi_i h = sum_j c_ij(theta) Phi_j, so that the point
  alpha goes to beta_j = sum_i alpha_i c_ij.
 */
inline std::vector<std::vector<ModVector>>
ActionOnHom(IrreducibleType const &tp, std::vector<ModMatrix> const &theta_pow,
            ModMatrix const &h, T_small_prime p) {
  size_t m = tp.Phi.size();
  ModMatrix basis;
  for (size_t j = 0; j < m; j++) {
    for (int a = 0; a < tp.f; a++) {
      AppendModRow(basis,
                   FlattenModMatrix(ModMatrixProduct(theta_pow[a], tp.Phi[j], p)));
    }
  }
  std::vector<std::vector<ModVector>> c(m, std::vector<ModVector>(m));
  for (size_t i = 0; i < m; i++) {
    std::optional<ModVector> sol = SolveRowsModP(
        basis, FlattenModMatrix(ModMatrixProduct(tp.Phi[i], h, p)), p);
    if (!sol) {
      std::cerr << "GROUP_NEIGHBOR: h does not preserve Hom_G(S0, L/pL)\n";
      throw TerminalException{1};
    }
    for (size_t j = 0; j < m; j++) {
      c[i][j] = ModVector(sol->begin() + j * tp.f, sol->begin() + (j + 1) * tp.f);
    }
  }
  return c;
}

inline std::vector<ModVector>
ApplyActionOnHom(FiniteFieldModP const &F,
                 std::vector<std::vector<ModVector>> const &c,
                 std::vector<ModVector> const &alpha) {
  size_t m = alpha.size();
  std::vector<ModVector> beta(m, ModVector(F.f, 0));
  for (size_t i = 0; i < m; i++) {
    if (!FqIsZero(alpha[i])) {
      for (size_t j = 0; j < m; j++) {
        beta[j] = FqAdd(F, beta[j], FqMul(F, alpha[i], c[i][j]));
      }
    }
  }
  return beta;
}

/*
  The lift of a totally isotropic G-submodule S of L/pL (rows of S_rows)
  to a G-stable submodule of (Z/p^2)^n, free of rank k, totally isotropic
  modulo p^2 (see the header comment). A: the Gram matrix, modulo p^2 for
  p odd and modulo 8 for p = 2; ListElt2 and ListEltInv2: the elements of G
  and their inverses modulo p^2. Nothing when the lift fails, which does not
  happen for p not dividing |G|.

  At p = 2 (|G| odd, the lattice even of odd determinant) S is totally
  singular for Q(x) = b(x, x) / 2 mod 2, the vectors of S~ must have
  b(s~, s~) in 8 Z and their pairings in 4 Z. On S, beta(x, y) = b(x~, y~)
  / 2 mod 2 and gamma(x) = b(x~, x~) / 4 mod 2 are G-invariant, gamma a
  quadratic form with polar form beta. Moving s~ to s~ + 2 phi(s), phi : S
  -> T G-equivariant, adds psi + psi^T to beta and psi(x, x) to gamma, psi(x,
  y) = b(x, phi(y)) mod 2. So psi must satisfy psi + psi^T = beta and
  psi(x, x) = gamma(x): the solutions are psi0 (upper triangular) plus the
  alternating forms, a G-stable affine space, and the sum of the g . psi0
  over G is a G-invariant solution, |G| being odd.
 */
inline std::optional<ModMatrix>
LiftIsotropicSubmodule(ModMatrix const &S_rows, ModMatrix const &A,
                       std::vector<ModMatrix> const &ListElt2,
                       std::vector<ModMatrix> const &ListEltInv2,
                       T_small_prime p) {
  T_residue p2 = p * p;
  int n = A.rows();
  int k = S_rows.rows();
  auto reduce = [&](ModMatrix X, T_residue mod) -> ModMatrix {
    for (int i = 0; i < X.rows(); i++) {
      for (int j = 0; j < X.cols(); j++) {
        X(i, j) %= mod;
      }
    }
    return X;
  };
  std::vector<ModMatrix> ListElt, ListEltInv;
  for (size_t i = 0; i < ListElt2.size(); i++) {
    ListElt.push_back(reduce(ListElt2[i], p));
    ListEltInv.push_back(reduce(ListEltInv2[i], p));
  }
  RowSpaceModP S = ComputeRowSpaceModP(S_rows, p);
  ModMatrix ES = EquivariantProjector(S, ListElt, ListEltInv, p);
  // the average modulo p^2 of the integer lift, then made idempotent
  ModMatrix E1 = ZeroMatrix<T_residue>(n, n);
  for (size_t ig = 0; ig < ListElt2.size(); ig++) {
    ModMatrix X = ModMatrixProduct(ModMatrixProduct(ListEltInv2[ig], ES, p2),
                                   ListElt2[ig], p2);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E1(i, j) = ResAdd(E1(i, j), X(i, j), p2);
      }
    }
  }
  T_residue inv = ResInv(ListElt2.size() % p2, p2);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      E1(i, j) = ResMul(E1(i, j), inv, p2);
    }
  }
  ModMatrix E1sq = ModMatrixProduct(E1, E1, p2);
  ModMatrix E1cu = ModMatrixProduct(E1sq, E1, p2);
  ModMatrix E2(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // 3 E^2 - 2 E^3
      E2(i, j) = ResSub(ResMul(T_residue(3), E1sq(i, j), p2),
                        ResMul(T_residue(2), E1cu(i, j), p2),
                        p2);
    }
  }
  ModMatrix St = ModMatrixProduct(S.basis, E2, p2);
  // the norms modulo 8 at p = 2: b(s~, s~) mod 8 depends on s~ mod 4 only,
  // the lattice being even
  T_residue pnorm = (p == 2) ? 8 : p2;
  // b(x, y) modulo mod for the rows x, y of matrices of residues
  auto bil = [&](ModMatrix const &X, int ix, ModMatrix const &Y, int iy,
                 T_residue mod) -> T_residue {
    T_residue sum = 0;
    for (int a = 0; a < n; a++) {
      T_residue x_a = X(ix, a) % mod;
      if (x_a != 0) {
        T_residue t = 0;
        for (int b = 0; b < n; b++) {
          t = ResAdd(t, ResMul(A(a, b) % mod, Y(iy, b) % mod, mod), mod);
        }
        sum = ResAdd(sum, ResMul(x_a, t, mod), mod);
      }
    }
    return sum;
  };
  // beta = b(s~_i, s~_j) / p modulo p
  ModMatrix beta(k, k);
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      T_residue v = bil(St, i, St, j, p2);
      if (v % p != 0) {
        std::cerr << "GROUP_NEIGHBOR: the submodule is not isotropic\n";
        throw TerminalException{1};
      }
      beta(i, j) = v / p;
    }
  }
  // p = 2: gamma = b(s~_i, s~_i) / 4 modulo 2, S being totally singular
  ModVector gamma(k, 0);
  if (p == 2) {
    for (int i = 0; i < k; i++) {
      T_residue v = bil(St, i, St, i, 8);
      if (v % 4 != 0) {
        std::cerr << "GROUP_NEIGHBOR: the submodule is not singular\n";
        throw TerminalException{1};
      }
      gamma[i] = v / 4;
    }
  }
  // T: a G-stable complement of S^perp, paired perfectly with S
  ModMatrix AS = ModMatrixProduct(reduce(A, p), ModMatrix(S.basis.transpose()), p);
  RowSpaceModP Sperp = ComputeRowSpaceModP(LeftNullspaceModP(AS, p), p);
  ModMatrix Eperp = EquivariantProjector(Sperp, ListElt, ListEltInv, p);
  RowSpaceModP T = ComplementOfProjector(Eperp, p);
  if (T.basis.rows() != k) {
    std::cerr << "GROUP_NEIGHBOR: the complement of S^perp has dimension "
              << T.basis.rows() << " instead of " << k << "\n";
    throw TerminalException{1};
  }
  // M(i, l) = b(s_i, tau_l), invertible
  ModMatrix M(k, k);
  for (int i = 0; i < k; i++) {
    for (int l = 0; l < k; l++) {
      M(i, l) = bil(S.basis, i, T.basis, l, p);
    }
  }
  /*
    psi(s_i, s_j) = b(s_i, phi(s_j)) mod p, phi : S -> T the correction:
    -beta / 2 for p odd; at p = 2 the G-invariant solution of psi + psi^T =
    beta, psi(s_i, s_i) = gamma_i (see above), with the action R_g of the
    elements on the basis of S.
   */
  ModMatrix psi = ZeroMatrix<T_residue>(k, k);
  if (p != 2) {
    T_residue half = ResInv(T_residue(2), p);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        psi(i, j) = ResNeg(ResMul(beta(i, j), half, p), p);
      }
    }
  } else {
    ModMatrix psi0 = ZeroMatrix<T_residue>(k, k);
    for (int i = 0; i < k; i++) {
      psi0(i, i) = gamma[i] % 2;
      for (int j = i + 1; j < k; j++) {
        psi0(i, j) = beta(i, j) % 2;
      }
    }
    for (auto const &g : ListElt) {
      ModMatrix Sg = ModMatrixProduct(S.basis, g, 2);
      ModMatrix R(k, k);
      for (int i = 0; i < k; i++) {
        ModVector row(n);
        for (int a = 0; a < n; a++) {
          row[a] = Sg(i, a);
        }
        std::optional<ModVector> r = SolveRowsModP(S.basis, row, 2);
        if (!r) {
          std::cerr << "GROUP_NEIGHBOR: S is not G-stable\n";
          throw TerminalException{1};
        }
        for (int l = 0; l < k; l++) {
          R(i, l) = (*r)[l];
        }
      }
      // psi += R psi0 R^T over F_2
      ModMatrix Rpsi = ModMatrixProduct(R, psi0, 2);
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          T_residue v = 0;
          for (int l = 0; l < k; l++) {
            v ^= (Rpsi(i, l) & R(j, l));
          }
          psi(i, j) ^= v;
        }
      }
    }
  }
  // t_j = sum_l x_l tau_l with sum_l M(i, l) x_l = psi(i, j)
  ModMatrix Mt(M.transpose());
  for (int j = 0; j < k; j++) {
    ModVector rhs(k);
    for (int i = 0; i < k; i++) {
      rhs[i] = psi(i, j);
    }
    // x Mt = rhs, i.e. sum_l x_l M(i, l) = rhs_i
    std::optional<ModVector> x = SolveRowsModP(Mt, rhs, p);
    if (!x) {
      std::cerr << "GROUP_NEIGHBOR: S and T are not in duality\n";
      throw TerminalException{1};
    }
    for (int l = 0; l < k; l++) {
      for (int c = 0; c < n; c++) {
        T_residue add = p * ResMul((*x)[l], T.basis(l, c), p);
        St(j, c) = ResAdd(St(j, c), add, p2);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    if (bil(St, i, St, i, pnorm) != 0) {
      return {};
    }
    for (int j = 0; j < k; j++) {
      if (bil(St, i, St, j, p2) != 0) {
        return {};
      }
    }
  }
#ifdef SANITY_CHECK_GROUP_NEIGHBOR_MOD_P
  // G-stable modulo p^2: every s~ g is in the span of the s~ modulo p^2,
  // tested on the reduction of the coordinates
  for (auto const &g : ListElt2) {
    ModMatrix img = reduce(ModMatrixProduct(St, g, p2), p);
    RowSpaceModP rs = ComputeRowSpaceModP(reduce(St, p), p);
    for (int r = 0; r < img.rows(); r++) {
      ModMatrix X = rs.basis;
      AppendModRow(X, FlattenModMatrix(img.row(r)));
      if (RankModPrime(X, p) != static_cast<size_t>(k)) {
        std::cerr << "GROUP_NEIGHBOR: the lift is not G-stable\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return St;
}

/*
  Generators of p N for the neighbor N = L_S + (1/p) S~: p^2 e_i, p times
  lifts of a basis of the kernel of x -> x A S~^T modulo p, and the rows of
  S~.
 */
inline ModMatrix GroupNeighborGenerators(ModMatrix const &St,
                                         ModMatrix const &A, T_small_prime p) {
  int n = A.rows();
  int k = St.rows();
  ModMatrix A_p(n, n), St_p(k, n);
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) {
      A_p(a, b) = A(a, b) % p;
    }
  }
  for (int j = 0; j < k; j++) {
    for (int b = 0; b < n; b++) {
      St_p(j, b) = St(j, b) % p;
    }
  }
  ModMatrix ker =
      LeftNullspaceModP(ModMatrixProduct(A_p, ModMatrix(St_p.transpose()), p), p);
  ModMatrix gens = ZeroMatrix<T_residue>(n + ker.rows() + k, n);
  for (int i = 0; i < n; i++) {
    gens(i, i) = p * p;
  }
  for (int r = 0; r < ker.rows(); r++) {
    for (int i = 0; i < n; i++) {
      gens(n + r, i) = p * ker(r, i);
    }
  }
  gens.bottomRows(k) = St;
  return gens;
}

// clang-format off
#endif  // SRC_GENERA_GROUPNEIGHBORMODP_H_
// clang-format on
