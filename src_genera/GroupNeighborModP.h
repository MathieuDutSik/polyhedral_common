// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GROUPNEIGHBORMODP_H_
#define SRC_GENERA_GROUPNEIGHBORMODP_H_

// clang-format off
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
                                         int64_t p) {
  GroupSubmodule sm{space, {}};
  for (auto const &g : gens) {
    ModMatrix img = ModMatrixProduct(space.basis, g, p);
    ModMatrix a;
    for (auto const &row : img) {
      a.push_back(CoordinatesInRowSpaceModP(space, row));
    }
    sm.act.push_back(a);
  }
  return sm;
}

// The G-submodule spanned by the orbit of the rows of V.
inline RowSpaceModP SpinModule(ModMatrix const &V, std::vector<ModMatrix> const &gens,
                              int64_t p) {
  RowSpaceModP rs = ComputeRowSpaceModP(V, p);
  while (true) {
    ModMatrix rows = rs.basis;
    for (auto const &g : gens) {
      ModMatrix img = ModMatrixProduct(rs.basis, g, p);
      rows.insert(rows.end(), img.begin(), img.end());
    }
    RowSpaceModP rs2 = ComputeRowSpaceModP(rows, p);
    if (rs2.basis.size() == rs.basis.size()) {
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
                                           GroupSubmodule const &T, int64_t p) {
  size_t k = S.space.basis.size();
  size_t l = T.space.basis.size();
  size_t n_var = k * l;
  // the equations, as columns of a matrix with n_var rows
  ModMatrix M(n_var);
  for (size_t ig = 0; ig < S.act.size(); ig++) {
    ModMatrix const &a = S.act[ig];
    ModMatrix const &b = T.act[ig];
    for (size_t i = 0; i < k; i++) {
      for (size_t j = 0; j < l; j++) {
        // (a X - X b)(i, j) = sum_u a(i,u) X(u,j) - sum_v X(i,v) b(v,j)
        std::vector<int64_t> col(n_var, 0);
        for (size_t u = 0; u < k; u++) {
          col[u * l + j] = ResMod(col[u * l + j] + a[i][u], p);
        }
        for (size_t v = 0; v < l; v++) {
          col[i * l + v] = ResMod(col[i * l + v] - b[v][j], p);
        }
        for (size_t w = 0; w < n_var; w++) {
          M[w].push_back(col[w]);
        }
      }
    }
  }
  ModMatrix ker = LeftNullspaceModP(M, p);
  std::vector<ModMatrix> res;
  for (auto const &x : ker) {
    ModMatrix X(k, std::vector<int64_t>(l));
    for (size_t i = 0; i < k; i++) {
      for (size_t j = 0; j < l; j++) {
        X[i][j] = x[i * l + j];
      }
    }
    res.push_back(X);
  }
  return res;
}

// The minimal polynomial of a square matrix, by the first dependence of
// its powers.
inline PolyQ ModMinimalPolynomial(ModMatrix const &X, int64_t p) {
  size_t k = X.size();
  ModMatrix flat;
  ModMatrix pw = ModMatrixIdentity(k);
  for (size_t d = 0; d <= k; d++) {
    std::vector<int64_t> v;
    for (auto const &row : pw) {
      v.insert(v.end(), row.begin(), row.end());
    }
    if (!flat.empty()) {
      std::optional<std::vector<int64_t>> sol = SolveRowsModP(flat, v, p);
      if (sol) {
        // X^d = sum_i sol_i X^i
        PolyQ mu(d + 1, 0);
        for (size_t i = 0; i < d; i++) {
          mu[i] = ResMod(-(*sol)[i], p);
        }
        mu[d] = 1;
        return mu;
      }
    }
    flat.push_back(v);
    pw = ModMatrixProduct(pw, X, p);
  }
  std::cerr << "GROUP_NEIGHBOR: no minimal polynomial found\n";
  throw TerminalException{1};
}

// Rabin's test of irreducibility over F_p.
inline bool IsIrreducibleModP(PolyQ const &f, int64_t p) {
  int d = static_cast<int>(f.size()) - 1;
  if (d <= 0) {
    return false;
  }
  if (d == 1) {
    return true;
  }
  PolyQ x{0, 1};
  auto frob_power = [&](int e) -> PolyQ {
    PolyQ y = x;
    for (int i = 0; i < e; i++) {
      y = PolyPowModQ(y, static_cast<uint64_t>(p), f, p);
    }
    return y;
  };
  // f divides x^{p^d} - x
  PolyQ y = PolySubQ(frob_power(d), x, p);
  PolyTrim(y);
  if (!y.empty() && !(y.size() == 1 && y[0] == 0)) {
    std::pair<PolyQ, PolyQ> qr = PolyDivModQ(y, f, p);
    PolyQ rem = qr.second;
    PolyTrim(rem);
    if (!rem.empty() && !(rem.size() == 1 && rem[0] == 0)) {
      return false;
    }
  }
  // gcd(x^{p^{d/r}} - x, f) = 1 for the primes r dividing d
  for (int r = 2; r <= d; r++) {
    if (d % r != 0) {
      continue;
    }
    bool is_prime = true;
    for (int s = 2; s * s <= r; s++) {
      if (r % s == 0) {
        is_prime = false;
      }
    }
    if (!is_prime) {
      continue;
    }
    PolyQ z = PolySubQ(frob_power(d / r), x, p);
    PolyQ g = PolyGcdQ(z, f, p);
    PolyTrim(g);
    if (g.size() > 1) {
      return false;
    }
  }
  return true;
}

inline ModMatrix PolyOfMatrix(PolyQ const &f, ModMatrix const &X, int64_t p) {
  size_t k = X.size();
  ModMatrix res(k, std::vector<int64_t>(k, 0));
  for (size_t i = f.size(); i-- > 0;) {
    res = ModMatrixProduct(res, X, p);
    for (size_t j = 0; j < k; j++) {
      res[j][j] = ResMod(res[j][j] + f[i], p);
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
                     int64_t p, std::mt19937_64 &rng) {
  size_t n = W.basis[0].size();
  std::uniform_int_distribution<int64_t> unif(0, p - 1);
  auto random_combination = [&](ModMatrix const &B) -> std::vector<int64_t> {
    std::vector<int64_t> v(B.empty() ? n : B[0].size(), 0);
    for (auto const &row : B) {
      int64_t c = unif(rng);
      for (size_t j = 0; j < v.size(); j++) {
        v[j] = ResMod(v[j] + ResMul(c, row[j], p), p);
      }
    }
    return v;
  };
  std::vector<int64_t> v;
  do {
    v = random_combination(W.basis);
  } while (std::all_of(v.begin(), v.end(), [](int64_t x) { return x == 0; }));
  RowSpaceModP S = SpinModule({v}, gens, p);
  for (int iter = 0; iter < 10000; iter++) {
    GroupSubmodule sm = MakeGroupSubmodule(S, gens, p);
    std::vector<ModMatrix> End = HomGroupModules(sm, sm, p);
    size_t d = End.size();
    if (d == 1) {
      size_t k = S.basis.size();
      return {sm, {ModMatrixIdentity(k), PolyQ{ResMod(-1, p), 1}}};
    }
    // a random element of End_G(S)
    size_t k = S.basis.size();
    ModMatrix x(k, std::vector<int64_t>(k, 0));
    for (auto const &E : End) {
      int64_t c = unif(rng);
      for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < k; j++) {
          x[i][j] = ResMod(x[i][j] + ResMul(c, E[i][j], p), p);
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
        y = PolyPowModQ(y, static_cast<uint64_t>(p), mu, p);
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
          y = PolyPowModQ(y, static_cast<uint64_t>(p), mu, p);
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
    if (ker.empty() || ker.size() == k) {
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
                                   int64_t p) {
  size_t n = ListElt[0].size();
  ModMatrix E0(n, std::vector<int64_t>(n, 0));
  for (size_t i = 0; i < U.pivots.size(); i++) {
    E0[U.pivots[i]] = U.basis[i];
  }
  ModMatrix E(n, std::vector<int64_t>(n, 0));
  for (size_t ig = 0; ig < ListElt.size(); ig++) {
    ModMatrix X = ModMatrixProduct(ModMatrixProduct(ListEltInv[ig], E0, p), ListElt[ig], p);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        E[i][j] = (E[i][j] + X[i][j]) % p;
      }
    }
  }
  int64_t inv = ResInv(ResMod(static_cast<int64_t>(ListElt.size()), p), p);
  for (auto &row : E) {
    for (auto &x : row) {
      x = ResMul(x, inv, p);
    }
  }
  return E;
}

// The row space of I - E.
inline RowSpaceModP ComplementOfProjector(ModMatrix const &E, int64_t p) {
  size_t n = E.size();
  ModMatrix C = ModMatrixIdentity(n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      C[i][j] = ResMod(C[i][j] - E[i][j], p);
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
IrreducibleTypes(size_t n, std::vector<ModMatrix> const &gens,
                 std::vector<ModMatrix> const &ListElt,
                 std::vector<ModMatrix> const &ListEltInv, ModMatrix const &A,
                 int64_t p, std::mt19937_64 &rng) {
  std::vector<IrreducibleType> res;
  RowSpaceModP Full = ComputeRowSpaceModP(ModMatrixIdentity(n), p);
  GroupSubmodule FullMod = MakeGroupSubmodule(Full, gens, p);
  ModMatrix found;
  RowSpaceModP Rest = Full;
  while (!Rest.basis.empty()) {
    std::pair<GroupSubmodule, std::pair<ModMatrix, PolyQ>> irr =
        IrreducibleSubmodule(Rest, gens, p, rng);
    IrreducibleType tp;
    tp.S0 = irr.first;
    tp.theta = irr.second.first;
    tp.mu = irr.second.second;
    tp.f = static_cast<int>(tp.mu.size()) - 1;
    size_t k = tp.S0.space.basis.size();
    std::vector<ModMatrix> Hom = HomGroupModules(tp.S0, FullMod, p);
    // the k x n matrices of the maps (the full module has the standard
    // basis)
    std::vector<ModMatrix> HomMat = Hom;
    // an F_q-basis
    std::vector<ModMatrix> theta_pow{ModMatrixIdentity(k)};
    for (int a = 1; a < tp.f; a++) {
      theta_pow.push_back(ModMatrixProduct(theta_pow.back(), tp.theta, p));
    }
    auto flatten = [&](ModMatrix const &X) -> std::vector<int64_t> {
      std::vector<int64_t> v;
      for (auto const &row : X) {
        v.insert(v.end(), row.begin(), row.end());
      }
      return v;
    };
    ModMatrix span;
    for (auto const &X : HomMat) {
      ModMatrix test = span;
      test.push_back(flatten(X));
      if (RankModPrime(test, p) == span.size()) {
        continue;
      }
      tp.Phi.push_back(X);
      for (auto const &tp_a : theta_pow) {
        span.push_back(flatten(ModMatrixProduct(tp_a, X, p)));
      }
      span = ComputeRowSpaceModP(span, p).basis;
    }
    // the isotypic part
    ModMatrix iso_rows;
    for (auto const &X : tp.Phi) {
      for (auto const &tp_a : theta_pow) {
        ModMatrix Y = ModMatrixProduct(tp_a, X, p);
        iso_rows.insert(iso_rows.end(), Y.begin(), Y.end());
      }
    }
    RowSpaceModP Iso = ComputeRowSpaceModP(iso_rows, p);
    // b on the isotypic part
    ModMatrix B = ModMatrixProduct(ModMatrixProduct(Iso.basis, A, p),
                         [&] {
                           ModMatrix T(n, std::vector<int64_t>(Iso.basis.size()));
                           for (size_t i = 0; i < n; i++) {
                             for (size_t j = 0; j < Iso.basis.size(); j++) {
                               T[i][j] = Iso.basis[j][i];
                             }
                           }
                           return T;
                         }(),
                         p);
    tp.isotypic_isotropic = true;
    for (auto const &row : B) {
      for (auto x : row) {
        if (x != 0) {
          tp.isotypic_isotropic = false;
        }
      }
    }
#ifdef DEBUG_GROUP_NEIGHBOR_MOD_P
    std::cerr << "GROUP_NEIGHBOR: p=" << p << " type dim=" << k
              << " f=" << tp.f << " m=" << tp.Phi.size()
              << " isotypic dim=" << Iso.basis.size()
              << " isotropic=" << tp.isotypic_isotropic << "\n";
#endif
    if (Iso.basis.size() != k * tp.Phi.size()) {
      std::cerr << "GROUP_NEIGHBOR: the isotypic part has dimension "
                << Iso.basis.size() << " instead of " << k * tp.Phi.size()
                << "\n";
      throw TerminalException{1};
    }
    found.insert(found.end(), Iso.basis.begin(), Iso.basis.end());
    res.push_back(tp);
    RowSpaceModP Found = ComputeRowSpaceModP(found, p);
    if (Found.basis.size() == n) {
      break;
    }
    ModMatrix E = EquivariantProjector(Found, ListElt, ListEltInv, p);
    Rest = ComplementOfProjector(E, p);
  }
  return res;
}

/*
  F_q = F_p[theta] / mu: elements as coefficient vectors of length f, and
  an encoding as integers in [0, q).
 */
struct FiniteFieldModP {
  int64_t p;
  int f;
  PolyQ mu;
  int64_t q;
};

inline std::vector<int64_t> FqMul(FiniteFieldModP const &F,
                                  std::vector<int64_t> const &a,
                                  std::vector<int64_t> const &b) {
  PolyQ prod = PolyMulQ(a, b, F.p);
  std::pair<PolyQ, PolyQ> qr = PolyDivModQ(prod, F.mu, F.p);
  std::vector<int64_t> r = qr.second;
  r.resize(F.f, 0);
  return r;
}

inline std::vector<int64_t> FqAdd(FiniteFieldModP const &F,
                                  std::vector<int64_t> const &a,
                                  std::vector<int64_t> const &b) {
  std::vector<int64_t> r(F.f);
  for (int i = 0; i < F.f; i++) {
    r[i] = (a[i] + b[i]) % F.p;
  }
  return r;
}

inline bool FqIsZero(std::vector<int64_t> const &a) {
  return std::all_of(a.begin(), a.end(), [](int64_t x) { return x == 0; });
}

inline std::vector<int64_t> FqInverse(FiniteFieldModP const &F,
                                      std::vector<int64_t> const &a) {
  PolyQ a2 = a;
  PolyTrim(a2);
  PolyQ inv = PolyInverseModQ(a2, F.mu, F.p);
  inv.resize(F.f, 0);
  return inv;
}

inline int64_t FqEncode(FiniteFieldModP const &F,
                        std::vector<int64_t> const &a) {
  int64_t c = 0;
  for (int i = F.f; i-- > 0;) {
    c = c * F.p + a[i];
  }
  return c;
}

inline std::vector<int64_t> FqDecode(FiniteFieldModP const &F, int64_t c) {
  std::vector<int64_t> a(F.f);
  for (int i = 0; i < F.f; i++) {
    a[i] = c % F.p;
    c /= F.p;
  }
  return a;
}

/*
  The points of P^{m-1}(F_q), normalized with the first nonzero coordinate
  equal to 1, indexed by the position i of that coordinate (q^{m-1-i}
  points each) and the remaining coordinates in base q.
 */
inline size_t NumberProjectivePointsFq(int m, int64_t q) {
  // the points are indexed by size_t: a count beyond it is refused, as in
  // NumberProjectiveLines, rather than wrapped around
  size_t const max_val = std::numeric_limits<size_t>::max();
  size_t qs = static_cast<size_t>(q);
  size_t tot = 0, pw = 1;
  for (int i = 0; i < m; i++) {
    if (tot > max_val - pw) {
      std::cerr << "GENUS: |P^" << (m - 1) << "(F_" << q
                << ")| overflows size_t\n";
      throw TerminalException{1};
    }
    tot += pw;
    if (i + 1 < m) {
      if (pw > max_val / qs) {
        std::cerr << "GENUS: " << q << "^" << (i + 1) << " overflows size_t\n";
        throw TerminalException{1};
      }
      pw *= qs;
    }
  }
  return tot;
}

inline std::vector<std::vector<int64_t>>
ProjectivePointFqOfIndex(FiniteFieldModP const &F, int m, size_t idx) {
  // the blocks in the order i = m-1 (one point) up to i = 0
  size_t q = static_cast<size_t>(F.q);
  std::vector<size_t> block(m);
  size_t pw = 1;
  for (int i = m - 1; i >= 0; i--) {
    block[i] = pw;
    pw *= q;
  }
  int lead = m - 1;
  for (int i = m - 1; i >= 0; i--) {
    if (idx < block[i]) {
      lead = i;
      break;
    }
    idx -= block[i];
  }
  std::vector<std::vector<int64_t>> x(m, std::vector<int64_t>(F.f, 0));
  x[lead][0] = 1;
  for (int j = m - 1; j > lead; j--) {
    x[j] = FqDecode(F, static_cast<int64_t>(idx % q));
    idx /= q;
  }
  return x;
}

inline size_t IndexOfProjectivePointFq(FiniteFieldModP const &F, int m,
                                       std::vector<std::vector<int64_t>> x) {
  size_t q = static_cast<size_t>(F.q);
  int lead = 0;
  while (lead < m && FqIsZero(x[lead])) {
    lead++;
  }
  std::vector<int64_t> inv = FqInverse(F, x[lead]);
  for (int j = lead; j < m; j++) {
    x[j] = FqMul(F, x[j], inv);
  }
  size_t offset = 0, pw = 1;
  for (int i = m - 1; i > lead; i--) {
    offset += pw;
    pw *= q;
  }
  size_t idx = 0;
  for (int j = lead + 1; j < m; j++) {
    idx = idx * q + static_cast<size_t>(FqEncode(F, x[j]));
  }
  return offset + idx;
}

/*
  F_q for a small q as tables on the codes of FqEncode, and the projective
  points and the action on Hom of the orbit enumeration in those codes: the
  same indices as ProjectivePointFqOfIndex and IndexOfProjectivePointFq,
  without the polynomial arithmetic and the allocations of every product
  (the orbits of P^11(F_4), 5.6 10^6 points, under the centralizer of
  (A2^12, w) took about ten minutes per class).
 */
struct FqTables {
  int64_t q;
  std::vector<int> add;
  std::vector<int> mul;
  std::vector<int> inv;
};

inline size_t const fq_tables_max_q = 4096;

inline FqTables BuildFqTables(FiniteFieldModP const &F) {
  int64_t q = F.q;
  FqTables tb{q, std::vector<int>(q * q), std::vector<int>(q * q),
              std::vector<int>(q, 0)};
  for (int64_t a = 0; a < q; a++) {
    std::vector<int64_t> fa = FqDecode(F, a);
    for (int64_t b = 0; b < q; b++) {
      std::vector<int64_t> fb = FqDecode(F, b);
      tb.add[a * q + b] = static_cast<int>(FqEncode(F, FqAdd(F, fa, fb)));
      int c = static_cast<int>(FqEncode(F, FqMul(F, fa, fb)));
      tb.mul[a * q + b] = c;
      if (c == 1) {
        tb.inv[a] = static_cast<int>(b);
      }
    }
  }
  return tb;
}

// The point of index idx as codes, x of size m.
inline void ProjectivePointCodesOfIndex(FqTables const &tb, int m, size_t idx,
                                        std::vector<int> &x) {
  size_t q = static_cast<size_t>(tb.q);
  size_t pw = 1;
  // the block of the lead i holds q^{m-1-i} points, in the order i = m-1
  // (one point) up to i = 0
  int lead = m - 1;
  for (int i = m - 1; i >= 0; i--) {
    if (idx < pw) {
      lead = i;
      break;
    }
    idx -= pw;
    pw *= q;
  }
  for (int j = 0; j < lead; j++) {
    x[j] = 0;
  }
  x[lead] = 1;
  for (int j = m - 1; j > lead; j--) {
    x[j] = static_cast<int>(idx % q);
    idx /= q;
  }
}

// The index of the point of codes x, which is normalized in place.
inline size_t IndexOfProjectivePointCodes(FqTables const &tb, int m,
                                          std::vector<int> &x) {
  size_t q = static_cast<size_t>(tb.q);
  int lead = 0;
  while (lead < m && x[lead] == 0) {
    lead++;
  }
  int inv = tb.inv[x[lead]];
  for (int j = lead; j < m; j++) {
    x[j] = tb.mul[x[j] * tb.q + inv];
  }
  size_t offset = 0, pw = 1;
  for (int i = m - 1; i > lead; i--) {
    offset += pw;
    pw *= q;
  }
  size_t idx = 0;
  for (int j = lead + 1; j < m; j++) {
    idx = idx * q + static_cast<size_t>(x[j]);
  }
  return offset + idx;
}

// The action of ApplyActionOnHom in codes: cc[i * m + j] the code of c[i][j].
inline std::vector<int>
ActionOnHomCodes(FiniteFieldModP const &F,
                 std::vector<std::vector<std::vector<int64_t>>> const &c) {
  size_t m = c.size();
  std::vector<int> cc(m * m);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < m; j++) {
      cc[i * m + j] = static_cast<int>(FqEncode(F, c[i][j]));
    }
  }
  return cc;
}

inline void ApplyActionOnHomCodes(FqTables const &tb, int m,
                                  std::vector<int> const &cc,
                                  std::vector<int> const &alpha,
                                  std::vector<int> &beta) {
  for (int j = 0; j < m; j++) {
    beta[j] = 0;
  }
  for (int i = 0; i < m; i++) {
    if (alpha[i] == 0) {
      continue;
    }
    int const *row = cc.data() + i * m;
    int const *mul_a = tb.mul.data() + alpha[i] * tb.q;
    for (int j = 0; j < m; j++) {
      beta[j] = tb.add[beta[j] * tb.q + mul_a[row[j]]];
    }
  }
}

inline std::vector<int64_t> FlattenModMat(ModMatrix const &X) {
  std::vector<int64_t> v;
  for (auto const &row : X) {
    v.insert(v.end(), row.begin(), row.end());
  }
  return v;
}

inline std::vector<ModMatrix> ThetaPowers(IrreducibleType const &tp, int64_t p) {
  size_t k = tp.S0.space.basis.size();
  std::vector<ModMatrix> res{ModMatrixIdentity(k)};
  for (int a = 1; a < tp.f; a++) {
    res.push_back(ModMatrixProduct(res.back(), tp.theta, p));
  }
  return res;
}

// sum_i alpha_i(theta) Phi_i, a k x n matrix whose rows span the submodule.
inline ModMatrix HomOfPoint(IrreducibleType const &tp,
                         std::vector<ModMatrix> const &theta_pow,
                         std::vector<std::vector<int64_t>> const &alpha,
                         int64_t p) {
  size_t k = tp.S0.space.basis.size();
  size_t n = tp.Phi[0][0].size();
  ModMatrix res(k, std::vector<int64_t>(n, 0));
  for (size_t i = 0; i < tp.Phi.size(); i++) {
    for (int a = 0; a < tp.f; a++) {
      int64_t c = alpha[i][a];
      if (c == 0) {
        continue;
      }
      ModMatrix Y = ModMatrixProduct(theta_pow[a], tp.Phi[i], p);
      for (size_t u = 0; u < k; u++) {
        for (size_t v = 0; v < n; v++) {
          res[u][v] = ResMod(res[u][v] + ResMul(c, Y[u][v], p), p);
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
inline std::vector<std::vector<std::vector<int64_t>>>
ActionOnHom(IrreducibleType const &tp, std::vector<ModMatrix> const &theta_pow,
            ModMatrix const &h, int64_t p) {
  size_t m = tp.Phi.size();
  ModMatrix basis;
  for (size_t j = 0; j < m; j++) {
    for (int a = 0; a < tp.f; a++) {
      basis.push_back(FlattenModMat(ModMatrixProduct(theta_pow[a], tp.Phi[j], p)));
    }
  }
  std::vector<std::vector<std::vector<int64_t>>> c(
      m, std::vector<std::vector<int64_t>>(m));
  for (size_t i = 0; i < m; i++) {
    std::optional<std::vector<int64_t>> sol = SolveRowsModP(
        basis, FlattenModMat(ModMatrixProduct(tp.Phi[i], h, p)), p);
    if (!sol) {
      std::cerr << "GROUP_NEIGHBOR: h does not preserve Hom_G(S0, L/pL)\n";
      throw TerminalException{1};
    }
    for (size_t j = 0; j < m; j++) {
      c[i][j] = std::vector<int64_t>(sol->begin() + j * tp.f,
                                     sol->begin() + (j + 1) * tp.f);
    }
  }
  return c;
}

inline std::vector<std::vector<int64_t>>
ApplyActionOnHom(FiniteFieldModP const &F,
                 std::vector<std::vector<std::vector<int64_t>>> const &c,
                 std::vector<std::vector<int64_t>> const &alpha) {
  size_t m = alpha.size();
  std::vector<std::vector<int64_t>> beta(m, std::vector<int64_t>(F.f, 0));
  for (size_t i = 0; i < m; i++) {
    if (FqIsZero(alpha[i])) {
      continue;
    }
    for (size_t j = 0; j < m; j++) {
      beta[j] = FqAdd(F, beta[j], FqMul(F, alpha[i], c[i][j]));
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
                       std::vector<ModMatrix> const &ListEltInv2, int64_t p) {
  int64_t p2 = p * p;
  size_t n = A.size();
  size_t k = S_rows.size();
  auto reduce = [&](ModMatrix X, int64_t mod) -> ModMatrix {
    for (auto &row : X) {
      for (auto &x : row) {
        x = ResMod(x, mod);
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
  ModMatrix E1(n, std::vector<int64_t>(n, 0));
  for (size_t ig = 0; ig < ListElt2.size(); ig++) {
    ModMatrix X =
        ModMatrixProduct(ModMatrixProduct(ListEltInv2[ig], ES, p2), ListElt2[ig], p2);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        E1[i][j] = (E1[i][j] + X[i][j]) % p2;
      }
    }
  }
  int64_t inv = ResInv(ResMod(static_cast<int64_t>(ListElt2.size()), p2), p2);
  for (auto &row : E1) {
    for (auto &x : row) {
      x = ResMul(x, inv, p2);
    }
  }
  ModMatrix E1sq = ModMatrixProduct(E1, E1, p2);
  ModMatrix E1cu = ModMatrixProduct(E1sq, E1, p2);
  ModMatrix E2(n, std::vector<int64_t>(n));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      E2[i][j] = ResMod(3 * E1sq[i][j] - 2 * E1cu[i][j], p2);
    }
  }
  ModMatrix St = ModMatrixProduct(S.basis, E2, p2);
  // the norms modulo 8 at p = 2: b(s~, s~) mod 8 depends on s~ mod 4 only,
  // the lattice being even
  int64_t pnorm = (p == 2) ? 8 : p2;
  auto bil = [&](std::vector<int64_t> const &x, std::vector<int64_t> const &y,
                 int64_t mod) -> int64_t {
    int64_t sum = 0;
    for (size_t a = 0; a < n; a++) {
      if (x[a] == 0) {
        continue;
      }
      int64_t t = 0;
      for (size_t b = 0; b < n; b++) {
        t = (t + ResMul(ResMod(A[a][b], mod), y[b], mod)) % mod;
      }
      sum = (sum + ResMul(x[a], t, mod)) % mod;
    }
    return sum;
  };
  // beta = b(s~_i, s~_j) / p modulo p
  ModMatrix beta(k, std::vector<int64_t>(k));
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < k; j++) {
      int64_t v = bil(St[i], St[j], p2);
      if (v % p != 0) {
        std::cerr << "GROUP_NEIGHBOR: the submodule is not isotropic\n";
        throw TerminalException{1};
      }
      beta[i][j] = v / p;
    }
  }
  // p = 2: gamma = b(s~_i, s~_i) / 4 modulo 2, S being totally singular
  std::vector<int64_t> gamma(k, 0);
  if (p == 2) {
    for (size_t i = 0; i < k; i++) {
      int64_t v = bil(St[i], St[i], 8);
      if (v % 4 != 0) {
        std::cerr << "GROUP_NEIGHBOR: the submodule is not singular\n";
        throw TerminalException{1};
      }
      gamma[i] = v / 4;
    }
  }
  // T: a G-stable complement of S^perp, paired perfectly with S
  ModMatrix AS(n, std::vector<int64_t>(k));
  for (size_t a = 0; a < n; a++) {
    for (size_t j = 0; j < k; j++) {
      int64_t t = 0;
      for (size_t b = 0; b < n; b++) {
        t = (t + ResMul(ResMod(A[a][b], p), S.basis[j][b], p)) % p;
      }
      AS[a][j] = t;
    }
  }
  RowSpaceModP Sperp = ComputeRowSpaceModP(LeftNullspaceModP(AS, p), p);
  ModMatrix Eperp = EquivariantProjector(Sperp, ListElt, ListEltInv, p);
  RowSpaceModP T = ComplementOfProjector(Eperp, p);
  if (T.basis.size() != k) {
    std::cerr << "GROUP_NEIGHBOR: the complement of S^perp has dimension "
              << T.basis.size() << " instead of " << k << "\n";
    throw TerminalException{1};
  }
  // M(i, l) = b(s_i, tau_l), invertible
  ModMatrix M(k, std::vector<int64_t>(k));
  for (size_t i = 0; i < k; i++) {
    for (size_t l = 0; l < k; l++) {
      M[i][l] = bil(S.basis[i], T.basis[l], p);
    }
  }
  /*
    psi(s_i, s_j) = b(s_i, phi(s_j)) mod p, phi : S -> T the correction:
    -beta / 2 for p odd; at p = 2 the G-invariant solution of psi + psi^T =
    beta, psi(s_i, s_i) = gamma_i (see above), with the action R_g of the
    elements on the basis of S.
   */
  ModMatrix psi(k, std::vector<int64_t>(k, 0));
  if (p != 2) {
    int64_t half = ResInv(2, p);
    for (size_t i = 0; i < k; i++) {
      for (size_t j = 0; j < k; j++) {
        psi[i][j] = ResMod(-ResMul(beta[i][j], half, p), p);
      }
    }
  } else {
    ModMatrix psi0(k, std::vector<int64_t>(k, 0));
    for (size_t i = 0; i < k; i++) {
      psi0[i][i] = ResMod(gamma[i], 2);
      for (size_t j = i + 1; j < k; j++) {
        psi0[i][j] = ResMod(beta[i][j], 2);
      }
    }
    // S^T, to read the coordinates of a vector of S in its basis
    ModMatrix SbT(n, std::vector<int64_t>(k));
    for (size_t a = 0; a < n; a++) {
      for (size_t i = 0; i < k; i++) {
        SbT[a][i] = S.basis[i][a];
      }
    }
    for (auto const &g : ListElt) {
      ModMatrix Sg = ModMatrixProduct(S.basis, g, 2);
      ModMatrix R(k, std::vector<int64_t>(k));
      for (size_t i = 0; i < k; i++) {
        std::optional<std::vector<int64_t>> r = SolveRowsModP(
            [&] {
              ModMatrix Sb(k, std::vector<int64_t>(n));
              for (size_t a = 0; a < k; a++) {
                Sb[a] = S.basis[a];
              }
              return Sb;
            }(),
            Sg[i], 2);
        if (!r) {
          std::cerr << "GROUP_NEIGHBOR: S is not G-stable\n";
          throw TerminalException{1};
        }
        R[i] = *r;
      }
      // psi += R psi0 R^T
      ModMatrix Rpsi = ModMatrixProduct(R, psi0, 2);
      for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < k; j++) {
          int64_t v = 0;
          for (size_t l = 0; l < k; l++) {
            v ^= (Rpsi[i][l] & R[j][l]);
          }
          psi[i][j] ^= v;
        }
      }
    }
    (void)SbT;
  }
  // t_j = sum_l x_l tau_l with sum_l M(i, l) x_l = psi(i, j)
  ModMatrix Mt(k, std::vector<int64_t>(k));
  for (size_t i = 0; i < k; i++) {
    for (size_t l = 0; l < k; l++) {
      Mt[l][i] = M[i][l];
    }
  }
  for (size_t j = 0; j < k; j++) {
    std::vector<int64_t> rhs(k);
    for (size_t i = 0; i < k; i++) {
      rhs[i] = psi[i][j];
    }
    // x Mt = rhs, i.e. sum_l x_l M(i, l) = rhs_i
    std::optional<std::vector<int64_t>> x = SolveRowsModP(Mt, rhs, p);
    if (!x) {
      std::cerr << "GROUP_NEIGHBOR: S and T are not in duality\n";
      throw TerminalException{1};
    }
    for (size_t l = 0; l < k; l++) {
      for (size_t c = 0; c < n; c++) {
        St[j][c] = ResMod(St[j][c] + p * ResMul((*x)[l], T.basis[l][c], p), p2);
      }
    }
  }
  for (size_t i = 0; i < k; i++) {
    if (bil(St[i], St[i], pnorm) != 0) {
      return {};
    }
    for (size_t j = 0; j < k; j++) {
      if (bil(St[i], St[j], p2) != 0) {
        return {};
      }
    }
  }
#ifdef SANITY_CHECK_GROUP_NEIGHBOR_MOD_P
  // G-stable modulo p^2: every s~ g is in the span of the s~ modulo p^2,
  // tested on the reduction of the coordinates
  for (auto const &g : ListElt2) {
    ModMatrix img = ModMatrixProduct(St, g, p2);
    ModMatrix Smod = reduce(St, p);
    RowSpaceModP rs = ComputeRowSpaceModP(Smod, p);
    for (auto const &row : img) {
      std::vector<int64_t> r = row;
      for (auto &x : r) {
        x = ResMod(x, p);
      }
      if (RankModPrime([&] { ModMatrix X = rs.basis; X.push_back(r); return X; }(),
                  p) != k) {
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
inline ModMatrix GroupNeighborGenerators(ModMatrix const &St, ModMatrix const &A,
                                      int64_t p) {
  size_t n = A.size();
  size_t k = St.size();
  ModMatrix AS(n, std::vector<int64_t>(k));
  for (size_t a = 0; a < n; a++) {
    for (size_t j = 0; j < k; j++) {
      int64_t t = 0;
      for (size_t b = 0; b < n; b++) {
        t = (t + ResMul(ResMod(A[a][b], p), ResMod(St[j][b], p), p)) % p;
      }
      AS[a][j] = t;
    }
  }
  ModMatrix ker = LeftNullspaceModP(AS, p);
  ModMatrix gens;
  for (size_t i = 0; i < n; i++) {
    std::vector<int64_t> v(n, 0);
    v[i] = p * p;
    gens.push_back(v);
  }
  for (auto const &x : ker) {
    std::vector<int64_t> v(n);
    for (size_t i = 0; i < n; i++) {
      v[i] = p * x[i];
    }
    gens.push_back(v);
  }
  for (auto const &row : St) {
    gens.push_back(row);
  }
  return gens;
}

// clang-format off
#endif  // SRC_GENERA_GROUPNEIGHBORMODP_H_
// clang-format on
