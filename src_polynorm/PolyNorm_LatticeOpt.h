// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_LATTICEOPT_H_
#define SRC_POLYNORM_POLYNORM_LATTICEOPT_H_

// clang-format off
#include "PolyNorm_Covering.h"
#include "PolyNorm_Packing.h"
#include "Basic_random.h"
#include "ClassicLLL.h"
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_LATTICE_OPT
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_LATTICE_OPT
#endif

#ifdef TIMINGS
#define TIMINGS_POLYNORM_LATTICE_OPT
#endif

/*
  Optimal lattices for packing and covering with a polytope.

  The lattice is L = Z^n W, the basis vectors being the rows of the n x n
  matrix W, and P = {x : ell_j(x) <= 1} is the recentered polytope of
  PolyNorm_Basic.h. A lattice vector is x = z W (z a row of integers) and
  ell_j(zW) = z W ell_j^T, which is linear in W. Equivalently Z^n is the
  lattice and the polytope has the linear forms Lmat W^T, which is how the
  packing and covering code of Z^n is reused.

  Packing. The translates P + L form a packing exactly when no nonzero
  lattice vector is in the interior of D = P - P, so the densest lattice
  packing is
      min |det W|  subject to  ||zW||_D >= 1 for all z != 0,
  and its density is vol(P) / |det W|. Each constraint is a disjunction over
  the facets d_j of D (z W d_j^T >= 1 for some j), each disjunct is linear
  in W, and det is a cubic (for n = 3) form in W.

  Covering. A hole of the lattice is an empty translate p - mu P pinned by
  n+1 pairs (z_a, i_a), meaning ell_{i_a}(p - z_a W) = mu. For fixed pairs
  the system in (p, mu) has the constant matrix with rows (ell_{i_a}, -1),
  so the radius of the hole
      mu_S(W) = sum_a c_a z_a W ell_{i_a}^T = <G_S, W>,
      G_S = sum_a c_a z_a^T ell_{i_a}
  is LINEAR in W, c being the last row of the inverse of that matrix. The
  thinnest lattice covering is max |det W| subject to mu(W) <= 1 and its
  density vol(P) / |det W|.

  The search is a sequential linear programming method with a trust region:
  at the current W the constraints that are close to active are linearized
  (they already are linear; only the choice of the facet of each contact
  and of the pairs of each hole is frozen), det is replaced by its first
  order expansion det(W + X) = det(W) (1 + tr(W^{-1} X)) and the linear
  program gives the step. Every trial point is evaluated exactly with the
  packing and covering code, so the linear programs only guide the search
  and need not be exact models.

  At the end of the search the active constraints are imposed as
  equalities, the critical point of det on that affine space is computed
  by Newton's method in floating point and recognized as a rational
  matrix, and then checked exactly:
  --- the lattice has the expected packing or covering radius;
  --- first order: s W^{-T} = sum_k lambda_k h_k with lambda >= 0 (a
      Karush-Kuhn-Tucker point of min s det(W) s.t. <h_k, W> >= b_k; for
      packing s = 1, h_k = z^T d_j; for covering s = -1, h_k = -G_S);
  --- second order: on T = {X : <h_k, X> = 0 for lambda_k > 0} the Hessian
      of s det is s det(W) ((tr Y)^2 - tr(Y^2)) with Y = W^{-1} X, and
      tr Y = 0 on T, so the requirement is that -s tr(Y^2) is positive
      definite on T. For packing (s = 1) the matrices Y of T must be
      "mostly antisymmetric", which needs dim T <= n(n-1)/2, that is at
      least n(n+1)/2 independent contacts (Swinnerton-Dyer's bound).
  These two conditions are the classical sufficient conditions for a
  strict local minimum of the program. For packing the feasible region
  near W is exactly the union over the choices of one active facet per
  contact of polyhedra, so checking each choice ("branch") certifies a
  strict local optimum of the packing problem. For covering the region
  {mu(W) <= 1} is contained near W in the polyhedron of the constraints
  mu_S(W) <= 1 for the holes having exactly n+1 lattice points on their
  boundary (choosing one active facet per lattice point gives a lower
  bound on mu), which certifies the covering case when every deepest hole
  has exactly n+1 lattice points. Holes with more lattice points are
  reported and the certificate is not claimed.

  Nothing here proves global optimality, except in the cases where the
  density is 1 (a tiling). The 2-dimensional packing problem has an exact
  global algorithm in PolyNorm_Packing2D.h.
 */

// The flattening of an n x n matrix into a vector of length n^2.
template <typename T> MyVector<T> PolyNorm_FlattenMat(MyMatrix<T> const &M) {
  int n = M.rows();
  MyVector<T> V(n * n);
  for (int k = 0; k < n; k++) {
    for (int b = 0; b < n; b++) {
      V(k * n + b) = M(k, b);
    }
  }
  return V;
}

template <typename T>
MyMatrix<T> PolyNorm_UnflattenVect(MyVector<T> const &V, int n) {
  MyMatrix<T> M(n, n);
  for (int k = 0; k < n; k++) {
    for (int b = 0; b < n; b++) {
      M(k, b) = V(k * n + b);
    }
  }
  return M;
}

// <A, B> = sum_{k,b} A(k,b) B(k,b).
template <typename T>
T PolyNorm_MatScalar(MyMatrix<T> const &A, MyMatrix<T> const &B) {
  int n = A.rows();
  int p = A.cols();
  T sum(0);
  for (int k = 0; k < n; k++) {
    for (int b = 0; b < p; b++) {
      sum += A(k, b) * B(k, b);
    }
  }
  return sum;
}

// The entries as a std::vector, usable as a key of an ordered container.
template <typename T> std::vector<T> PolyNorm_Key(MyVector<T> const &V) {
  std::vector<T> key(V.size());
  for (int k = 0; k < V.size(); k++) {
    key[k] = V(k);
  }
  return key;
}

template <typename T> double PolyNorm_ToDouble(T const &x) {
  return UniversalScalarConversion<double, T>(x);
}

template <typename T>
MyMatrix<double> PolyNorm_ToDoubleMat(MyMatrix<T> const &M) {
  return UniversalMatrixConversion<double, T>(M);
}

// The rational m / 2^e closest to x, for a grid fine enough for the trust
// region but keeping the denominators of the iterates bounded.
template <typename T> T PolyNorm_RoundToGrid(double x, int e) {
  double scaled = std::ldexp(x, e);
  long num = std::lround(scaled);
  T val(num);
  T den(1);
  for (int i = 0; i < e; i++) {
    den *= 2;
  }
  return val / den;
}

// The simplest rational within tol of x: the first convergent of the
// continued fraction of x that is close enough.
template <typename T>
std::optional<T> PolyNorm_RationalApprox(double x, double tol) {
  double max_den = 1.0e12;
  double h_prev = 1, h_prev2 = 0, k_prev = 0, k_prev2 = 1;
  double y = x;
  for (int iter = 0; iter < 64; iter++) {
    double a = std::floor(y);
    double h = a * h_prev + h_prev2;
    double k = a * k_prev + k_prev2;
    if (k > max_den || std::fabs(h) > 1.0e15) {
      return {};
    }
    if (std::fabs(x - h / k) <= tol) {
      T hT(static_cast<long>(h));
      T kT(static_cast<long>(k));
      return hT / kT;
    }
    double frac = y - a;
    if (frac == 0) {
      return {};
    }
    y = 1 / frac;
    h_prev2 = h_prev;
    h_prev = h;
    k_prev2 = k_prev;
    k_prev = k;
  }
  return {};
}

// The data of the polytope for the lattice Z^n W, expressed for Z^n: the
// vertices v W^{-1}, the linear forms Lmat W^T. The symmetry group is
// taken trivial: a generic lattice breaks the symmetries of P, and the
// covering search then only identifies the nodes up to translation.
template <typename T, typename Tint, typename Tgroup>
PolyNormData<T, Tint, Tgroup>
PolyNorm_TransformedData(PolyNormData<T, Tint, Tgroup> const &data,
                         MyMatrix<T> const &W) {
  using Tidx = typename Tgroup::Telt::Tidx;
  int n = data.n;
  PolyNormData<T, Tint, Tgroup> res;
  res.n = n;
  res.center = ZeroVector<T>(n);
  MyMatrix<T> Winv = Inverse(W);
  int n_vert = data.EXT.rows();
  res.EXT = MyMatrix<T>(n_vert, n + 1);
  for (int i = 0; i < n_vert; i++) {
    res.EXT(i, 0) = 1;
    for (int b = 0; b < n; b++) {
      T sum(0);
      for (int k = 0; k < n; k++) {
        sum += data.EXT(i, k + 1) * Winv(k, b);
      }
      res.EXT(i, b + 1) = sum;
    }
  }
  res.Lmat = data.Lmat * W.transpose();
  // vectface is not copyable.
  vectface ListIncd(n_vert);
  for (auto &f : data.ListIncd) {
    ListIncd.push_back(f);
  }
  res.ListIncd = std::move(ListIncd);
  res.LmatDiff = data.LmatDiff * W.transpose();
  Tgroup GRP(static_cast<Tidx>(n_vert));
  res.GRP = GRP;
  int m = data.Lmat.rows();
  std::vector<int> id_perm(m);
  for (int j = 0; j < m; j++) {
    id_perm[j] = j;
  }
  std::vector<PolyNormSymm<Tint>> ListElt{{IdentityMat<Tint>(n), id_perm}};
  res.ListElt = std::move(ListElt);
  return res;
}

/*
  Packing evaluation.
 */

template <typename T, typename Tint> struct PolyNormLatticePackingEval {
  // The first minimum of D for the lattice: the lattice alpha^{-1} Z^n W
  // is a packing lattice, so the packing density is
  // vol(P) alpha^n / |det W|.
  T alpha;
  T det;
  // The contact vectors z, one per pair +-z.
  std::vector<MyVector<Tint>> ListContact;
};

// The canonical sign of a vector: first nonzero coordinate positive.
template <typename Tint>
bool PolyNorm_IsCanonicalSign(MyVector<Tint> const &z) {
  for (int k = 0; k < z.size(); k++) {
    if (z(k) > 0) {
      return true;
    }
    if (z(k) < 0) {
      return false;
    }
  }
  return false;
}

template <typename T, typename Tint>
PolyNormLatticePackingEval<T, Tint>
PolyNorm_EvaluatePacking(MyMatrix<T> const &LmatDiff, MyMatrix<T> const &W,
                         std::ostream &os) {
  MyMatrix<T> LmatDiffW = LmatDiff * W.transpose();
  PolyNormPacking<T, Tint> pack =
      ComputePolyNormPacking<T, Tint>(LmatDiffW, os);
  std::vector<MyVector<Tint>> ListContact;
  for (auto &z : pack.ListContact) {
    if (PolyNorm_IsCanonicalSign(z)) {
      ListContact.push_back(z);
    }
  }
  return {pack.alpha, T_abs(DeterminantMat(W)), std::move(ListContact)};
}

/*
  Covering evaluation.
 */

// A hole of the lattice in the original coordinates: p - mu P is an empty
// translate whose boundary contains the lattice points z W of the pairs.
template <typename T, typename Tint> struct PolyNormHole {
  MyVector<T> p;
  T mu;
  std::vector<std::pair<MyVector<Tint>, int>> ListPair;
  // The number of distinct lattice points among the pairs.
  int n_point;
};

template <typename T, typename Tint> struct PolyNormLatticeCoveringEval {
  T mu;
  T det;
  // The holes whose radius is at least the collecting ratio times mu,
  // empty when the holes are not collected.
  std::vector<PolyNormHole<T, Tint>> ListHole;
};

// The tight pairs of the point p' of the transformed problem (lattice
// Z^n, linear forms LmatW) at the radius mu: the lattice points z with
// ||p - z|| = mu and the facets i with ell_i(p - z) = mu. With a relative
// tolerance tol > 0, the pairs that are tight up to tol: near an optimum
// the holes are often degenerate (their lattice points on lower
// dimensional faces of the translates) and the search only sees the
// degeneracy appear through such nearly tight pairs.
template <typename T, typename Tint>
std::vector<std::pair<MyVector<Tint>, int>>
PolyNorm_TightPairs(MyMatrix<T> const &LmatW, MyVector<T> const &p, T const &mu,
                    T const &tol, std::ostream &os) {
  int m = LmatW.rows();
  int n = LmatW.cols();
  T mu_low = mu * (1 - tol);
  T mu_upp = mu * (1 + tol);
  std::vector<MyVector<Tint>> ListZ =
      PolyNorm_BlockingPoints<T, Tint>(LmatW, p, mu_upp, os);
  std::vector<std::pair<MyVector<Tint>, int>> ListPair;
  MyVector<T> diff(n);
  for (auto &z : ListZ) {
    diff = p - UniversalVectorConversion<T, Tint>(z);
    T gauge = PolyNorm_Gauge(LmatW, diff);
    if (gauge < mu_low) {
      continue;
    }
    for (int i = 0; i < m; i++) {
      T val(0);
      for (int k = 0; k < n; k++) {
        val += LmatW(i, k) * diff(k);
      }
      if (val >= mu_low) {
        ListPair.push_back({z, i});
      }
    }
  }
  return ListPair;
}

template <typename T, typename Tint, typename Tgroup>
PolyNormLatticeCoveringEval<T, Tint>
PolyNorm_EvaluateCovering(PolyNormData<T, Tint, Tgroup> const &data,
                          MyMatrix<T> const &W,
                          std::optional<T> const &collect_ratio,
                          std::ostream &os, T const &pair_tol = T(0)) {
  PolyNormData<T, Tint, Tgroup> dataW =
      PolyNorm_TransformedData<T, Tint, Tgroup>(data, W);
  PolyNormCoveringSolver<T, Tint, Tgroup> solver(dataW, os);
  solver.collect_ratio = collect_ratio;
  PolyNormCovering<T, Tint> cov = solver.Run();
  T det = T_abs(DeterminantMat(W));
  std::vector<PolyNormHole<T, Tint>> ListHole;
  if (collect_ratio) {
    T thr = *collect_ratio * cov.mu;
    MyMatrix<T> WT = W.transpose();
    // The holes are identified up to translation by their pairs, translated
    // so that the smallest lattice point is 0.
    std::set<std::vector<std::pair<std::vector<Tint>, int>>> set_key;
    for (auto &ePair : solver.ListSolved) {
      MyVector<T> const &p = ePair.first;
      T const &mu = ePair.second;
      if (mu < thr) {
        continue;
      }
      std::vector<std::pair<MyVector<Tint>, int>> ListPair =
          PolyNorm_TightPairs<T, Tint>(dataW.Lmat, p, mu, pair_tol, os);
      std::set<std::vector<Tint>> set_z;
      for (auto &eP : ListPair) {
        set_z.insert(PolyNorm_Key(eP.first));
      }
      std::vector<Tint> const &z_min = *set_z.begin();
      std::vector<std::pair<std::vector<Tint>, int>> key;
      for (auto &eP : ListPair) {
        std::vector<Tint> V = PolyNorm_Key(eP.first);
        for (size_t k = 0; k < V.size(); k++) {
          V[k] -= z_min[k];
        }
        key.push_back({V, eP.second});
      }
      std::sort(key.begin(), key.end());
      if (set_key.count(key) > 0) {
        continue;
      }
      set_key.insert(key);
      PolyNormHole<T, Tint> hole{WT * p, mu, std::move(ListPair),
                                 static_cast<int>(set_z.size())};
      ListHole.push_back(std::move(hole));
    }
  }
  return {cov.mu, det, std::move(ListHole)};
}

/*
  The linear constraints.

  A constraint is <h, W> >= b with h an n x n matrix. For packing it is
  z W d_j^T >= alpha, that is h = z^T d_j. For covering it is mu_S(W) <=
  mu, that is h = -G_S and b = -mu.
 */

template <typename T> struct PolyNormLinCons {
  MyMatrix<T> h;
  T b;
};

// The constraints of the contacts. For each contact z with
// ||zW||_D <= ratio alpha, the facets attaining the gauge are the active
// ones; the linear program uses one of them (chosen at random among ties),
// the certification all of them.
template <typename T, typename Tint> struct PolyNormContactFacets {
  MyVector<Tint> z;
  T val;
  std::vector<int> ListFacet;
};

template <typename T, typename Tint>
std::vector<PolyNormContactFacets<T, Tint>>
PolyNorm_PackingNearContacts(MyMatrix<T> const &LmatDiff, MyMatrix<T> const &W,
                             T const &radius, std::ostream &os) {
  MyMatrix<T> LmatDiffW = LmatDiff * W.transpose();
  std::vector<MyVector<Tint>> ListPt =
      PolyNorm_LatticePoints<T, Tint>(LmatDiffW, radius, os);
  int m = LmatDiffW.rows();
  int n = LmatDiffW.cols();
  std::vector<PolyNormContactFacets<T, Tint>> ListRet;
  for (auto &z : ListPt) {
    if (!PolyNorm_IsCanonicalSign(z)) {
      continue;
    }
    std::optional<T> val;
    std::vector<int> ListFacet;
    for (int j = 0; j < m; j++) {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += LmatDiffW(j, k) * z(k);
      }
      if (!val || scal > *val) {
        val = scal;
        ListFacet.clear();
      }
      if (scal == *val) {
        ListFacet.push_back(j);
      }
    }
    ListRet.push_back({z, *val, std::move(ListFacet)});
  }
  return ListRet;
}

// h = z^T d_j for the contact z and the facet j of D.
template <typename T, typename Tint>
MyMatrix<T> PolyNorm_PackingGradient(MyMatrix<T> const &LmatDiff,
                                     MyVector<Tint> const &z, int j) {
  int n = LmatDiff.cols();
  MyMatrix<T> h(n, n);
  for (int k = 0; k < n; k++) {
    T zk = UniversalScalarConversion<T, Tint>(z(k));
    for (int b = 0; b < n; b++) {
      h(k, b) = zk * LmatDiff(j, b);
    }
  }
  return h;
}

// For n+1 pairs (z_a, i_a), the matrix G_S with mu_S(W) = <G_S, W>, or
// nothing when the normals are affinely dependent.
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
PolyNorm_HoleGradient(MyMatrix<T> const &Lmat,
                      std::vector<std::pair<MyVector<Tint>, int>> const &S) {
  int n = Lmat.cols();
  MyMatrix<T> M(n + 1, n + 1);
  for (int a = 0; a <= n; a++) {
    int i = S[a].second;
    for (int k = 0; k < n; k++) {
      M(a, k) = Lmat(i, k);
    }
    M(a, n) = -1;
  }
  if (DeterminantMat(M) == 0) {
    return {};
  }
  MyMatrix<T> Minv = Inverse(M);
  MyMatrix<T> G = ZeroMatrix<T>(n, n);
  for (int a = 0; a <= n; a++) {
    T c = Minv(n, a);
    int i = S[a].second;
    for (int k = 0; k < n; k++) {
      T zk = UniversalScalarConversion<T, Tint>(S[a].first(k));
      for (int b = 0; b < n; b++) {
        G(k, b) += c * zk * Lmat(i, b);
      }
    }
  }
  return G;
}

// The matrices G_S of the systems of a hole. When the hole has exactly n+1
// lattice points, these are the selections of one facet per lattice
// point: each gives a lower bound mu_S(W) <= mu(W) near the current W,
// which is what the certification needs. Otherwise all the choices of n+1
// distinct lattice points with one facet each are taken, which only
// guides the search. At most max_system are generated.
template <typename T, typename Tint>
std::vector<MyMatrix<T>> PolyNorm_HoleSystems(MyMatrix<T> const &Lmat,
                                              PolyNormHole<T, Tint> const &hole,
                                              size_t max_system) {
  int n = Lmat.cols();
  // Group the pairs by lattice point.
  std::map<std::vector<Tint>, std::pair<MyVector<Tint>, std::vector<int>>>
      map_z;
  for (auto &eP : hole.ListPair) {
    auto &entry = map_z[PolyNorm_Key(eP.first)];
    entry.first = eP.first;
    entry.second.push_back(eP.second);
  }
  std::vector<std::pair<MyVector<Tint>, std::vector<int>>> ListZ;
  for (auto &kv : map_z) {
    ListZ.push_back(kv.second);
  }
  int n_z = ListZ.size();
  std::vector<MyMatrix<T>> ListG;
  if (n_z < n + 1) {
    return ListG;
  }
  // Enumerate the subsets of n+1 lattice points and for each the
  // selections of facets.
  std::vector<int> sub(n + 1);
  for (int a = 0; a <= n; a++) {
    sub[a] = a;
  }
  while (true) {
    std::vector<int> sel(n + 1, 0);
    while (true) {
      std::vector<std::pair<MyVector<Tint>, int>> S;
      for (int a = 0; a <= n; a++) {
        auto const &eZ = ListZ[sub[a]];
        S.push_back({eZ.first, eZ.second[sel[a]]});
      }
      std::optional<MyMatrix<T>> optG = PolyNorm_HoleGradient<T, Tint>(Lmat, S);
      if (optG) {
        ListG.push_back(*optG);
        if (ListG.size() >= max_system) {
          return ListG;
        }
      }
      int a = 0;
      while (a <= n) {
        sel[a]++;
        if (sel[a] < static_cast<int>(ListZ[sub[a]].second.size())) {
          break;
        }
        sel[a] = 0;
        a++;
      }
      if (a > n) {
        break;
      }
    }
    // Next subset in lexicographic order.
    int a = n;
    while (a >= 0 && sub[a] == n_z - (n + 1) + a) {
      a--;
    }
    if (a < 0) {
      break;
    }
    sub[a]++;
    for (int b = a + 1; b <= n; b++) {
      sub[b] = sub[b - 1] + 1;
    }
  }
  return ListG;
}

// The covering constraints mu_S(W') <= max(mu, mu_S(W)) from the holes:
// h = -G_S, b = -max(mu, mu_S(W)). For exactly tight pairs mu_S(W) = mu;
// the systems of nearly tight pairs can exceed mu at W and the maximum
// keeps W feasible. Duplicates (the same hole found from several of its
// lattice points, or a hole and its opposite for a symmetric P) are
// removed.
template <typename T, typename Tint>
std::vector<PolyNormLinCons<T>>
PolyNorm_CoveringConstraints(MyMatrix<T> const &Lmat,
                             std::vector<PolyNormHole<T, Tint>> const &ListHole,
                             T const &mu, MyMatrix<T> const &W,
                             bool only_nondegenerate, size_t max_system) {
  int n = Lmat.cols();
  std::set<std::vector<T>> set_G;
  std::vector<PolyNormLinCons<T>> ListCons;
  for (auto &hole : ListHole) {
    if (only_nondegenerate && hole.n_point != n + 1) {
      continue;
    }
    for (auto &G : PolyNorm_HoleSystems<T, Tint>(Lmat, hole, max_system)) {
      std::vector<T> V = PolyNorm_Key(PolyNorm_FlattenMat(G));
      if (set_G.count(V) > 0) {
        continue;
      }
      set_G.insert(V);
      T muS = PolyNorm_MatScalar(G, W);
      T b = muS > mu ? -muS : -mu;
      ListCons.push_back({-G, b});
    }
  }
  return ListCons;
}

/*
  The step: the linear program
     minimize s tr(W^{-1} X)
     subject to <h_k, W + X> >= b_k, |X(k,b)| <= r.
 */
template <typename T>
std::optional<MyMatrix<T>>
PolyNorm_LinearStep(MyMatrix<T> const &W,
                    std::vector<PolyNormLinCons<T>> const &ListCons, int s,
                    T const &r, std::ostream &os) {
  int n = W.rows();
  int nn = n * n;
  int n_cons = ListCons.size();
  MyMatrix<T> ListIneq = ZeroMatrix<T>(n_cons + 2 * nn, nn + 1);
  int pos = 0;
  for (auto &eCons : ListCons) {
    ListIneq(pos, 0) = PolyNorm_MatScalar(eCons.h, W) - eCons.b;
    MyVector<T> V = PolyNorm_FlattenMat(eCons.h);
    for (int u = 0; u < nn; u++) {
      ListIneq(pos, u + 1) = V(u);
    }
    pos++;
  }
  for (int u = 0; u < nn; u++) {
    ListIneq(pos, 0) = r;
    ListIneq(pos, u + 1) = 1;
    pos++;
    ListIneq(pos, 0) = r;
    ListIneq(pos, u + 1) = -1;
    pos++;
  }
  // tr(W^{-1} X) = sum_{k,b} Winv(b,k) X(k,b) = <Winv^T, X>.
  MyMatrix<T> WinvT = Inverse(W).transpose();
  MyVector<T> ToBeMinimized = ZeroVector<T>(nn + 1);
  MyVector<T> V = PolyNorm_FlattenMat(WinvT);
  for (int u = 0; u < nn; u++) {
    ToBeMinimized(u + 1) = T(s) * V(u);
  }
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
  if (!eSol.DirectSolution || !eSol.DualSolution) {
    return {};
  }
  MyVector<T> X(nn);
  for (int u = 0; u < nn; u++) {
    X(u) = (*eSol.DirectSolution)(u);
  }
  return PolyNorm_UnflattenVect(X, n);
}

/*
  The polishing: the critical point of det on the affine space
  {W : <h_k, W> = b_k}, recognized as a rational matrix.
 */

template <typename T>
std::optional<MyMatrix<T>>
PolyNorm_Polish(MyMatrix<T> const &W_approx,
                std::vector<PolyNormLinCons<T>> const &ListEq, double tol,
                double &det_numeric, [[maybe_unused]] std::ostream &os) {
  // |det| of the floating point critical point when Newton's method
  // converged, -1 otherwise. It is the useful output when the critical
  // point is irrational.
  det_numeric = -1;
  int n = W_approx.rows();
  int nn = n * n;
  int k = ListEq.size();
  // The affine space: C vec(W) = b.
  MyMatrix<T> C(k, nn);
  MyVector<T> bvec(k);
  for (int u = 0; u < k; u++) {
    MyVector<T> V = PolyNorm_FlattenMat(ListEq[u].h);
    for (int v = 0; v < nn; v++) {
      C(u, v) = V(v);
    }
    bvec(u) = ListEq[u].b;
  }
  // A particular solution and the directions.
  MyMatrix<T> CT = C.transpose();
  std::optional<MyVector<T>> opt = SolutionMat(CT, bvec);
  if (!opt) {
#ifdef DEBUG_POLYNORM_LATTICE_OPT
    os << "POLYNORM_LATTICE_OPT: Polish, the active constraints are "
       << "inconsistent\n";
#endif
    return {};
  }
  MyVector<T> w0 = *opt;
  // NullspaceTrMat(C) gives the vectors x orthogonal to the rows of C,
  // that is the kernel of C, as rows.
  MyMatrix<T> NSP = NullspaceTrMat(C);
  int d = NSP.rows();
  MyMatrix<T> W0 = PolyNorm_UnflattenVect(w0, n);
  std::vector<MyMatrix<T>> ListM;
  for (int i = 0; i < d; i++) {
    MyVector<T> V = GetMatrixRow(NSP, i);
    ListM.push_back(PolyNorm_UnflattenVect(V, n));
  }
  if (d == 0) {
    if (DeterminantMat(W0) == 0) {
      return {};
    }
    return W0;
  }
  // Starting point: the least squares projection of W_approx.
  MyMatrix<double> W0_d = PolyNorm_ToDoubleMat(W0);
  std::vector<MyMatrix<double>> ListM_d;
  for (auto &M : ListM) {
    ListM_d.push_back(PolyNorm_ToDoubleMat(M));
  }
  MyMatrix<double> Wa_d = PolyNorm_ToDoubleMat(W_approx);
  MyMatrix<double> Gram(d, d);
  MyVector<double> rhs(d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      Gram(i, j) = (ListM_d[i].array() * ListM_d[j].array()).sum();
    }
    rhs(i) = (ListM_d[i].array() * (Wa_d - W0_d).array()).sum();
  }
  MyVector<double> lambda = Gram.fullPivLu().solve(rhs);
  auto get_W = [&](MyVector<double> const &lam) -> MyMatrix<double> {
    MyMatrix<double> Wd = W0_d;
    for (int i = 0; i < d; i++) {
      Wd += lam(i) * ListM_d[i];
    }
    return Wd;
  };
  // Newton's method on the gradient of log|det|: g_i = tr(W^{-1} M_i),
  // H_ij = -tr(W^{-1} M_i W^{-1} M_j).
  double norm_g = 0;
  for (int iter = 0; iter < 200; iter++) {
    MyMatrix<double> Wd = get_W(lambda);
    MyMatrix<double> Winv = Wd.inverse();
    std::vector<MyMatrix<double>> ListP;
    MyVector<double> g(d);
    for (int i = 0; i < d; i++) {
      ListP.push_back(Winv * ListM_d[i]);
      g(i) = ListP[i].trace();
    }
    MyMatrix<double> H(d, d);
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        H(i, j) = -(ListP[i] * ListP[j]).trace();
      }
    }
    norm_g = g.norm();
    if (norm_g < 1.0e-15) {
      break;
    }
    MyVector<double> delta = H.completeOrthogonalDecomposition().solve(-g);
    lambda += delta;
    if (delta.norm() < 1.0e-16 * (1 + lambda.norm())) {
      break;
    }
  }
#ifdef DEBUG_POLYNORM_LATTICE_OPT
  os << "POLYNORM_LATTICE_OPT: Polish, d=" << d << " |grad|=" << norm_g << "\n";
#endif
  if (norm_g < 1.0e-12) {
    det_numeric = std::fabs(get_W(lambda).determinant());
  }
  // Recognition of lambda as a rational vector.
  MyMatrix<T> W = W0;
  for (int i = 0; i < d; i++) {
    std::optional<T> opt_q = PolyNorm_RationalApprox<T>(lambda(i), tol);
    if (!opt_q) {
      return {};
    }
    W += *opt_q * ListM[i];
  }
  if (DeterminantMat(W) == 0) {
    return {};
  }
  // The gradient must vanish exactly.
  MyMatrix<T> Winv = Inverse(W);
  for (int i = 0; i < d; i++) {
    MyMatrix<T> P = Winv * ListM[i];
    T tr(0);
    for (int u = 0; u < n; u++) {
      tr += P(u, u);
    }
    if (tr != 0) {
#ifdef DEBUG_POLYNORM_LATTICE_OPT
      os << "POLYNORM_LATTICE_OPT: Polish, the rational point is not "
         << "critical, i=" << i << " tr=" << tr << "\n";
#endif
      return {};
    }
  }
  return W;
}

/*
  The certification: Karush-Kuhn-Tucker multipliers and second order
  condition for min s det(W) subject to <h_k, W> >= b_k, all of which are
  active at W.
 */

// Whether the quadratic form t^T Q t is positive on the standard simplex
// {t >= 0, sum t = 1}, exactly. The minimum is attained in the relative
// interior of some face, at a critical point of the form restricted to the
// affine hull of the face: the solution of [Q_S 1; 1^T 0] (t, -nu) = (0, 1).
// A face whose system is singular either has no critical point or a
// family of them along which the form is constant, and then the minimum is
// also attained on a smaller face. So it suffices to check the nonsingular
// faces, of which there are at most 2^p - 1.
template <typename T>
std::optional<bool> PolyNorm_IsStrictlyCopositive(MyMatrix<T> const &Q,
                                                  int max_size) {
  int p = Q.rows();
  if (p > max_size) {
    return {};
  }
  for (uint64_t mask = 1; mask < (uint64_t(1) << p); mask++) {
    std::vector<int> S;
    for (int i = 0; i < p; i++) {
      if (mask & (uint64_t(1) << i)) {
        S.push_back(i);
      }
    }
    int k = S.size();
    MyMatrix<T> B = ZeroMatrix<T>(k + 1, k + 1);
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) {
        B(a, b) = Q(S[a], S[b]);
      }
      B(a, k) = 1;
      B(k, a) = 1;
    }
    if (DeterminantMat(B) == 0) {
      continue;
    }
    MyVector<T> rhs = ZeroVector<T>(k + 1);
    rhs(k) = 1;
    MyVector<T> sol = Inverse(B) * rhs;
    bool interior = true;
    for (int a = 0; a < k; a++) {
      if (sol(a) <= 0) {
        interior = false;
      }
    }
    if (!interior) {
      continue;
    }
    T val(0);
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) {
        val += sol(a) * Q(S[a], S[b]) * sol(b);
      }
    }
    if (val <= 0) {
      return false;
    }
  }
  return true;
}

// The extreme rays of the pointed cone {y : F y >= 0}, F of rank d = the
// number of columns: each is the kernel of d-1 independent rows, with the
// sign making it feasible.
template <typename T>
std::vector<MyVector<T>> PolyNorm_ConeExtremeRays(MyMatrix<T> const &F) {
  int k = F.rows();
  int d = F.cols();
  std::vector<MyVector<T>> ListRay;
  std::set<std::vector<T>> set_ray;
  if (d == 1) {
    for (int sign : {1, -1}) {
      MyVector<T> y(1);
      y(0) = sign;
      if ((F * y).minCoeff() >= 0) {
        ListRay.push_back(y);
      }
    }
    return ListRay;
  }
  std::vector<int> sub(d - 1);
  for (int a = 0; a < d - 1; a++) {
    sub[a] = a;
  }
  if (k < d - 1) {
    return ListRay;
  }
  while (true) {
    MyMatrix<T> Fsub = SelectRow(F, sub);
    MyMatrix<T> NSP = NullspaceTrMat(Fsub);
    if (NSP.rows() == 1) {
      MyVector<T> y = GetMatrixRow(NSP, 0);
      for (int sign : {1, -1}) {
        MyVector<T> ys = T(sign) * y;
        MyVector<T> Fy = F * ys;
        bool feasible = true;
        for (int u = 0; u < k; u++) {
          if (Fy(u) < 0) {
            feasible = false;
          }
        }
        if (feasible) {
          // Normalization: first nonzero entry of absolute value 1.
          T scal(0);
          for (int u = 0; u < d; u++) {
            if (ys(u) != 0) {
              scal = T_abs(ys(u));
              break;
            }
          }
          MyVector<T> yn = ys / scal;
          std::vector<T> key = PolyNorm_Key(yn);
          if (set_ray.count(key) == 0) {
            set_ray.insert(key);
            ListRay.push_back(yn);
          }
        }
      }
    }
    int a = d - 2;
    while (a >= 0 && sub[a] == k - (d - 1) + a) {
      a--;
    }
    if (a < 0) {
      break;
    }
    sub[a]++;
    for (int b = a + 1; b < d - 1; b++) {
      sub[b] = sub[b - 1] + 1;
    }
  }
  return ListRay;
}

struct PolyNormKKT {
  // Whether there are multipliers lambda >= 0.
  bool first_order;
  // The number of positive multipliers in a solution of maximal support.
  int n_positive;
  // The dimension of the space T on which the second order condition is
  // checked, and whether the form is positive on the critical cone.
  int dim_tangent;
  bool second_order;
  // The number of extreme rays of the critical cone, when the form is not
  // positive definite on T and the cone has to be used.
  int n_ray;
};

template <typename T>
PolyNormKKT PolyNorm_CheckKKT(MyMatrix<T> const &W,
                              std::vector<MyMatrix<T>> const &ListH, int s,
                              std::ostream &os) {
  int n = W.rows();
  int nn = n * n;
  int k = ListH.size();
  MyMatrix<T> Winv = Inverse(W);
  T sign_det = DeterminantMat(W) > 0 ? T(1) : T(-1);
  // The gradient of s det(W) divided by |det W|: s sign(det) W^{-T}.
  MyMatrix<T> grad = T(s) * sign_det * Winv.transpose();
  MyVector<T> grad_V = PolyNorm_FlattenMat(grad);
  // Variables (lambda_1..k, y_1..k). Constraints: sum lambda h = grad
  // (two inequalities per coordinate), lambda >= 0, y <= lambda, y <= 1,
  // and maximize sum y (y >= 0 is not needed).
  int n_var = 2 * k;
  MyMatrix<T> ListIneq = ZeroMatrix<T>(2 * nn + 3 * k, n_var + 1);
  int pos = 0;
  std::vector<MyVector<T>> ListH_V;
  for (auto &h : ListH) {
    ListH_V.push_back(PolyNorm_FlattenMat(h));
  }
  for (int u = 0; u < nn; u++) {
    ListIneq(pos, 0) = -grad_V(u);
    ListIneq(pos + 1, 0) = grad_V(u);
    for (int i = 0; i < k; i++) {
      ListIneq(pos, 1 + i) = ListH_V[i](u);
      ListIneq(pos + 1, 1 + i) = -ListH_V[i](u);
    }
    pos += 2;
  }
  for (int i = 0; i < k; i++) {
    ListIneq(pos, 1 + i) = 1;
    pos++;
    ListIneq(pos, 1 + i) = 1;
    ListIneq(pos, 1 + k + i) = -1;
    pos++;
    ListIneq(pos, 0) = 1;
    ListIneq(pos, 1 + k + i) = -1;
    pos++;
  }
  MyVector<T> ToBeMinimized = ZeroVector<T>(n_var + 1);
  for (int i = 0; i < k; i++) {
    ToBeMinimized(1 + k + i) = -1;
  }
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
  PolyNormKKT res{false, 0, -1, false, 0};
  if (!eSol.DirectSolution || !eSol.DualSolution) {
    return res;
  }
  res.first_order = true;
  MyVector<T> const &x = *eSol.DirectSolution;
  std::vector<int> ListPos;
  for (int i = 0; i < k; i++) {
    if (x(i) > 0) {
      ListPos.push_back(i);
    }
  }
  res.n_positive = ListPos.size();
  std::vector<int> ListZero;
  for (int i = 0; i < k; i++) {
    if (x(i) == 0) {
      ListZero.push_back(i);
    }
  }
  // The tangent space T = {X : <h_i, X> = 0, lambda_i > 0}.
  MyMatrix<T> Cmat(ListPos.size(), nn);
  for (size_t u = 0; u < ListPos.size(); u++) {
    for (int v = 0; v < nn; v++) {
      Cmat(u, v) = ListH_V[ListPos[u]](v);
    }
  }
  MyMatrix<T> NSP = ListPos.empty() ? IdentityMat<T>(nn) : NullspaceTrMat(Cmat);
  int d = NSP.rows();
  res.dim_tangent = d;
  if (d == 0) {
    res.second_order = true;
    return res;
  }
  // The Hessian of s det divided by |det W| on T is
  // -s sign(det) tr(Y_a Y_b), Y_a = W^{-1} X_a.
  std::vector<MyMatrix<T>> ListY;
  for (int a = 0; a < d; a++) {
    MyVector<T> V = GetMatrixRow(NSP, a);
    ListY.push_back(Winv * PolyNorm_UnflattenVect(V, n));
  }
  MyMatrix<T> Q(d, d);
  for (int a = 0; a < d; a++) {
    for (int b = a; b < d; b++) {
      MyMatrix<T> P = ListY[a] * ListY[b];
      T tr(0);
      for (int u = 0; u < n; u++) {
        tr += P(u, u);
      }
      Q(a, b) = -T(s) * sign_det * tr;
      Q(b, a) = Q(a, b);
    }
  }
  if (IsPositiveDefinite(Q, os)) {
    res.second_order = true;
    return res;
  }
  // The critical cone C = {X in T : <h_i, X> >= 0 for lambda_i = 0}. The
  // second order sufficient condition is that the form is positive on
  // C minus 0. When C is pointed, that is strict copositivity on its
  // extreme rays.
  MyMatrix<T> F(ListZero.size(), d);
  for (size_t u = 0; u < ListZero.size(); u++) {
    for (int a = 0; a < d; a++) {
      MyVector<T> V = GetMatrixRow(NSP, a);
      T scal(0);
      for (int v = 0; v < nn; v++) {
        scal += ListH_V[ListZero[u]](v) * V(v);
      }
      F(u, a) = scal;
    }
  }
  if (ListZero.empty()) {
    // C = T, on which the form is not positive definite.
    return res;
  }
  // C = L + C' with L = ker F the lineality space and C' = C cap M for
  // M = rowspace(F), a complement of L on which F is injective, so that C'
  // is pointed. With the coordinates (l, m), the form is positive on
  // C minus 0 if it is positive definite on L and its Schur complement
  // S = Q_MM - Q_ML Q_LL^{-1} Q_LM (the minimum over l) is strictly
  // copositive on C'.
  MyMatrix<T> Lbasis = NullspaceTrMat(F);
  MyMatrix<T> Mbasis = RowReduction(F);
  int dL = Lbasis.rows();
  MyMatrix<T> S;
  if (dL > 0) {
    MyMatrix<T> Q_LL = Lbasis * Q * Lbasis.transpose();
    if (!IsPositiveDefinite(Q_LL, os)) {
      return res;
    }
    MyMatrix<T> Q_LM = Lbasis * Q * Mbasis.transpose();
    MyMatrix<T> Q_MM = Mbasis * Q * Mbasis.transpose();
    S = Q_MM - Q_LM.transpose() * Inverse(Q_LL) * Q_LM;
  } else {
    S = Mbasis * Q * Mbasis.transpose();
  }
  MyMatrix<T> F_M = F * Mbasis.transpose();
  std::vector<MyVector<T>> ListRay = PolyNorm_ConeExtremeRays(F_M);
  res.n_ray = ListRay.size();
  if (ListRay.empty()) {
    // C' = {0}, so C = L.
    res.second_order = true;
    return res;
  }
  int p = ListRay.size();
  MyMatrix<T> Qr(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      Qr(i, j) = ListRay[i].dot(S * ListRay[j]);
    }
  }
  std::optional<bool> opt = PolyNorm_IsStrictlyCopositive(Qr, 20);
  res.second_order = opt && *opt;
#ifdef DEBUG_POLYNORM_LATTICE_OPT
  os << "POLYNORM_LATTICE_OPT: KKT dL=" << dL << " dM=" << Mbasis.rows()
     << " p=" << p << " copositive=" << res.second_order << "\n";
#endif
  return res;
}

/*
  The certification of a lattice basis W*, for packing and for covering.
 */

struct PolyNormCertificate {
  bool exact_point;
  // Packing: every branch (choice of an active facet per contact) is a
  // strict local minimum. Covering: the polyhedral outer approximation
  // has a strict local maximum of det.
  bool local_optimum;
  // Covering: whether every deepest hole has n+1 lattice points, which the
  // validity of the certificate needs.
  bool holes_nondegenerate;
  // For packing, the number of branches checked.
  size_t n_branch;
  int n_contact;
  int n_constraint;
  int n_positive;
  int dim_tangent;
  std::string message;
};

template <typename T, typename Tint>
PolyNormCertificate
PolyNorm_CertifyPacking(MyMatrix<T> const &LmatDiff, MyMatrix<T> const &W,
                        size_t max_branch, std::ostream &os) {
  PolyNormCertificate cert{true, false, true, 0, 0, 0, 0, 0, ""};
  PolyNormLatticePackingEval<T, Tint> eval =
      PolyNorm_EvaluatePacking<T, Tint>(LmatDiff, W, os);
  // The contacts with all their active facets.
  std::vector<PolyNormContactFacets<T, Tint>> ListCF =
      PolyNorm_PackingNearContacts<T, Tint>(LmatDiff, W, eval.alpha, os);
  std::vector<PolyNormContactFacets<T, Tint>> ListContact;
  for (auto &eCF : ListCF) {
    if (eCF.val == eval.alpha) {
      ListContact.push_back(eCF);
    }
  }
  cert.n_contact = ListContact.size();
  size_t n_branch = 1;
  for (auto &eCF : ListContact) {
    n_branch *= eCF.ListFacet.size();
    if (n_branch > max_branch) {
      cert.message = "too many branches";
      return cert;
    }
  }
  cert.n_branch = n_branch;
  int n_c = ListContact.size();
  std::vector<size_t> sel(n_c, 0);
  bool all_ok = true;
  int min_positive = std::numeric_limits<int>::max();
  int max_tangent = 0;
  while (true) {
    std::vector<MyMatrix<T>> ListH;
    for (int u = 0; u < n_c; u++) {
      ListH.push_back(PolyNorm_PackingGradient<T, Tint>(
          LmatDiff, ListContact[u].z, ListContact[u].ListFacet[sel[u]]));
    }
    PolyNormKKT kkt = PolyNorm_CheckKKT(W, ListH, 1, os);
    cert.n_constraint = ListH.size();
    min_positive = std::min(min_positive, kkt.n_positive);
    max_tangent = std::max(max_tangent, kkt.dim_tangent);
    if (!kkt.first_order) {
      cert.message = "first order condition fails on a branch";
      all_ok = false;
      break;
    }
    if (!kkt.second_order) {
      cert.message = "second order condition fails or is degenerate on a "
                     "branch";
      all_ok = false;
      break;
    }
    int u = 0;
    while (u < n_c) {
      sel[u]++;
      if (sel[u] < ListContact[u].ListFacet.size()) {
        break;
      }
      sel[u] = 0;
      u++;
    }
    if (u == n_c) {
      break;
    }
  }
  cert.n_positive = min_positive;
  cert.dim_tangent = max_tangent;
  cert.local_optimum = all_ok;
  if (all_ok) {
    cert.message = "strict local optimum on every branch";
  }
  return cert;
}

template <typename T, typename Tint, typename Tgroup>
PolyNormCertificate
PolyNorm_CertifyCovering(PolyNormData<T, Tint, Tgroup> const &data,
                         MyMatrix<T> const &W, size_t max_system,
                         std::ostream &os) {
  int n = data.n;
  PolyNormCertificate cert{true, false, true, 1, 0, 0, 0, 0, ""};
  PolyNormLatticeCoveringEval<T, Tint> eval =
      PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, W, T(1), os);
  for (auto &hole : eval.ListHole) {
    if (hole.n_point != n + 1) {
      cert.holes_nondegenerate = false;
    }
  }
  cert.n_contact = eval.ListHole.size();
  std::vector<PolyNormLinCons<T>> ListCons =
      PolyNorm_CoveringConstraints<T, Tint>(data.Lmat, eval.ListHole, eval.mu,
                                            W, false, max_system);
  std::vector<MyMatrix<T>> ListH;
  for (auto &eCons : ListCons) {
    ListH.push_back(eCons.h);
  }
  cert.n_constraint = ListH.size();
  PolyNormKKT kkt = PolyNorm_CheckKKT(W, ListH, -1, os);
  cert.n_positive = kkt.n_positive;
  cert.dim_tangent = kkt.dim_tangent;
  if (!kkt.first_order) {
    cert.message = "first order condition fails";
    return cert;
  }
  if (!kkt.second_order) {
    cert.message = "second order condition fails or is degenerate";
    return cert;
  }
  cert.local_optimum = true;
  if (cert.holes_nondegenerate) {
    cert.message = "strict local optimum";
  } else {
    cert.message = "strict local optimum of the polyhedral model only, some "
                   "deepest hole has more than n+1 lattice points";
  }
  return cert;
}

/*
  The search.
 */

struct PolyNormLatticeOptParam {
  // "packing" or "covering".
  std::string problem;
  int n_start;
  uint64_t seed;
  int max_iter;
  // The trust region radius relative to the largest entry of W, at the
  // start and at which the search stops.
  double r_init;
  double r_min;
  // The contacts and holes within this ratio of the extreme value are
  // given to the linear program.
  double near_ratio;
  // The bound on the number of systems per hole.
  size_t max_system;
  // The bound on the number of branches of the packing certification.
  size_t max_branch;
  // When the linear programming steps stall, a pattern search along random
  // directions with exact evaluations is tried, from the radius r_init
  // down to r_pattern_min; the search then resumes. At most max_round such
  // rounds.
  double r_pattern_min;
  int max_round;
};

template <typename T> struct PolyNormLatticeResult {
  // The basis (rows) of the best lattice, normalized to a packing radius
  // or covering radius 1.
  MyMatrix<T> W;
  // The density, exact when the polished point is used.
  T density;
  double density_d;
  bool polished;
  // When the critical point could not be recognized as rational (it is
  // then typically algebraic irrational), the density at the floating
  // point critical point of the active constraints, if one is consistent
  // with the search.
  std::optional<double> density_refined;
  PolyNormCertificate cert;
};

// The cost to minimize: |det W| / alpha^n for packing, mu^n / |det W| for
// covering, from which the density is vol(P) / cost and vol(P) cost.
template <typename T, typename Tint, typename Tgroup>
T PolyNorm_LatticeCost(PolyNormData<T, Tint, Tgroup> const &data,
                       MyMatrix<T> const &W, bool is_packing,
                       std::ostream &os) {
  int n = data.n;
  if (is_packing) {
    PolyNormLatticePackingEval<T, Tint> eval =
        PolyNorm_EvaluatePacking<T, Tint>(data.LmatDiff, W, os);
    T pw(1);
    for (int i = 0; i < n; i++) {
      pw *= eval.alpha;
    }
    return eval.det / pw;
  }
  PolyNormLatticeCoveringEval<T, Tint> eval =
      PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, W, {}, os);
  T pw(1);
  for (int i = 0; i < n; i++) {
    pw *= eval.mu;
  }
  return pw / eval.det;
}

template <typename T>
T PolyNorm_DensityFromCost(T const &vol, T const &cost, bool is_packing) {
  if (is_packing) {
    return vol / cost;
  }
  return vol * cost;
}

// Reduce the basis with LLL for the Euclidean form and make det positive.
template <typename T, typename Tint>
MyMatrix<T> PolyNorm_ReduceBasis(MyMatrix<T> const &W, std::ostream &os) {
  LLLbasis<T, Tint> red = LLLbasisReduction<T, Tint>(W, os);
  MyMatrix<T> Wred = red.LattRed;
  if (DeterminantMat(Wred) < 0) {
    for (int b = 0; b < Wred.cols(); b++) {
      Wred(0, b) = -Wred(0, b);
    }
  }
  return Wred;
}

template <typename T>
MyMatrix<T> PolyNorm_RoundMatrix(MyMatrix<double> const &Wd, int e) {
  int n = Wd.rows();
  MyMatrix<T> W(n, n);
  for (int k = 0; k < n; k++) {
    for (int b = 0; b < n; b++) {
      W(k, b) = PolyNorm_RoundToGrid<T>(Wd(k, b), e);
    }
  }
  return W;
}

// The constraints of the linear program at W, with the scale value
// (alpha or mu) of W.
template <typename T, typename Tint, typename Tgroup>
std::pair<std::vector<PolyNormLinCons<T>>, T>
PolyNorm_SearchConstraints(PolyNormData<T, Tint, Tgroup> const &data,
                           MyMatrix<T> const &W, bool is_packing,
                           T const &near_ratio, T const &pair_tol,
                           size_t max_system, std::ostream &os) {
  std::vector<PolyNormLinCons<T>> ListCons;
  if (is_packing) {
    PolyNormLatticePackingEval<T, Tint> eval =
        PolyNorm_EvaluatePacking<T, Tint>(data.LmatDiff, W, os);
    T radius = eval.alpha / near_ratio;
    std::vector<PolyNormContactFacets<T, Tint>> ListCF =
        PolyNorm_PackingNearContacts<T, Tint>(data.LmatDiff, W, radius, os);
    for (auto &eCF : ListCF) {
      size_t idx = random_index(eCF.ListFacet.size());
      int j = eCF.ListFacet[idx];
      ListCons.push_back(
          {PolyNorm_PackingGradient<T, Tint>(data.LmatDiff, eCF.z, j),
           eval.alpha});
    }
    return {ListCons, eval.alpha};
  }
  PolyNormLatticeCoveringEval<T, Tint> eval =
      PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, W, near_ratio, os,
                                                 pair_tol);
  ListCons = PolyNorm_CoveringConstraints<T, Tint>(
      data.Lmat, eval.ListHole, eval.mu, W, false, max_system);
  return {ListCons, eval.mu};
}

// The equalities for the polishing: the constraints within the relative
// tolerance tol of being active, for the normalized lattice (radius 1).
template <typename T>
std::vector<PolyNormLinCons<T>>
PolyNorm_ActiveEqualities(MyMatrix<T> const &W,
                          std::vector<PolyNormLinCons<T>> const &ListCons,
                          T const &scale, double tol) {
  std::vector<PolyNormLinCons<T>> ListEq;
  std::set<std::vector<T>> set_h;
  for (auto &eCons : ListCons) {
    // The constraints are <h, W> >= alpha (packing, b > 0) and
    // <h, W> >= -mu (covering, b < 0, possibly below -mu for the nearly
    // tight systems). At the optimum, normalized to radius 1, the active
    // ones are the equalities <h, W> = +-1.
    T b_norm = eCons.b > 0 ? T(1) : T(-1);
    T gap_T = PolyNorm_MatScalar(eCons.h, W) - b_norm * scale;
    double gap = PolyNorm_ToDouble(gap_T);
    double scale_d = std::fabs(PolyNorm_ToDouble(scale));
    if (std::fabs(gap) <= tol * scale_d) {
      std::vector<T> V = PolyNorm_Key(PolyNorm_FlattenMat(eCons.h));
      if (set_h.count(V) > 0) {
        continue;
      }
      set_h.insert(V);
      ListEq.push_back({eCons.h, b_norm});
    }
  }
  return ListEq;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<PolyNormLatticeResult<T>>
PolyNorm_LatticeOptimizeFrom(PolyNormData<T, Tint, Tgroup> const &data,
                             T const &vol, MyMatrix<T> const &W_start,
                             PolyNormLatticeOptParam const &param,
                             std::ostream &os) {
  bool is_packing = param.problem == "packing";
  int s = is_packing ? 1 : -1;
  // Rational form of the near ratio, e.g. 0.95 -> 19/20.
  T near_ratio = *PolyNorm_RationalApprox<T>(param.near_ratio, 1.0e-9);
  MyMatrix<T> W = PolyNorm_ReduceBasis<T, Tint>(W_start, os);
  T cost = PolyNorm_LatticeCost(data, W, is_packing, os);
  double r = param.r_init;
  int n_accept = 0;
  int nn = data.n * data.n;
  int iter = 0;
  for (int i_round = 0; i_round < param.max_round; i_round++) {
    for (; iter < param.max_iter; iter++) {
      if (r < param.r_min) {
        break;
      }
      // The pairs nearly tight within a multiple of the trust region.
      double tol_d = std::min(1.0e-2, std::max(1.0e-7, 20 * r));
      T pair_tol = *PolyNorm_RationalApprox<T>(tol_d, tol_d * 1.0e-3);
      std::pair<std::vector<PolyNormLinCons<T>>, T> pair =
          PolyNorm_SearchConstraints(data, W, is_packing, near_ratio, pair_tol,
                                     param.max_system, os);
      MyMatrix<double> Wd = PolyNorm_ToDoubleMat(W);
      double scale = Wd.cwiseAbs().maxCoeff();
      double r_abs = r * scale;
      // The grid is about 10^{-4} of the trust region.
      int e = static_cast<int>(std::ceil(std::log2(1.0e4 / r_abs)));
      if (e < 0) {
        e = 0;
      }
      T r_T = PolyNorm_RoundToGrid<T>(r_abs, e + 4);
      std::optional<MyMatrix<T>> optX =
          PolyNorm_LinearStep(W, pair.first, s, r_T, os);
      if (!optX) {
#ifdef DEBUG_POLYNORM_LATTICE_OPT
        os << "POLYNORM_LATTICE_OPT: iter=" << iter << " the LP failed\n";
#endif
        r /= 4;
        continue;
      }
      MyMatrix<double> Wtry_d = Wd + PolyNorm_ToDoubleMat(*optX);
      MyMatrix<T> W_try = PolyNorm_RoundMatrix<T>(Wtry_d, e);
      if (DeterminantMat(W_try) <= 0) {
        r /= 4;
        continue;
      }
      T cost_try = PolyNorm_LatticeCost(data, W_try, is_packing, os);
      if (cost_try < cost) {
        W = W_try;
        cost = cost_try;
        n_accept++;
        r = std::min(2 * r, param.r_init);
        if (n_accept % 10 == 0) {
          W = PolyNorm_ReduceBasis<T, Tint>(W, os);
        }
      } else {
        r /= 4;
      }
#ifdef DEBUG_POLYNORM_LATTICE_OPT
      os << "POLYNORM_LATTICE_OPT: iter=" << iter << " r=" << r
         << " |cons|=" << pair.first.size() << " density="
         << PolyNorm_ToDouble(PolyNorm_DensityFromCost(vol, cost, is_packing))
         << "\n";
#endif
    }
    // The pattern search. The linear model can be wrong near holes whose
    // pinning points lie on lower dimensional faces of the translates, and
    // then every step is rejected although W is not a local optimum.
    bool improved = false;
    double rp = param.r_init;
    while (rp >= param.r_pattern_min && !improved) {
      MyMatrix<double> Wd = PolyNorm_ToDoubleMat(W);
      double scale = Wd.cwiseAbs().maxCoeff();
      double r_abs = rp * scale;
      int e = static_cast<int>(std::ceil(std::log2(1.0e4 / r_abs)));
      if (e < 0) {
        e = 0;
      }
      for (int i_try = 0; i_try < 2 * nn && !improved; i_try++) {
        MyMatrix<double> D(data.n, data.n);
        for (int k = 0; k < data.n; k++) {
          for (int b = 0; b < data.n; b++) {
            D(k, b) = random_int(-1000, 1000) / 1000.0;
          }
        }
        for (double sign : {1.0, -1.0}) {
          MyMatrix<T> W_try = PolyNorm_RoundMatrix<T>(Wd + sign * r_abs * D, e);
          if (DeterminantMat(W_try) <= 0) {
            continue;
          }
          T cost_try = PolyNorm_LatticeCost(data, W_try, is_packing, os);
          if (cost_try < cost) {
            W = W_try;
            cost = cost_try;
            improved = true;
            break;
          }
        }
      }
      rp /= 2;
    }
#ifdef DEBUG_POLYNORM_LATTICE_OPT
    os << "POLYNORM_LATTICE_OPT: round=" << i_round << " pattern search "
       << (improved ? "improved" : "failed") << " density="
       << PolyNorm_ToDouble(PolyNorm_DensityFromCost(vol, cost, is_packing))
       << "\n";
#endif
    if (!improved) {
      break;
    }
    W = PolyNorm_ReduceBasis<T, Tint>(W, os);
    r = std::max(2 * rp, param.r_min * 4);
  }
  // Normalization to radius 1.
  T scale_val;
  if (is_packing) {
    scale_val = PolyNorm_EvaluatePacking<T, Tint>(data.LmatDiff, W, os).alpha;
  } else {
    scale_val = PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, W, {}, os).mu;
  }
  PolyNormLatticeResult<T> res;
  res.W = W / scale_val;
  res.density = PolyNorm_DensityFromCost(vol, cost, is_packing);
  res.density_d = PolyNorm_ToDouble(res.density);
  res.polished = false;
  res.cert = {false, false, true, 0, 0, 0, 0, 0, "not polished"};
  // Polishing, with decreasing tolerances on the active set. The same
  // tolerance is used for the nearly tight pairs of the holes.
  for (double tol : {1.0e-4, 1.0e-5, 1.0e-6, 1.0e-3, 1.0e-7}) {
    T pair_tol = *PolyNorm_RationalApprox<T>(tol, tol * 1.0e-3);
    std::pair<std::vector<PolyNormLinCons<T>>, T> pair =
        PolyNorm_SearchConstraints(data, W, is_packing, near_ratio, pair_tol,
                                   param.max_system, os);
    std::vector<PolyNormLinCons<T>> ListEq =
        PolyNorm_ActiveEqualities(W, pair.first, pair.second, tol);
    for (double tol_rat : {1.0e-9, 1.0e-11, 1.0e-7}) {
      double det_numeric;
      std::optional<MyMatrix<T>> optW =
          PolyNorm_Polish(res.W, ListEq, tol_rat, det_numeric, os);
      if (det_numeric > 0) {
        // Both problems are normalized to radius 1, so the density is
        // vol(P) / |det W|. The floating point critical point is kept
        // only when it is consistent with the search and not worse.
        double dens_num = PolyNorm_ToDouble(vol) / det_numeric;
        double rel = std::fabs(dens_num - res.density_d) / res.density_d;
        bool not_worse = is_packing ? dens_num >= res.density_d * (1 - 1.0e-12)
                                    : dens_num <= res.density_d * (1 + 1.0e-12);
        if (rel < 1.0e-5 && not_worse) {
          bool better = !res.density_refined ||
                        (is_packing ? dens_num > *res.density_refined
                                    : dens_num < *res.density_refined);
          if (better) {
            res.density_refined = dens_num;
          }
        }
      }
      if (!optW) {
        continue;
      }
      MyMatrix<T> Wp = *optW;
      if (DeterminantMat(Wp) < 0) {
        continue;
      }
      T cost_p = PolyNorm_LatticeCost(data, Wp, is_packing, os);
      if (cost_p > cost) {
        // The polished point must be at least as good as the search.
        continue;
      }
      res.W = Wp;
      res.density = PolyNorm_DensityFromCost(vol, cost_p, is_packing);
      res.density_d = PolyNorm_ToDouble(res.density);
      res.polished = true;
      res.density_refined.reset();
      // Normalize to radius 1 (the equalities already do that, but the
      // radius is recomputed for safety).
      T sc;
      if (is_packing) {
        sc = PolyNorm_EvaluatePacking<T, Tint>(data.LmatDiff, Wp, os).alpha;
      } else {
        sc = PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, Wp, {}, os).mu;
      }
      res.W = Wp / sc;
      if (is_packing) {
        res.cert = PolyNorm_CertifyPacking<T, Tint>(data.LmatDiff, res.W,
                                                    param.max_branch, os);
      } else {
        res.cert = PolyNorm_CertifyCovering<T, Tint, Tgroup>(
            data, res.W, param.max_system, os);
      }
      return res;
    }
  }
  return res;
}

// A random starting basis: integer entries in [-8, 8] divided by 8, with
// nonzero determinant.
template <typename T> MyMatrix<T> PolyNorm_RandomBasis(int n) {
  while (true) {
    MyMatrix<T> W(n, n);
    for (int k = 0; k < n; k++) {
      for (int b = 0; b < n; b++) {
        W(k, b) = T(random_int(-8, 8)) / T(8);
      }
    }
    T det = DeterminantMat(W);
    if (det != 0) {
      if (det < 0) {
        for (int b = 0; b < n; b++) {
          W(0, b) = -W(0, b);
        }
      }
      return W;
    }
  }
}

template <typename T, typename Tint, typename Tgroup>
std::vector<PolyNormLatticeResult<T>>
PolyNorm_LatticeOptimize(PolyNormData<T, Tint, Tgroup> const &data,
                         PolyNormLatticeOptParam const &param,
                         std::vector<MyMatrix<T>> const &ListStart,
                         std::ostream &os) {
  T vol = rev_search::Kernel_VolumePolytope(data.EXT);
  set_random_seed(param.seed);
  std::vector<PolyNormLatticeResult<T>> ListRes;
  std::vector<MyMatrix<T>> ListW = ListStart;
  for (int i = 0; i < param.n_start; i++) {
    ListW.push_back(PolyNorm_RandomBasis<T>(data.n));
  }
  for (auto &W0 : ListW) {
#ifdef TIMINGS_POLYNORM_LATTICE_OPT
    MicrosecondTime time;
#endif
    std::optional<PolyNormLatticeResult<T>> opt =
        PolyNorm_LatticeOptimizeFrom(data, vol, W0, param, os);
#ifdef TIMINGS_POLYNORM_LATTICE_OPT
    os << "|POLYNORM_LATTICE_OPT: one start|=" << time << "\n";
#endif
    if (opt) {
      ListRes.push_back(*opt);
    }
  }
  return ListRes;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_LATTICEOPT_H_
// clang-format on
