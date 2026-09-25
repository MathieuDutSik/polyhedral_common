// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_COVERING_H_
#define SRC_POLYNORM_POLYNORM_COVERING_H_

// clang-format off
#include "PolyNorm_Basic.h"
#include "POLY_LinearProgramming.h"
#include <algorithm>
#include <optional>
#include <map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_COVERING
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_COVERING
#endif

#ifdef TIMINGS
#define TIMINGS_POLYNORM_COVERING
#endif

/*
  The covering problem.

  With g_z(p) = ||p - z||_P and f(p) = min_{z in Z^n} g_z(p), the covering
  radius is mu(P) = max_p f(p): the point p is covered by mu P + z exactly
  when g_z(p) <= mu. A point attaining the maximum is a last-covered point
  (Cslovjecsek, Malikiosis, Naszodi, Schymura, Combinatorica 2022) and the
  maximal empty translate p - mu P around it is the analogue of a Delaunay
  polytope: its interior contains no lattice point and the lattice points on
  its boundary pin it.

  The function f is the minimum of the convex piecewise linear functions
  g_z, so the feasible region {(p, mu) : mu <= f(p)} is an intersection of
  unions of halfspaces, one union per lattice point z (mu <= ell_i(p - z)
  for at least one facet i). That is a disjunctive program and it is solved
  exactly by branch and bound:
  --- A node is a set C of pairs (z, i) meaning that p lies in the cone of
      the facet i at z (ell_i(p - z) >= ell_j(p - z) for all j) and that
      mu <= ell_i(p - z). Its linear program maximizes mu under these
      constraints, which is an upper bound for the node.
  --- At the optimum (p*, mu*), the value f(p*) is a lower bound on mu(P)
      and if f(p*) = mu* the node is solved. Otherwise some lattice point z'
      has g_{z'}(p*) < mu* and the node is split into the m children
      C + (z', i), i in [m], one per facet cone of z'.
  --- The translation symmetry is used by requiring 0 to be a nearest
      lattice point of p, which is no restriction: the root nodes are
      {(0, i)} together with the constraint p in mu P.
  --- The linear symmetries {A in GL_n(Z) : PA = P} reduce the root nodes to
      one per orbit of facets and, when the group is listed, prune every
      node equivalent to an already processed one under the affine group,
      through a canonical form of C.
  The lattice points that can block a point p with g_0(p) = mu are the ones
  of the translate p - mu P, which are enumerated at each node; they all lie
  in mu0 (P - P) for an upper bound mu0 of mu(P), so the search is finite.
  The linear programs are solved exactly over the rationals.

  For cross-checking, PolyNorm_CoveringBruteForce enumerates the systems of
  the paper: n+1 pairs (z_j, i_j) with affinely independent normals, one of
  them at z = 0, solved for (p, mu) and kept when the translate is empty.
 */

// A lattice point and a facet of P: the constraint that p lies in the cone
// of the facet i at z, and that the translate mu P + z is blocked by that
// facet at p.
template <typename T, typename Tint> struct PolyNormPair {
  MyVector<Tint> z;
  MyVector<T> z_T;
  int i;
};

template <typename T, typename Tint> struct PolyNormCovering {
  T mu;
  // A last-covered point, in the recentered coordinates.
  MyVector<T> p;
  // The pairs (z, i) with ell_i(p - z) = mu: the lattice points z on the
  // boundary of the empty translate p - mu P, with the facets they touch.
  std::vector<PolyNormPair<T, Tint>> ListTight;
  // The upper bound that bounds the root linear programs.
  T mu0;
  // Search statistics.
  size_t n_node;
  size_t n_lp;
  size_t n_skip_canonical;
};

// The smallest t such that a translate of the polytope with vertices ListQ
// fits in tP, that is min over s of max_k ||q_k - s||_P, by the linear
// program in (s, t) with the constraints t >= ell_j(q_k - s).
template <typename T>
T PolyNorm_TranslatedFit(MyMatrix<T> const &Lmat,
                         std::vector<MyVector<T>> const &ListQ,
                         std::ostream &os) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  int n_q = ListQ.size();
  MyMatrix<T> ListIneq = ZeroMatrix<T>(m * n_q, n + 2);
  int pos = 0;
  for (auto &q : ListQ) {
    for (int j = 0; j < m; j++) {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += Lmat(j, k) * q(k);
        ListIneq(pos, k + 1) = Lmat(j, k);
      }
      ListIneq(pos, 0) = -scal;
      ListIneq(pos, n + 1) = 1;
      pos++;
    }
  }
  MyVector<T> ToBeMinimized = ZeroVector<T>(n + 2);
  ToBeMinimized(n + 1) = 1;
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
#ifdef SANITY_CHECK_POLYNORM_COVERING
  if (!eSol.DirectSolution || !eSol.DualSolution) {
    std::cerr << "POLYNORM_COVERING: The translated fit LP should have an "
              << "optimal solution\n";
    throw TerminalException{1};
  }
#endif
  return eSol.OptimalValue;
}

// An upper bound on the covering radius from lattice polytopes Q of known
// covering radius: if a translate of Q is contained in tP then
// mu(P) <= t mu(Q), since the covering radius is monotone under inclusion
// and invariant under translation. The polytopes used are the simplices
// conv(0, +-e_1, ..., +-e_n) (mu = n), the cross polytope conv(+-e_j)
// (mu = n/2) and the cube [-1,1]^n (mu = 1/2). The fitting factor t is
// the translated fit of the vertices.
template <typename T>
T PolyNorm_CoveringUpperBound(MyMatrix<T> const &Lmat, std::ostream &os) {
  int n = Lmat.cols();
  std::optional<T> bound;
  auto update = [&](T const &val) -> void {
    if (!bound || val < *bound) {
      bound = val;
    }
  };
  // The cross polytope.
  {
    std::vector<MyVector<T>> ListQ;
    for (int j = 0; j < n; j++) {
      MyVector<T> e = ZeroVector<T>(n);
      e(j) = 1;
      ListQ.push_back(e);
      e(j) = -1;
      ListQ.push_back(e);
    }
    T t = PolyNorm_TranslatedFit(Lmat, ListQ, os);
    update(t * T(n) / T(2));
  }
  // The cube, with 2^n vertices.
  if (n <= 12) {
    std::vector<MyVector<T>> ListQ;
    size_t n_sign = size_t(1) << n;
    for (size_t s = 0; s < n_sign; s++) {
      MyVector<T> v(n);
      for (int j = 0; j < n; j++) {
        if ((s >> j) & 1) {
          v(j) = 1;
        } else {
          v(j) = -1;
        }
      }
      ListQ.push_back(v);
    }
    T t = PolyNorm_TranslatedFit(Lmat, ListQ, os);
    update(t / T(2));
  }
  // The simplices, one per choice of signs (all of them for small n,
  // the two constant ones otherwise).
  {
    std::vector<size_t> ListSign;
    if (n <= 8) {
      size_t n_sign = size_t(1) << n;
      for (size_t s = 0; s < n_sign; s++) {
        ListSign.push_back(s);
      }
    } else {
      ListSign.push_back(0);
      ListSign.push_back((size_t(1) << n) - 1);
    }
    for (auto &s : ListSign) {
      std::vector<MyVector<T>> ListQ;
      ListQ.push_back(ZeroVector<T>(n));
      for (int j = 0; j < n; j++) {
        MyVector<T> e = ZeroVector<T>(n);
        if ((s >> j) & 1) {
          e(j) = 1;
        } else {
          e(j) = -1;
        }
        ListQ.push_back(e);
      }
      T t = PolyNorm_TranslatedFit(Lmat, ListQ, os);
      update(t * T(n));
    }
  }
#ifdef DEBUG_POLYNORM_COVERING
  os << "POLYNORM_COVERING: upper bound mu0=" << *bound << "\n";
#endif
  return *bound;
}

// The lattice points z with g_z(p) <= mu, that is the lattice points of the
// closed translate p - mu P: the inequalities are mu - ell_j(p - z) >= 0.
template <typename T, typename Tint>
std::vector<MyVector<Tint>> PolyNorm_BlockingPoints(MyMatrix<T> const &Lmat,
                                                    MyVector<T> const &p,
                                                    T const &mu,
                                                    std::ostream &os) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  MyMatrix<T> FAC(m, n + 1);
  for (int j = 0; j < m; j++) {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += Lmat(j, k) * p(k);
      FAC(j, k + 1) = Lmat(j, k);
    }
    FAC(j, 0) = mu - scal;
  }
  return GetListIntegralPoint<T, Tint>(FAC, os);
}

// The linear program of a node: maximize mu over (p, mu) in R^{n+1} with
// mu <= mu0, p in mu P and for each (z, i) in C the cone and blocking
// constraints. The rows follow the convention of the simplex code: the
// vector (b, a_p, a_mu) means b + a_p . p + a_mu mu >= 0. Returns the
// optimal (p, mu) when feasible.
template <typename T, typename Tint>
std::optional<std::pair<MyVector<T>, T>>
PolyNorm_SolveNodeLP(MyMatrix<T> const &Lmat, T const &mu0,
                     std::vector<PolyNormPair<T, Tint>> const &C,
                     std::ostream &os) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  int n_row = 1 + m + m * C.size();
  MyMatrix<T> ListIneq = ZeroMatrix<T>(n_row, n + 2);
  int pos = 0;
  // mu0 - mu >= 0
  ListIneq(pos, 0) = mu0;
  ListIneq(pos, n + 1) = -1;
  pos++;
  // mu - ell_j(p) >= 0
  for (int j = 0; j < m; j++) {
    for (int k = 0; k < n; k++) {
      ListIneq(pos, k + 1) = -Lmat(j, k);
    }
    ListIneq(pos, n + 1) = 1;
    pos++;
  }
  for (auto &ePair : C) {
    int i = ePair.i;
    std::vector<T> ell_z(m);
    for (int j = 0; j < m; j++) {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += Lmat(j, k) * ePair.z_T(k);
      }
      ell_z[j] = scal;
    }
    // ell_i(p) - ell_i(z) - mu >= 0
    ListIneq(pos, 0) = -ell_z[i];
    for (int k = 0; k < n; k++) {
      ListIneq(pos, k + 1) = Lmat(i, k);
    }
    ListIneq(pos, n + 1) = -1;
    pos++;
    // ell_i(p - z) - ell_j(p - z) >= 0
    for (int j = 0; j < m; j++) {
      if (j != i) {
        ListIneq(pos, 0) = ell_z[j] - ell_z[i];
        for (int k = 0; k < n; k++) {
          ListIneq(pos, k + 1) = Lmat(i, k) - Lmat(j, k);
        }
        pos++;
      }
    }
  }
  MyVector<T> ToBeMinimized = ZeroVector<T>(n + 2);
  ToBeMinimized(n + 1) = -1;
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
  if (!eSol.DirectSolution) {
    // Infeasible
    return {};
  }
#ifdef SANITY_CHECK_POLYNORM_COVERING
  if (!eSol.DualSolution) {
    std::cerr << "POLYNORM_COVERING: The node LP is unbounded, which the "
              << "constraint mu <= mu0 forbids\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> const &x = *eSol.DirectSolution;
  MyVector<T> p(n);
  for (int k = 0; k < n; k++) {
    p(k) = x(k);
  }
  T mu = x(n);
  return std::pair<MyVector<T>, T>{std::move(p), std::move(mu)};
}

template <typename T, typename Tint, typename Tgroup>
struct PolyNormCoveringSolver {
  PolyNormData<T, Tint, Tgroup> const &data;
  std::ostream &os;
  int n;
  int m;
  T mu0;
  T mu_best;
  MyVector<T> p_best;
  size_t n_node;
  size_t n_lp;
  size_t n_skip_canonical;
  size_t n_canonical;
  // The orbit of each facet under the symmetry group, for the invariants.
  std::vector<int> facet_orbit;
  // The relative position of two facets, facet_pos(i, i') = ell_i(c_i')
  // with c_i' the isobarycenter of the vertices of the facet i'. It is
  // invariant under the symmetries since they map facet isobarycenters
  // to facet isobarycenters.
  MyMatrix<T> facet_pos;
  // When the group elements are listed: for each pair of facets (i, i'),
  // the indices of the elements whose facet permutation maps i to i'. An
  // equivalence between two nodes must map a chosen pair (z, i) of the
  // first to some pair (z', i') of the second, which leaves only these
  // elements to try.
  std::vector<std::vector<std::vector<size_t>>> elt_by_facet_image;
  // The nodes processed so far, keyed by their invariant. Equivalence is
  // tested only among nodes with the same invariant.
  std::map<std::vector<T>, std::vector<std::vector<PolyNormPair<T, Tint>>>>
      map_visited;

  PolyNormCoveringSolver(PolyNormData<T, Tint, Tgroup> const &_data,
                         std::ostream &_os)
      : data(_data), os(_os), n(_data.n), m(_data.Lmat.rows()),
        mu_best(0), n_node(0), n_lp(0), n_skip_canonical(0),
        n_canonical(0) {
    mu0 = PolyNorm_CoveringUpperBound(data.Lmat, os);
    p_best = ZeroVector<T>(n);
    facet_orbit = PolyNorm_FacetOrbitIds(data.ListGen, m);
    facet_pos = MyMatrix<T>(m, m);
    for (int i2 = 0; i2 < m; i2++) {
      Face const &f = data.ListIncd[i2];
      MyVector<T> cent = ZeroVector<T>(n);
      boost::dynamic_bitset<>::size_type k = f.find_first();
      while (k != boost::dynamic_bitset<>::npos) {
        for (int j = 0; j < n; j++) {
          cent(j) += data.EXT(k, j + 1);
        }
        k = f.find_next(k);
      }
      cent /= T(int(f.count()));
      for (int i1 = 0; i1 < m; i1++) {
        T scal(0);
        for (int j = 0; j < n; j++) {
          scal += data.Lmat(i1, j) * cent(j);
        }
        facet_pos(i1, i2) = scal;
      }
    }
    if (data.ListElt) {
      elt_by_facet_image.resize(m);
      for (int i = 0; i < m; i++) {
        elt_by_facet_image[i].resize(m);
      }
      std::vector<PolyNormSymm<Tint>> const &ListElt = *data.ListElt;
      for (size_t u = 0; u < ListElt.size(); u++) {
        for (int i = 0; i < m; i++) {
          elt_by_facet_image[i][ListElt[u].facet_perm[i]].push_back(u);
        }
      }
    }
#ifdef DEBUG_POLYNORM_COVERING
    os << "POLYNORM_COVERING: mu0=" << mu0 << "\n";
#endif
  }

  // The value f(p) = min_z g_z(p) at a point p with g_0(p) = mu, together
  // with a lattice point attaining it. Only the lattice points of the
  // translate p - mu P can do better than 0, so they are the ones
  // enumerated.
  std::pair<T, MyVector<Tint>> EvaluateF(MyVector<T> const &p,
                                         T const &mu) const {
    std::vector<MyVector<Tint>> ListZ =
        PolyNorm_BlockingPoints<T, Tint>(data.Lmat, p, mu, os);
    T f_val = mu;
    MyVector<Tint> z_min = ZeroVector<Tint>(n);
    MyVector<T> diff(n);
    for (auto &z : ListZ) {
      diff = p - UniversalVectorConversion<T, Tint>(z);
      T val = PolyNorm_Gauge(data.Lmat, diff);
      if (val < f_val) {
        f_val = val;
        z_min = z;
      }
    }
    return {f_val, z_min};
  }

  // The image of C under the symmetry of index u followed by the
  // translation making the pair a land on the lattice point z_target, as
  // a sorted list of (z, i) flattened for comparison.
  std::vector<Tint> ImageSorted(std::vector<PolyNormPair<T, Tint>> const &C,
                                size_t u, size_t a,
                                MyVector<Tint> const &z_target) const {
    PolyNormSymm<Tint> const &eElt = (*data.ListElt)[u];
    size_t len = C.size();
    MyVector<Tint> shift = z_target - eElt.A.transpose() * C[a].z;
    std::vector<std::vector<Tint>> ListImg(len);
    for (size_t b = 0; b < len; b++) {
      MyVector<Tint> zA = eElt.A.transpose() * C[b].z + shift;
      std::vector<Tint> V(n + 1);
      for (int k = 0; k < n; k++) {
        V[k] = zA(k);
      }
      V[n] = eElt.facet_perm[C[b].i];
      ListImg[b] = std::move(V);
    }
    std::sort(ListImg.begin(), ListImg.end());
    std::vector<Tint> V;
    V.reserve(len * (n + 1));
    for (auto &eImg : ListImg) {
      V.insert(V.end(), eImg.begin(), eImg.end());
    }
    return V;
  }

  // Whether C1 and C2 are equivalent under the affine group generated by
  // the symmetries and the translations. The pair 0 of C1 must be mapped
  // to some pair b of C2 in the same facet orbit, by an element mapping
  // the facet of the one to the facet of the other; the translation is
  // then determined.
  bool IsEquivalent(std::vector<PolyNormPair<T, Tint>> const &C1,
                    std::vector<PolyNormPair<T, Tint>> const &C2) const {
    if (C1.size() != C2.size()) {
      return false;
    }
    MyVector<Tint> zero = ZeroVector<Tint>(n);
    std::vector<Tint> V2 = ImageSorted(C2, 0, 0, zero);
    // The identity is at some index; the image by it of C2 normalized at
    // its pair 0 is what the images of C1 are compared to, so C2 is sent
    // to its own normalization by the identity. Instead of locating the
    // identity, C2 is normalized directly.
    std::vector<std::vector<Tint>> ListImg2;
    for (size_t b = 0; b < C2.size(); b++) {
      std::vector<Tint> V(n + 1);
      for (int k = 0; k < n; k++) {
        V[k] = C2[b].z(k) - C2[0].z(k);
      }
      V[n] = C2[b].i;
      ListImg2.push_back(std::move(V));
    }
    std::sort(ListImg2.begin(), ListImg2.end());
    V2.clear();
    for (auto &eImg : ListImg2) {
      V2.insert(V2.end(), eImg.begin(), eImg.end());
    }
    int i1 = C1[0].i;
    for (size_t b = 0; b < C2.size(); b++) {
      int i2 = C2[b].i;
      if (facet_orbit[i1] != facet_orbit[i2]) {
        continue;
      }
      MyVector<Tint> z_target = C2[b].z - C2[0].z;
      for (auto &u : elt_by_facet_image[i1][i2]) {
        std::vector<Tint> V1 = ImageSorted(C1, u, 0, z_target);
        if (V1 == V2) {
          return true;
        }
      }
    }
    return false;
  }

  // An invariant of C under the symmetries and the translations: for each
  // pair a, the orbit of its facet and the sorted list over the other pairs
  // b of the relative positions of the facets i_a and i_b and of the values
  // of ell_{i_a}, ell_{i_b} and of the gauge on z_b - z_a and on z_a - z_b;
  // the descriptors are sorted. The values are invariant
  // since ell_{sigma(i)}(xA) = ell_i(x) for a symmetry A with facet
  // permutation sigma, and differences of lattice points do not see the
  // translations.
  std::vector<T> Invariant(std::vector<PolyNormPair<T, Tint>> const &C) const {
    size_t len = C.size();
    auto ell = [&](int i, MyVector<T> const &x) -> T {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += data.Lmat(i, k) * x(k);
      }
      return scal;
    };
    std::vector<std::vector<T>> ListDesc(len);
    MyVector<T> diff(n);
    for (size_t a = 0; a < len; a++) {
      std::vector<std::vector<T>> ListVal;
      for (size_t b = 0; b < len; b++) {
        if (b != a) {
          diff = C[b].z_T - C[a].z_T;
          std::vector<T> val{facet_pos(C[a].i, C[b].i),
                             facet_pos(C[b].i, C[a].i),
                             ell(C[a].i, diff), ell(C[b].i, diff),
                             PolyNorm_Gauge(data.Lmat, diff)};
          diff = -diff;
          val.push_back(ell(C[a].i, diff));
          val.push_back(ell(C[b].i, diff));
          val.push_back(PolyNorm_Gauge(data.Lmat, diff));
          ListVal.push_back(std::move(val));
        }
      }
      std::sort(ListVal.begin(), ListVal.end());
      std::vector<T> desc;
      desc.push_back(T(facet_orbit[C[a].i]));
      for (auto &val : ListVal) {
        desc.insert(desc.end(), val.begin(), val.end());
      }
      ListDesc[a] = std::move(desc);
    }
    std::sort(ListDesc.begin(), ListDesc.end());
    std::vector<T> V;
    for (auto &desc : ListDesc) {
      V.insert(V.end(), desc.begin(), desc.end());
    }
    return V;
  }

  // Whether the node is new, that is not equivalent to a processed one.
  // Records it when it is. The equivalence is tested only against the
  // nodes with the same invariant.
  bool IsNewNode(std::vector<PolyNormPair<T, Tint>> const &C) {
    if (!data.ListElt) {
      return true;
    }
    std::vector<T> inv = Invariant(C);
    std::vector<std::vector<PolyNormPair<T, Tint>>> &ListEntry =
        map_visited[std::move(inv)];
    for (auto &eEntry : ListEntry) {
      n_canonical++;
      if (IsEquivalent(C, eEntry)) {
        n_skip_canonical++;
        return false;
      }
    }
    ListEntry.push_back(C);
    return true;
  }

  void Process(std::vector<PolyNormPair<T, Tint>> const &C,
               MyVector<T> const &p_star, T const &mu_star) {
    n_node++;
#ifdef DEBUG_POLYNORM_COVERING
    os << "POLYNORM_COVERING: node " << n_node << " |C|=" << C.size()
       << " mu*=" << mu_star << " mu_best=" << mu_best << "\n";
#endif
    if (mu_star <= mu_best) {
      return;
    }
    std::pair<T, MyVector<Tint>> pair = EvaluateF(p_star, mu_star);
    T const &f_val = pair.first;
    if (f_val > mu_best) {
      mu_best = f_val;
      p_best = p_star;
#ifdef DEBUG_POLYNORM_COVERING
      os << "POLYNORM_COVERING: new incumbent mu_best=" << mu_best << "\n";
#endif
    }
    if (f_val == mu_star) {
      // The bound is attained: the node is solved.
      return;
    }
    // Branching on the most violated lattice point.
    MyVector<Tint> const &z = pair.second;
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    struct Child {
      std::vector<PolyNormPair<T, Tint>> C;
      MyVector<T> p;
      T mu;
    };
    std::vector<Child> ListChild;
    for (int i = 0; i < m; i++) {
      std::vector<PolyNormPair<T, Tint>> Cnew = C;
      Cnew.push_back({z, z_T, i});
      if (!IsNewNode(Cnew)) {
        continue;
      }
      n_lp++;
      std::optional<std::pair<MyVector<T>, T>> opt =
          PolyNorm_SolveNodeLP(data.Lmat, mu0, Cnew, os);
      if (!opt) {
        continue;
      }
      if (opt->second <= mu_best) {
        continue;
      }
      ListChild.push_back({std::move(Cnew), std::move(opt->first),
                           std::move(opt->second)});
    }
    // Best bound first, so that good incumbents are found early.
    std::stable_sort(ListChild.begin(), ListChild.end(),
                     [](Child const &a, Child const &b) -> bool {
                       return a.mu > b.mu;
                     });
    for (auto &eChild : ListChild) {
      Process(eChild.C, eChild.p, eChild.mu);
    }
  }

  PolyNormCovering<T, Tint> Run() {
#ifdef TIMINGS_POLYNORM_COVERING
    MicrosecondTime time;
#endif
    std::vector<int> ListRep = PolyNorm_FacetOrbitRepresentatives(data);
#ifdef DEBUG_POLYNORM_COVERING
    os << "POLYNORM_COVERING: m=" << m << " facet orbits=" << ListRep.size()
       << "\n";
#endif
    MyVector<Tint> zero = ZeroVector<Tint>(n);
    MyVector<T> zero_T = ZeroVector<T>(n);
    for (auto &i : ListRep) {
      std::vector<PolyNormPair<T, Tint>> C{{zero, zero_T, i}};
      if (!IsNewNode(C)) {
        continue;
      }
      n_lp++;
      std::optional<std::pair<MyVector<T>, T>> opt =
          PolyNorm_SolveNodeLP(data.Lmat, mu0, C, os);
#ifdef SANITY_CHECK_POLYNORM_COVERING
      if (!opt) {
        std::cerr << "POLYNORM_COVERING: The root node of facet " << i
                  << " is infeasible\n";
        throw TerminalException{1};
      }
#endif
      Process(C, opt->first, opt->second);
    }
#ifdef TIMINGS_POLYNORM_COVERING
    os << "|POLYNORM_COVERING: branch and bound|=" << time << "\n";
#endif
#ifdef DEBUG_POLYNORM_COVERING
    os << "POLYNORM_COVERING: n_node=" << n_node << " n_lp=" << n_lp
       << " n_skip_canonical=" << n_skip_canonical
       << " n_equivalence_test=" << n_canonical << "\n";
#endif
    // The lattice points on the boundary of the empty translate, with the
    // facets they touch: the facets attaining the gauge, for the lattice
    // points whose gauge is the covering radius.
    std::vector<PolyNormPair<T, Tint>> ListTight;
    std::vector<MyVector<Tint>> ListZ =
        PolyNorm_BlockingPoints<T, Tint>(data.Lmat, p_best, mu_best, os);
    MyVector<T> diff(n);
    for (auto &z : ListZ) {
      MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
      diff = p_best - z_T;
      if (PolyNorm_Gauge(data.Lmat, diff) == mu_best) {
        for (int i = 0; i < m; i++) {
          T val(0);
          for (int k = 0; k < n; k++) {
            val += data.Lmat(i, k) * diff(k);
          }
          if (val == mu_best) {
            ListTight.push_back({z, z_T, i});
          }
        }
      }
    }
    return {mu_best, p_best, std::move(ListTight), mu0,
            n_node,  n_lp,   n_skip_canonical};
  }
};

template <typename T, typename Tint, typename Tgroup>
PolyNormCovering<T, Tint>
ComputePolyNormCovering(PolyNormData<T, Tint, Tgroup> const &data,
                        std::ostream &os) {
  PolyNormCoveringSolver<T, Tint, Tgroup> solver(data, os);
  return solver.Run();
}

// The enumeration of the paper, in the setting of this file: the last
// covered point p is translated so that 0 is a nearest lattice point, and
// the n+1 pairs (z_j, i_j) determining it have z_1 = 0 and the other z_j
// among the lattice points of mu0 (P - P). Every such system with affinely
// independent normals is solved and the solution kept when p is in mu P and
// no translate mu P + z contains p in its interior. The number of systems is
// m C(|L| m - 1, n), so this is only for small cases.
template <typename T, typename Tint, typename Tgroup>
std::pair<T, MyVector<T>>
PolyNorm_CoveringBruteForce(PolyNormData<T, Tint, Tgroup> const &data,
                            std::ostream &os) {
  int n = data.n;
  int m = data.Lmat.rows();
  T mu0 = PolyNorm_CoveringUpperBound(data.Lmat, os);
  std::vector<MyVector<Tint>> ListL =
      PolyNorm_LatticePoints<T, Tint>(data.LmatDiff, mu0, os);
  std::vector<MyVector<T>> ListL_T;
  for (auto &z : ListL) {
    ListL_T.push_back(UniversalVectorConversion<T, Tint>(z));
  }
  size_t n_latt = ListL.size();
  // The pairs, q = u * m + i for the lattice point u and the facet i, with
  // the value ell_i(z_u).
  size_t n_pair = n_latt * m;
  std::vector<T> ell_pair(n_pair);
  for (size_t u = 0; u < n_latt; u++) {
    for (int i = 0; i < m; i++) {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += data.Lmat(i, k) * ListL_T[u](k);
      }
      ell_pair[u * m + i] = scal;
    }
  }
  size_t u_zero = std::numeric_limits<size_t>::max();
  for (size_t u = 0; u < n_latt; u++) {
    if (IsZeroVector(ListL[u])) {
      u_zero = u;
    }
  }
#ifdef SANITY_CHECK_POLYNORM_COVERING
  if (u_zero == std::numeric_limits<size_t>::max()) {
    std::cerr << "POLYNORM_COVERING: The origin is missing from L\n";
    throw TerminalException{1};
  }
#endif
  auto is_valid = [&](MyVector<T> const &p, T const &mu) -> bool {
    if (mu > mu0) {
      return false;
    }
    for (int j = 0; j < m; j++) {
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += data.Lmat(j, k) * p(k);
      }
      if (scal > mu) {
        return false;
      }
    }
    MyVector<T> diff(n);
    for (size_t u = 0; u < n_latt; u++) {
      diff = p - ListL_T[u];
      if (PolyNorm_Gauge(data.Lmat, diff) < mu) {
        return false;
      }
    }
    return true;
  };
  T mu_max(0);
  MyVector<T> p_max = ZeroVector<T>(n);
  [[maybe_unused]] size_t n_system = 0;
  MyMatrix<T> M(n + 1, n + 1);
  MyVector<T> b(n + 1);
  auto set_row = [&](int row, size_t q) -> void {
    int i = q % m;
    for (int k = 0; k < n; k++) {
      M(row, k) = data.Lmat(i, k);
    }
    M(row, n) = -1;
    b(row) = ell_pair[q];
  };
  for (int i1 = 0; i1 < m; i1++) {
    size_t q0 = u_zero * m + i1;
    set_row(0, q0);
    // The other pairs, from which n are chosen.
    std::vector<size_t> others;
    for (size_t q = 0; q < n_pair; q++) {
      if (q != q0) {
        others.push_back(q);
      }
    }
    size_t n_oth = others.size();
    if (n_oth < size_t(n)) {
      continue;
    }
    std::vector<size_t> comb(n);
    for (int r = 0; r < n; r++) {
      comb[r] = r;
    }
    while (true) {
      for (int r = 0; r < n; r++) {
        set_row(r + 1, others[comb[r]]);
      }
      n_system++;
      T det = DeterminantMat(M);
      if (det != 0) {
        MyMatrix<T> Minv = Inverse(M);
        MyVector<T> x = Minv * b;
        T mu = x(n);
        if (mu > mu_max) {
          MyVector<T> p(n);
          for (int k = 0; k < n; k++) {
            p(k) = x(k);
          }
          if (is_valid(p, mu)) {
            mu_max = mu;
            p_max = p;
          }
        }
      }
      // Next combination in lexicographic order.
      int r = n - 1;
      while (r >= 0 && comb[r] == n_oth - n + r) {
        r--;
      }
      if (r < 0) {
        break;
      }
      comb[r]++;
      for (int r2 = r + 1; r2 < n; r2++) {
        comb[r2] = comb[r2 - 1] + 1;
      }
    }
  }
#ifdef DEBUG_POLYNORM_COVERING
  os << "POLYNORM_COVERING: brute force |L|=" << n_latt
     << " n_system=" << n_system << " mu=" << mu_max << "\n";
#endif
  return {mu_max, p_max};
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_COVERING_H_
// clang-format on
