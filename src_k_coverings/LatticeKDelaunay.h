// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_K_COVERINGS_LATTICEKDELAUNAY_H_
#define SRC_K_COVERINGS_LATTICEKDELAUNAY_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "LatticeStabEquiCan.h"
#include "POLY_DualDesc_reverse_search.h"
#include <algorithm>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_K_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_K_DELAUNAY
#endif

#ifdef TIMINGS
#define TIMINGS_K_DELAUNAY
#endif

/*
  Order-k Delaunay tiling of a lattice.

  Let L = Z^n with the positive definite Gram matrix Q and let k >= 1. The
  order-k Voronoi diagram k-V(L) is the subdivision of R^n into the cells

      cell(S) = { x : |x - s|_Q <= |x - p|_Q for s in S, p in L \ S },

  S running over the k-subsets of L with a nonempty cell. The k-covering
  radius of (L, Q) is the largest distance from a point of R^n to its k-th
  nearest lattice point; that distance function is convex on each cell, so
  the maximum is reached at a vertex of k-V(L).

  A vertex v of k-V(L) is characterized by its sphere: the lattice points at
  distance exactly r from v form the set P_0, the ones strictly inside form
  the set P_- of cardinality i, and v is a vertex if and only if P_0 spans
  R^n affinely and i < k < i + |P_0|. The k-subsets whose cell contains v
  are the sets P_- cup T with T a (k-i)-subset of P_0.

  Following Aurenhammer (A new duality result concerning Voronoi diagrams,
  1990), the lifted point set

      Q_k = { (sum(S), sum_{p in S} Q[p]) : S a k-subset of L }

  has its lower convex hull dual to k-V(L): the facet dual to the vertex v
  projects to the polytope

      D(v) = sum(P_-) + Delta_j(P_0),   j = k - i,

  where Delta_j(P_0) is the hypersimplex of P_0, the convex hull of the
  sums of the j-subsets of P_0. Those tiles form a tiling of R^n, which is
  the order-k Delaunay tiling of (L, Q) (order-k Delaunay mosaic in the
  terminology of Edelsbrunner and Nikitenko). For k = 1 it is the Delaunay
  tessellation. The coordinates used here are the sum coordinates, so the
  tiles are lattice polytopes and a translation by w in L acts on them as a
  translation by k w.

  A tile is stored by its sphere data: the matrix EXT of the points of P_0
  and the matrix INT of the points of P_-, both with the homogenizing 1 in
  the first column. P_- is determined by P_0 (it is the set of lattice
  points strictly inside the circumsphere of P_0) but is kept since every
  operation on a tile needs it.

  The facets of a tile D(v) correspond to the hyperplanes H spanned by
  subsets B = P_0 cap H, with A the points of P_0 on one side of H and C
  the points on the other side, subject to |A| < j < |A| + |B|: the facet
  is sum(P_-) + sum(A) + Delta_{j - |A|}(B). The tile adjacent across that
  facet is dual to the other endpoint of the edge of k-V(L) starting at v
  along which the sphere keeps B on its boundary and swallows A: its center
  moves on the line v - t Q^{-1} l for t > 0, l being the normal covector
  of H oriented so that A is on the negative side. Along that sweep the
  points of C leave the sphere at once, and the next vertex is met at the
  first t > 0 where either an outside lattice point enters the sphere or a
  point of P_- leaves it. The exit events are read off P_-, the entering
  event is found by a local descent followed by an exact enumeration of the
  lattice points of the candidate ball, in the spirit of
  FindAdjacentDelaunayPolytope.

  The initial tile comes from the vertex of the cell of a k-nearest set S
  found by linear programming, generalizing FindDelaunayPolytope_direction:
  the cell of S is described by the inequalities |x - s|^2 <= |x - p|^2 and
  the constraint set is grown until a vertex of the current approximation
  is certified by the enumeration of the lattice points of its ball.

  The enumeration of the tiles up to the group Aut(L, Q) of affine lattice
  isometries runs through the generic adjacency scheme, with the stabilizer
  and equivalence tests of the Delaunay code applied to P_0: an affine
  isometry maps a sphere to a sphere, so it maps P_- along with P_0.
 */

template <typename Tint> struct KDelaunayTile {
  // The lattice points on the sphere (P_0), homogeneous coordinates.
  MyMatrix<Tint> EXT;
  // The lattice points strictly inside the sphere (P_-), homogeneous
  // coordinates.
  MyMatrix<Tint> INT;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, KDelaunayTile<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", eRec.EXT);
  ar &make_nvp("INT", eRec.INT);
}
} // namespace boost::serialization

// Whether the sphere data (EXT, INT) is that of a vertex of k-V(L): P_0
// spans R^n affinely and |P_-| < k < |P_-| + |P_0|.
template <typename Tint>
bool IsValidKDelaunayTile(KDelaunayTile<Tint> const &tile, int const &k) {
  int n = tile.EXT.cols() - 1;
  int i = tile.INT.rows();
  int m = tile.EXT.rows();
  if (i >= k || k >= i + m) {
    return false;
  }
  if (RankMat(tile.EXT) != n + 1) {
    return false;
  }
  return true;
}

// The sphere of a tile: circumcenter and squared radius of P_0.
template <typename T, typename Tint>
CP<T> GetKDelaunayTileSphere(MyMatrix<T> const &GramMat,
                             KDelaunayTile<Tint> const &tile) {
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(tile.EXT);
  return CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT_T);
}

// The tile of a sphere given by its center (non-homogeneous) and squared
// radius: the lattice points of the closed ball are enumerated and split
// into the boundary and the interior.
template <typename T, typename Tint>
KDelaunayTile<Tint> GetKDelaunayTileFromSphere(CVPSolver<T, Tint> const &solver,
                                               MyVector<T> const &center,
                                               T const &SquareRadius) {
  int n = center.size();
  std::vector<MyVector<Tint>> ball =
      solver.at_most_dist_vectors(center, SquareRadius);
  std::vector<MyVector<Tint>> l_ext, l_int;
  for (auto &q : ball) {
    T dist = solver.comp_norm_diff(q, center);
    if (dist == SquareRadius) {
      l_ext.push_back(q);
    } else {
      l_int.push_back(q);
    }
  }
  auto to_matrix = [&](std::vector<MyVector<Tint>> const &l) -> MyMatrix<Tint> {
    int n_row = l.size();
    MyMatrix<Tint> M(n_row, n + 1);
    for (int i_row = 0; i_row < n_row; i_row++) {
      M(i_row, 0) = 1;
      for (int i = 0; i < n; i++) {
        M(i_row, i + 1) = l[i_row](i);
      }
    }
    return M;
  };
  return {to_matrix(l_ext), to_matrix(l_int)};
}

// The homogeneous rows of a matrix as a hashable set of non-homogeneous
// integral vectors.
template <typename Tint>
std::unordered_set<MyVector<Tint>>
GetRowSetNonHomogeneous(MyMatrix<Tint> const &M) {
  std::unordered_set<MyVector<Tint>> set;
  int n = M.cols() - 1;
  int n_row = M.rows();
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<Tint> V(n);
    for (int i = 0; i < n; i++) {
      V(i) = M(i_row, i + 1);
    }
    set.insert(V);
  }
  return set;
}

/*
  The vertices of a tile in sum coordinates: sum(P_-) + sum(T) for T
  running over the j-subsets of P_0, j = k - |P_-|. Distinct subsets can
  have the same sum when P_0 is not affinely independent, so the sums are
  made unique; for each vertex one subset realizing it is kept, which the
  induced action of the stabilizer of P_0 on the vertices needs.
 */
template <typename Tint> struct KDelaunayTileVertices {
  // The vertices, homogeneous coordinates.
  MyMatrix<Tint> EXTsum;
  // For each vertex, the rows of EXT of one subset realizing it.
  std::vector<std::vector<int>> ListSubset;
};

template <typename Tint>
KDelaunayTileVertices<Tint>
GetKDelaunayTileVertices(KDelaunayTile<Tint> const &tile, int const &k) {
  int n = tile.EXT.cols() - 1;
  int m = tile.EXT.rows();
  int i = tile.INT.rows();
  int j = k - i;
#ifdef SANITY_CHECK_K_DELAUNAY
  if (j <= 0 || j >= m) {
    std::cerr << "K_DELAUNAY: GetKDelaunayTileVertices: j=" << j
              << " must satisfy 0 < j < m=" << m << "\n";
    throw TerminalException{1};
  }
#endif
  MyVector<Tint> base = ZeroVector<Tint>(n);
  for (int i_row = 0; i_row < i; i_row++) {
    for (int u = 0; u < n; u++) {
      base(u) += tile.INT(i_row, u + 1);
    }
  }
  std::unordered_map<MyVector<Tint>, size_t> map;
  std::vector<MyVector<Tint>> l_vert;
  std::vector<std::vector<int>> ListSubset;
  // Enumeration of the j-subsets of {0, ..., m-1} in lexicographic order.
  std::vector<int> subset(j);
  for (int u = 0; u < j; u++) {
    subset[u] = u;
  }
  while (true) {
    MyVector<Tint> V = base;
    for (int u = 0; u < j; u++) {
      for (int w = 0; w < n; w++) {
        V(w) += tile.EXT(subset[u], w + 1);
      }
    }
    if (map.find(V) == map.end()) {
      map[V] = l_vert.size();
      l_vert.push_back(V);
      ListSubset.push_back(subset);
    }
    // Next subset.
    int pos = j - 1;
    while (pos >= 0 && subset[pos] == m - j + pos) {
      pos--;
    }
    if (pos < 0) {
      break;
    }
    subset[pos]++;
    for (int u = pos + 1; u < j; u++) {
      subset[u] = subset[u - 1] + 1;
    }
  }
  int n_vert = l_vert.size();
  MyMatrix<Tint> EXTsum(n_vert, n + 1);
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    EXTsum(i_vert, 0) = 1;
    for (int u = 0; u < n; u++) {
      EXTsum(i_vert, u + 1) = l_vert[i_vert](u);
    }
  }
  return {std::move(EXTsum), std::move(ListSubset)};
}

/*
  The permutation group induced on the vertices of the tile by a group of
  permutations of the rows of P_0 (the stabilizer of the tile). A
  permutation g of P_0 extends to an affine isometry, which maps the sum of
  a subset T to the sum of g(T); equal sums therefore have equal images and
  the induced permutation is well defined.
 */
template <typename Tint, typename Tgroup>
Tgroup GetInducedGroupOnVertices(KDelaunayTile<Tint> const &tile,
                                 KDelaunayTileVertices<Tint> const &vert,
                                 Tgroup const &GRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = tile.EXT.cols() - 1;
  int n_vert = vert.EXTsum.rows();
  std::unordered_map<MyVector<Tint>, size_t> map;
  for (int i_vert = 0; i_vert < n_vert; i_vert++) {
    MyVector<Tint> V(n);
    for (int u = 0; u < n; u++) {
      V(u) = vert.EXTsum(i_vert, u + 1);
    }
    map[V] = i_vert;
  }
  MyVector<Tint> base = ZeroVector<Tint>(n);
  int i = tile.INT.rows();
  for (int i_row = 0; i_row < i; i_row++) {
    for (int u = 0; u < n; u++) {
      base(u) += tile.INT(i_row, u + 1);
    }
  }
  std::vector<Telt> ListGen;
  for (auto &eGen : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> l_pos(n_vert);
    for (int i_vert = 0; i_vert < n_vert; i_vert++) {
      MyVector<Tint> V = base;
      for (auto &i_row : vert.ListSubset[i_vert]) {
        int i_img = OnPoints(i_row, eGen);
        for (int u = 0; u < n; u++) {
          V(u) += tile.EXT(i_img, u + 1);
        }
      }
      auto iter = map.find(V);
#ifdef SANITY_CHECK_K_DELAUNAY
      if (iter == map.end()) {
        std::cerr << "K_DELAUNAY: GetInducedGroupOnVertices: the image of a "
                     "vertex is not a vertex\n";
        throw TerminalException{1};
      }
#endif
      l_pos[i_vert] = iter->second;
    }
    ListGen.push_back(Telt(l_pos));
  }
  return Tgroup(ListGen, n_vert);
}

/*
  The data of a tile over the field: the sphere and the affine coordinates
  of P_- and P_0, computed once per tile and shared by its facets.
 */
template <typename T> struct KDelaunayTileGeometry {
  int n;
  int i;
  int m;
  int j;
  MyMatrix<T> EXT_T;
  MyMatrix<T> INT_T;
  // The center (non-homogeneous) and the squared radius.
  MyVector<T> cent;
  T SquareRadius;
};

template <typename T, typename Tint>
KDelaunayTileGeometry<T>
GetKDelaunayTileGeometry(MyMatrix<T> const &GramMat,
                         KDelaunayTile<Tint> const &tile, int const &k) {
  int n = GramMat.rows();
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(tile.EXT);
  MyMatrix<T> INT_T = UniversalMatrixConversion<T, Tint>(tile.INT);
  CP<T> cp = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT_T);
  MyVector<T> cent(n);
  for (int u = 0; u < n; u++) {
    cent(u) = cp.eCent(u + 1);
  }
  int i = tile.INT.rows();
  int m = tile.EXT.rows();
  return {n, i, m, k - i, std::move(EXT_T), std::move(INT_T),
          std::move(cent), cp.SquareRadius};
}

// |p - v|_Q^2 - r^2 for the sphere of the tile: negative inside, zero on
// the sphere, positive outside.
template <typename T>
T KDelaunaySphereExcess(MyMatrix<T> const &GramMat,
                        KDelaunayTileGeometry<T> const &geom,
                        MyVector<T> const &p) {
  MyVector<T> diff = p - geom.cent;
  return EvaluationQuadForm<T, T>(GramMat, diff) - geom.SquareRadius;
}

/*
  The adjacent tile across the facet of the tile with facet functional
  ell_hom (homogeneous, nonnegative on the tile and zero on the facet), by
  the sweep described at the top of the file.
 */
template <typename T, typename Tint>
KDelaunayTile<Tint>
FindAdjacentKDelaunayTile(CVPSolver<T, Tint> const &solver,
                          MyMatrix<Tint> const &ShvGraverBasis,
                          MyMatrix<T> const &Qinv,
                          KDelaunayTile<Tint> const &tile,
                          KDelaunayTileGeometry<T> const &geom,
                          MyVector<T> const &ell_hom, int const &k,
                          [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_K_DELAUNAY
  MicrosecondTime time;
#endif
  MyMatrix<T> const &GramMat = solver.GramMat;
  int n = geom.n;
  int j = geom.j;
  MyVector<T> ell(n);
  for (int u = 0; u < n; u++) {
    ell(u) = ell_hom(u + 1);
  }
  // The values of the linear part on P_0 and the threshold theta: the facet
  // subsets T are the ones minimizing ell . sum(T), i.e. taking the j
  // smallest values; A is the set of points strictly below the j-th
  // smallest value, B the points at that value, C the rest.
  std::vector<T> vals(geom.m);
  for (int i_row = 0; i_row < geom.m; i_row++) {
    T val(0);
    for (int u = 0; u < n; u++) {
      AddMul(val, ell(u), geom.EXT_T(i_row, u + 1));
    }
    vals[i_row] = val;
  }
  std::vector<T> sorted = vals;
  std::sort(sorted.begin(), sorted.end());
  T theta = sorted[j - 1];
  std::vector<int> A, B, C;
  for (int i_row = 0; i_row < geom.m; i_row++) {
    if (vals[i_row] < theta) {
      A.push_back(i_row);
    } else if (vals[i_row] == theta) {
      B.push_back(i_row);
    } else {
      C.push_back(i_row);
    }
  }
#ifdef SANITY_CHECK_K_DELAUNAY
  {
    int sizA = A.size();
    int sizB = B.size();
    if (!(sizA < j && j < sizA + sizB)) {
      std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: |A|=" << sizA
                << " |B|=" << sizB << " j=" << j
                << " do not satisfy |A| < j < |A| + |B|\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Bmat = SelectRow(geom.EXT_T, B);
    if (RankMat(Bmat) != n) {
      std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: B does not span a "
                   "hyperplane, rank=" << RankMat(Bmat) << " n=" << n << "\n";
      throw TerminalException{1};
    }
  }
#endif
  // The sweep: center c(t) = v - t u with u = Q^{-1} ell. For a point p the
  // signed distance to the sphere is |p - c(t)|^2 - r(t)^2 =
  // phi(p) - 2 t (theta - ell . p) with phi the excess at t = 0, so p is on
  // the moving sphere at t_p = phi(p) / (2 (theta - ell . p)).
  MyVector<T> u = Qinv * ell;
  MyVector<T> b(n);
  for (int w = 0; w < n; w++) {
    b(w) = geom.EXT_T(B[0], w + 1);
  }
  auto get_phi = [&](MyVector<T> const &p) -> T {
    return KDelaunaySphereExcess(GramMat, geom, p);
  };
  auto get_gam = [&](MyVector<T> const &p) -> T {
    T val = theta;
    for (int w = 0; w < n; w++) {
      SubMul(val, ell(w), p(w));
    }
    return val;
  };
  // Exit events: points of P_- on the positive side of H.
  std::optional<T> t_exit;
  for (int i_row = 0; i_row < geom.i; i_row++) {
    MyVector<T> p(n);
    for (int w = 0; w < n; w++) {
      p(w) = geom.INT_T(i_row, w + 1);
    }
    T gam = get_gam(p);
    if (gam < 0) {
      T phi = get_phi(p);
      T t = phi / (2 * gam);
      if (!t_exit || t < *t_exit) {
        t_exit = t;
      }
    }
  }
  // Entering candidate: a lattice point outside the sphere on the negative
  // side of H, starting from a point of B moved by a step of the lattice.
  auto get_move = [&]() -> MyVector<T> {
    for (auto &eMove : solver.get_seed_differences()) {
      T eScal(0);
      for (int w = 0; w < n; w++) {
        AddMul(eScal, ell(w), eMove(w));
      }
      if (eScal < 0) {
        return eMove;
      }
    }
    std::cerr << "K_DELAUNAY: no move of the lattice goes to the negative "
                 "side of the facet\n";
    throw TerminalException{1};
  };
  MyVector<T> eMove = get_move();
  MyVector<T> cand = b;
  while (true) {
    cand += eMove;
    if (get_phi(cand) > 0) {
      break;
    }
  }
  auto get_t = [&](MyVector<T> const &p) -> T {
    return get_phi(p) / (2 * get_gam(p));
  };
  T t_cand = get_t(cand);
  auto graver_descent = [&]() -> void {
    int n_graver = ShvGraverBasis.rows();
    MyVector<T> q(n);
    while (true) {
      bool IsImprovement = false;
      for (int i_graver = 0; i_graver < n_graver; i_graver++) {
        for (int w = 0; w < n; w++) {
          q(w) = cand(w) +
                 UniversalScalarConversion<T, Tint>(ShvGraverBasis(i_graver, w));
        }
        if (get_gam(q) > 0 && get_phi(q) > 0) {
          T t_q = get_t(q);
          if (t_q < t_cand) {
            cand = q;
            t_cand = t_q;
            IsImprovement = true;
          }
        }
      }
      if (!IsImprovement) {
        return;
      }
    }
  };
  graver_descent();
  std::unordered_set<MyVector<Tint>> setIn = GetRowSetNonHomogeneous(tile.INT);
  for (auto &i_row : A) {
    MyVector<Tint> V(n);
    for (int w = 0; w < n; w++) {
      V(w) = tile.EXT(i_row, w + 1);
    }
    setIn.insert(V);
  }
  // Verification: the lattice points of the ball at the candidate time are
  // enumerated; a point strictly inside that is neither in P_- nor in A is
  // an outside point that entered earlier, and becomes the candidate.
  while (true) {
    T t_star = t_cand;
    if (t_exit && *t_exit < t_star) {
      t_star = *t_exit;
    }
    MyVector<T> cent_new = geom.cent - t_star * u;
    MyVector<T> diff = b - cent_new;
    T SquareRadius_new = EvaluationQuadForm<T, T>(GramMat, diff);
    std::vector<MyVector<Tint>> ball =
        solver.at_most_dist_vectors(cent_new, SquareRadius_new);
    std::optional<MyVector<Tint>> bad;
    for (auto &q : ball) {
      T dist = solver.comp_norm_diff(q, cent_new);
      if (dist < SquareRadius_new) {
        if (setIn.find(q) == setIn.end()) {
          bad = q;
          break;
        }
      }
    }
    if (!bad) {
      KDelaunayTile<Tint> tile_new =
          GetKDelaunayTileFromSphere<T, Tint>(solver, cent_new, SquareRadius_new);
#ifdef SANITY_CHECK_K_DELAUNAY
      if (!IsValidKDelaunayTile(tile_new, k)) {
        std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: the adjacent "
                     "tile is not valid |INT|=" << tile_new.INT.rows()
                  << " |EXT|=" << tile_new.EXT.rows() << " k=" << k << "\n";
        throw TerminalException{1};
      }
      std::unordered_set<MyVector<Tint>> setExt =
          GetRowSetNonHomogeneous(tile_new.EXT);
      for (auto &i_row : B) {
        MyVector<Tint> V(n);
        for (int w = 0; w < n; w++) {
          V(w) = tile.EXT(i_row, w + 1);
        }
        if (setExt.find(V) == setExt.end()) {
          std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: a point of B "
                       "is not on the adjacent sphere\n";
          throw TerminalException{1};
        }
      }
#endif
#ifdef TIMINGS_K_DELAUNAY
      os << "|K_DELAUNAY: FindAdjacentKDelaunayTile|=" << time << "\n";
#endif
      return tile_new;
    }
    MyVector<T> q_T = UniversalVectorConversion<T, Tint>(*bad);
#ifdef SANITY_CHECK_K_DELAUNAY
    if (!(get_gam(q_T) > 0 && get_phi(q_T) > 0)) {
      std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: the unexpected "
                   "interior point is not an entering point\n";
      throw TerminalException{1};
    }
    if (!(get_t(q_T) < t_star)) {
      std::cerr << "K_DELAUNAY: FindAdjacentKDelaunayTile: the entering time "
                   "of the unexpected point is not smaller\n";
      throw TerminalException{1};
    }
#endif
    cand = q_T;
    t_cand = get_t(cand);
    graver_descent();
  }
}

/*
  A vertex of the cell of a k-nearest set, by linear programming. The
  random point x0 provides the k-nearest set S (the k-th and (k+1)-th
  distances have to differ, otherwise nothing is returned and the caller
  changes x0). The polyhedron of x with |x - s|^2 <= |x - p|^2 (s in S, p
  in a growing finite set) contains cell(S); the random objective is
  minimized over it and the optimal vertex x is certified when the lattice
  points strictly inside the ball of radius max_s |x - s| all belong to S;
  the offending points are otherwise added to the constraints. A vertex of
  the enclosing polyhedron in cell(S) is a vertex of cell(S), hence of
  k-V(L), and the ball of x gives the tile.
 */
template <typename T, typename Tint>
std::optional<KDelaunayTile<Tint>>
FindInitialKDelaunayTile_direction(CVPSolver<T, Tint> const &solver,
                                   int const &k, MyVector<T> const &x0,
                                   MyVector<T> const &TheRandomDirection,
                                   [[maybe_unused]] std::ostream &os) {
  MyMatrix<T> const &GramMat = solver.GramMat;
  int n = GramMat.rows();
  // The k nearest lattice points of x0, from a ball grown until it holds at
  // least k + 1 points.
  resultCVP<T, Tint> near = solver.nearest_vectors(x0);
  T R = near.TheNorm;
  std::vector<MyVector<Tint>> ball;
  while (true) {
    ball = solver.at_most_dist_vectors(x0, R);
    if (static_cast<int>(ball.size()) >= k + 1) {
      break;
    }
    R = 2 * R + 1;
  }
  std::vector<std::pair<T, MyVector<Tint>>> l_pair;
  for (auto &q : ball) {
    l_pair.push_back({solver.comp_norm_diff(q, x0), q});
  }
  std::sort(l_pair.begin(), l_pair.end(),
            [](std::pair<T, MyVector<Tint>> const &a,
               std::pair<T, MyVector<Tint>> const &b) -> bool {
              return a.first < b.first;
            });
  if (l_pair[k - 1].first == l_pair[k].first) {
#ifdef DEBUG_K_DELAUNAY
    os << "K_DELAUNAY: FindInitialKDelaunayTile_direction: tie at the k-th "
          "distance, retrying\n";
#endif
    return {};
  }
  std::vector<MyVector<Tint>> S;
  std::unordered_set<MyVector<Tint>> setS;
  for (int u = 0; u < k; u++) {
    S.push_back(l_pair[u].second);
    setS.insert(l_pair[u].second);
  }
  std::vector<MyVector<T>> S_T;
  std::vector<T> S_norm;
  for (auto &s : S) {
    MyVector<T> s_T = UniversalVectorConversion<T, Tint>(s);
    S_norm.push_back(EvaluationQuadForm<T, T>(GramMat, s_T));
    S_T.push_back(s_T);
  }
  // The inequality |x - s|^2 <= |x - p|^2, i.e.
  // Q[p] - Q[s] - 2 (p - s)^T Q x >= 0, as (constant, linear part).
  std::vector<MyVector<T>> ListIneq_vect;
  auto insert_ineq = [&](int i_s, MyVector<T> const &p_T) -> void {
    MyVector<T> eIneq(n + 1);
    eIneq(0) = EvaluationQuadForm<T, T>(GramMat, p_T) - S_norm[i_s];
    MyVector<T> diff = p_T - S_T[i_s];
    MyVector<T> Qdiff = GramMat * diff;
    for (int u = 0; u < n; u++) {
      eIneq(u + 1) = -2 * Qdiff(u);
    }
    ListIneq_vect.push_back(eIneq);
  };
  int n_S = S.size();
  // For each s and each lattice step the first point along that step
  // outside S: those bound the polyhedron.
  for (int i_s = 0; i_s < n_S; i_s++) {
    for (auto &eMove : solver.get_seed_differences()) {
      MyVector<Tint> eMove_i = UniversalVectorConversion<Tint, T>(eMove);
      MyVector<Tint> p = S[i_s];
      while (true) {
        p += eMove_i;
        if (setS.find(p) == setS.end()) {
          break;
        }
      }
      insert_ineq(i_s, UniversalVectorConversion<T, Tint>(p));
    }
  }
  for (auto &ePair : l_pair) {
    if (setS.find(ePair.second) == setS.end()) {
      MyVector<T> p_T = UniversalVectorConversion<T, Tint>(ePair.second);
      for (int i_s = 0; i_s < n_S; i_s++) {
        insert_ineq(i_s, p_T);
      }
    }
  }
  while (true) {
    MyMatrix<T> ListIneq = MatrixFromVectorFamily(ListIneq_vect);
    LpSolution<T> eSol =
        SIMPLEX_LinearProgramming(ListIneq, TheRandomDirection, os);
    if (!eSol.DirectSolution) {
#ifdef DEBUG_K_DELAUNAY
      os << "K_DELAUNAY: FindInitialKDelaunayTile_direction: no primal "
            "solution, retrying\n";
#endif
      return {};
    }
    MyVector<T> const &x = *eSol.DirectSolution;
    if (!IsVertexOfPolyhedron(ListIneq, x)) {
#ifdef DEBUG_K_DELAUNAY
      os << "K_DELAUNAY: FindInitialKDelaunayTile_direction: the optimum is "
            "not a vertex, retrying\n";
#endif
      return {};
    }
    T R2(0);
    for (int i_s = 0; i_s < n_S; i_s++) {
      MyVector<T> diff = x - S_T[i_s];
      T dist = EvaluationQuadForm<T, T>(GramMat, diff);
      if (dist > R2) {
        R2 = dist;
      }
    }
    std::vector<MyVector<Tint>> ball_x = solver.at_most_dist_vectors(x, R2);
    std::vector<MyVector<Tint>> bad;
    for (auto &q : ball_x) {
      T dist = solver.comp_norm_diff(q, x);
      if (dist < R2 && setS.find(q) == setS.end()) {
        bad.push_back(q);
      }
    }
    if (bad.empty()) {
      KDelaunayTile<Tint> tile =
          GetKDelaunayTileFromSphere<T, Tint>(solver, x, R2);
#ifdef SANITY_CHECK_K_DELAUNAY
      if (!IsValidKDelaunayTile(tile, k)) {
        std::cerr << "K_DELAUNAY: FindInitialKDelaunayTile_direction: the "
                     "tile is not valid |INT|=" << tile.INT.rows()
                  << " |EXT|=" << tile.EXT.rows() << " k=" << k << "\n";
        throw TerminalException{1};
      }
#endif
#ifdef DEBUG_K_DELAUNAY
      os << "K_DELAUNAY: FindInitialKDelaunayTile_direction: found a tile "
            "with |INT|=" << tile.INT.rows() << " |EXT|=" << tile.EXT.rows()
         << "\n";
#endif
      return tile;
    }
    for (auto &q : bad) {
      MyVector<T> q_T = UniversalVectorConversion<T, Tint>(q);
      for (int i_s = 0; i_s < n_S; i_s++) {
        insert_ineq(i_s, q_T);
      }
    }
  }
}

template <typename T, typename Tint>
KDelaunayTile<Tint> FindInitialKDelaunayTile(CVPSolver<T, Tint> const &solver,
                                             int const &k, std::ostream &os) {
  int n = solver.GramMat.rows();
  int N = 3;
  while (true) {
    // A random rational point of the cube [-1, 1]^n with a denominator
    // growing along the attempts, which keeps the ties at the k-th distance
    // exceptional.
    int denom = 101 + 2 * N;
    MyVector<T> x0(n);
    for (int u = 0; u < n; u++) {
      int num = random() % (2 * denom + 1);
      x0(u) = T(num - denom) / T(denom);
    }
    MyVector<T> TheRandomDirection = FuncRandomDirection<T>(n + 1, N);
    std::optional<KDelaunayTile<Tint>> opt =
        FindInitialKDelaunayTile_direction<T, Tint>(solver, k, x0,
                                                    TheRandomDirection, os);
    if (opt) {
      return *opt;
    }
    N += 1;
  }
}

/*
  Enumeration data. The objects are tiles, the adjacencies are facets (as
  faces of the vertex set of the tile) with the adjacent tile, and the
  equivalence is that of the affine lattice isometries applied to P_0. The
  transformation eBigMat of an adjacency maps the orbit representative onto
  the actual adjacent tile: rep.EXT * eBigMat = adjacent.EXT as sets.
 */
template <typename Tint> struct KDelaunay_AdjI {
  Face eInc;
  KDelaunayTile<Tint> tile;
};

template <typename Tint> struct KDelaunay_AdjO_spec {
  Face eInc;
  MyMatrix<Tint> eBigMat;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, KDelaunay_AdjI<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("tile", eRec.tile);
}
template <class Archive, typename Tint>
inline void serialize(Archive &ar, KDelaunay_AdjO_spec<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
}
} // namespace boost::serialization

template <typename Tint> struct KDelaunay_AdjO {
  Face eInc;
  MyMatrix<Tint> eBigMat;
  int iOrb;
};

template <typename Tint, typename Tgroup> struct KDelaunay_Entry {
  KDelaunayTile<Tint> tile;
  Tgroup GRP;
  std::vector<KDelaunay_AdjO<Tint>> ListAdj;
};

template <typename Tint, typename Tgroup> struct KDelaunayTesselation {
  int k;
  std::vector<KDelaunay_Entry<Tint, Tgroup>> l_tiles;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, KDelaunay_AdjO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
  ar &make_nvp("iOrb", eRec.iOrb);
}
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, KDelaunay_Entry<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("tile", eRec.tile);
  ar &make_nvp("GRP", eRec.GRP);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, KDelaunayTesselation<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("k", eRec.k);
  ar &make_nvp("l_tiles", eRec.l_tiles);
}
} // namespace boost::serialization

template <typename Tint, typename Tgroup> struct KDelaunay_Obj {
  KDelaunayTile<Tint> tile;
  Tgroup GRP;
};

namespace boost::serialization {
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, KDelaunay_Obj<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("tile", eRec.tile);
  ar &make_nvp("GRP", eRec.GRP);
}
} // namespace boost::serialization

// The stabilizer of the tile and the adjacent tile across each orbit
// representative of facets.
template <typename T, typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<KDelaunay_AdjI<Tint>>>
ComputeKDelaunayGroupAndAdjacencies(DataLattice<T, Tint, Tgroup> &eData,
                                    MyMatrix<T> const &Qinv, int const &k,
                                    KDelaunayTile<Tint> const &tile) {
  std::ostream &os = eData.rddo.os;
  MyMatrix<T> const &GramMat = eData.solver.GramMat;
#ifdef TIMINGS_K_DELAUNAY
  MicrosecondTime time;
#endif
  KDelaunayTileGeometry<T> geom = GetKDelaunayTileGeometry(GramMat, tile, k);
  Tgroup GRP = Delaunay_Stabilizer(eData, geom.EXT_T);
#ifdef TIMINGS_K_DELAUNAY
  os << "|K_DELAUNAY: Delaunay_Stabilizer|=" << time << "\n";
#endif
  KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(tile, k);
  Tgroup GRPsum = GetInducedGroupOnVertices(tile, vert, GRP);
  MyMatrix<T> EXTsum_T = UniversalMatrixConversion<T, Tint>(vert.EXTsum);
#ifdef DEBUG_K_DELAUNAY
  os << "K_DELAUNAY: tile |INT|=" << tile.INT.rows() << " |EXT|="
     << tile.EXT.rows() << " |GRP|=" << GRP.size() << " |vertices|="
     << vert.EXTsum.rows() << "\n";
#endif
  vectface TheOutput = DualDescriptionRecordFullDim(EXTsum_T, GRPsum, eData.rddo);
#ifdef TIMINGS_K_DELAUNAY
  os << "|K_DELAUNAY: DualDescriptionRecordFullDim|=" << time << "\n";
#endif
  SubsetRankOneSolver<T> ext_solver(EXTsum_T);
  std::vector<KDelaunay_AdjI<Tint>> ListAdj;
  for (auto &eOrbB : TheOutput) {
    MyVector<T> ell_hom = ext_solver.GetPositiveKernelVector(eOrbB);
    KDelaunayTile<Tint> tile_adj = FindAdjacentKDelaunayTile<T, Tint>(
        eData.solver, eData.ShvGraverBasis, Qinv, tile, geom, ell_hom, k, os);
    ListAdj.push_back({eOrbB, std::move(tile_adj)});
  }
#ifdef DEBUG_K_DELAUNAY
  os << "K_DELAUNAY: |ListAdj|=" << ListAdj.size() << "\n";
#endif
  return {std::move(GRP), std::move(ListAdj)};
}

template <typename T, typename Tint, typename Tgroup> struct DataLatticeKFunc {
  DataLattice<T, Tint, Tgroup> &data;
  int k;
  MyMatrix<T> Qinv;
  using Tobj = KDelaunay_Obj<Tint, Tgroup>;
  using TadjI = KDelaunay_AdjI<Tint>;
  using TadjO = KDelaunay_AdjO_spec<Tint>;
  DataLatticeKFunc(DataLattice<T, Tint, Tgroup> &_data, int const &_k)
      : data(_data), k(_k), Qinv(Inverse(_data.solver.GramMat)) {}
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
    KDelaunayTile<Tint> tile =
        FindInitialKDelaunayTile<T, Tint>(data.solver, k, data.rddo.os);
    return {std::move(tile), {}};
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    size_t hash = ComputeInvariantDelaunay(data, seed, x.tile.EXT, data.rddo.os);
    size_t hash_int = std::hash<int>()(x.tile.INT.rows());
    hash ^= hash_int + 0x9e3779b8 + (hash << 6) + (hash >> 2);
    return hash;
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    MyMatrix<T> EXT1_T = UniversalMatrixConversion<T, Tint>(x.tile.EXT);
    MyMatrix<T> EXT2_T = UniversalMatrixConversion<T, Tint>(y.tile.EXT);
    std::optional<MyMatrix<T>> opt =
        Delaunay_TestEquivalence<T, Tint, Tgroup>(data, EXT1_T, EXT2_T);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> eBigMat = UniversalMatrixConversion<Tint, T>(*opt);
#ifdef SANITY_CHECK_K_DELAUNAY
    // The isometry maps the sphere of x to the sphere of y, hence P_- to P_-.
    MyMatrix<Tint> INTimg = x.tile.INT * eBigMat;
    std::unordered_set<MyVector<Tint>> setInt =
        GetRowSetNonHomogeneous(y.tile.INT);
    std::unordered_set<MyVector<Tint>> setImg = GetRowSetNonHomogeneous(INTimg);
    if (setInt != setImg) {
      std::cerr << "K_DELAUNAY: f_repr: the equivalence does not map P_- to "
                   "P_-\n";
      throw TerminalException{1};
    }
#endif
    return TadjO{y.eInc, std::move(eBigMat)};
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    Tobj x_ret{x.tile, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n + 1);
    TadjO ret{x.eInc, eBigMat};
    return {std::move(x_ret), std::move(ret)};
  }
  std::optional<std::vector<TadjI>> f_adj(Tobj &x) {
    std::pair<Tgroup, std::vector<TadjI>> pair =
        ComputeKDelaunayGroupAndAdjacencies<T, Tint, Tgroup>(data, Qinv, k,
                                                             x.tile);
    x.GRP = pair.first;
    return pair.second;
  }
  Tobj f_adji_obj(TadjI const &x) { return {x.tile, {}}; }
  size_t f_complexity(Tobj const &x) { return x.tile.EXT.rows(); }
};

template <typename T, typename Tint, typename Tgroup>
KDelaunayTesselation<Tint, Tgroup>
KDelaunayTesselation_From_DatabaseEntries_Serial(
    int const &k,
    std::vector<DatabaseEntry_Serial<
        typename DataLatticeKFunc<T, Tint, Tgroup>::Tobj,
        typename DataLatticeKFunc<T, Tint, Tgroup>::TadjO>> const &l_ent) {
  std::vector<KDelaunay_Entry<Tint, Tgroup>> l_tiles;
  for (auto &eEnt : l_ent) {
    std::vector<KDelaunay_AdjO<Tint>> ListAdj;
    for (auto &eAdj : eEnt.ListAdj) {
      ListAdj.push_back({eAdj.x.eInc, eAdj.x.eBigMat, eAdj.iOrb});
    }
    l_tiles.push_back({eEnt.x.tile, eEnt.x.GRP, std::move(ListAdj)});
  }
  return {k, std::move(l_tiles)};
}

template <typename T, typename Tint, typename Tgroup, typename Fincorrect>
std::optional<KDelaunayTesselation<Tint, Tgroup>>
EnumerationKDelaunayTiles(DataLattice<T, Tint, Tgroup> &data, int const &k,
                          Fincorrect f_incorrect,
                          int const &max_runtime_second) {
  using Tdata = DataLatticeKFunc<T, Tint, Tgroup>;
  Tdata data_func(data, k);
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = std::vector<DatabaseEntry_Serial<Tobj, TadjO>>;
  std::optional<Tout> opt_result = EnumerateAndStore_Serial<Tdata>(
      data_func, f_incorrect, max_runtime_second);
  if (!opt_result) {
    return {};
  }
  return KDelaunayTesselation_From_DatabaseEntries_Serial<T, Tint, Tgroup>(
      k, *opt_result);
}

// The order-k Delaunay tiling of a Gram matrix with the standard
// heuristics, the plain entry point for the other programs.
template <typename T, typename Tint, typename Tgroup>
KDelaunayTesselation<Tint, Tgroup>
ComputeKDelaunayTesselation(MyMatrix<T> const &GramMat, int const &k,
                            PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
                            std::ostream &os) {
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  auto f_incorrect =
      [&]([[maybe_unused]] KDelaunay_Obj<Tint, Tgroup> const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::optional<KDelaunayTesselation<Tint, Tgroup>> opt =
      EnumerationKDelaunayTiles<T, Tint, Tgroup, decltype(f_incorrect)>(
          data, k, f_incorrect, max_runtime_second);
  return unfold_opt(opt, "The order-k Delaunay tiling");
}

// Consistency of the recorded adjacencies: the adjacent tile, mapped by
// eBigMat, has to share the facet eInc with the tile (the vertices of the
// tile on the facet are vertices of the adjacent tile, the others are not).
template <typename Tint, typename Tgroup>
void check_k_delaunay_tessellation(KDelaunayTesselation<Tint, Tgroup> const &DT,
                                   [[maybe_unused]] std::ostream &os) {
  int k = DT.k;
  for (auto &eEnt : DT.l_tiles) {
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(eEnt.tile, k);
    int n_vert = vert.EXTsum.rows();
    for (auto &eAdj : eEnt.ListAdj) {
      KDelaunayTile<Tint> const &tile2 = DT.l_tiles[eAdj.iOrb].tile;
      KDelaunayTile<Tint> tile_adj{tile2.EXT * eAdj.eBigMat,
                                   tile2.INT * eAdj.eBigMat};
      KDelaunayTileVertices<Tint> vert_adj =
          GetKDelaunayTileVertices(tile_adj, k);
      ContainerMatrix<Tint> cont(vert_adj.EXTsum);
      Face eIncEff(n_vert);
      for (int i_vert = 0; i_vert < n_vert; i_vert++) {
        MyVector<Tint> V = GetMatrixRow(vert.EXTsum, i_vert);
        std::optional<size_t> opt = cont.GetIdx_v(V);
        if (opt) {
          eIncEff[i_vert] = 1;
        }
      }
      if (eIncEff != eAdj.eInc) {
        std::cerr << "K_DELAUNAY: Inconsistency in the adjacency of the "
                     "order-k tiling\n";
        throw TerminalException{1};
      }
    }
  }
}

// The squared k-covering radius: the maximum of the squared radius of the
// spheres of the tiles.
template <typename T, typename Tint, typename Tgroup>
T MaximumKSphereRadiusSquared(MyMatrix<T> const &GramMat,
                              KDelaunayTesselation<Tint, Tgroup> const &DT) {
  T TheCovSqr(0);
  for (auto &eEnt : DT.l_tiles) {
    CP<T> cp = GetKDelaunayTileSphere<T, Tint>(GramMat, eEnt.tile);
    if (cp.SquareRadius > TheCovSqr) {
      TheCovSqr = cp.SquareRadius;
    }
  }
  return TheCovSqr;
}

/*
  The volume identity of the tiling. In sum coordinates the lattice
  translations act on the tiles as the translations by k Z^n, and the
  number of tiles of the orbit of a tile D per fundamental domain of k Z^n
  is |Aut(Q)| / |Stab(D)|, so

      sum_D |Aut(Q)| / |Stab(D)| vol(D) = k^n

  over the orbit representatives. The order of Aut(Q) is obtained from its
  permutation action on an invariant vector family.
 */
template <typename T, typename Tint, typename Tgroup> struct KTilingVolumeCheck {
  T total;
  T expected;
  bool correct;
  typename Tgroup::Tint order_aut;
};

template <typename T, typename Tint, typename Tgroup>
KTilingVolumeCheck<T, Tint, Tgroup>
CheckKDelaunayTesselationVolume(MyMatrix<T> const &GramMat,
                                KDelaunayTesselation<Tint, Tgroup> const &DT,
                                std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using TintGroup = typename Tgroup::Tint;
  int n = GramMat.rows();
  int k = DT.k;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<T, Tint>(GramMat, os);
  int n_shv = SHV.rows();
  std::unordered_map<MyVector<Tint>, size_t> map;
  for (int i_row = 0; i_row < n_shv; i_row++) {
    map[GetMatrixRow(SHV, i_row)] = i_row;
  }
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
  std::vector<Telt> ListPermGen;
  for (auto &eGen : ListGen) {
    std::vector<Tidx> l_pos(n_shv);
    for (int i_row = 0; i_row < n_shv; i_row++) {
      MyVector<Tint> V = GetMatrixRow(SHV, i_row);
      MyVector<Tint> Vimg = eGen.transpose() * V;
      auto iter = map.find(Vimg);
      if (iter == map.end()) {
        std::cerr << "K_DELAUNAY: CheckKDelaunayTesselationVolume: the "
                     "automorphism does not preserve the invariant family\n";
        throw TerminalException{1};
      }
      l_pos[i_row] = iter->second;
    }
    ListPermGen.push_back(Telt(l_pos));
  }
  Tgroup GRPaut(ListPermGen, n_shv);
  TintGroup order_aut = GRPaut.size();
  T total(0);
  for (auto &eEnt : DT.l_tiles) {
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(eEnt.tile, k);
    MyMatrix<T> EXTsum_T = UniversalMatrixConversion<T, Tint>(vert.EXTsum);
    T vol = rev_search::Kernel_VolumePolytope<T>(EXTsum_T);
    TintGroup order_stab = eEnt.GRP.size();
    TintGroup n_orbit = order_aut / order_stab;
    T n_orbit_T = UniversalScalarConversion<T, TintGroup>(n_orbit);
    total += n_orbit_T * vol;
  }
  T expected(1);
  for (int u = 0; u < n; u++) {
    expected *= T(k);
  }
  return {total, expected, total == expected, order_aut};
}

template <typename T, typename Tint, typename Tgroup>
void WriteKTilingVolumeCheckGAP(std::ostream &os_out,
                                KTilingVolumeCheck<T, Tint, Tgroup> const &x) {
  os_out << "return rec(total:=" << x.total << ", expected:=" << x.expected
         << ", order_aut:=" << x.order_aut
         << ", correct:=" << (x.correct ? "true" : "false") << ");\n";
}

// The k-covering density: with r the k-covering radius, Theta_k = kappa_n
// r^n / sqrt(det Q) is the average number of balls of radius r centered at
// the lattice points covering a point of space, at least k. The normalized
// value Theta_k / k is at least 1.
template <typename T> struct ResultKCov {
  int k;
  ResultCov<T> res;
  double KCovDensityNormalized;
};

template <typename T>
ResultKCov<T> ComputeKCoveringDensityFromDimDetCov(int TheDim, int k, T TheDet,
                                                   T TheCov) {
  ResultCov<T> res = ComputeCoveringDensityFromDimDetCov<T>(TheDim, TheDet, TheCov);
  double norm = res.CovDensity / static_cast<double>(k);
  return {k, res, norm};
}

template <typename T> std::string to_stringGAP(ResultKCov<T> const &x) {
  return "rec(k:=" + std::to_string(x.k) +
         ", KCovDensityNormalized:=" + std::to_string(x.KCovDensityNormalized) +
         ", KCovDensity:=" + std::to_string(x.res.CovDensity) +
         ", KCoveringRadius:=" + std::to_string(x.res.CoveringRadius) +
         ", DeterminantLattice:=" + std::to_string(x.res.DeterminantLattice) +
         ", VolumeBall:=" + std::to_string(x.res.VolumeBall) +
         ", TheDim:=" + std::to_string(x.res.TheDim) +
         ", TheDet:=" + std::format("{}", x.res.TheDet) +
         ", TheCov:=" + std::format("{}", x.res.TheCov) + ")";
}

template <typename Tint, typename Tgroup>
void WriteKDelaunayTesselationGAP(std::ostream &os_out,
                                  KDelaunayTesselation<Tint, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  int k = DT.k;
  os_out << "rec(k:=" << k << ", ListTile:=[";
  size_t n_tile = DT.l_tiles.size();
  for (size_t i_tile = 0; i_tile < n_tile; i_tile++) {
    KDelaunay_Entry<Tint, Tgroup> const &eEnt = DT.l_tiles[i_tile];
    MyMatrix<Tint> const &EXT = eEnt.tile.EXT;
    if (i_tile > 0) {
      os_out << ",";
    }
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(eEnt.tile, k);
    os_out << "rec(EXT:=" << StringMatrixGAP(EXT) << ",\n";
    os_out << "INT:=" << StringMatrixGAP(eEnt.tile.INT) << ",\n";
    os_out << "Vertices:=" << StringMatrixGAP(vert.EXTsum) << ",\n";
    std::vector<Telt> LGen = eEnt.GRP.SmallGeneratingSet();
    if (LGen.empty()) {
      LGen.push_back(eEnt.GRP.get_identity());
    }
    std::string str_perm = "[", str_matr = "[";
    bool IsFirst = true;
    for (auto &eElt : LGen) {
      if (!IsFirst) {
        str_perm += ",";
        str_matr += ",";
      }
      IsFirst = false;
      MyMatrix<Tint> M = RepresentVertexPermutation(EXT, EXT, eElt);
      str_perm += GapStyleString(eElt);
      str_matr += StringMatrixGAP(M);
    }
    str_perm += "]";
    str_matr += "]";
    os_out << "TheStab:=rec(PermutationStabilizer:=Group(" << str_perm
           << "), ListMatrGens:=" << str_matr << "), ";
    os_out << "Adjacencies:=[";
    IsFirst = true;
    for (auto &eAdj : eEnt.ListAdj) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "rec(iTile:=" << (eAdj.iOrb + 1) << ", eInc:=[";
      std::vector<int> V = FaceToVector<int>(eAdj.eInc);
      for (size_t u = 0; u < V.size(); u++) {
        if (u > 0) {
          os_out << ",";
        }
        os_out << (V[u] + 1);
      }
      os_out << "],\n";
      os_out << "eBigMat:=" << StringMatrixGAP(eAdj.eBigMat) << ")";
    }
    os_out << "])";
  }
  os_out << "])";
}

template <typename Tint, typename Tgroup>
void WriteKDelaunayTesselationPYTHON(
    std::ostream &os_out, KDelaunayTesselation<Tint, Tgroup> const &DT) {
  using Telt = typename Tgroup::Telt;
  int k = DT.k;
  os_out << "{\"k\":" << k << ", \"ListTile\":[";
  size_t n_tile = DT.l_tiles.size();
  for (size_t i_tile = 0; i_tile < n_tile; i_tile++) {
    KDelaunay_Entry<Tint, Tgroup> const &eEnt = DT.l_tiles[i_tile];
    MyMatrix<Tint> const &EXT = eEnt.tile.EXT;
    int n_vert = EXT.rows();
    if (i_tile > 0) {
      os_out << ",";
    }
    KDelaunayTileVertices<Tint> vert = GetKDelaunayTileVertices(eEnt.tile, k);
    os_out << "{\"EXT\":" << StringMatrixPYTHON(EXT);
    os_out << ", \"INT\":" << StringMatrixPYTHON(eEnt.tile.INT);
    os_out << ", \"Vertices\":" << StringMatrixPYTHON(vert.EXTsum);
    std::vector<Telt> LGen = eEnt.GRP.SmallGeneratingSet();
    os_out << ", \"TheStab\":[";
    bool IsFirst = true;
    for (auto &eElt : LGen) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "[";
      for (int i_vert = 0; i_vert < n_vert; i_vert++) {
        if (i_vert > 0) {
          os_out << ",";
        }
        os_out << OnPoints(i_vert, eElt);
      }
      os_out << "]";
    }
    os_out << "], \"Adjacencies\":[";
    IsFirst = true;
    for (auto &eAdj : eEnt.ListAdj) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "{\"iTile\":" << eAdj.iOrb << ", \"eInc\":[";
      std::vector<int> V = FaceToVector<int>(eAdj.eInc);
      for (size_t u = 0; u < V.size(); u++) {
        if (u > 0) {
          os_out << ",";
        }
        os_out << V[u];
      }
      os_out << "], \"eBigMat\":" << StringMatrixPYTHON(eAdj.eBigMat) << "}";
    }
    os_out << "]}";
  }
  os_out << "]}";
}

template <typename T, typename Tint, typename Tgroup>
void WriteKDelaunayTesselation(std::string const &OutFormat,
                               std::ostream &os_out, MyMatrix<T> const &GramMat,
                               KDelaunayTesselation<Tint, Tgroup> const &DT) {
  if (OutFormat == "nothing") {
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteKDelaunayTesselationGAP(os_out, DT);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "PYTHON") {
    return WriteKDelaunayTesselationPYTHON(os_out, DT);
  }
  if (OutFormat == "GAP_Covering") {
    T TheCovSqr = MaximumKSphereRadiusSquared<T, Tint, Tgroup>(GramMat, DT);
    T TheDet = DeterminantMat(GramMat);
    int TheDim = GramMat.rows();
    ResultKCov<T> x =
        ComputeKCoveringDensityFromDimDetCov<T>(TheDim, DT.k, TheDet, TheCovSqr);
    os_out << "return " << to_stringGAP(x) << ";\n";
    return;
  }
  std::cerr << "K_DELAUNAY: WriteKDelaunayTesselation failed for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

template <typename Tint, typename Tgroup>
void write_k_delaunay_to_file(std::string const &file,
                              KDelaunayTesselation<Tint, Tgroup> const &DT) {
  std::ofstream ofs(file);
  boost::archive::text_oarchive oa(ofs);
  oa << DT;
}

template <typename Tint, typename Tgroup>
KDelaunayTesselation<Tint, Tgroup>
read_k_delaunay_from_file(std::string const &file) {
  KDelaunayTesselation<Tint, Tgroup> DT;
  std::ifstream ifs(file);
  boost::archive::text_iarchive ia(ifs);
  ia >> DT;
  return DT;
}

// The order-k tiling, computed or read from CacheFile ("none" disables the
// caching). A cached tiling of another order is not used.
template <typename T, typename Tint, typename Tgroup>
KDelaunayTesselation<Tint, Tgroup>
get_k_delaunay_tessellation_serial(DataLattice<T, Tint, Tgroup> &data,
                                   int const &k, std::string const &CacheFile,
                                   int const &max_runtime_second,
                                   [[maybe_unused]] std::ostream &os) {
  auto compute = [&]() -> KDelaunayTesselation<Tint, Tgroup> {
    auto f_incorrect =
        [&]([[maybe_unused]] KDelaunay_Obj<Tint, Tgroup> const &x) -> bool {
      return false;
    };
    std::optional<KDelaunayTesselation<Tint, Tgroup>> opt =
        EnumerationKDelaunayTiles<T, Tint, Tgroup, decltype(f_incorrect)>(
            data, k, f_incorrect, max_runtime_second);
    KDelaunayTesselation<Tint, Tgroup> DT =
        unfold_opt(opt, "The order-k Delaunay tiling");
#ifdef SANITY_CHECK_K_DELAUNAY
    check_k_delaunay_tessellation(DT, os);
#endif
    return DT;
  };
  if (CacheFile == "none") {
    return compute();
  }
  if (FILE_IsExistingFile(CacheFile)) {
    KDelaunayTesselation<Tint, Tgroup> DT =
        read_k_delaunay_from_file<Tint, Tgroup>(CacheFile);
    if (DT.k == k) {
      return DT;
    }
#ifdef DEBUG_K_DELAUNAY
    os << "K_DELAUNAY: the cache file holds a tiling of order " << DT.k
       << " instead of " << k << ", recomputing\n";
#endif
  }
  KDelaunayTesselation<Tint, Tgroup> DT = compute();
  write_k_delaunay_to_file<Tint, Tgroup>(CacheFile, DT);
  return DT;
}

// clang-format off
#endif  // SRC_K_COVERINGS_LATTICEKDELAUNAY_H_
// clang-format on
