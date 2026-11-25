// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_MINIMALCOVERING_H_
#define SRC_POLYNORM_POLYNORM_MINIMALCOVERING_H_

// clang-format off
#include "PolyNorm_Covering.h"
#include "POLY_DirectTriangulation.h"
#include <algorithm>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
# define DEBUG_POLYNORM_MINCOV
#endif

#ifdef SANITY_CHECK
# define SANITY_CHECK_POLYNORM_MINCOV
#endif

#ifdef TIMINGS
# define TIMINGS_POLYNORM_MINCOV
#endif

/*
  Minimal covering polytopes.

  A convex body K is covering when K + Z^n = R^n, and minimal covering when
  no convex body L strictly contained in K is covering (Codenotti, Freyer,
  Krivokuca, "Minimal covering bodies and Brunn-Minkowski type inequalities
  for the covering radius"). Every covering body contains a minimal one
  (their Proposition 2.1) and a minimal covering body is a polytope
  (Theorem 2.4). Given a covering polytope P, the code below returns a
  minimal covering polytope Q contained in P.

  Removing vertices (Q = conv of a subset of the vertices of P) is not
  enough: the rectangle [0,3/2] x [0,1] contains the unit square, but the
  convex hull of three of its vertices has area 3/4 and does not cover. The
  vertices have to be cut. Lemma 2.2 of the paper says that a covering
  polytope K is minimal exactly when for every vertex v, every u exposing v
  and every eps > 0, the truncation
      K_{u,eps} = K cap {x : <x,u> <= <v,u> - eps}
  is not covering. Its proof also shows that whether some truncation of v
  is covering does not depend on u. So the algorithm is:
  --- Shrink P to mu(P) P around its isobarycenter, which is covering of
      covering radius 1 and contained in P.
  --- Reduce K by the candidate of smallest volume among the following,
      and restart on the new polytope:
      (1) remove a vertex, K -> conv(other vertices), if that is covering
          (exact covering radius of PolyNorm_Covering.h);
      (2) push a facet inward, to the deepest covering translate;
      (3) cut a vertex v: with u the sum of the facet normals at v (it
          exposes v), t* = min{t : K cap {<x,u> <= t} is covering}, and v
          is cut at t* if t* < <v,u>.
  --- When no vertex can be cut, K is minimal covering. A covering K of
      volume 1 is a lattice tile and is minimal without further check.
  The cuts (3) alone are correct but need not terminate: on the rectangle
  [0,3/2] x [0,1] they produce a vertex that gets flatter and flatter on
  an edge of the square [1/2,3/2] x [0,1], with rationals of hundreds of
  digits. Steps (1) and (2) add no new facet direction, and with them the
  rectangle takes one push. Taking the first candidate that applies
  instead of the smallest one also goes on forever, on the octahedron.
  Termination is not proved, hence the limit on the number of steps.

  The threshold t*. With h(x) the minimum of <x - z, u> over the lattice
  points z with x - z in K, the points that K_t = K cap {<x,u> <= t} does
  not cover are those with h(x) > t, so t* = sup_x h(x). The set of
  uncovered points is open, so the supremum can be taken over the generic
  points, those for which no x - z is on the boundary of K. Since K is
  covering, a generic x has some z with x - z in K, and after a translation
  z = 0: x is in K with <x,u> >= s. For every other lattice point z one of
  the following holds:
  --- (a) x - z is in K and <x - z, u> >= s;
  --- (b_j) x - z is outside K through the facet j: it lies in the cone of
      the facet j (ell_j(x - z) >= ell_k(x - z) for all k) with
      ell_j(x - z) >= 1.
  The cones of the facets partition R^n, so the options (a), (b_1), ...,
  (b_m) cover all the cases. This is a disjunctive program, solved by branch
  and bound as the covering radius in PolyNorm_Covering.h:
  --- A node is a set of choices (z, option), its linear program maximizes s
      over (x, s), which bounds h on the generic points of the node.
  --- A node whose region has an empty interior contains no generic point
      and is dropped; this is decided by a second linear program that
      maximizes a common slack of the inequalities. Keeping such nodes would
      be wrong, not only slow: for the rectangle [0,3/2] x [0,1] and
      u = (1,0), the point x = 0 lies on the boundary of every translate,
      and the closed relaxation would let s go up to 3/2.
  --- At the optimum (x*, s*), the lattice points z with x* - z in K are
      enumerated (the lattice points of x* - K). If one of them, not yet in
      the node, has <x* - z, u> < s*, the node is split on it into m + 1
      children. Otherwise the node is solved and h approaches s* near x*:
      the segment from x* to a point satisfying the inequalities of the node
      strictly is made of points satisfying them strictly, and for the
      lattice points outside the node the inequalities are either strict at
      x* (x* - z outside K) or have value at least s*.
  The z of a node are lattice points of K - K, so the search is finite. The
  search stops as soon as a solved node reaches <v,u>: the vertex then
  cannot be cut.

  All the linear programs are exact (SIMPLEX_LinearProgramming), and so is
  t*. The polytopes are handled by their vertices and recentered at their
  isobarycenter in each round, so that the facets are linear forms
  ell_j(x) <= 1.
 */

// The facet inequalities of the polytope with homogeneous vertices EXT, as
// rows (b, a) meaning b + a.x >= 0.
template <typename T>
MyMatrix<T> PolyNorm_MinCovFacets(MyMatrix<T> const &EXT, std::ostream &os) {
  return DirectFacetComputationInequalities_General(EXT, "cdd", os);
}

// The vertices of the bounded polytope {x : b + a.x >= 0}, normalized to
// homogeneous coordinates (1, v). The row (1, 0, ..., 0) is added so that
// the cone is the homogenization of the polytope.
template <typename T>
MyMatrix<T> PolyNorm_MinCovVertices(MyMatrix<T> const &FAC, std::ostream &os) {
  int n_fac = FAC.rows();
  int n1 = FAC.cols();
  MyMatrix<T> FACext(n_fac + 1, n1);
  for (int i = 0; i < n_fac; i++) {
    for (int k = 0; k < n1; k++) {
      FACext(i, k) = FAC(i, k);
    }
  }
  FACext(n_fac, 0) = 1;
  for (int k = 1; k < n1; k++) {
    FACext(n_fac, k) = 0;
  }
  MyMatrix<T> EXT = DirectFacetComputationInequalities_General(FACext, "cdd",
                                                               os);
  int n_vert = EXT.rows();
  for (int i = 0; i < n_vert; i++) {
    T x0 = EXT(i, 0);
#ifdef SANITY_CHECK_POLYNORM_MINCOV
    if (x0 <= 0) {
      std::cerr << "POLYNORM_MINCOV: The polytope is not bounded\n";
      throw TerminalException{1};
    }
#endif
    for (int k = 0; k < n1; k++) {
      EXT(i, k) /= x0;
    }
  }
  return EXT;
}

// The isobarycenter of homogeneous vertices.
template <typename T> MyVector<T> PolyNorm_MinCovCenter(MyMatrix<T> const &EXT) {
  int n_vert = EXT.rows();
  int n = EXT.cols() - 1;
  MyVector<T> center = ZeroVector<T>(n);
  for (int i = 0; i < n_vert; i++) {
    for (int k = 0; k < n; k++) {
      center(k) += EXT(i, k + 1);
    }
  }
  center /= T(n_vert);
  return center;
}

// The polytope data used by the covering radius code, without computing
// the symmetry group (it is taken trivial). The covering radius does not
// need the difference body, which is left empty.
template <typename T, typename Tint, typename Tgroup>
PolyNormData<T, Tint, Tgroup>
PolyNorm_BuildDataTrivialGroup(MyMatrix<T> const &EXT_input,
                               std::ostream &os) {
  int n_vert = EXT_input.rows();
  int n = EXT_input.cols() - 1;
  PolyNormData<T, Tint, Tgroup> data;
  data.n = n;
  data.center = PolyNorm_MinCovCenter(EXT_input);
  MyMatrix<T> EXT(n_vert, n + 1);
  for (int i = 0; i < n_vert; i++) {
    EXT(i, 0) = 1;
    for (int k = 0; k < n; k++) {
      EXT(i, k + 1) = EXT_input(i, k + 1) - data.center(k);
    }
  }
  data.EXT = EXT;
  data.ListIncd = DirectFacetComputationIncidence(EXT, "cdd", os);
  data.Lmat = PolyNorm_LinearForms(EXT, data.ListIncd);
  data.LmatDiff = MyMatrix<T>(0, n);
  data.GRP = Tgroup(n_vert);
  return data;
}

// The covering radius of the polytope with homogeneous vertices EXT.
template <typename T, typename Tint, typename Tgroup>
T PolyNorm_MinCovCoveringRadius(MyMatrix<T> const &EXT, std::ostream &os) {
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildDataTrivialGroup<T, Tint, Tgroup>(EXT, os);
  return ComputePolyNormCovering<T, Tint, Tgroup>(data, os).mu;
}

// One choice of the branch and bound: the lattice point z and the option,
// -1 for (a) (x - z in K with <x - z, u> >= s) and j >= 0 for (b_j).
template <typename T, typename Tint> struct PolyNormCutChoice {
  MyVector<Tint> z;
  MyVector<T> z_T;
  int option;
};

template <typename T, typename Tint> struct PolyNormVertexCut {
  // The threshold t* in the recentered coordinates: K cap {<x,u> <= t} is
  // covering exactly when t >= t*.
  T t_star;
  // The value <v,u> at the vertex; the vertex can be cut iff t_star < s_max.
  T s_max;
  // A point of K - z around which the uncovered points of K_t accumulate
  // as t goes to t*, for the solved node attaining t*.
  MyVector<T> x_star;
  size_t n_node;
  size_t n_lp;
};

// The threshold t* of the cut of the polytope K = {x : Lmat x <= 1} in the
// direction u, the search stopping as soon as s_max is reached.
template <typename T, typename Tint> struct PolyNormVertexCutSolver {
  MyMatrix<T> const &Lmat;
  MyVector<T> const &u;
  T s_max;
  std::ostream &os;
  int n;
  int m;
  std::optional<T> s_best;
  MyVector<T> x_best;
  size_t n_node;
  size_t n_lp;

  PolyNormVertexCutSolver(MyMatrix<T> const &_Lmat, MyVector<T> const &_u,
                          T const &_s_max, std::ostream &_os)
      : Lmat(_Lmat), u(_u), s_max(_s_max), os(_os), n(_Lmat.cols()),
        m(_Lmat.rows()), n_node(0), n_lp(0) {
    x_best = ZeroVector<T>(n);
  }

  T ell(int j, MyVector<T> const &x) const {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += Lmat(j, k) * x(k);
    }
    return scal;
  }

  T scal_u(MyVector<T> const &x) const {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += u(k) * x(k);
    }
    return scal;
  }

  // The inequalities of a node over (x, s), as rows (b, a_x, a_s) meaning
  // b + a_x . x + a_s s >= 0.
  MyMatrix<T>
  NodeInequalities(std::vector<PolyNormCutChoice<T, Tint>> const &C) const {
    int n_row = 1 + m + 1;
    for (auto &eChoice : C) {
      n_row += (eChoice.option == -1) ? m + 1 : m;
    }
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_row, n + 2);
    int pos = 0;
    // s_max - s >= 0, implied by the other rows but keeps the program
    // visibly bounded.
    ListIneq(pos, 0) = s_max;
    ListIneq(pos, n + 1) = -1;
    pos++;
    // The root: 1 - ell_j(x) >= 0 and <x,u> - s >= 0.
    for (int j = 0; j < m; j++) {
      ListIneq(pos, 0) = 1;
      for (int k = 0; k < n; k++) {
        ListIneq(pos, k + 1) = -Lmat(j, k);
      }
      pos++;
    }
    for (int k = 0; k < n; k++) {
      ListIneq(pos, k + 1) = u(k);
    }
    ListIneq(pos, n + 1) = -1;
    pos++;
    for (auto &eChoice : C) {
      std::vector<T> ell_z(m);
      for (int j = 0; j < m; j++) {
        ell_z[j] = ell(j, eChoice.z_T);
      }
      if (eChoice.option == -1) {
        // 1 - ell_j(x - z) >= 0 and <x - z, u> - s >= 0
        for (int j = 0; j < m; j++) {
          ListIneq(pos, 0) = 1 + ell_z[j];
          for (int k = 0; k < n; k++) {
            ListIneq(pos, k + 1) = -Lmat(j, k);
          }
          pos++;
        }
        ListIneq(pos, 0) = -scal_u(eChoice.z_T);
        for (int k = 0; k < n; k++) {
          ListIneq(pos, k + 1) = u(k);
        }
        ListIneq(pos, n + 1) = -1;
        pos++;
      } else {
        int i = eChoice.option;
        // ell_i(x - z) - 1 >= 0
        ListIneq(pos, 0) = -ell_z[i] - 1;
        for (int k = 0; k < n; k++) {
          ListIneq(pos, k + 1) = Lmat(i, k);
        }
        pos++;
        // ell_i(x - z) - ell_j(x - z) >= 0
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
    }
    return ListIneq;
  }

  // The optimum (x*, s*) of the linear program of the node, if feasible.
  std::optional<std::pair<MyVector<T>, T>>
  SolveNodeLP(MyMatrix<T> const &ListIneq) {
    n_lp++;
    MyVector<T> ToBeMinimized = ZeroVector<T>(n + 2);
    ToBeMinimized(n + 1) = -1;
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.DirectSolution) {
      return {};
    }
#ifdef SANITY_CHECK_POLYNORM_MINCOV
    if (!eSol.DualSolution) {
      std::cerr << "POLYNORM_MINCOV: The node LP is unbounded, which the "
                << "constraint s <= s_max forbids\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> const &sol = *eSol.DirectSolution;
    MyVector<T> x(n);
    for (int k = 0; k < n; k++) {
      x(k) = sol(k);
    }
    T s = sol(n);
    return std::pair<MyVector<T>, T>{std::move(x), std::move(s)};
  }

  // Whether the region of the node has a nonempty interior: some (x, s)
  // satisfies all its inequalities strictly. The linear program maximizes
  // the common slack delta <= 1.
  bool HasInterior(MyMatrix<T> const &ListIneq) {
    n_lp++;
    int n_row = ListIneq.rows();
    MyMatrix<T> ListIneqSlack(n_row + 1, n + 3);
    for (int r = 0; r < n_row; r++) {
      for (int k = 0; k < n + 2; k++) {
        ListIneqSlack(r, k) = ListIneq(r, k);
      }
      ListIneqSlack(r, n + 2) = -1;
    }
    for (int k = 0; k < n + 2; k++) {
      ListIneqSlack(n_row, k) = 0;
    }
    ListIneqSlack(n_row, 0) = 1;
    ListIneqSlack(n_row, n + 2) = -1;
    MyVector<T> ToBeMinimized = ZeroVector<T>(n + 3);
    ToBeMinimized(n + 2) = -1;
    LpSolution<T> eSol =
        SIMPLEX_LinearProgramming(ListIneqSlack, ToBeMinimized, os);
    if (!eSol.DirectSolution) {
      return false;
    }
    return (*eSol.DirectSolution)(n + 1) > 0;
  }

  bool IsPruned(T const &s_star) const {
    if (!s_best) {
      return false;
    }
    return s_star <= *s_best || *s_best == s_max;
  }

  void Process(std::vector<PolyNormCutChoice<T, Tint>> const &C,
               MyVector<T> const &x_star, T const &s_star) {
    n_node++;
#ifdef DEBUG_POLYNORM_MINCOV
    os << "POLYNORM_MINCOV: node " << n_node << " |C|=" << C.size()
       << " s*=" << s_star << "\n";
#endif
    if (IsPruned(s_star)) {
      return;
    }
    // The lattice points z with x* - z in K: 1 - ell_j(x*) + ell_j(z) >= 0.
    MyMatrix<T> FAC(m, n + 1);
    for (int j = 0; j < m; j++) {
      FAC(j, 0) = 1 - ell(j, x_star);
      for (int k = 0; k < n; k++) {
        FAC(j, k + 1) = Lmat(j, k);
      }
    }
    std::vector<MyVector<Tint>> ListZ = GetListIntegralPoint<T, Tint>(FAC, os);
    std::optional<std::pair<T, MyVector<Tint>>> branch;
    T u_x = scal_u(x_star);
    for (auto &z : ListZ) {
      if (IsZeroVector(z)) {
        continue;
      }
      bool in_node = false;
      for (auto &eChoice : C) {
        if (eChoice.z == z) {
          in_node = true;
          break;
        }
      }
      if (in_node) {
        continue;
      }
      MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
      T val = u_x - scal_u(z_T);
      if (val < s_star) {
        if (!branch || val < branch->first) {
          branch = std::pair<T, MyVector<Tint>>{val, z};
        }
      }
    }
    if (!branch) {
      // The node is solved: h approaches s* near x*.
      s_best = s_star;
      x_best = x_star;
#ifdef DEBUG_POLYNORM_MINCOV
      os << "POLYNORM_MINCOV: solved node, s_best=" << s_star << "\n";
#endif
      return;
    }
    MyVector<Tint> const &z = branch->second;
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    struct Child {
      std::vector<PolyNormCutChoice<T, Tint>> C;
      MyVector<T> x;
      T s;
    };
    std::vector<Child> ListChild;
    for (int option = -1; option < m; option++) {
      std::vector<PolyNormCutChoice<T, Tint>> Cnew = C;
      Cnew.push_back({z, z_T, option});
      MyMatrix<T> ListIneq = NodeInequalities(Cnew);
      std::optional<std::pair<MyVector<T>, T>> opt = SolveNodeLP(ListIneq);
      if (!opt) {
        continue;
      }
      if (IsPruned(opt->second)) {
        continue;
      }
      if (!HasInterior(ListIneq)) {
        continue;
      }
      ListChild.push_back(
          {std::move(Cnew), std::move(opt->first), std::move(opt->second)});
    }
    // Best bound first, so that good incumbents are found early.
    std::stable_sort(ListChild.begin(), ListChild.end(),
                     [](Child const &a, Child const &b) -> bool {
                       return a.s > b.s;
                     });
    for (auto &eChild : ListChild) {
      Process(eChild.C, eChild.x, eChild.s);
    }
  }

  PolyNormVertexCut<T, Tint> Run() {
    std::vector<PolyNormCutChoice<T, Tint>> C;
    MyMatrix<T> ListIneq = NodeInequalities(C);
    std::optional<std::pair<MyVector<T>, T>> opt = SolveNodeLP(ListIneq);
#ifdef SANITY_CHECK_POLYNORM_MINCOV
    if (!opt) {
      std::cerr << "POLYNORM_MINCOV: The root node is infeasible\n";
      throw TerminalException{1};
    }
#endif
    Process(C, opt->first, opt->second);
#ifdef SANITY_CHECK_POLYNORM_MINCOV
    if (!s_best) {
      std::cerr << "POLYNORM_MINCOV: No solved node, the polytope is not "
                << "covering\n";
      throw TerminalException{1};
    }
#endif
    return {*s_best, s_max, x_best, n_node, n_lp};
  }
};

template <typename T> struct PolyNormMinimalCovering {
  // The covering radius of the input; the input is first shrunk by it.
  T mu_input;
  // Whether the input is covering (mu_input <= 1); the other fields are
  // set only when it is.
  bool is_covering;
  // The minimal covering polytope, homogeneous vertices and facet
  // inequalities (b, a) meaning b + a.x >= 0, in the input coordinates.
  MyMatrix<T> EXT;
  MyMatrix<T> FAC;
  // The reductions done after the shrinking: vertices removed, facets
  // pushed inward and vertices cut.
  size_t n_remove;
  size_t n_push;
  size_t n_cut;
  // Whether every vertex was checked to be uncuttable; false when the
  // number of reductions reached the limit.
  bool is_minimal;
  size_t n_node;
  size_t n_lp;
};

// The polytope with homogeneous vertices EXT and facets FAC (in the input
// coordinates) recentered at the isobarycenter of its vertices, with its
// facets as linear forms.
template <typename T> struct PolyNormMinCovCentered {
  MyVector<T> center;
  // In the coordinates y = x - center the facet b + a.x >= 0 is
  // ell(y) = -a.y / (b + a.center) <= 1.
  MyMatrix<T> Lmat;
};

template <typename T>
PolyNormMinCovCentered<T> PolyNorm_MinCovCentered(MyMatrix<T> const &EXT,
                                                  MyMatrix<T> const &FAC) {
  int n = EXT.cols() - 1;
  int m = FAC.rows();
  MyVector<T> center = PolyNorm_MinCovCenter(EXT);
  MyMatrix<T> Lmat(m, n);
  for (int j = 0; j < m; j++) {
    T val = FAC(j, 0);
    for (int k = 0; k < n; k++) {
      val += FAC(j, k + 1) * center(k);
    }
    for (int k = 0; k < n; k++) {
      Lmat(j, k) = -FAC(j, k + 1) / val;
    }
  }
  return {std::move(center), std::move(Lmat)};
}

// The deepest covering cut of the polytope in the direction u: with
// s_max = max_{x in K} <x - center, u>, the threshold t* of
// K cap {<x - center, u> <= t}. Returns the cut as an inequality (b, a) in
// the input coordinates when t* < s_max, nothing otherwise.
template <typename T, typename Tint>
std::optional<MyVector<T>>
PolyNorm_MinCovDeepestCut(MyMatrix<T> const &EXT,
                          PolyNormMinCovCentered<T> const &cent,
                          MyVector<T> const &u, size_t &n_node, size_t &n_lp,
                          std::ostream &os) {
  int n = EXT.cols() - 1;
  T s_max(0);
  for (int i = 0; i < EXT.rows(); i++) {
    T val(0);
    for (int k = 0; k < n; k++) {
      val += u(k) * (EXT(i, k + 1) - cent.center(k));
    }
    if (i == 0 || val > s_max) {
      s_max = val;
    }
  }
  PolyNormVertexCutSolver<T, Tint> solver(cent.Lmat, u, s_max, os);
  PolyNormVertexCut<T, Tint> cut = solver.Run();
  n_node += cut.n_node;
  n_lp += cut.n_lp;
#ifdef DEBUG_POLYNORM_MINCOV
  os << "POLYNORM_MINCOV: cut u=" << StringVector(u) << " s_max=" << s_max
     << " t*=" << cut.t_star << " n_node=" << cut.n_node << "\n";
#endif
  if (cut.t_star == s_max) {
    return {};
  }
  // <x - center, u> <= t*: t* + <center,u> - <x,u> >= 0.
  MyVector<T> eCut(n + 1);
  eCut(0) = cut.t_star;
  for (int k = 0; k < n; k++) {
    eCut(0) += cent.center(k) * u(k);
    eCut(k + 1) = -u(k);
  }
  return eCut;
}

// The direction exposing the vertex iv: the sum of the normals of the
// facets containing it, in the centered normalization.
template <typename T>
MyVector<T> PolyNorm_MinCovVertexDirection(MyMatrix<T> const &EXT,
                                           PolyNormMinCovCentered<T> const &cent,
                                           int iv) {
  int n = EXT.cols() - 1;
  int m = cent.Lmat.rows();
  MyVector<T> u = ZeroVector<T>(n);
  for (int j = 0; j < m; j++) {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += cent.Lmat(j, k) * (EXT(iv, k + 1) - cent.center(k));
    }
    if (scal == 1) {
      for (int k = 0; k < n; k++) {
        u(k) += cent.Lmat(j, k);
      }
    }
  }
  return u;
}

// The polytope K cap {b + a.x >= 0} by its vertices.
template <typename T>
MyMatrix<T> PolyNorm_MinCovApplyCut(MyMatrix<T> const &FAC,
                                    MyVector<T> const &eCut,
                                    std::ostream &os) {
  int n_fac = FAC.rows();
  int n1 = FAC.cols();
  MyMatrix<T> FACcut(n_fac + 1, n1);
  for (int j = 0; j < n_fac; j++) {
    for (int k = 0; k < n1; k++) {
      FACcut(j, k) = FAC(j, k);
    }
  }
  for (int k = 0; k < n1; k++) {
    FACcut(n_fac, k) = eCut(k);
  }
  return PolyNorm_MinCovVertices(FACcut, os);
}

// The kind of a reduction of the polytope.
enum class PolyNormMinCovStep { remove, push, cut };

// The volume of the polytope with homogeneous vertices EXT.
template <typename T>
T PolyNorm_MinCovVolume(MyMatrix<T> const &EXT, std::ostream &os) {
  return DirectVolumePolytope(EXT, TriangulationProgram::lrs, os);
}

// One reduction of the covering polytope K, or nothing when K is minimal
// covering. The candidates are:
// --- the removal of a vertex, K -> conv(other vertices), when that is still
//     full dimensional and covering;
// --- the deepest covering translate of each facet;
// --- the deepest covering cut of each vertex in the direction of the sum of
//     its facet normals.
// and the one of smallest volume is chosen. Taking instead the first one
// that applies leads to endless sequences of shallower and shallower
// reductions (on the rectangle [0,3/2] x [0,1], the cuts alone produce a
// vertex that gets flatter and flatter on an edge of the square
// [1/2,3/2] x [0,1]; on the octahedron the pushes of the facets behave
// alike). The removals and the pushes add no new facet direction. The vertex
// cuts decide minimality: by Lemma 2.2 of the paper, K is minimal covering
// when no vertex can be cut.
template <typename T, typename Tint, typename Tgroup>
std::optional<std::pair<PolyNormMinCovStep, MyMatrix<T>>>
PolyNorm_MinCovReduce(MyMatrix<T> const &EXT, MyMatrix<T> const &FAC,
                      size_t &n_node, size_t &n_lp, std::ostream &os) {
  int n = EXT.cols() - 1;
  int n_vert = EXT.rows();
  std::optional<std::pair<PolyNormMinCovStep, MyMatrix<T>>> best;
  T best_vol(0);
  auto insert = [&](PolyNormMinCovStep step, MyMatrix<T> const &EXTnew)
      -> void {
    T vol = PolyNorm_MinCovVolume(EXTnew, os);
#ifdef DEBUG_POLYNORM_MINCOV
    os << "POLYNORM_MINCOV: candidate of kind " << static_cast<int>(step)
       << " vol=" << vol << "\n";
#endif
    if (!best || vol < best_vol) {
      best = std::pair<PolyNormMinCovStep, MyMatrix<T>>{step, EXTnew};
      best_vol = vol;
    }
  };
  for (int iv = 0; iv < n_vert; iv++) {
    std::vector<int> ListKeep;
    for (int i = 0; i < n_vert; i++) {
      if (i != iv) {
        ListKeep.push_back(i);
      }
    }
    MyMatrix<T> EXTred = SelectRow(EXT, ListKeep);
    if (RankMat(EXTred) < n + 1) {
      continue;
    }
    T mu = PolyNorm_MinCovCoveringRadius<T, Tint, Tgroup>(EXTred, os);
#ifdef DEBUG_POLYNORM_MINCOV
    os << "POLYNORM_MINCOV: removal of vertex " << iv << " mu=" << mu
       << "\n";
#endif
    if (mu <= 1) {
      MyMatrix<T> FACred = PolyNorm_MinCovFacets(EXTred, os);
      insert(PolyNormMinCovStep::remove, PolyNorm_MinCovVertices(FACred, os));
    }
  }
  PolyNormMinCovCentered<T> cent = PolyNorm_MinCovCentered(EXT, FAC);
  int m = FAC.rows();
  for (int j = 0; j < m; j++) {
    MyVector<T> u = GetMatrixRow(cent.Lmat, j);
    std::optional<MyVector<T>> opt =
        PolyNorm_MinCovDeepestCut<T, Tint>(EXT, cent, u, n_node, n_lp, os);
    if (opt) {
      insert(PolyNormMinCovStep::push, PolyNorm_MinCovApplyCut(FAC, *opt, os));
    }
  }
  for (int iv = 0; iv < n_vert; iv++) {
    MyVector<T> u = PolyNorm_MinCovVertexDirection(EXT, cent, iv);
    std::optional<MyVector<T>> opt =
        PolyNorm_MinCovDeepestCut<T, Tint>(EXT, cent, u, n_node, n_lp, os);
    if (opt) {
      MyMatrix<T> EXTcut = PolyNorm_MinCovApplyCut(FAC, *opt, os);
#ifdef SANITY_CHECK_POLYNORM_MINCOV
      // The cut polytope has covering radius exactly 1, and pushing the cut
      // further by a little makes it not covering.
      T mu_cut = PolyNorm_MinCovCoveringRadius<T, Tint, Tgroup>(EXTcut, os);
      if (mu_cut != 1) {
        std::cerr << "POLYNORM_MINCOV: The cut polytope has covering radius "
                  << mu_cut << " instead of 1\n";
        throw TerminalException{1};
      }
      // The largest slack of the cut on the vertices of the cut polytope:
      // the cut can be moved by a fraction of it without emptying it.
      T slack(0);
      for (int i = 0; i < EXTcut.rows(); i++) {
        T val(0);
        for (int k = 0; k <= n; k++) {
          val += (*opt)(k) * EXTcut(i, k);
        }
        if (i == 0 || val > slack) {
          slack = val;
        }
      }
      MyVector<T> eCutEps = *opt;
      eCutEps(0) -= slack / T(1024);
      MyMatrix<T> EXTeps = PolyNorm_MinCovApplyCut(FAC, eCutEps, os);
      T mu_eps = PolyNorm_MinCovCoveringRadius<T, Tint, Tgroup>(EXTeps, os);
      if (mu_eps <= 1) {
        std::cerr << "POLYNORM_MINCOV: The cut is not the deepest one, "
                  << "mu_eps=" << mu_eps << "\n";
        throw TerminalException{1};
      }
#endif
      insert(PolyNormMinCovStep::cut, EXTcut);
    }
  }
  return best;
}

// A minimal covering polytope contained in the polytope with homogeneous
// vertices EXT_input, when that polytope is covering (is_covering of the
// result). At most max_step reductions are done.
template <typename T, typename Tint, typename Tgroup>
PolyNormMinimalCovering<T>
PolyNorm_MinimalCoveringPolytope(MyMatrix<T> const &EXT_input,
                                 size_t max_step, std::ostream &os) {
#ifdef TIMINGS_POLYNORM_MINCOV
  MicrosecondTime time;
#endif
  int n = EXT_input.cols() - 1;
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT_input, os);
  T mu = ComputePolyNormCovering<T, Tint, Tgroup>(data, os).mu;
#ifdef DEBUG_POLYNORM_MINCOV
  os << "POLYNORM_MINCOV: covering radius of the input mu=" << mu << "\n";
#endif
  PolyNormMinimalCovering<T> res{mu, mu <= 1, {}, {}, 0, 0, 0, false, 0, 0};
  if (!res.is_covering) {
    return res;
  }
  // The shrinking: center + mu (P - center) is in P since center is.
  int n_vert = data.EXT.rows();
  MyMatrix<T> EXT(n_vert, n + 1);
  for (int i = 0; i < n_vert; i++) {
    EXT(i, 0) = 1;
    for (int k = 0; k < n; k++) {
      EXT(i, k + 1) = data.center(k) + mu * data.EXT(i, k + 1);
    }
  }
  MyMatrix<T> FAC = PolyNorm_MinCovFacets(EXT, os);
  size_t n_step = 0;
  while (true) {
    // A covering polytope of volume 1 is a lattice tile, and is minimal:
    // a convex body strictly contained in it has a smaller volume.
    if (PolyNorm_MinCovVolume(EXT, os) == 1) {
      res.is_minimal = true;
      break;
    }
    std::optional<std::pair<PolyNormMinCovStep, MyMatrix<T>>> opt =
        PolyNorm_MinCovReduce<T, Tint, Tgroup>(EXT, FAC, res.n_node, res.n_lp,
                                               os);
    if (!opt) {
      res.is_minimal = true;
      break;
    }
    if (n_step == max_step) {
      break;
    }
    n_step++;
    switch (opt->first) {
    case PolyNormMinCovStep::remove:
      res.n_remove++;
      break;
    case PolyNormMinCovStep::push:
      res.n_push++;
      break;
    case PolyNormMinCovStep::cut:
      res.n_cut++;
      break;
    }
    EXT = std::move(opt->second);
    FAC = PolyNorm_MinCovFacets(EXT, os);
#ifdef DEBUG_POLYNORM_MINCOV
    os << "POLYNORM_MINCOV: step " << n_step << " (remove=" << res.n_remove
       << " push=" << res.n_push << " cut=" << res.n_cut << "), now "
       << EXT.rows() << " vertices and " << FAC.rows() << " facets\n";
    for (int i = 0; i < EXT.rows(); i++) {
      os << "POLYNORM_MINCOV:   vertex";
      for (int k = 0; k < n; k++) {
        os << " " << UniversalScalarConversion<double, T>(EXT(i, k + 1));
      }
      os << "\n";
    }
#endif
  }
  res.EXT = EXT;
  res.FAC = FAC;
#ifdef TIMINGS_POLYNORM_MINCOV
  os << "|POLYNORM_MINCOV: minimal covering polytope|=" << time << "\n";
#endif
  return res;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_MINIMALCOVERING_H_
// clang-format on
