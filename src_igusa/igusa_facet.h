// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_FACET_H_
#define SRC_IGUSA_IGUSA_FACET_H_

// clang-format off
#include "igusa.h"
#include <functional>
#include <optional>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA_FACET
#endif

#ifdef TIMINGS
#define TIMINGS_IGUSA_FACET
#endif

/*
  The facets of the Igusa polyhedron P = conv(I) directly, instead of as a
  side product of the local cones of the vertices.

  A facet is an inequality tr(F X) >= c valid on I with equality on a set
  of rank dim(T) - 1. Its incidence is

      Inc(F) = { X in I : tr(F X) = c }.

  If F is positive definite, tr(F X) = c bounds X (the eigenvalues of
  F^{1/2} X F^{1/2} are in [0, c]) and Inc(F) is finite. If F is only
  positive semidefinite, for example F = v v^T, Inc(F) is infinite and
  invariant under unipotent elements of the stabilizer of F.

  The code here is for the finite case. The incident points are found by
  the theory of Voronoi:
  * The minimum of tr(F X) over the Ryshkov polyhedron {X : X[v] >= 1} is
    attained at a perfect form P* (in general, a vertex of the Ryshkov
    polyhedron of the T-space), and linear programming duality gives
    F = sum_v lambda_v v v^T (projected on T) with lambda_v >= 0 over the
    minimal vectors v of P* and sum_v lambda_v = min. This minimum is at
    most c; the difference c - min is the slack of the facet.
  * For X in Inc(F) the values m_v = X[v] are integers >= 1 with
    sum_v lambda_v m_v = c. So m_v - 1 <= slack / lambda_v and there are few
    possible vectors (m_v). Each one defines the affine subspace X[v] = m_v,
    which determines X when the v with lambda_v > 0 span the T-space; in
    general the integral points of the face of P cut by these equalities
    are enumerated by splitting on the coordinates (integer programming).
 */

template <typename T, typename Tint>
MyVector<T> igusa_trace_function(IgusaSpace<T, Tint> const &space,
                                 MyMatrix<T> const &F) {
  int dim = space.dim;
  MyVector<T> f = ZeroVector<T>(dim + 1);
  for (int j = 0; j < dim; j++) {
    f(j + 1) = frobenius_inner(F, space.ListMatInt[j]);
  }
  return f;
}

/*
  The integral points x of the face of P defined by ListEqua.(1,x) = 0,
  which has to be bounded. The equalities x_k = t are added one coordinate
  at a time.
 */
template <typename T, typename Tint>
std::vector<MyVector<T>>
igusa_face_integral_points(IgusaSpace<T, Tint> &space,
                           MyMatrix<T> const &ListEqua, std::ostream &os) {
  int dim = space.dim;
  MyMatrix<T> ListExtraIneq(0, dim + 1);
  std::vector<MyVector<T>> l_point;
#ifdef DEBUG_IGUSA_FACET
  size_t n_ilp = 0;
#endif
  std::function<void(MyMatrix<T> const &, int)> f_rec =
      [&](MyMatrix<T> const &Equa, int k) -> void {
    if (k == dim) {
      // all the coordinates are fixed: the point is the solution
      std::optional<IgusaMinimum<T>> opt = igusa_integral_minimization(
          space, ZeroVector<T>(dim + 1), ListExtraIneq, Equa, os);
      if (opt) {
        l_point.push_back(opt->x);
      }
      return;
    }
    MyVector<T> g = ZeroVector<T>(dim + 1);
    g(k + 1) = T(1);
    std::optional<IgusaMinimum<T>> opt_min =
        igusa_integral_minimization(space, g, ListExtraIneq, Equa, os);
#ifdef DEBUG_IGUSA_FACET
    n_ilp++;
#endif
    if (!opt_min) {
      return;
    }
    MyVector<T> g_neg = -g;
    std::optional<IgusaMinimum<T>> opt_max =
        igusa_integral_minimization(space, g_neg, ListExtraIneq, Equa, os);
#ifdef DEBUG_IGUSA_FACET
    n_ilp++;
#endif
    T x_min = opt_min->value;
    T x_max = -unfold_opt(opt_max, "The face is not empty").value;
    int n_row = Equa.rows();
    MyMatrix<T> EquaNew(n_row + 1, dim + 1);
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j <= dim; j++) {
        EquaNew(i, j) = Equa(i, j);
      }
    }
    for (int j = 0; j <= dim; j++) {
      EquaNew(n_row, j) = T(0);
    }
    EquaNew(n_row, k + 1) = T(1);
    for (T t = x_min; t <= x_max; t += T(1)) {
      EquaNew(n_row, 0) = -t;
      f_rec(EquaNew, k + 1);
    }
  };
  f_rec(ListEqua, 0);
#ifdef DEBUG_IGUSA_FACET
  os << "IGUSA_FACET: face integral points, |points|=" << l_point.size()
     << " n_ilp=" << n_ilp << "\n";
#endif
  return l_point;
}

template <typename T> struct IgusaFacetIncidence {
  MyMatrix<T> F;
  T rhs;
  // The minimum of tr(F X) over the Ryshkov polyhedron and a minimizer
  T ryshkov_min;
  MyMatrix<T> ryshkov_minimizer;
  // The minimum of tr(F X) over I, which is rhs for a valid facet
  T igusa_min;
  std::vector<MyMatrix<T>> l_incident;
  // The rank of the family of incident points as affine points
  int rank;
  T slack;
};

/*
  The weights lambda_v >= 0 with f = sum_v lambda_v e_v, e_v(j) =
  ListMatInt[j][v], over the minimal vectors v of the minimizer P* of f over
  the Ryshkov polyhedron, from the dual of the linear program restricted to
  the cuts tight at P*.
 */
template <typename T, typename Tint>
std::vector<std::pair<MyVector<Tint>, T>>
igusa_voronoi_weights(IgusaSpace<T, Tint> const &space, MyVector<T> const &f,
                      MyMatrix<T> const &Pstar, std::ostream &os) {
  int dim = space.dim;
  std::vector<MyVector<Tint>> l_v = computeLevel_GramMat<T, Tint>(Pstar, T(1), os);
  int n_v = l_v.size();
  MyMatrix<T> Ineq(n_v, dim + 1);
  for (int i = 0; i < n_v; i++) {
    Ineq(i, 0) = T(-1);
    for (int j = 0; j < dim; j++) {
      Ineq(i, j + 1) = EvaluationQuadForm<T, Tint>(space.ListMatInt[j], l_v[i]);
    }
  }
  LpSolution<T> eSol = SIMPLEX_LinearProgramming(Ineq, f, os);
  MyVector<T> lambda = unfold_opt(eSol.DualSolution, "The dual weights");
  // The sign convention of the dual solution: the weights are of one sign
  bool all_nonpositive = true;
  for (int i = 0; i < n_v; i++) {
    if (lambda(i) > 0) {
      all_nonpositive = false;
    }
  }
  if (all_nonpositive) {
    lambda = -lambda;
  }
  std::vector<std::pair<MyVector<Tint>, T>> l_weight;
  MyVector<T> check = ZeroVector<T>(dim);
  for (int i = 0; i < n_v; i++) {
    if (lambda(i) < 0) {
      std::cerr << "IGUSA_FACET: negative dual weight\n";
      throw TerminalException{1};
    }
    if (lambda(i) > 0) {
      l_weight.push_back({l_v[i], lambda(i)});
      for (int j = 0; j < dim; j++) {
        check(j) += lambda(i) * Ineq(i, j + 1);
      }
    }
  }
  for (int j = 0; j < dim; j++) {
    if (check(j) != f(j + 1)) {
      std::cerr << "IGUSA_FACET: the dual weights do not decompose f\n";
      throw TerminalException{1};
    }
  }
  return l_weight;
}

template <typename T, typename Tint>
IgusaFacetIncidence<T> igusa_facet_incidence(IgusaSpace<T, Tint> &space,
                                             MyMatrix<T> const &F,
                                             T const &rhs, std::ostream &os) {
  int dim = space.dim;
  if (!IsPositiveDefinite(F, os)) {
    std::cerr << "IGUSA_FACET: F has to be positive definite, the incidence "
                 "is infinite otherwise\n";
    throw TerminalException{1};
  }
  MyVector<T> f = igusa_trace_function(space, F);
  MyMatrix<T> ListExtraIneq(0, dim + 1);
  MyMatrix<T> NoEqua(0, dim + 1);
  // The minimizer over the Ryshkov polyhedron has its minimal vectors among
  // the short vectors of F^{-1}: seeding the cuts with them avoids a long
  // sequence of cutting planes when F is far from the identity.
  MyMatrix<T> Finv = Inverse(F);
  T min_diag = Finv(0, 0);
  for (int i = 1; i < Finv.rows(); i++) {
    if (Finv(i, i) < min_diag) {
      min_diag = Finv(i, i);
    }
  }
  for (auto &v : computeLevel_GramMat<T, Tint>(Finv, T(3) * min_diag, os)) {
    (void)igusa_insert_cut(space, v);
  }
  std::optional<IgusaMinimum<T>> opt_lp =
      igusa_lp_minimization(space, f, ListExtraIneq, NoEqua, os);
  IgusaMinimum<T> res_lp =
      unfold_opt(opt_lp, "The minimum over the Ryshkov polyhedron");
  MyMatrix<T> Pstar = igusa_matrix(space, res_lp.x);
  // The validity of the facet: the minimum of tr(F X) over I is rhs
  std::optional<IgusaMinimum<T>> opt_ilp =
      igusa_integral_minimization(space, f, ListExtraIneq, NoEqua, os);
  IgusaMinimum<T> res_ilp = unfold_opt(opt_ilp, "The minimum over I");
  if (res_ilp.value != rhs) {
    std::cerr << "IGUSA_FACET: the minimum of tr(F X) over I is "
              << res_ilp.value << " and not rhs=" << rhs << "\n";
    throw TerminalException{1};
  }
  std::vector<std::pair<MyVector<Tint>, T>> l_weight =
      igusa_voronoi_weights(space, f, Pstar, os);
  int n_w = l_weight.size();
  T slack = rhs - res_lp.value;
#ifdef DEBUG_IGUSA_FACET
  os << "IGUSA_FACET: min over Ryshkov=" << res_lp.value << " slack=" << slack
     << " |support|=" << n_w << "\n";
#endif
  if (slack < 0) {
    std::cerr << "IGUSA_FACET: rhs is below the minimum over the Ryshkov "
                 "polyhedron, this is not a valid inequality\n";
    throw TerminalException{1};
  }
  // The rows e_v of the support, and their rank
  MyMatrix<T> Esupp(n_w, dim);
  for (int i = 0; i < n_w; i++) {
    for (int j = 0; j < dim; j++) {
      Esupp(i, j) =
          EvaluationQuadForm<T, Tint>(space.ListMatInt[j], l_weight[i].first);
    }
  }
  int rank_supp = RankMat(Esupp);
  MyMatrix<T> EsuppT = Esupp.transpose();
  // The vectors m with m_v >= 1 integral and sum lambda_v (m_v - 1) = slack
  std::vector<MyVector<T>> l_x;
  std::set<MyVector<T>> set_x;
  MyVector<T> m(n_w);
  size_t n_m = 0;
  auto f_treat = [&]() -> void {
    n_m++;
    MyMatrix<T> Equa(n_w + 1, dim + 1);
    for (int i = 0; i < n_w; i++) {
      Equa(i, 0) = -m(i);
      for (int j = 0; j < dim; j++) {
        Equa(i, j + 1) = Esupp(i, j);
      }
    }
    for (int j = 0; j <= dim; j++) {
      Equa(n_w, j) = f(j);
    }
    Equa(n_w, 0) -= rhs;
    std::vector<MyVector<T>> l_sol;
    if (rank_supp == dim) {
      // X is determined by the values m_v
      std::optional<MyVector<T>> opt = SolutionMat(EsuppT, m);
      if (opt && IsIntegralVector(*opt) &&
          IsPositiveDefinite(igusa_matrix(space, *opt), os)) {
        l_sol.push_back(*opt);
      }
    } else {
      l_sol = igusa_face_integral_points(space, Equa, os);
    }
    for (auto &x : l_sol) {
      if (set_x.count(x) == 0) {
        set_x.insert(x);
        l_x.push_back(x);
      }
    }
  };
  std::function<void(int, T const &)> f_rec = [&](int i, T const &rem) -> void {
    if (i == n_w) {
      if (rem == 0) {
        f_treat();
      }
      return;
    }
    T const &lambda = l_weight[i].second;
    for (Tint k = 0; T(k) * lambda <= rem; k++) {
      m(i) = T(1) + T(k);
      T rem2 = rem - T(k) * lambda;
      f_rec(i + 1, rem2);
    }
  };
  f_rec(0, slack);
#ifdef DEBUG_IGUSA_FACET
  os << "IGUSA_FACET: number of vectors m=" << n_m << " |Inc|=" << l_x.size()
     << " rank of the support=" << rank_supp << "\n";
#endif
  std::vector<MyMatrix<T>> l_incident;
  MyMatrix<T> AffPoints(l_x.size(), dim + 1);
  for (size_t i = 0; i < l_x.size(); i++) {
    l_incident.push_back(igusa_matrix(space, l_x[i]));
    AffPoints(i, 0) = T(1);
    for (int j = 0; j < dim; j++) {
      AffPoints(i, j + 1) = l_x[i](j);
    }
  }
  int rank = l_x.empty() ? 0 : RankMat(AffPoints);
  return {F,     rhs,
          res_lp.value, Pstar,
          res_ilp.value, std::move(l_incident),
          rank,  slack};
}

inline FullNamelist NAMELIST_GetStandard_IGUSA_FACET_INCIDENCE() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    std::map<std::string, bool> ListBoolValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["IlpMethod"] = "default";
    ListStringValues["OutFile"] = "stderr";
    // The matrix F of the facet tr(F X) >= rhs
    ListStringValues["FileFacet"] = "unset";
    // The right hand side, a rational number
    ListStringValues["FacetRhs"] = "unset";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
    BlockDATA.setListBoolValues(ListBoolValues);
    ListBlock["DATA"] = BlockDATA;
  }
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_IGUSA_IGUSA_FACET_H_
// clang-format on
