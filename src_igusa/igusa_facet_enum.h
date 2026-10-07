// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_FACET_ENUM_H_
#define SRC_IGUSA_IGUSA_FACET_ENUM_H_

// clang-format off
#include "igusa_facet.h"
#include "LatticeStabEquiCan.h"
#include "SignatureSymmetric.h"
#include "Tspace_Generation.h"
#include "fractions.h"
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA_FACET_ENUM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_IGUSA_FACET_ENUM
#endif

/*
  Enumeration of the facets of P = conv(I_n) centered on the facets, for the
  full space of symmetric matrices (T-space "Classic").

  Structure of the facets (lifting lemma). A facet tr(F X) >= c has F
  positive semidefinite. Let W be the saturation of the image of F, of rank
  r, with basis B (r x n). Then F = B^T F0 B with F0 positive definite and
  tr(F X) = tr(F0 B X B^T) depends only on the restriction of X to W.
  Every form of I_r is such a restriction, so tr(F0 Y) >= c is valid on
  I_r, and it is a facet of conv(I_r). Conversely every facet of conv(I_r)
  with F0 positive definite lifts to a facet of conv(I_n). So the orbits of
  facets of conv(I_n) are the orbits of the full rank facets of conv(I_r),
  r <= n, and the full rank facets are the bounded ones. A facet is thus
  identified by (r, c, canonical form of F0).

  The adjacency decomposition is done on the full rank facets:
  * Inc(F) by igusa_facet_incidence (finite since F is positive definite);
  * the orbits under Stab(F) = Aut(F) of the facets of the polytope
    conv(Inc(F)), the ridges of P in the facet;
  * for each ridge, with g >= 0 on Inc(F) and g = 0 on the ridge, the
    adjacent facet g + s f >= 0, f = tr(F X) - c, for the smallest s making
    it valid on I (igusa_facet_flip);
  * the type of the adjacent facet: full rank (a new facet of conv(I_n) to
    process, if not equivalent to a known one) or of rank r < n (a lift of
    a full rank facet of conv(I_r)).
  The adjacent facets of rank < n are recorded but not processed: their
  incidence is infinite.
 */

// Scale (F, rhs) by a positive factor so that F is integral with coprime
// entries.
template <typename T>
std::pair<MyMatrix<T>, T> igusa_normalize_facet(MyMatrix<T> const &F,
                                                T const &rhs) {
  FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(F);
  return {fr.TheMat, fr.TheMult * rhs};
}

template <typename T> struct IgusaFacetType {
  int rank;
  T rhs;
  // F0 on a basis of W, integral with coprime entries
  MyMatrix<T> F0;
  // rank, rhs and canonical form of F0: the key of the orbit
  std::string key;
};

template <typename T, typename Tint, typename Tgroup>
IgusaFacetType<T> igusa_facet_type(MyMatrix<T> const &F, T const &rhs,
                                   std::ostream &os) {
  MyMatrix<T> B = IntegralSpaceSaturation(RowReduction(F));
  int r = B.rows();
  MyMatrix<T> BBt = B * B.transpose();
  MyMatrix<T> BBt_inv = Inverse(BBt);
  MyMatrix<T> F0 = BBt_inv * B * F * B.transpose() * BBt_inv;
#ifdef SANITY_CHECK_IGUSA_FACET_ENUM
  if (B.transpose() * F0 * B != F) {
    std::cerr << "IGUSA_FACET_ENUM: F is not B^T F0 B\n";
    throw TerminalException{1};
  }
#endif
  std::pair<MyMatrix<T>, T> pair = igusa_normalize_facet(F0, rhs);
  // ComputeCanonicalForm gives the canonical basis B; the canonical form
  // is B F0 B^T
  MyMatrix<Tint> B_can =
      ComputeCanonicalForm<T, Tint, Tgroup>(pair.first, os);
  MyMatrix<T> B_can_T = UniversalMatrixConversion<T, Tint>(B_can);
  MyMatrix<T> can = B_can_T * pair.first * B_can_T.transpose();
  std::ostringstream os_key;
  os_key << r << " " << pair.second << " " << StringMatrixGAP(can);
  return {r, pair.second, pair.first, os_key.str()};
}

/*
  The generators g of the common stabilizer of the matrices M of ListM, which
  are functionals on the forms (they transform as M -> g^T M g), acting on
  the forms by X -> g X g^T. The first matrix has to be positive definite.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
igusa_functional_stabilizer(std::vector<MyMatrix<T>> const &ListM,
                            std::ostream &os) {
  std::vector<MyMatrix<Tint>> l_h =
      ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListM, os);
  std::vector<MyMatrix<Tint>> l_g;
  for (auto &h : l_h) {
    MyMatrix<T> h_T = UniversalMatrixConversion<T, Tint>(h);
    bool is_left = true, is_right = true;
    for (auto &M : ListM) {
      if (h_T.transpose() * M * h_T != M) {
        is_left = false;
      }
      if (h_T * M * h_T.transpose() != M) {
        is_right = false;
      }
    }
    if (is_left) {
      l_g.push_back(h);
    } else if (is_right) {
      l_g.push_back(h.transpose());
    } else {
      std::cerr << "IGUSA_FACET_ENUM: not a common automorphism\n";
      throw TerminalException{1};
    }
  }
  return l_g;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> igusa_facet_stabilizer(MyMatrix<T> const &F,
                                                   std::ostream &os) {
  return igusa_functional_stabilizer<T, Tint, Tgroup>({F}, os);
}

template <typename T> struct IgusaRidge {
  Face face;
  // g(0) + sum_j g(j+1) x_j, >= 0 on Inc(F) and = 0 on the ridge
  MyVector<T> g;
};

/*
  The orbits under the stabilizer of the facets of conv(Inc(F)). The
  points x of Inc(F) lie in the hyperplane f(x) = 0: the dual description is
  done on independent columns of the rows (1, x), and the inequality is
  extended by zero on the other columns, which gives the same values on
  the points.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<IgusaRidge<T>>
igusa_facet_ridges(IgusaSpace<T, Tint> const &space,
                   std::vector<MyVector<T>> const &l_x,
                   std::vector<MyMatrix<Tint>> const &l_g, Tgroup &GRP_out,
                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int dim = space.dim;
  int n_pt = l_x.size();
  MyMatrix<T> EXT(n_pt, dim + 1);
  std::map<MyVector<T>, int> map_pt;
  for (int i = 0; i < n_pt; i++) {
    EXT(i, 0) = T(1);
    for (int j = 0; j < dim; j++) {
      EXT(i, j + 1) = l_x[i](j);
    }
    map_pt[l_x[i]] = i;
  }
  std::vector<Telt> l_gens;
  for (auto &g : l_g) {
    MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
    std::vector<Tidx> eList(n_pt);
    for (int i = 0; i < n_pt; i++) {
      MyMatrix<T> X = igusa_matrix(space, l_x[i]);
      MyMatrix<T> Y = g_T * X * g_T.transpose();
      MyVector<T> y = igusa_coordinates(space, Y);
      eList[i] = map_pt.at(y);
    }
    l_gens.push_back(Telt(eList));
  }
  Tgroup GRP(l_gens, n_pt);
  GRP_out = GRP;
  vectface vf = DualDescriptionStandard(EXT, GRP, os);
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT);
  std::vector<int> const &cols = eSelect.ListColSelect;
  MyMatrix<T> EXTsel = SelectColumn(EXT, cols);
  std::vector<IgusaRidge<T>> l_ridge;
  for (auto &face : vf) {
    MyVector<T> ineq_sel = FindFacetInequality(EXTsel, face);
    MyVector<T> g = ZeroVector<T>(dim + 1);
    for (size_t k = 0; k < cols.size(); k++) {
      g(cols[k]) = ineq_sel(k);
    }
    l_ridge.push_back({face, g});
  }
  return l_ridge;
}

/*
  An X in I_n with the restriction B X B^T = Y, for B (r x n) primitive and
  Y in I_r.
 */
template <typename T>
MyMatrix<T> igusa_extend_form(MyMatrix<T> const &B, MyMatrix<T> const &Y) {
  int r = B.rows();
  int n = B.cols();
  MyMatrix<T> Compl = SubspaceCompletionInt(B, n);
  MyMatrix<T> U(n, n);
  for (int i = 0; i < r; i++) {
    U.row(i) = B.row(i);
  }
  for (int i = 0; i < n - r; i++) {
    U.row(r + i) = Compl.row(i);
  }
  MyMatrix<T> Z = IdentityMat<T>(n);
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < r; j++) {
      Z(i, j) = Y(i, j);
    }
  }
  MyMatrix<T> Uinv = Inverse(U);
  return Uinv * Z * Uinv.transpose();
}

/*
  A face of P given by the equality tr(G X) = cG, where tr(G X) >= cG is a
  facet of P (G positive semidefinite). Without it, the whole of P.
 */
template <typename T> struct IgusaFaceRestriction {
  bool active;
  MyMatrix<T> G;
  T cG;
};

// The points X and X + t v v^T for the smallest power of two t with h < 0.
template <typename T, typename Tint>
MyMatrix<T> igusa_ray_violation(IgusaSpace<T, Tint> const &space,
                                MyVector<T> const &h, MyMatrix<T> const &M,
                                MyVector<Tint> const &v,
                                MyMatrix<T> const &Ridge0) {
  MyMatrix<T> vvT = ZeroMatrix<T>(space.n, space.n);
  for (int i = 0; i < space.n; i++) {
    for (int j = 0; j < space.n; j++) {
      vvT(i, j) = UniversalScalarConversion<T, Tint>(v(i) * v(j));
    }
  }
  T Mv = EvaluationQuadForm<T, Tint>(M, v);
  T h0 = ILP_EvaluateRow(h, igusa_coordinates(space, Ridge0));
  T t = T(1);
  while (h0 + t * Mv >= 0) {
    t *= T(2);
  }
  return Ridge0 + t * vvT;
}

/*
  If h(x) = h(0) + sum_j h(j+1) x_j >= 0 is valid on I_n (on the points of
  the face if face.active), returns none. Otherwise returns such a point X
  with h(X) < 0. Ridge0 is a point of the face with h(Ridge0) = 0.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
igusa_violating_point(IgusaSpace<T, Tint> &space, MyVector<T> const &h,
                      MyMatrix<T> const &Ridge0,
                      IgusaFaceRestriction<T> const &face, std::ostream &os) {
  int dim = space.dim;
  MyMatrix<T> M = igusa_functional_matrix(space, h);
  MyMatrix<T> NoIneq(0, dim + 1);
  MyMatrix<T> NoEqua(0, dim + 1);
  if (face.active) {
    // On the face, h can be changed by multiples of tr(G X) - cG.
    MyVector<T> fG = igusa_trace_function(space, face.G);
    fG(0) = -face.cG;
    MyMatrix<T> Equa(1, dim + 1);
    for (int j = 0; j <= dim; j++) {
      Equa(0, j) = fG(j);
    }
    T lambda(0);
    for (int iter = 0; iter < 60; iter++) {
      MyMatrix<T> Ml = M + lambda * face.G;
      if (IsPositiveDefinite(Ml, os)) {
        MyVector<T> hl = h + lambda * fG;
        igusa_insert_seed_cuts(space, Ml, os);
        std::optional<IgusaMinimum<T>> opt =
            igusa_integral_minimization(space, hl, NoIneq, Equa, os);
        IgusaMinimum<T> res = unfold_opt(opt, "The face is not empty");
        if (res.value < 0) {
          return igusa_matrix(space, res.x);
        }
        return {};
      }
      lambda = (lambda == 0) ? T(1) : T(2) * lambda;
    }
    // M is not positive definite on the kernel of G: a vector v of that
    // kernel with M[v] < 0 gives the recession direction v v^T of the face.
    MyMatrix<T> K = NullspaceIntMat(face.G);
    MyMatrix<T> MK = K * M * K.transpose();
    if (IsPositiveSemiDefinite(MK, os)) {
      std::cerr << "IGUSA_FACET_ENUM: degenerate functional on the face\n";
      throw TerminalException{1};
    }
    MyVector<Tint> vK = igusa_non_positive_vector<T, Tint>(MK, T(0), os);
    MyVector<T> vK_T = UniversalVectorConversion<T, Tint>(vK);
    MyVector<T> v_T = K.transpose() * vK_T;
    MyVector<Tint> v = UniversalVectorConversion<Tint, T>(v_T);
    return igusa_ray_violation(space, h, M, v, Ridge0);
  }
  if (IsPositiveDefinite(M, os)) {
    igusa_insert_seed_cuts(space, M, os);
    std::optional<IgusaMinimum<T>> opt =
        igusa_integral_minimization(space, h, NoIneq, NoEqua, os);
    IgusaMinimum<T> res = unfold_opt(opt, "The minimum over I");
    if (res.value < 0) {
      return igusa_matrix(space, res.x);
    }
    return {};
  }
  if (!IsPositiveSemiDefinite(M, os)) {
    MyVector<Tint> v = igusa_non_positive_vector<T, Tint>(M, T(0), os);
    return igusa_ray_violation(space, h, M, v, Ridge0);
  }
  // M positive semidefinite and singular: restrict to the image W of M.
  MyMatrix<T> B = IntegralSpaceSaturation(RowReduction(M));
  int r = B.rows();
  MyMatrix<T> BBt = B * B.transpose();
  MyMatrix<T> BBt_inv = Inverse(BBt);
  MyMatrix<T> M0 = BBt_inv * B * M * B.transpose() * BBt_inv;
  LinSpaceMatrix<T> LinSpa_r = ComputeCanonicalSpace<T>(r);
  IgusaSpace<T, Tint> space_r =
      build_igusa_space<T, Tint>(LinSpa_r, space.params, os);
  igusa_insert_initial_cuts(space_r);
  igusa_insert_seed_cuts(space_r, M0, os);
  MyVector<T> h_r = igusa_trace_function(space_r, M0);
  h_r(0) = h(0);
  std::optional<IgusaMinimum<T>> opt =
      igusa_integral_minimization(space_r, h_r, NoIneq, NoEqua, os);
  IgusaMinimum<T> res = unfold_opt(opt, "The minimum over I_r");
  if (res.value < 0) {
    MyMatrix<T> Y = igusa_matrix(space_r, res.x);
    return igusa_extend_form(B, Y);
  }
  return {};
}

/*
  The threshold a0 = inf { a : M + a F is positive definite on K }, with K
  the rows of a basis of the subspace (the identity if K has no row), F
  positive semidefinite and positive definite on K: returns a_lo < a0 <= a_hi
  with a_hi - a_lo <= 2^-40 (exact rationals, by bisection).
 */
template <typename T>
std::pair<T, T> igusa_pd_threshold(MyMatrix<T> const &M, MyMatrix<T> const &F,
                                   MyMatrix<T> const &K, std::ostream &os) {
  auto is_pd = [&](T const &a) -> bool {
    MyMatrix<T> Ma = M + a * F;
    if (K.rows() == 0) {
      return IsPositiveDefinite(Ma, os);
    }
    MyMatrix<T> MaK = K * Ma * K.transpose();
    return IsPositiveDefinite(MaK, os);
  };
  T a_lo(-1), a_hi(1);
  while (!is_pd(a_hi)) {
    a_hi *= T(2);
  }
  while (is_pd(a_lo)) {
    a_lo *= T(2);
  }
  T eps = T(1);
  for (int i = 0; i < 40; i++) {
    eps /= T(2);
  }
  while (a_hi - a_lo > eps) {
    T a_mid = (a_lo + a_hi) / T(2);
    if (is_pd(a_mid)) {
      a_hi = a_mid;
    } else {
      a_lo = a_mid;
    }
  }
  return {a_lo, a_hi};
}

// The simplest rational (first continued fraction approximant) in [a, b].
template <typename T> T igusa_simplest_rational(T const &a, T const &b) {
  for (auto &q : get_sequence_continuous_fraction_approximant(b)) {
    if (a <= q && q <= b) {
      return q;
    }
  }
  return b;
}

template <typename T> struct IgusaFlipResult {
  MyMatrix<T> F;
  T rhs;
  // true if no point beyond the threshold was found: the new inequality has
  // a singular matrix (on the kernel of G for a rotation inside a face)
  bool degenerate;
};

/*
  The rotation of the hyperplane tr(F X) = rhs around the face where g = 0:
  the inequality g + s f >= 0, f = tr(F X) - rhs, for the smallest s >= 0
  such that it is valid on I, or on the points of the face if face.active.
  This gives the facet of P adjacent along a ridge, or, on a facet of P,
  the facet of that facet adjacent along a face of codimension 3. The s is
  increased by the points X violating g + s f >= 0 (s := -g(X) / f(X)).
  If F is positive definite they are first looked for in the bounded region
  f <= K; otherwise g has to be positive definite, so that g + s f is.
 */
template <typename T, typename Tint>
IgusaFlipResult<T>
igusa_facet_flip(IgusaSpace<T, Tint> &space, MyMatrix<T> const &F,
                 T const &rhs, MyVector<T> const &g,
                 MyMatrix<T> const &Ridge0,
                 IgusaFaceRestriction<T> const &face, std::ostream &os) {
  int dim = space.dim;
  MyVector<T> f = igusa_trace_function(space, F);
  f(0) = -rhs;
  bool F_pd = IsPositiveDefinite(F, os);
  T K = rhs;
  // The matrix of the result is positive semidefinite (on the kernel of G
  // for a face): s is at least the threshold a0 of g + a f, and starting
  // just above it keeps all the functionals g + s f coercive.
  MyMatrix<T> Mg = igusa_functional_matrix(space, g);
  MyMatrix<T> Kface(0, space.n);
  if (face.active) {
    Kface = NullspaceIntMat(face.G);
  }
  std::pair<T, T> thr = igusa_pd_threshold(Mg, F, Kface, os);
  // Start at a distance delta = 1 above the threshold, where the problems
  // are well conditioned; if no point violates there, s* is below and the
  // distance is halved, down to 2^-10, below which the result is taken at
  // the threshold (a matrix with a kernel).
  T delta(1);
  T delta_min(1);
  for (int i = 0; i < 10; i++) {
    delta_min /= T(2);
  }
  T s = thr.second + delta;
  bool updated = false;
  MyMatrix<T> Equa(0, dim + 1);
  if (face.active) {
    MyVector<T> fG = igusa_trace_function(space, face.G);
    Equa = MyMatrix<T>(1, dim + 1);
    for (int j = 0; j <= dim; j++) {
      Equa(0, j) = fG(j);
    }
    Equa(0, 0) = -face.cG;
  }
#ifdef DEBUG_IGUSA_FACET_ENUM
  size_t iter = 0;
#endif
  auto update = [&](MyMatrix<T> const &X) -> void {
    MyVector<T> x = igusa_coordinates(space, X);
    T fX = ILP_EvaluateRow(f, x);
    T gX = ILP_EvaluateRow(g, x);
    if (fX <= 0) {
      std::cerr << "IGUSA_FACET_ENUM: a violating point on the facet\n";
      throw TerminalException{1};
    }
    s = -gX / fX;
    updated = true;
  };
  while (true) {
#ifdef DEBUG_IGUSA_FACET_ENUM
    iter++;
#endif
    MyVector<T> h = g + s * f;
    if (F_pd) {
      // the bounded region f <= K
      MyMatrix<T> Ineq(1, dim + 1);
      Ineq(0, 0) = K - f(0);
      for (int j = 0; j < dim; j++) {
        Ineq(0, j + 1) = -f(j + 1);
      }
      std::optional<IgusaMinimum<T>> opt =
          igusa_integral_minimization(space, h, Ineq, Equa, os);
      IgusaMinimum<T> res = unfold_opt(opt, "The bounded region is not empty");
      if (res.value < 0) {
        update(igusa_matrix(space, res.x));
        continue;
      }
    }
    std::optional<MyMatrix<T>> opt_X =
        igusa_violating_point(space, h, Ridge0, face, os);
    if (opt_X) {
      update(*opt_X);
      K *= T(2);
      continue;
    }
    if (!updated && delta > delta_min) {
      delta /= T(2);
      s = thr.second + delta;
      continue;
    }
    if (!updated) {
      // No point beyond the threshold: the result is at the threshold a0,
      // with a singular matrix. For a facet of P, a0 is rational and found
      // as the simplest rational of the interval, then checked.
      T q = igusa_simplest_rational(thr.first, thr.second);
      MyVector<T> hq = g + q * f;
      MyMatrix<T> Mq = igusa_functional_matrix(space, hq);
      bool is_singular_psd =
          !face.active && IsPositiveSemiDefinite(Mq, os) &&
          !IsPositiveDefinite(Mq, os);
      if (is_singular_psd &&
          !igusa_violating_point(space, hq, Ridge0, face, os)) {
        T rhs_q = -hq(0);
        std::pair<MyMatrix<T>, T> pair = igusa_normalize_facet(Mq, rhs_q);
#ifdef DEBUG_IGUSA_FACET_ENUM
        os << "IGUSA_FACET_ENUM: flip at the threshold s=" << q << "\n";
#endif
        return {pair.first, pair.second, false};
      }
#ifdef DEBUG_IGUSA_FACET_ENUM
      os << "IGUSA_FACET_ENUM: flip degenerate\n";
#endif
      return {Mq, -hq(0), true};
    }
    MyMatrix<T> G = igusa_functional_matrix(space, h);
#ifdef DEBUG_IGUSA_FACET_ENUM
    os << "IGUSA_FACET_ENUM: flip done, iter=" << iter << " s=" << s << "\n";
#endif
    T rhs_new = -h(0);
    std::pair<MyMatrix<T>, T> pair = igusa_normalize_facet(G, rhs_new);
    return {pair.first, pair.second, false};
  }
}

template <typename T> struct IgusaFacetOrbitEntry {
  MyMatrix<T> F;
  T rhs;
  std::string key;
  int n_inc;
  std::string stab_size;
  int n_ridge_orbit;
  // for each orbit of ridges: the size of the ridge and the key and rank
  // of the adjacent facet
  std::vector<int> l_ridge_size;
  std::vector<std::string> l_adj_key;
  std::vector<int> l_adj_rank;
  std::vector<T> l_adj_rhs;
};

/*
  A bounded ridge R of P inside a facet tr(G X) >= cG of rank < n:
  Inc(R) = { X in I : tr(G X) = cG, tr(H X) = cH } with H positive definite
  and tr(H X) >= cH valid on the points of the facet G. The orbit of the
  pair (R, G) under GL_n(Z) is identified by the canonical form of the pair
  (S^{-1}, G), S the sum of the forms of Inc(R), whose barycenter is in the
  relative interior of R.
 */
template <typename T> struct IgusaRidgeEntry {
  MyMatrix<T> G;
  T cG;
  MyMatrix<T> H;
  T cH;
  std::vector<MyMatrix<T>> l_inc;
  std::string key;
  int n_subridge_orbit;
  int n_unbounded;
};

template <typename T, typename Tint, typename Tgroup>
std::string igusa_ridge_key(MyMatrix<T> const &G,
                            std::vector<MyMatrix<T>> const &l_inc,
                            std::ostream &os) {
  int n = G.rows();
  MyMatrix<T> S = ZeroMatrix<T>(n, n);
  for (auto &X : l_inc) {
    S += X;
  }
  MyMatrix<T> Sinv = Inverse(S);
  std::vector<MyMatrix<T>> ListM{Sinv, G};
  MyMatrix<Tint> B = ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListM, os);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  std::ostringstream os_key;
  MyMatrix<T> Sinv_can = B_T * Sinv * B_T.transpose();
  MyMatrix<T> G_can = B_T * G * B_T.transpose();
  os_key << l_inc.size() << " " << StringMatrixGAP(Sinv_can) << " "
         << StringMatrixGAP(G_can);
  return os_key.str();
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
igusa_ridge_stabilizer(MyMatrix<T> const &G,
                       std::vector<MyMatrix<T>> const &l_inc,
                       std::ostream &os) {
  int n = G.rows();
  MyMatrix<T> S = ZeroMatrix<T>(n, n);
  for (auto &X : l_inc) {
    S += X;
  }
  std::vector<MyMatrix<T>> ListM{Inverse(S), G};
  return igusa_functional_stabilizer<T, Tint, Tgroup>(ListM, os);
}

/*
  The incidences of the full rank facets already computed, by orbit. For a
  facet equivalent to a known one, the incident forms are mapped by an
  isometry instead of being recomputed: if h^T F h = F_rep, then
  X -> h X h^T maps Inc(F_rep) onto Inc(F).
 */
template <typename T, typename Tint, typename Tgroup> struct IgusaIncidenceCache {
  std::map<std::string, std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> map;

  std::vector<MyMatrix<T>> get(IgusaSpace<T, Tint> &space, MyMatrix<T> const &F,
                               T const &rhs, std::ostream &os) {
    std::pair<MyMatrix<T>, T> pair = igusa_normalize_facet(F, rhs);
    IgusaFacetType<T> type =
        igusa_facet_type<T, Tint, Tgroup>(pair.first, pair.second, os);
    auto iter = map.find(type.key);
    if (iter == map.end()) {
      IgusaFacetIncidence<T> inc =
          igusa_facet_incidence(space, pair.first, pair.second, os);
      map[type.key] = {pair.first, inc.l_incident};
      return inc.l_incident;
    }
    MyMatrix<T> const &Frep = iter->second.first;
    std::optional<MyMatrix<Tint>> opt =
        ArithmeticEquivalence<T, Tint, Tgroup>(pair.first, Frep, os);
    MyMatrix<Tint> P = unfold_opt(opt, "The facets are equivalent");
    MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
    MyMatrix<T> h;
    if (P_T * pair.first * P_T.transpose() == Frep) {
      h = P_T.transpose();
    } else if (P_T.transpose() * pair.first * P_T == Frep) {
      h = P_T;
    } else {
      std::cerr << "IGUSA_FACET_ENUM: the equivalence does not map F\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<T>> l_inc;
    for (auto &X : iter->second.second) {
      l_inc.push_back(h * X * h.transpose());
    }
#ifdef SANITY_CHECK_IGUSA_FACET_ENUM
    for (auto &Y : l_inc) {
      if (frobenius_inner(pair.first, Y) != pair.second) {
        std::cerr << "IGUSA_FACET_ENUM: a mapped point is not incident\n";
        throw TerminalException{1};
      }
    }
#endif
    return l_inc;
  }
};

template <typename T> struct IgusaFacetEnumResult {
  std::vector<IgusaFacetOrbitEntry<T>> l_facet;
  std::vector<IgusaRidgeEntry<T>> l_ridge;
};

/*
  The adjacency decomposition on the full rank facets and on the bounded
  ridges inside the facets of lower rank, starting from ListInit, up to
  max_obj processed objects (0: no limit).
  * A full rank facet: Inc(F), the orbits of ridges, and the adjacent facet
    along each of them. A full rank one is inserted; along a ridge R with a
    facet G of lower rank, the bounded ridge (R, G) is inserted.
  * A bounded ridge (R, G): the orbits of facets S of the polytope R, and for
    each one the other ridge R' of G containing S, by a rotation inside the
    hyperplane of G. If R' is bounded, it is inserted, and so is the facet
    of P on the other side of R' (or the bounded ridge (R', Q) if that
    facet Q is of lower rank).
 */
template <typename T, typename Tint, typename Tgroup>
IgusaFacetEnumResult<T>
igusa_facet_enumeration(IgusaSpace<T, Tint> &space,
                        std::vector<std::pair<MyMatrix<T>, T>> const &ListInit,
                        int max_obj, std::ostream &os) {
  int n = space.n;
  IgusaFacetEnumResult<T> result;
  std::vector<IgusaFacetOrbitEntry<T>> &l_orbit = result.l_facet;
  std::vector<IgusaRidgeEntry<T>> &l_ridge_obj = result.l_ridge;
  std::map<std::string, int> map_key;
  std::map<std::string, int> map_ridge_key;
  IgusaFaceRestriction<T> no_face{false, {}, T(0)};
  IgusaIncidenceCache<T, Tint, Tgroup> cache;
  auto f_insert = [&](MyMatrix<T> const &F, T const &rhs) -> void {
    IgusaFacetType<T> type = igusa_facet_type<T, Tint, Tgroup>(F, rhs, os);
    if (type.rank < n || map_key.count(type.key) > 0) {
      return;
    }
    map_key[type.key] = l_orbit.size();
    std::pair<MyMatrix<T>, T> pair = igusa_normalize_facet(F, rhs);
    IgusaFacetOrbitEntry<T> entry;
    entry.F = pair.first;
    entry.rhs = pair.second;
    entry.key = type.key;
    l_orbit.push_back(entry);
  };
  auto f_insert_ridge = [&](MyMatrix<T> const &G, T const &cG,
                            MyMatrix<T> const &H, T const &cH,
                            std::vector<MyMatrix<T>> const &l_inc) -> void {
    std::pair<MyMatrix<T>, T> pairG = igusa_normalize_facet(G, cG);
    std::string key =
        igusa_ridge_key<T, Tint, Tgroup>(pairG.first, l_inc, os);
    if (map_ridge_key.count(key) > 0) {
      return;
    }
    map_ridge_key[key] = l_ridge_obj.size();
    l_ridge_obj.push_back({pairG.first, pairG.second, H, cH, l_inc, key, 0, 0});
  };
  // The points of Inc(F) where g = 0
  auto get_face = [&](std::vector<MyMatrix<T>> const &l_inc,
                      Face const &face) -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> l_ret;
    for (size_t i = 0; i < l_inc.size(); i++) {
      if (face[i] == 1) {
        l_ret.push_back(l_inc[i]);
      }
    }
    return l_ret;
  };
  auto get_coords = [&](std::vector<MyMatrix<T>> const &l_inc)
      -> std::vector<MyVector<T>> {
    std::vector<MyVector<T>> l_x;
    for (auto &X : l_inc) {
      l_x.push_back(igusa_coordinates(space, X));
    }
    return l_x;
  };
  for (auto &pair : ListInit) {
    f_insert(pair.first, pair.second);
  }
  size_t pos_facet = 0, pos_ridge = 0;
  int n_obj = 0;
  while (pos_facet < l_orbit.size() || pos_ridge < l_ridge_obj.size()) {
    if (max_obj > 0 && n_obj >= max_obj) {
      break;
    }
    n_obj++;
    if (pos_facet < l_orbit.size()) {
      // A full rank facet
      MyMatrix<T> F = l_orbit[pos_facet].F;
      T rhs = l_orbit[pos_facet].rhs;
      IgusaFacetIncidence<T> inc;
      inc.l_incident = cache.get(space, F, rhs, os);
      std::vector<MyMatrix<Tint>> l_g =
          igusa_facet_stabilizer<T, Tint, Tgroup>(F, os);
      Tgroup GRP;
      std::vector<IgusaRidge<T>> l_ridge = igusa_facet_ridges<T, Tint, Tgroup>(
          space, get_coords(inc.l_incident), l_g, GRP, os);
      std::ostringstream os_size;
      os_size << GRP.size();
      l_orbit[pos_facet].n_inc = inc.l_incident.size();
      l_orbit[pos_facet].stab_size = os_size.str();
      l_orbit[pos_facet].n_ridge_orbit = l_ridge.size();
      for (auto &ridge : l_ridge) {
        MyMatrix<T> Ridge0 = inc.l_incident[ridge.face.find_first()];
        IgusaFlipResult<T> flip =
            igusa_facet_flip(space, F, rhs, ridge.g, Ridge0, no_face, os);
        if (flip.degenerate) {
          std::cerr << "IGUSA_FACET_ENUM: degenerate flip of a facet\n";
          throw TerminalException{1};
        }
        std::pair<MyMatrix<T>, T> adj{flip.F, flip.rhs};
        IgusaFacetType<T> type =
            igusa_facet_type<T, Tint, Tgroup>(adj.first, adj.second, os);
        l_orbit[pos_facet].l_ridge_size.push_back(ridge.face.count());
        l_orbit[pos_facet].l_adj_key.push_back(type.key);
        l_orbit[pos_facet].l_adj_rank.push_back(type.rank);
        l_orbit[pos_facet].l_adj_rhs.push_back(type.rhs);
        if (type.rank == n) {
          f_insert(adj.first, adj.second);
        } else {
          f_insert_ridge(adj.first, adj.second, F, rhs,
                         get_face(inc.l_incident, ridge.face));
        }
      }
      os << "IGUSA_FACET_ENUM: facet " << pos_facet << " rhs=" << rhs
         << " |Inc|=" << l_orbit[pos_facet].n_inc
         << " |Stab on Inc|=" << l_orbit[pos_facet].stab_size
         << " ridge orbits=" << l_orbit[pos_facet].n_ridge_orbit
         << " |facets|=" << l_orbit.size() << " |ridges|=" << l_ridge_obj.size()
         << "\n";
      pos_facet++;
      continue;
    }
    // A bounded ridge (R, G) inside the facet G of lower rank
    IgusaRidgeEntry<T> obj = l_ridge_obj[pos_ridge];
    IgusaFaceRestriction<T> face{true, obj.G, obj.cG};
    std::vector<MyMatrix<Tint>> l_g =
        igusa_ridge_stabilizer<T, Tint, Tgroup>(obj.G, obj.l_inc, os);
    Tgroup GRP;
    std::vector<IgusaRidge<T>> l_sub = igusa_facet_ridges<T, Tint, Tgroup>(
        space, get_coords(obj.l_inc), l_g, GRP, os);
    int n_unbounded = 0;
    for (auto &sub : l_sub) {
      MyMatrix<T> Ridge0 = obj.l_inc[sub.face.find_first()];
      // the other ridge R' of G through the subridge
      IgusaFlipResult<T> flip_rot =
          igusa_facet_flip(space, obj.H, obj.cH, sub.g, Ridge0, face, os);
      if (flip_rot.degenerate) {
        // R' has a singular functional on the kernel of G: it is unbounded
        n_unbounded++;
        continue;
      }
      std::pair<MyMatrix<T>, T> rot{flip_rot.F, flip_rot.rhs};
      // R' is bounded (the rotation is not degenerate): its functional is
      // positive definite modulo G, and rot + lambda G is positive definite
      // above a threshold. A small lambda keeps the slack of the incidence
      // computation small.
#ifdef DEBUG_IGUSA_FACET_ENUM
      os << "IGUSA_FACET_ENUM: step bounded test\n";
#endif
      MyMatrix<T> NoK(0, n);
      std::pair<T, T> thr_l = igusa_pd_threshold(rot.first, obj.G, NoK, os);
      // an integer at distance at least 1 above the threshold, so that Hp is
      // well conditioned
      T lambda_int = UniversalFloorScalarInteger<T, T>(thr_l.second) + T(2);
      std::optional<T> opt_lambda = lambda_int;
      MyMatrix<T> Hp = rot.first + (*opt_lambda) * obj.G;
      T cHp = rot.second + (*opt_lambda) * obj.cG;
      // The facet of P on the other side of R' first: it needs only the
      // inequality of R' and a point of R' (the points of the subridge are
      // in R'). Then Inc(R') is Inc(H'') on the hyperplane of G if H'' has
      // full rank (fast, since H'' is a facet), or the incidence of the
      // valid inequality tr((G + Q) X) >= cG + cQ, positive definite and
      // tight exactly on R', if the facet Q on the other side has lower
      // rank.
#ifdef DEBUG_IGUSA_FACET_ENUM
      os << "IGUSA_FACET_ENUM: step flip across R'\n";
#endif
      MyVector<T> gp = igusa_trace_function(space, Hp);
      gp(0) = -cHp;
      IgusaFlipResult<T> flip_adj =
          igusa_facet_flip(space, obj.G, obj.cG, gp, Ridge0, no_face, os);
      if (flip_adj.degenerate) {
        std::cerr << "IGUSA_FACET_ENUM: degenerate flip from a ridge\n";
        throw TerminalException{1};
      }
      std::pair<MyMatrix<T>, T> adj{flip_adj.F, flip_adj.rhs};
      IgusaFacetType<T> type =
          igusa_facet_type<T, Tint, Tgroup>(adj.first, adj.second, os);
#ifdef DEBUG_IGUSA_FACET_ENUM
      os << "IGUSA_FACET_ENUM: step incidence of R', adjacent rank="
         << type.rank << "\n";
#endif
      std::vector<MyMatrix<T>> l_incp;
      if (type.rank == n) {
        std::vector<MyMatrix<T>> l_inc_adj =
            cache.get(space, adj.first, adj.second, os);
        for (auto &X : l_inc_adj) {
          if (frobenius_inner(obj.G, X) == obj.cG) {
            l_incp.push_back(X);
          }
        }
      } else {
        MyMatrix<T> GQ = obj.G + adj.first;
        T cGQ = obj.cG + adj.second;
        IgusaFacetIncidence<T> inc_gq =
            igusa_facet_incidence(space, GQ, cGQ, os, false);
        l_incp = inc_gq.l_incident;
      }
#ifdef DEBUG_IGUSA_FACET_ENUM
      os << "IGUSA_FACET_ENUM: step insert R', |Inc(R')|=" << l_incp.size()
         << "\n";
#endif
      f_insert_ridge(obj.G, obj.cG, Hp, cHp, l_incp);
      if (type.rank == n) {
        f_insert(adj.first, adj.second);
      } else {
        f_insert_ridge(adj.first, adj.second, MyMatrix<T>(obj.G + adj.first),
                       obj.cG + adj.second, l_incp);
      }
    }
    l_ridge_obj[pos_ridge].n_subridge_orbit = l_sub.size();
    l_ridge_obj[pos_ridge].n_unbounded = n_unbounded;
    os << "IGUSA_FACET_ENUM: ridge " << pos_ridge
       << " |Inc|=" << obj.l_inc.size() << " |Stab|=" << GRP.size()
       << " subridge orbits=" << l_sub.size()
       << " unbounded=" << n_unbounded << " |facets|=" << l_orbit.size()
       << " |ridges|=" << l_ridge_obj.size() << "\n";
    pos_ridge++;
  }
  return result;
}

template <typename T>
void WriteFacetOrbitsGAP(std::ostream &os_out,
                         IgusaFacetEnumResult<T> const &result) {
  os_out << "return rec(ListFacet:=[";
  for (size_t i = 0; i < result.l_facet.size(); i++) {
    auto const &e = result.l_facet[i];
    if (i > 0) {
      os_out << ",\n";
    }
    os_out << "rec(F:=";
    WriteMatrixGAP(os_out, e.F);
    os_out << ", rhs:=" << e.rhs << ", nInc:=" << e.n_inc
           << ", StabSizeOnInc:=" << e.stab_size
           << ", nRidgeOrbit:=" << e.n_ridge_orbit << ", ListAdjacent:=[";
    for (size_t k = 0; k < e.l_adj_key.size(); k++) {
      if (k > 0) {
        os_out << ", ";
      }
      os_out << "rec(RidgeSize:=" << e.l_ridge_size[k]
             << ", rank:=" << e.l_adj_rank[k] << ", rhs:=" << e.l_adj_rhs[k]
             << ")";
    }
    os_out << "])";
  }
  os_out << "],\n ListRidgeInLowerRankFacet:=[";
  for (size_t i = 0; i < result.l_ridge.size(); i++) {
    auto const &e = result.l_ridge[i];
    if (i > 0) {
      os_out << ",\n";
    }
    os_out << "rec(G:=";
    WriteMatrixGAP(os_out, e.G);
    os_out << ", cG:=" << e.cG << ", nInc:=" << e.l_inc.size()
           << ", nSubridgeOrbit:=" << e.n_subridge_orbit
           << ", nUnbounded:=" << e.n_unbounded << ")";
  }
  os_out << "]);\n";
}

inline FullNamelist NAMELIST_GetStandard_IGUSA_FACET_ENUMERATION() {
  std::map<std::string, SingleBlock> ListBlock;
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    std::map<std::string, bool> ListBoolValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["IlpMethod"] = "default";
    ListStringValues["OutFile"] = "stderr";
    // The initial full rank facets: a list of matrices F (ListMatrix
    // format) and, in FileInitialRhs, the list of the right hand sides
    ListStringValues["FileInitialFacets"] = "unset";
    ListStringValues["FileInitialRhs"] = "unset";
    // Stop after processing that many objects, facets of full rank and
    // bounded ridges in the facets of lower rank (0: no limit)
    ListIntValues["MaxOrbit"] = 0;
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
#endif  // SRC_IGUSA_IGUSA_FACET_ENUM_H_
// clang-format on
