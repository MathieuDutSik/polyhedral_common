// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_COVERINGMAXDET_H_
#define SRC_DELAUNAY_COVERINGMAXDET_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  Optimization of the sphere covering density over one iso-Delaunay domain,
  by the determinant maximization ("MAXDET") formulation of Dutour Sikiric,
  Schuermann and Vallentin, "A generalization of Voronoi's reduction theory
  and its application", Duke Math. J. 142 (2008), Section 6.2.

  --- The problem ---

  Let D be a Delone (iso-Delaunay) subdivision and T(D) its secondary cone
  inside a T-space T, described by the linear inequalities a . x >= 0 with
  Q(x) = sum_u x_u A_u for a basis (A_1, ..., A_dim) of T. The covering
  density

      Theta(Q) = kappa_n mu(Q)^n / sqrt(det Q)

  is invariant under scaling of Q, so minimizing it over T(D) amounts to
  normalizing mu(Q) <= 1 and maximizing det Q.

  The normalization is a linear matrix inequality. By Proposition 6.1 of the
  reference (due to Delone, Dolbilin, Ryshkov and Stogrin), an n-dimensional
  simplex conv{0, v_1, ..., v_n} has circumradius at most 1 for the inner
  product (y, z) = y^T Q z if and only if

                 / 4          (v_1,v_1)  ...  (v_n,v_n) \
                 | (v_1,v_1)  (v_1,v_1)  ...  (v_1,v_n) |
      BR(Q)  =   |    .           .       .        .    |  >= 0.
                 |    .           .        .       .    |
                 \ (v_n,v_n)  (v_n,v_1)  ...  (v_n,v_n) /

  Every entry is a linear function of Q, hence of x. All the vertices of a
  Delone polytope lie on its (empty) circumsphere, so one full-dimensional
  simplex spanned by vertices of the polytope carries the same circumradius
  as the polytope itself. Taking one such simplex per orbit of Delone
  polytopes (Proposition 6.2) gives mu(Q) <= 1 as a single block-diagonal
  LMI, and the covering optimization over the domain becomes

      minimize   -log det Q(x)
      subject to a . x >= 0 for the facets, BR_j(x) >= 0 for the orbits.

  This is a determinant maximization problem: convex, and with a strictly
  feasible starting point readily available (an interior point of the cone,
  scaled down until every circumradius is below 1).

  --- Periodic point sets ---

  Nothing above uses that the point set is a lattice: only that the vertices
  of the Delone polytopes are known vectors and that Q is the only unknown.
  For the periodic point sets of PeriodicStructures.h the cosets are fixed
  rational vectors, so this stays true, and the vertex matrices EXT that the
  periodic Delone geometry produces carry the coordinates scaled by the coset
  denominator N. Those scaled coordinates are used as they are, which is the
  same optimization problem: the point set is then the one scaled by N, and
  the L-type cone, the LMI and the optimizer are unchanged.

  The covering density is not unchanged, and this is what point_density
  accounts for. With the true squared circumradius R_true^2 the computation
  returns R^2 = N^2 R_true^2, and a periodic set with m cosets has m points
  per fundamental cell of Z^n, of covolume sqrt(det Q). Hence

      Theta = kappa_n R_true^n m / sqrt(det Q)
            = [ kappa_n R^n / sqrt(det Q) ] * m / N^n,

  so point_density = m / N^n is the factor by which the lattice formula has
  to be multiplied. It is 1 for a lattice, and being a constant it does not
  move the optimizer, only the density reported at it.

  --- What is exact and what is not ---

  The problem data (the simplices, the LMI blocks, the cone inequalities) is
  built in exact arithmetic over T. The optimization itself is a numerical
  interior point method in floating point: the reported optimum is a
  floating point Gram matrix and the covering density is evaluated from it.
  It is therefore an estimate, not a certificate. Certification would require
  rounding the optimum to a rational Gram matrix and re-evaluating the
  covering radius and the cone membership exactly; MaxSquaredCircumRadius and
  IsInCone below are already type-generic for that purpose.

  --- Solver ---

  The solver is the barrier path-following method below. It depends on
  nothing outside this repository, so the covering optimization is available
  in a plain build.
 */

#ifdef DEBUG
#define DEBUG_COVERING_MAXDET
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_COVERING_MAXDET
#endif

#ifdef TIMINGS
#define TIMINGS_COVERING_MAXDET
#endif

namespace covering_maxdet {

// ---------------------------------------------------------------------------
// Problem data
// ---------------------------------------------------------------------------

/*
  The covering optimization problem attached to one iso-Delaunay domain.
  Built exactly over T; converted to a floating point type for the solve.
 */
template <typename T> struct CoveringData {
  // The dimension n of the ambient space.
  int n;
  // The number of optimization variables, that is dim T.
  int dim;
  // The basis (A_1, ..., A_dim) of the T-space: Q(x) = sum_u x_u A_u.
  std::vector<MyMatrix<T>> ListA;
  // The defining inequalities of the L-type cone in T-space coordinates,
  // one row a per inequality, the constraint being a . x >= 0.
  MyMatrix<T> Alin;
  // One full-dimensional simplex per orbit of Delone polytopes. Entry j is
  // the n x n matrix whose rows are v_1, ..., v_n, the differences from a
  // chosen vertex of the polytope to n further affinely independent ones.
  std::vector<MyMatrix<T>> ListSimplex;
  // The factor by which the lattice covering density formula has to be
  // multiplied, see the header comment: 1 for a lattice, m / N^n for a
  // periodic point set of m cosets with common denominator N.
  T point_density;
};

/*
  Convert the problem data to another (typically floating point) scalar type.
 */
template <typename Tout, typename T>
CoveringData<Tout> ConvertCoveringData(CoveringData<T> const &cd) {
  std::vector<MyMatrix<Tout>> ListA;
  ListA.reserve(cd.ListA.size());
  for (auto &eMat : cd.ListA) {
    ListA.push_back(UniversalMatrixConversion<Tout, T>(eMat));
  }
  std::vector<MyMatrix<Tout>> ListSimplex;
  ListSimplex.reserve(cd.ListSimplex.size());
  for (auto &eMat : cd.ListSimplex) {
    ListSimplex.push_back(UniversalMatrixConversion<Tout, T>(eMat));
  }
  MyMatrix<Tout> Alin = UniversalMatrixConversion<Tout, T>(cd.Alin);
  Tout point_density = UniversalScalarConversion<Tout, T>(cd.point_density);
  return {cd.n,          cd.dim, std::move(ListA), std::move(Alin),
          std::move(ListSimplex), point_density};
}

/*
  A full-dimensional simplex spanned by vertices of a Delone polytope, as the
  n differences from the first vertex to n affinely independent others.

  EXT is the homogeneous vertex matrix of the Delone geometry: a leading
  column of 1 followed by the (possibly scaled, see the header comment)
  coordinates. Delone polytopes are n-dimensional, so n such differences
  always exist.
 */
template <typename T, typename Tint>
MyMatrix<T> GetSimplexFromDelaunay(MyMatrix<Tint> const &EXT) {
  int n = EXT.cols() - 1;
  int nbVert = EXT.rows();
  MyMatrix<T> Sel(n, n);
  int n_sel = 0;
  for (int iVert = 1; iVert < nbVert && n_sel < n; iVert++) {
    for (int i = 0; i < n; i++) {
      Sel(n_sel, i) =
          UniversalScalarConversion<T, Tint>(EXT(iVert, i + 1) - EXT(0, i + 1));
    }
    MyMatrix<T> Test = Sel.topRows(n_sel + 1);
    if (RankMat(Test) == n_sel + 1) {
      n_sel++;
    }
  }
  if (n_sel < n) {
    std::cerr << "COVERING_MAXDET: GetSimplexFromDelaunay: the Delone polytope "
                 "spans only "
              << n_sel << " dimensions instead of " << n
              << ", so it is not a Delone polytope of a positive definite "
                 "form\n";
    throw TerminalException{1};
  }
  return Sel;
}

/*
  Assemble the covering optimization problem of an iso-Delaunay domain. FAC
  holds the defining inequalities of its L-type cone in T-space coordinates
  (redundant ones are harmless, they only add 1 x 1 barrier blocks).
 */
template <typename T, typename Tint, typename Tgroup>
CoveringData<T> BuildCoveringData(IsoDelaunayDomain<T, Tint, Tgroup> const &x,
                                  LinSpaceMatrix<T> const &LinSpa,
                                  MyMatrix<T> const &FAC,
                                  T const &point_density,
                                  [[maybe_unused]] std::ostream &os) {
  int n = LinSpa.n;
  int dim = LinSpa.ListMat.size();
  if (FAC.cols() != dim) {
    std::cerr << "COVERING_MAXDET: BuildCoveringData: the inequalities have "
              << FAC.cols() << " columns while the T-space has dimension "
              << dim
              << ". The T-space used here has to be the one the iso-Delaunay "
                 "domain was computed with\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListSimplex;
  ListSimplex.reserve(x.DT.l_dels.size());
  for (auto &eDel : x.DT.l_dels) {
    ListSimplex.push_back(GetSimplexFromDelaunay<T, Tint>(eDel.EXT));
  }
#ifdef DEBUG_COVERING_MAXDET
  os << "COVERING_MAXDET: BuildCoveringData n=" << n << " dim=" << dim
     << " |ListSimplex|=" << ListSimplex.size() << " |FAC|=" << FAC.rows()
     << "\n";
#endif
  return {n, dim, LinSpa.ListMat, FAC, std::move(ListSimplex), point_density};
}

// ---------------------------------------------------------------------------
// Geometry of a fixed Gram matrix
// ---------------------------------------------------------------------------

/*
  The squared circumradius of the simplex conv{0, v_1, ..., v_n}, the rows of
  V being the v_i, with respect to the inner product given by Q.

  The circumcenter c = sum_j a_j v_j is determined by 2 (v_i, c) = (v_i, v_i),
  that is G a = q / 2 with G_ij = (v_i, v_j) and q_i = (v_i, v_i). Hence
  R^2 = (c, c) = a^T G a = q^T G^{-1} q / 4.
 */
template <typename T>
T SquaredCircumRadiusSimplex(MyMatrix<T> const &V, MyMatrix<T> const &Q) {
  int n = V.rows();
  MyMatrix<T> G = V * Q * V.transpose();
  MyVector<T> q(n);
  for (int i = 0; i < n; i++) {
    q(i) = G(i, i);
  }
  MyMatrix<T> Ginv = Inverse(G);
  MyVector<T> a = Ginv * q;
  T sum(0);
  for (int i = 0; i < n; i++) {
    sum += q(i) * a(i);
  }
  return sum / T(4);
}

/*
  The squared covering radius mu(Q)^2 of the point set, that is the largest
  squared circumradius over the orbits of Delone polytopes. Valid for Q in
  the closed L-type cone of the domain, where the stored simplices really are
  the Delone polytopes.
 */
template <typename T>
T MaxSquaredCircumRadius(CoveringData<T> const &cd, MyMatrix<T> const &Q) {
  bool is_first = true;
  T maxval(0);
  for (auto &V : cd.ListSimplex) {
    T val = SquaredCircumRadiusSimplex(V, Q);
    if (is_first || val > maxval) {
      maxval = val;
      is_first = false;
    }
  }
  return maxval;
}

// Whether the point x lies in the closed L-type cone.
template <typename T>
bool IsInCone(CoveringData<T> const &cd, MyVector<T> const &x) {
  int n_ineq = cd.Alin.rows();
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    T sum(0);
    for (int u = 0; u < cd.dim; u++) {
      sum += cd.Alin(i_ineq, u) * x(u);
    }
    if (sum < 0) {
      return false;
    }
  }
  return true;
}

// The Gram matrix Q(x) = sum_u x_u A_u.
template <typename T>
MyMatrix<T> GetGramMatrix(CoveringData<T> const &cd, MyVector<T> const &x) {
  MyMatrix<T> Q = ZeroMatrix<T>(cd.n, cd.n);
  for (int u = 0; u < cd.dim; u++) {
    Q += x(u) * cd.ListA[u];
  }
  return Q;
}

// ---------------------------------------------------------------------------
// The linear matrix inequalities
// ---------------------------------------------------------------------------

/*
  The circumradius block BR of the simplex of index i_simplex, as the dim + 1
  symmetric matrices of its affine expansion: entry 0 is the constant term and
  entry u + 1 is the coefficient of x_u.
 */
template <typename T>
std::vector<MyMatrix<T>> GetLmiBlock(CoveringData<T> const &cd, int i_simplex) {
  int n = cd.n;
  MyMatrix<T> const &V = cd.ListSimplex[i_simplex];
  std::vector<MyMatrix<T>> ListB;
  ListB.reserve(cd.dim + 1);
  // The constant term of the affine expansion. Named Bconst rather than B0
  // because B0 is a baud-rate macro of <termios.h>.
  MyMatrix<T> Bconst = ZeroMatrix<T>(n + 1, n + 1);
  Bconst(0, 0) = T(4);
  ListB.push_back(std::move(Bconst));
  for (int u = 0; u < cd.dim; u++) {
    MyMatrix<T> M = V * cd.ListA[u] * V.transpose();
    MyMatrix<T> B = ZeroMatrix<T>(n + 1, n + 1);
    for (int i = 0; i < n; i++) {
      B(0, i + 1) = M(i, i);
      B(i + 1, 0) = M(i, i);
      for (int j = 0; j < n; j++) {
        B(i + 1, j + 1) = M(i, j);
      }
    }
    ListB.push_back(std::move(B));
  }
  return ListB;
}

/*
  The full system handed to the solver: the affine map of the objective
  (the Gram matrix), the LMI blocks and the linear inequalities.
 */
template <typename Tfloat> struct MaxdetSystem {
  int n;
  int dim;
  // Q(x) = sum_u ListA[u] x_u, the argument of the log-determinant objective.
  std::vector<MyMatrix<Tfloat>> ListA;
  // ListLmi[j][0] + sum_u ListLmi[j][u + 1] x_u >= 0.
  std::vector<std::vector<MyMatrix<Tfloat>>> ListLmi;
  // Alin x >= 0, one row per inequality.
  MyMatrix<Tfloat> Alin;
  // The parameter of the barrier, which bounds the duality gap by nu / t.
  Tfloat nu;
};

template <typename Tfloat>
MaxdetSystem<Tfloat> GetMaxdetSystem(CoveringData<Tfloat> const &cd) {
  std::vector<std::vector<MyMatrix<Tfloat>>> ListLmi;
  int n_simplex = cd.ListSimplex.size();
  ListLmi.reserve(n_simplex);
  for (int i = 0; i < n_simplex; i++) {
    ListLmi.push_back(GetLmiBlock(cd, i));
  }
  // The barrier is -log det Q - sum_j log det BR_j - sum_k log(a_k . x),
  // of parameter n + n_simplex (n + 1) + n_ineq.
  Tfloat nu = Tfloat(cd.n) + Tfloat(n_simplex) * Tfloat(cd.n + 1) +
              Tfloat(cd.Alin.rows());
  return {cd.n, cd.dim, cd.ListA, std::move(ListLmi), cd.Alin, nu};
}

// ---------------------------------------------------------------------------
// Numerical primitives
// ---------------------------------------------------------------------------

template <typename Tfloat> struct CholeskyInfo {
  bool is_pd;
  Tfloat logdet;
  MyMatrix<Tfloat> inv;
};

/*
  Cholesky-based positive definiteness test, log-determinant and inverse.
  Eigen's LLT can succeed on a matrix that is only barely indefinite, so the
  diagonal of the factor is checked as well.
 */
template <typename Tfloat>
CholeskyInfo<Tfloat> GetCholeskyInfo(MyMatrix<Tfloat> const &M, bool need_inv) {
  int m = M.rows();
  Eigen::LLT<MyMatrix<Tfloat>> llt(M);
  if (llt.info() != Eigen::Success) {
    return {false, Tfloat(0), MyMatrix<Tfloat>(0, 0)};
  }
  MyMatrix<Tfloat> L = llt.matrixL();
  Tfloat logdet(0);
  for (int i = 0; i < m; i++) {
    if (!(L(i, i) > 0)) {
      return {false, Tfloat(0), MyMatrix<Tfloat>(0, 0)};
    }
    logdet += Tfloat(2) * std::log(L(i, i));
  }
  MyMatrix<Tfloat> inv(0, 0);
  if (need_inv) {
    inv = llt.solve(IdentityMat<Tfloat>(m));
  }
  return {true, logdet, std::move(inv)};
}

// The value at x of the affine matrix map ListM[0] + sum_u ListM[u + 1] x_u.
template <typename Tfloat>
MyMatrix<Tfloat> EvalAffine(std::vector<MyMatrix<Tfloat>> const &ListM,
                            MyVector<Tfloat> const &x) {
  MyMatrix<Tfloat> M = ListM[0];
  int dim = x.size();
  for (int u = 0; u < dim; u++) {
    M += x(u) * ListM[u + 1];
  }
  return M;
}

// The value at x of the linear matrix map sum_u ListM[u] x_u.
template <typename Tfloat>
MyMatrix<Tfloat> EvalLinear(std::vector<MyMatrix<Tfloat>> const &ListM,
                            MyVector<Tfloat> const &x) {
  int dim = x.size();
  MyMatrix<Tfloat> M = x(0) * ListM[0];
  for (int u = 1; u < dim; u++) {
    M += x(u) * ListM[u];
  }
  return M;
}

/*
  Add to g and H the gradient and the Hessian of w * (-log det M(x)), where
  M(x) is the affine map of ListM read from index shift on:

      d_u  (-log det M) = -tr(M^{-1} M_u)
      d_uv (-log det M) =  tr(M^{-1} M_u M^{-1} M_v).
 */
template <typename Tfloat>
void AccumulateLogDet(Tfloat const &w, MyMatrix<Tfloat> const &Minv,
                      std::vector<MyMatrix<Tfloat>> const &ListM, int shift,
                      MyVector<Tfloat> &g, MyMatrix<Tfloat> &H) {
  int dim = g.size();
  std::vector<MyMatrix<Tfloat>> W(dim);
  for (int u = 0; u < dim; u++) {
    W[u] = Minv * ListM[u + shift];
  }
  for (int u = 0; u < dim; u++) {
    g(u) -= w * W[u].trace();
    for (int v = u; v < dim; v++) {
      // tr(W_u W_v), computed without forming the product.
      Tfloat val = w * (W[u].array() * W[v].transpose().array()).sum();
      H(u, v) += val;
      if (v > u) {
        H(v, u) += val;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// The barrier method
// ---------------------------------------------------------------------------

template <typename Tfloat> struct MaxdetOptions {
  // Stop the outer loop when the bound nu / t on the duality gap of
  // -log det Q falls below this.
  Tfloat tol_gap;
  // Multiplicative increase of t between outer iterations.
  Tfloat mu;
  // Stop the Newton loop when half the squared Newton decrement is below
  // this.
  Tfloat newton_tol;
  // Armijo parameter of the backtracking line search.
  Tfloat alpha;
  int max_outer;
  int max_newton;
  int max_backtrack;
};

template <typename Tfloat> MaxdetOptions<Tfloat> GetDefaultMaxdetOptions() {
  return {Tfloat(1e-9), Tfloat(15), Tfloat(1e-13), Tfloat(0.01),
          200,          200,        80};
}

template <typename Tfloat> struct MaxdetResult {
  // Whether the solve reached the requested tolerance.
  bool success = false;
  // Whether the returned point is a strictly feasible point of the domain,
  // which is a weaker and more useful statement: the barrier method never
  // leaves the feasible set, so cov_density is then a covering density the
  // point set really attains, an upper bound on the optimum of the domain
  // and hence on the covering density of the point set. A solve that ran
  // out of iterations still delivers that.
  bool has_point = false;
  std::string message;
  // The optimizer in T-space coordinates and the Gram matrix it defines.
  MyVector<Tfloat> x;
  MyMatrix<Tfloat> Q;
  // The value -log det Q of the objective and the bound nu / t reached on
  // its distance to the optimum. A negative gap_bound means that no such
  // bound was reached, the solve having stopped before its first outer
  // iteration.
  Tfloat obj = Tfloat(0);
  Tfloat gap_bound = Tfloat(-1);
  // The geometry of Q: mu(Q)^2, det Q and the covering density, together
  // with the point density factor that entered the last.
  Tfloat cov_radius_sq = Tfloat(0);
  Tfloat det = Tfloat(0);
  Tfloat cov_density = Tfloat(0);
  Tfloat point_density = Tfloat(1);
  int n_newton = 0;
};

/*
  The value at x of the weighted barrier

      phi_t(x) = -(1 + 1/t) log det Q(x)
                 - (1/t) ( sum_j log det BR_j(x) + sum_k log(a_k . x) ),

  which is the objective -log det Q plus 1/t times the barrier of the
  feasible set. Dividing the usual t phi_0 + psi by t keeps the coefficients
  of order one, which matters once t reaches 1e12.

  Returns nothing when x is not strictly feasible.
 */
template <typename Tfloat>
std::optional<Tfloat> BarrierValue(MaxdetSystem<Tfloat> const &sys,
                                   MyVector<Tfloat> const &x,
                                   Tfloat const &t) {
  Tfloat w_bar = Tfloat(1) / t;
  MyMatrix<Tfloat> Q = EvalLinear(sys.ListA, x);
  CholeskyInfo<Tfloat> ci_Q = GetCholeskyInfo(Q, false);
  if (!ci_Q.is_pd) {
    return {};
  }
  Tfloat val = -(Tfloat(1) + w_bar) * ci_Q.logdet;
  for (auto &ListB : sys.ListLmi) {
    MyMatrix<Tfloat> F = EvalAffine(ListB, x);
    CholeskyInfo<Tfloat> ci_F = GetCholeskyInfo(F, false);
    if (!ci_F.is_pd) {
      return {};
    }
    val -= w_bar * ci_F.logdet;
  }
  int n_ineq = sys.Alin.rows();
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    Tfloat s = sys.Alin.row(i_ineq).dot(x);
    if (!(s > 0)) {
      return {};
    }
    val -= w_bar * std::log(s);
  }
  return val;
}

/*
  Which of the conditions of BarrierValue a point fails, as a sentence. Used
  to report why a point that is strictly feasible in exact arithmetic is no
  longer so once converted: the answer says whether the Gram matrix, one of
  the circumradius blocks or one of the facets is what the conversion lost,
  which are three quite different things.
 */
template <typename Tfloat>
std::string DescribeInfeasibility(MaxdetSystem<Tfloat> const &sys,
                                  MyVector<Tfloat> const &x) {
  MyMatrix<Tfloat> Q = EvalLinear(sys.ListA, x);
  CholeskyInfo<Tfloat> ci_Q = GetCholeskyInfo(Q, false);
  if (!ci_Q.is_pd) {
    return "the Gram matrix is not positive definite, of smallest diagonal "
           "entry " +
           std::to_string(Q.diagonal().minCoeff()) + " and largest " +
           std::to_string(Q.diagonal().maxCoeff());
  }
  int n_block = sys.ListLmi.size();
  for (int j = 0; j < n_block; j++) {
    MyMatrix<Tfloat> F = EvalAffine(sys.ListLmi[j], x);
    CholeskyInfo<Tfloat> ci_F = GetCholeskyInfo(F, false);
    if (!ci_F.is_pd) {
      std::ostringstream os_str;
      os_str << "the circumradius block of orbit " << j
             << " is not positive definite; its largest entry is "
             << F.cwiseAbs().maxCoeff()
             << " and the Gram matrix has diagonal entries between "
             << Q.diagonal().minCoeff() << " and " << Q.diagonal().maxCoeff();
      return os_str.str();
    }
  }
  int n_ineq = sys.Alin.rows();
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    Tfloat s = sys.Alin.row(i_ineq).dot(x);
    if (!(s > 0)) {
      return "facet " + std::to_string(i_ineq) + " evaluates to " +
             std::to_string(s) + ", of largest coefficient " +
             std::to_string(sys.Alin.row(i_ineq).cwiseAbs().maxCoeff());
    }
  }
  return "no condition fails on re-examination, which should not happen";
}

/*
  The value, gradient and Hessian of the same barrier at a strictly feasible
  x. The caller has already checked feasibility through BarrierValue.
 */
template <typename Tfloat> struct BarrierDerivatives {
  Tfloat val;
  MyVector<Tfloat> g;
  MyMatrix<Tfloat> H;
};

template <typename Tfloat>
BarrierDerivatives<Tfloat> BarrierEvalDerivatives(
    MaxdetSystem<Tfloat> const &sys, MyVector<Tfloat> const &x,
    Tfloat const &t) {
  int dim = sys.dim;
  Tfloat w_bar = Tfloat(1) / t;
  MyVector<Tfloat> g = ZeroVector<Tfloat>(dim);
  MyMatrix<Tfloat> H = ZeroMatrix<Tfloat>(dim, dim);
  //
  MyMatrix<Tfloat> Q = EvalLinear(sys.ListA, x);
  CholeskyInfo<Tfloat> ci_Q = GetCholeskyInfo(Q, true);
  Tfloat val = -(Tfloat(1) + w_bar) * ci_Q.logdet;
  AccumulateLogDet(Tfloat(1) + w_bar, ci_Q.inv, sys.ListA, 0, g, H);
  //
  for (auto &ListB : sys.ListLmi) {
    MyMatrix<Tfloat> F = EvalAffine(ListB, x);
    CholeskyInfo<Tfloat> ci_F = GetCholeskyInfo(F, true);
    val -= w_bar * ci_F.logdet;
    AccumulateLogDet(w_bar, ci_F.inv, ListB, 1, g, H);
  }
  //
  int n_ineq = sys.Alin.rows();
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    MyVector<Tfloat> a = sys.Alin.row(i_ineq).transpose();
    Tfloat s = a.dot(x);
    val -= w_bar * std::log(s);
    Tfloat inv_s = Tfloat(1) / s;
    Tfloat inv_s2 = inv_s * inv_s;
    for (int u = 0; u < dim; u++) {
      g(u) -= w_bar * a(u) * inv_s;
      for (int v = 0; v < dim; v++) {
        H(u, v) += w_bar * a(u) * a(v) * inv_s2;
      }
    }
  }
  return {val, std::move(g), std::move(H)};
}

/*
  Newton's method on phi_t from a strictly feasible x, with a backtracking
  line search that rejects any step leaving the feasible set. Returns the
  number of accepted steps, or -1 when the search stalls.
 */
template <typename Tfloat>
int CenteringStep(MaxdetSystem<Tfloat> const &sys, MyVector<Tfloat> &x,
                  Tfloat const &t, MaxdetOptions<Tfloat> const &opts,
                  [[maybe_unused]] std::ostream &os) {
  int dim = sys.dim;
  int n_step = 0;
  for (int iter = 0; iter < opts.max_newton; iter++) {
    BarrierDerivatives<Tfloat> bd = BarrierEvalDerivatives(sys, x, t);
    // The Newton direction. The Hessian is positive definite in exact
    // arithmetic; a small multiple of the identity restores that when the
    // conditioning at large t degrades it.
    Eigen::LLT<MyMatrix<Tfloat>> llt(bd.H);
    MyVector<Tfloat> delta;
    if (llt.info() == Eigen::Success) {
      delta = llt.solve(-bd.g);
    } else {
      MyMatrix<Tfloat> Hreg = bd.H;
      Tfloat eps = Tfloat(1e-12) * (Tfloat(1) + bd.H.diagonal().maxCoeff());
      for (int u = 0; u < dim; u++) {
        Hreg(u, u) += eps;
      }
      Eigen::LDLT<MyMatrix<Tfloat>> ldlt(Hreg);
      delta = ldlt.solve(-bd.g);
    }
    Tfloat gdd = bd.g.dot(delta);
    // The squared Newton decrement, which bounds the suboptimality on
    // phi_t by lambda^2 / 2 once the iterate is in the quadratic region.
    Tfloat lambda2 = -gdd;
    if (!(lambda2 > 0)) {
      // Not a descent direction: the numerical Hessian is no longer usable.
      return n_step;
    }
    if (lambda2 / Tfloat(2) < opts.newton_tol) {
      return n_step;
    }
    Tfloat step(1);
    bool found = false;
    for (int i_back = 0; i_back < opts.max_backtrack; i_back++) {
      MyVector<Tfloat> x_new = x + step * delta;
      std::optional<Tfloat> opt_val = BarrierValue(sys, x_new, t);
      if (opt_val && *opt_val <= bd.val + opts.alpha * step * gdd) {
        x = x_new;
        found = true;
        break;
      }
      step /= Tfloat(2);
    }
    if (!found) {
      return n_step == 0 ? -1 : n_step;
    }
    n_step++;
  }
  return n_step;
}

/*
  The covering density of a Gram matrix whose squared circumradius maximum is
  cov_radius_sq: the lattice value kappa_n mu^n / sqrt(det Q), corrected by
  the point density factor of the header comment.
 */
template <typename Tfloat>
Tfloat GetCoveringDensity(CoveringData<Tfloat> const &cd, Tfloat const &det,
                          Tfloat const &cov_radius_sq) {
  ResultCov<Tfloat> rc = ComputeCoveringDensityFromDimDetCov<Tfloat>(
      cd.n, det, cov_radius_sq);
  return Tfloat(rc.CovDensity) * cd.point_density;
}

/*
  The geometric quantities attached to a point x of the L-type cone: the Gram
  matrix, its determinant, the squared covering radius and the resulting
  covering density. The convergence fields are left to the caller.
 */
template <typename Tfloat>
MaxdetResult<Tfloat> GetResultFromPoint(CoveringData<Tfloat> const &cd,
                                        MyVector<Tfloat> const &x) {
  MaxdetResult<Tfloat> res;
  res.x = x;
  res.Q = GetGramMatrix(cd, x);
  CholeskyInfo<Tfloat> ci = GetCholeskyInfo(res.Q, false);
  if (!ci.is_pd) {
    res.success = false;
    res.message = "the point does not define a positive definite form";
    return res;
  }
  res.success = true;
  res.has_point = true;
  res.obj = -ci.logdet;
  res.det = std::exp(ci.logdet);
  res.cov_radius_sq = MaxSquaredCircumRadius(cd, res.Q);
  res.cov_density = GetCoveringDensity(cd, res.det, res.cov_radius_sq);
  res.point_density = cd.point_density;
  return res;
}

/*
  The barrier method proper: follow the central path by repeatedly centering
  and multiplying t by mu, until nu / t bounds the gap on -log det Q below
  the tolerance.

  x_start must be strictly feasible; GetStartingPoint below produces one.
 */
template <typename Tfloat>
MaxdetResult<Tfloat> SolveCoveringMaxdet(CoveringData<Tfloat> const &cd,
                                         MyVector<Tfloat> const &x_start,
                                         MaxdetOptions<Tfloat> const &opts,
                                         std::ostream &os) {
#ifdef TIMINGS_COVERING_MAXDET
  MicrosecondTime time;
#endif
  MaxdetSystem<Tfloat> sys = GetMaxdetSystem(cd);
  MyVector<Tfloat> x = x_start;
  std::optional<Tfloat> opt_start = BarrierValue(sys, x, Tfloat(1));
  if (!opt_start) {
    MaxdetResult<Tfloat> res;
    res.success = false;
    res.message = "the starting point is not strictly feasible";
    res.x = x;
    res.n_newton = 0;
    return res;
  }
  Tfloat t(1);
  int n_newton = 0;
  bool stalled = false;
  int i_outer = 0;
  for (; i_outer < opts.max_outer; i_outer++) {
    if (sys.nu / t < opts.tol_gap) {
      break;
    }
    int n_step = CenteringStep(sys, x, t, opts, os);
    if (n_step < 0) {
      stalled = true;
      break;
    }
    n_newton += n_step;
#ifdef DEBUG_COVERING_MAXDET
    os << "COVERING_MAXDET: outer=" << i_outer << " t=" << t
       << " n_step=" << n_step << " gap_bound=" << (sys.nu / t) << "\n";
#endif
    t *= opts.mu;
  }
#ifdef TIMINGS_COVERING_MAXDET
  os << "|COVERING_MAXDET: SolveCoveringMaxdet|=" << time << "\n";
#endif
  MaxdetResult<Tfloat> res = GetResultFromPoint(cd, x);
  res.gap_bound = sys.nu / t;
  res.n_newton = n_newton;
  if (!res.success) {
    return res;
  }
  // A stall of the line search still leaves a feasible iterate, so the
  // result is reported; only the gap bound stops improving.
  if (stalled) {
    res.success = false;
    res.message = "the line search stalled before reaching the tolerance, the "
                  "reported point is feasible but the gap bound is loose";
  } else if (i_outer == opts.max_outer) {
    res.success = false;
    res.message = "the outer loop hit max_outer before reaching the tolerance";
  } else {
    res.message = "converged";
  }
  return res;
}

/*
  The exact margins of a point of the T-space: how far it is inside the
  L-type cone and how far every squared circumradius is below 1. Both are
  positive exactly when the point is strictly feasible.

  Kept exact and reported because a domain can be so thin that the margins,
  while positive, are below what a double can resolve. That is a property of
  the domain, not a bug, and the caller has to be able to tell the two
  apart.
 */
template <typename T> struct FeasibilityMargins {
  // min over the facets of a . x, and the facet attaining it.
  T cone_margin;
  // 1 - max over the orbits of the squared circumradius.
  T radius_margin;
};

template <typename T>
FeasibilityMargins<T> GetFeasibilityMargins(CoveringData<T> const &cd,
                                            MyVector<T> const &x) {
  T cone_margin(0);
  bool is_first = true;
  int n_ineq = cd.Alin.rows();
  for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
    T sum(0);
    for (int u = 0; u < cd.dim; u++) {
      sum += cd.Alin(i_ineq, u) * x(u);
    }
    if (is_first || sum < cone_margin) {
      cone_margin = sum;
      is_first = false;
    }
  }
  MyMatrix<T> Q = GetGramMatrix(cd, x);
  T radius_margin = T(1) - MaxSquaredCircumRadius(cd, Q);
  return {cone_margin, radius_margin};
}

/*
  A strictly feasible starting point: an interior point of the L-type cone,
  scaled down until every circumradius is strictly below 1. The scaling is
  legitimate because Q -> Q / s divides every squared circumradius by s while
  leaving the cone conditions untouched, both being homogeneous in Q.

  Returns nothing when the point, strictly feasible in exact arithmetic, is
  no longer so once converted to floating point, DescribeInfeasibility
  saying which condition was lost. It is a normal outcome the caller has to
  handle, and it is distinguished from a broken construction, which is a
  programming error and still throws.

  What is observed is not a small margin -- the margins can be perfectly
  healthy, of the order of 1e-2 -- but a huge dynamic range. The interior
  point the linear program returns is then a very anisotropic form, with
  diagonal entries spanning 1e7 to 1e15 on the same domain, and the
  circumradius blocks of the different orbits live at scales too far apart
  for a double to hold at once: the block of an orbit whose simplex is short
  in that metric loses its positive definiteness while the block setting the
  scale is fine.

  The scaling by 1 / (2 max_sq) cannot help, being exactly scale invariant:
  R^2(cQ) = c R^2(Q), so Q / (2 R^2(Q)) does not depend on the size of Q,
  only on its shape.

  Such a form is not what a domain of a point set looks like; it is what a
  badly chosen representative of one looks like. The domains themselves have
  bounded coordinates, and the anisotropy comes from a caller that reached
  the domain by composing many flips without ever reducing -- see the drift
  discussed in CoveringRecordSearch.h, which the walk there handles by
  restarting rather than by carrying on.
 */
template <typename T, typename Tfloat>
std::optional<MyVector<Tfloat>>
GetStartingPoint(CoveringData<T> const &cd, CoveringData<Tfloat> const &cd_f,
                 std::ostream &os) {
  MyVector<T> x0 = GetGeometricallyUniqueInteriorPoint(cd.Alin, os);
  MyMatrix<T> Q0 = GetGramMatrix(cd, x0);
  T max_sq = MaxSquaredCircumRadius(cd, Q0);
  // s = 2 max_sq brings the largest squared circumradius to 1/2. Note that
  // this scaling cannot be followed by a normalization of the coordinates:
  // mu(Q) <= 1 is precisely the constraint that is not scale invariant, so
  // rescaling x afterwards would move the circumradii back out of range.
  // The coordinates need none anyway, s growing with x0 by homogeneity.
  T s = T(2) * max_sq;
  MyVector<T> x_scaled = x0 / s;
  FeasibilityMargins<T> margins = GetFeasibilityMargins(cd, x_scaled);
  // The exact construction has to be strictly feasible; that it is not is a
  // programming error rather than a property of the input.
  if (margins.cone_margin <= 0 || margins.radius_margin <= 0) {
    std::cerr << "COVERING_MAXDET: GetStartingPoint: the scaled interior point "
                 "is not strictly feasible in exact arithmetic, which "
                 "contradicts the construction. cone_margin="
              << margins.cone_margin
              << " radius_margin=" << margins.radius_margin << "\n";
    throw TerminalException{1};
  }
  MyVector<Tfloat> x_ret(cd.dim);
  for (int u = 0; u < cd.dim; u++) {
    x_ret(u) = UniversalScalarConversion<Tfloat, T>(x_scaled(u));
  }
  MaxdetSystem<Tfloat> sys = GetMaxdetSystem(cd_f);
  std::optional<Tfloat> opt = BarrierValue(sys, x_ret, Tfloat(1));
  if (!opt) {
    os << "COVERING_MAXDET: GetStartingPoint: the starting point is strictly "
          "feasible exactly, with cone_margin="
       << margins.cone_margin << " and radius_margin=" << margins.radius_margin
       << ", but not after conversion to floating point: "
       << DescribeInfeasibility(sys, x_ret)
       << ". The domain cannot be optimized at this precision\n";
    return {};
  }
  return x_ret;
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

/*
  The optimization result as a GAP record. The Gram matrix is floating point,
  as is everything the numerical solve produces.
 */
template <typename Tfloat>
void WriteCoveringOptimumGAP(std::string const &FileName,
                             MaxdetResult<Tfloat> const &res) {
  std::ofstream os_out(FileName);
  // showpoint so that a value that happens to be integral still carries a
  // decimal point and is read back as a float rather than as an integer.
  os_out << std::setprecision(17) << std::showpoint;
  os_out << "return rec(success:=" << (res.success ? "true" : "false");
  os_out << ", has_point:=" << (res.has_point ? "true" : "false");
  os_out << ", message:=\"" << res.message << "\"";
  os_out << ", covering_density:=" << res.cov_density;
  os_out << ", point_density:=" << res.point_density;
  os_out << ", covering_radius_sq:=" << res.cov_radius_sq;
  os_out << ", det:=" << res.det;
  os_out << ", minus_log_det:=" << res.obj;
  os_out << ", gap_bound:=" << res.gap_bound;
  os_out << ", n_newton:=" << res.n_newton;
  os_out << ", x:=[";
  for (int u = 0; u < res.x.size(); u++) {
    if (u > 0) {
      os_out << ", ";
    }
    os_out << res.x(u);
  }
  os_out << "], GramMat:=[";
  for (int i = 0; i < res.Q.rows(); i++) {
    if (i > 0) {
      os_out << ", ";
    }
    os_out << "[";
    for (int j = 0; j < res.Q.cols(); j++) {
      if (j > 0) {
        os_out << ", ";
      }
      os_out << res.Q(i, j);
    }
    os_out << "]";
  }
  os_out << "]);\n";
}

// clang-format off
}  // namespace covering_maxdet
// clang-format on

// clang-format off
#endif  // SRC_DELAUNAY_COVERINGMAXDET_H_
// clang-format on
