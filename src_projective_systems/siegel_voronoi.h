// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PROJECTIVE_SYSTEMS_SIEGEL_VORONOI_H_
#define SRC_PROJECTIVE_SYSTEMS_SIEGEL_VORONOI_H_

// clang-format off
#include "projective_complex.h"
#include "hecke_chains.h"
#include "Tspace_Generation.h"
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SIEGEL_VORONOI
#endif

/*
  The cells of the Siegel space of a symplectic form J on Z^n cut out by
  the Voronoi decomposition of the positive definite forms in n variables, as
  MacPherson and McConnell did for Sp_4(Z) (J standard), for any J such as
  the form J_N of the paramodular group K(N) = Sp(J_N, Z).

  The matrices act on the row vectors: a configuration X of vectors of Z^n
  gives the cone C(X) of the forms sum_{v in X} p_v v^T v (p_v > 0), and
  g in Sp(J, Z) maps it to C(X g) by A -> g^T A g. The Siegel space of J
  is the cone of the positive forms A with (A J)^2 = -lambda I (for
  n = 2g, the symmetric space of Sp_{2g}(R), of dimension g (g + 1), times
  R_{>0}). Each cell
  is C(X) cap (Siegel space), for X the vectors of a face of a perfect
  cone, and its dimension is that of the symmetric space.

  The exact parts (the perfect cones, their faces and adjacencies, the
  equivalence of the configurations under Sp(J, Z)) are the ones of
  src_perfect and of projective_complex.h. Whether C(X) meets the Siegel
  space, and the dimension of the intersection, is decided numerically
  (SiegelOracle): a point of the intersection with all weights p_v at
  least a margin is searched by Levenberg-Marquardt, and the dimension is
  that of the kernel of the differential of the equations on the span of
  C(X) at that point. The results are checked afterwards by the complex
  itself (every ridge in exactly two top cells, the boundary of the
  boundary, the orbifold Euler characteristic).
 */

namespace psys {

struct SiegelIntersection {
  bool nonempty;
  // The dimension of C(X) cap (Siegel space) in the symmetric space.
  int dim;
  // The smallest weight p_v (the weights summing to 1) of the point found.
  double margin;
  // A point was found only with weights below 0.01 / |X|: near a face of
  // C(X) where the Siegel space is tangent, so not counted as a point.
  bool doubtful;
  // The weights p_v of the point found (empty when there is none).
  std::vector<double> weights;
};

class SiegelOracle {
private:
  using MatD = Eigen::MatrixXd;
  int n;
  MatD J;
  std::mt19937 gen;
  int n_start;
  // Whether the margins 1e-3 / k and 1e-4 / k are tried, to report the
  // doubtful faces.
  bool check_doubtful = true;
  // The margins (times 1 / |X|) of the points searched, from the largest.
  // A point is accepted with a margin of at least 0.01, so the margin 0.01
  // alone gives the same cells, faster; the larger ones find points
  // further inside.
  std::vector<double> l_margin{0.2, 0.05, 0.01};
  // The residual (A J)^2 - tr((A J)^2) / n I, as a vector of n^2 entries.
  static Eigen::VectorXd residual(MatD const &A, MatD const &J) {
    int n = A.rows();
    MatD AJ = A * J;
    MatD C = AJ * AJ;
    C -= (C.trace() / n) * MatD::Identity(n, n);
    return Eigen::Map<Eigen::VectorXd>(C.data(), n * n);
  }
  // The differential of the residual at A in the direction H.
  static Eigen::VectorXd differential(MatD const &A, MatD const &H, MatD const &J) {
    int n = A.rows();
    MatD C = H * J * A * J + A * J * H * J;
    C -= (C.trace() / n) * MatD::Identity(n, n);
    return Eigen::Map<Eigen::VectorXd>(C.data(), n * n);
  }
  /*
    The weights of a point of the cone, in groups: for the indices i of a
    group, p_i = off + tot softmax(q)_i (the softmax over the group), so
    that the weights of the group are at least off and sum to
    |group| off + tot.
   */
  struct WeightGroup {
    std::vector<int> l_idx;
    double off;
    double tot;
  };
  // The weights and their derivatives dp_i / dq_j.
  static std::pair<std::vector<double>, Eigen::MatrixXd> weights(Eigen::VectorXd const &q,
                                                                 std::vector<WeightGroup> const &l_group) {
    int k = q.size();
    std::vector<double> p(k, 0.0);
    Eigen::MatrixXd Dq = Eigen::MatrixXd::Zero(k, k);
    for (auto &wg : l_group) {
      double q_max = -std::numeric_limits<double>::infinity();
      for (int i : wg.l_idx) {
        q_max = std::max(q_max, q(i));
      }
      double sum = 0;
      for (int i : wg.l_idx) {
        sum += std::exp(q(i) - q_max);
      }
      for (int i : wg.l_idx) {
        double s_i = std::exp(q(i) - q_max) / sum;
        p[i] = wg.off + wg.tot * s_i;
        for (int j : wg.l_idx) {
          double s_j = std::exp(q(j) - q_max) / sum;
          Dq(i, j) = wg.tot * s_i * ((i == j ? 1.0 : 0.0) - s_j);
        }
      }
    }
    return {p, Dq};
  }
  /*
    A point of the Siegel space in the cone of the rays with the weights of
    l_group, by Levenberg-Marquardt on the residual from n_start random
    starts, preceded by a start from the weights init when given (positive
    weights of a nearby point, for instance of a cell containing this one).
    Returns the weights.
   */
  std::optional<std::vector<double>> search_point(std::vector<MatD> const &l_ray,
                                                  std::vector<WeightGroup> const &l_group,
                                                  std::optional<std::vector<double>> const &init = {}) {
    int k = l_ray.size();
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    auto f_form = [&](std::vector<double> const &p) -> MatD {
      MatD A = MatD::Zero(n, n);
      for (int i = 0; i < k; i++) {
        A += p[i] * l_ray[i];
      }
      return A;
    };
    int n_try = n_start + (init ? 1 : 0);
    for (int i_start = 0; i_start < n_try; i_start++) {
      Eigen::VectorXd q(k);
      if (init && i_start == 0) {
        // The softmax of log(p) is p normalized in each group.
        for (int i = 0; i < k; i++) {
          q(i) = std::log(std::max((*init)[i], 1e-12));
        }
      } else {
        for (int i = 0; i < k; i++) {
          q(i) = dist(gen);
        }
      }
      double mu = 1e-3;
      std::pair<std::vector<double>, Eigen::MatrixXd> pw = weights(q, l_group);
      MatD A = f_form(pw.first);
      Eigen::VectorXd r = residual(A, J);
      double cost = r.squaredNorm();
      for (int iter = 0; iter < 400; iter++) {
        if (cost < 1e-28 || mu > 1e14) {
          break;
        }
        Eigen::MatrixXd Jp(n * n, k);
        for (int i = 0; i < k; i++) {
          Jp.col(i) = differential(A, l_ray[i], J);
        }
        Eigen::MatrixXd Jq = Jp * pw.second;
        Eigen::MatrixXd H = Jq.transpose() * Jq;
        Eigen::VectorXd g = Jq.transpose() * r;
        bool improved = false;
        while (mu <= 1e14) {
          Eigen::MatrixXd Hmu = H;
          for (int i = 0; i < k; i++) {
            Hmu(i, i) += mu * (1 + H(i, i));
          }
          Eigen::VectorXd delta = Hmu.ldlt().solve(-g);
          Eigen::VectorXd q_new = q + delta;
          std::pair<std::vector<double>, Eigen::MatrixXd> pw_new = weights(q_new, l_group);
          MatD A_new = f_form(pw_new.first);
          Eigen::VectorXd r_new = residual(A_new, J);
          double cost_new = r_new.squaredNorm();
          if (cost_new < cost) {
            q = q_new;
            pw = std::move(pw_new);
            A = A_new;
            r = r_new;
            cost = cost_new;
            mu = std::max(mu / 3.0, 1e-12);
            improved = true;
            break;
          }
          mu *= 4.0;
        }
        if (!improved) {
          break;
        }
      }
      double scale = A.squaredNorm();
      if (cost / (scale * scale) < 1e-24) {
        return pw.first;
      }
    }
    return {};
  }
  static std::vector<MatD> get_rays(MyMatrix<Tint> const &EXT) {
    int n = EXT.cols();
    std::vector<MatD> l_ray;
    for (int i = 0; i < EXT.rows(); i++) {
      Eigen::VectorXd v(n);
      for (int j = 0; j < n; j++) {
        v(j) = EXT(i, j).get_d();
      }
      l_ray.push_back(v * v.transpose());
    }
    return l_ray;
  }

public:
  void set_check_doubtful(bool const &val) { check_doubtful = val; }
  void set_margins(std::vector<double> const &val) { l_margin = val; }
  SiegelOracle(MyMatrix<Tint> const &J_in, int const &_n_start, unsigned const &seed)
      : n(J_in.rows()), J(J_in.rows(), J_in.rows()), gen(seed), n_start(_n_start) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        J(i, j) = J_in(i, j).get_d();
      }
    }
  }
  /*
    A point of the Siegel space in the open cone C(X) of a cell known to
    exist (a top cell reached across a certified ridge, for instance): the
    margins t / |X| for t = 1e-2, 1e-3, 1e-4 are tried in turn, small
    margins being accepted here. Returns the weights.
   */
  std::optional<std::vector<double>> point_of_existing_cell(MyMatrix<Tint> const &EXT) {
    int k = EXT.rows();
    std::vector<MatD> l_ray = get_rays(EXT);
    std::vector<int> l_all(k);
    for (int i = 0; i < k; i++) {
      l_all[i] = i;
    }
    for (double t_base : {1e-2, 1e-3, 1e-4}) {
      double t = t_base / k;
      std::optional<std::vector<double>> opt = search_point(l_ray, {{l_all, t, 1 - k * t}});
      if (opt) {
        return opt;
      }
    }
    return {};
  }
  /*
    Whether the relative interior of C(X) (X the rows of EXT, up to sign)
    meets the Siegel space, and the dimension of the intersection. The
    margins are tried from the largest: a point found with weights at
    least t is in the relative interior.
   */
  SiegelIntersection intersect(MyMatrix<Tint> const &EXT, std::optional<std::vector<double>> const &init = {}) {
    int k = EXT.rows();
    std::vector<MatD> l_ray = get_rays(EXT);
    MyMatrix<T> ScalMat(k, n * (n + 1) / 2);
    for (int i = 0; i < k; i++) {
      int pos = 0;
      for (int a = 0; a < n; a++) {
        for (int b = a; b < n; b++) {
          ScalMat(i, pos) = T(EXT(i, a) * EXT(i, b));
          pos++;
        }
      }
    }
    int dim_cone = RankMat(ScalMat);
    /*
      The points with a margin of 0.2 / k or 0.05 / k or 0.01 / k are kept.
      For Sp_4(Z) all the cells have points with the margin 0.2 / k, and
      the points found only with the margins 1e-3 / k or 1e-4 / k are
      artefacts near a lower face where the Siegel space is tangent (cells
      of dimension 6 on faces of dimension 9): they are reported as
      doubtful and not counted.
     */
    std::vector<double> l_t_base = l_margin;
    l_t_base.push_back(1e-3);
    l_t_base.push_back(1e-4);
    for (double t_base : l_t_base) {
      double t = t_base / k;
      std::vector<int> l_all(k);
      for (int i = 0; i < k; i++) {
        l_all[i] = i;
      }
      std::optional<std::vector<double>> opt = search_point(l_ray, {{l_all, t, 1 - k * t}}, init);
      if (!opt) {
        continue;
      }
      if (t_base < 0.01) {
        if (!check_doubtful) {
          break;
        }
        return {false, -1, t, true, {}};
      }
      std::vector<double> const &p = *opt;
      MatD A = MatD::Zero(n, n);
      for (int i = 0; i < k; i++) {
        A += p[i] * l_ray[i];
      }
      // The rank relative to the largest singular value: the columns are
      // the differentials in the directions of the rays, whose norms grow
      // with the vectors.
      Eigen::MatrixXd Jp(n * n, k);
      for (int i = 0; i < k; i++) {
        Jp.col(i) = differential(A, l_ray[i], J) / l_ray[i].norm();
      }
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jp);
      Eigen::VectorXd sv = svd.singularValues();
      int rnk = 0;
      for (int i = 0; i < sv.size(); i++) {
        if (sv(i) > 1e-8 * sv(0)) {
          rnk++;
        }
      }
#ifdef DEBUG_SIEGEL_VORONOI
      std::cerr << "SIEGEL: singular values";
      for (int i = 0; i < sv.size(); i++) {
        std::cerr << " " << sv(i);
      }
      std::cerr << "\n";
#endif
      double margin = *std::min_element(p.begin(), p.end());
#ifdef DEBUG_SIEGEL_VORONOI
      std::cerr << "SIEGEL: k=" << k << " dim_cone=" << dim_cone << " rnk=" << rnk << " t=" << t
                << " margin=" << margin << "\n";
#endif
      return {true, dim_cone - rnk - 1, margin, false, p};
    }
    return {false, -1, 0.0, false, {}};
  }
  /*
    Whether the cell of the subset Y (the rows of EXT with f_in[i] = 1) is
    in the closure of the cell of X (all the rows of EXT): points of the
    cell of X whose weights outside Y sum to eps are searched, the weights
    of Y being at least t (1 - eps), for eps = 1e-2, 1e-3, 1e-4, and the
    margins t = 0.2 / |Y|, 0.05 / |Y|, 0.01 / |Y| as in intersect. The
    cell of Y is in the closure if such points exist down to eps = 1e-4.
   */
  bool in_closure(MyMatrix<Tint> const &EXT, std::vector<int> const &f_in) {
    int k = EXT.rows();
    std::vector<MatD> l_ray = get_rays(EXT);
    std::vector<int> l_in, l_out;
    for (int i = 0; i < k; i++) {
      if (f_in[i] == 1) {
        l_in.push_back(i);
      } else {
        l_out.push_back(i);
      }
    }
    int k_in = l_in.size();
    for (double eps : {1e-2, 1e-3, 1e-4}) {
      bool found = false;
      for (double t_base : {0.2, 0.05, 0.01}) {
        double t = t_base / k_in;
        std::vector<WeightGroup> l_group{{l_in, (1 - eps) * t, (1 - eps) * (1 - k_in * t)}, {l_out, 0.0, eps}};
        if (search_point(l_ray, l_group)) {
          found = true;
          break;
        }
      }
#ifdef DEBUG_SIEGEL_VORONOI
      std::cerr << "SIEGEL: in_closure k=" << k << " k_in=" << k_in << " eps=" << eps << " found=" << found << "\n";
#endif
      if (!found) {
        return false;
      }
    }
    return true;
  }
};

/*
  The rows of EXT up to sign, without repetition (the configurations of
  the projective systems are taken up to the units).
 */
inline MyMatrix<Tint> rows_up_to_sign(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  std::set<MyVector<Tint>, std::function<bool(MyVector<Tint> const &, MyVector<Tint> const &)>> set_row(
      [](MyVector<Tint> const &a, MyVector<Tint> const &b) -> bool { return lex_less(a, b); });
  for (int i = 0; i < EXT.rows(); i++) {
    set_row.insert(unit_normalize(ring, GetMatrixRow(EXT, i)));
  }
  MyMatrix<Tint> M(set_row.size(), EXT.cols());
  int i = 0;
  for (auto &v : set_row) {
    for (int j = 0; j < EXT.cols(); j++) {
      M(i, j) = v(j);
    }
    i++;
  }
  return M;
}

// Whether the configurations EXT1 and EXT2 are equivalent under G.
inline bool are_equivalent_configurations(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT1,
                                          MyMatrix<Tint> const &EXT2) {
  if (EXT1.rows() != EXT2.rows() || RankMat(EXT1) != RankMat(EXT2)) {
    return false;
  }
  if (ring_rank(ring, EXT1) == ring.m) {
    return !list_embedding(ring, EXT1, EXT2).empty();
  }
  for (auto &l_pos : get_all_embeddings(ring, EXT1, EXT2)) {
    if (static_cast<int>(l_pos.size()) == EXT2.rows()) {
      return true;
    }
  }
  return false;
}

/*
  The faces of a polyhedral cone given by the incidences of its facets
  with its rays: all the non-empty intersections of facets, and the cone
  itself.
 */
inline std::vector<Face> all_faces_from_facets(int const &n_ray, std::vector<Face> const &l_facet) {
  Face f_all(n_ray);
  f_all.set();
  std::set<Face> set_face{f_all};
  std::vector<Face> l_face{f_all};
  size_t pos = 0;
  while (pos < l_face.size()) {
    Face f = l_face[pos];
    for (auto &f_facet : l_facet) {
      Face g = f & f_facet;
      if (g.count() > 0 && set_face.insert(g).second) {
        l_face.push_back(g);
      }
    }
    pos++;
  }
  return l_face;
}


/*
  Exact certificates.

  The Siegel space of J = [[0, D], [-D, 0]] (D diagonal, n = 2g) has the
  rational parametrization
    A(X, W) = S^{-1} [[W, W X], [X W, X W X + W^{-1}]] S^{-1},  S = diag(I, D),
  for X symmetric and W positive definite: the middle matrix is g^T g for
  g = [[W^{1/2}, W^{1/2} X], [0, W^{-1/2}]] in Sp_{2g}(R) for the standard
  form J_0, and J = S J_0 S. A(X, W) is rational for X, W rational, and
  every point of the Siegel space is, up to a positive scalar, an A(X, W).

  A top cell C(X) cap (Siegel space) is open in the Siegel space, so a
  numerical point of it, with X, W rounded to rationals, is an exact point
  of the open cone: that certifies the top cell.

  For a ridge, the facet F between the top cones sigma and sigma': with
  rational (X_0, W_0) on the side of sigma and (X_1, W_1) on the side of
  sigma', the curve (X, W)(t) = (1 - t) (X_0, W_0) + t (X_1, W_1) lies in the
  Siegel space and det(W(t)) A(t) is a polynomial in t. If the functional of
  F has exactly one root in (0, 1), a simple one, and the functionals of the
  other facets of sigma and of sigma' have no root in [0, 1], the point of
  the root is in the relative interior of F, and the functional of F has a
  nonzero differential on the Siegel space there: F meets the Siegel space
  in a manifold of dimension g (g + 1) - 1. The roots are counted exactly by
  Sturm sequences.
 */

// The functionals of the facets of the cone of the v^T v (v the rows of
// EXT_cone), as vectors on the coordinates A(a, b), a <= b, positive on
// the rows outside the facet.
inline std::vector<MyVector<T>> facet_functionals(MyMatrix<Tint> const &EXT_cone, std::vector<Face> const &l_facet) {
  int n = EXT_cone.cols();
  int dim_sym = n * (n + 1) / 2;
  auto f_coeff = [&](MyVector<Tint> const &v) -> MyVector<T> {
    MyVector<T> c(dim_sym);
    int pos = 0;
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        c(pos) = T(v(a) * v(b)) * (a == b ? 1 : 2);
        pos++;
      }
    }
    return c;
  };
  int n_row = EXT_cone.rows();
  std::vector<MyVector<T>> l_fct;
  for (auto &f : l_facet) {
    MyMatrix<T> Mrow(f.count(), dim_sym);
    int i_f = 0;
    for (int i = 0; i < n_row; i++) {
      if (f[i] == 1) {
        AssignMatrixRow(Mrow, i_f, f_coeff(GetMatrixRow(EXT_cone, i)));
        i_f++;
      }
    }
    MyMatrix<T> NSP = NullspaceTrMat(Mrow);
    if (NSP.rows() != 1) {
      std::cerr << "SIEGEL: the facet should have a nullspace of dimension 1\n";
      throw TerminalException{1};
    }
    MyVector<T> L = GetMatrixRow(NSP, 0);
    for (int i = 0; i < n_row; i++) {
      if (f[i] == 0) {
        if (L.dot(f_coeff(GetMatrixRow(EXT_cone, i))) < 0) {
          L = -L;
        }
        break;
      }
    }
    l_fct.push_back(L);
  }
  return l_fct;
}

// The value of a facet functional on the symmetric matrix A.
template <typename Tf> Tf eval_functional(MyVector<T> const &L, MyMatrix<Tf> const &A) {
  int n = A.rows();
  Tf val(0);
  int pos = 0;
  for (int a = 0; a < n; a++) {
    for (int b = a; b < n; b++) {
      val += UniversalScalarConversion<Tf, T>(L(pos)) * A(a, b) * (a == b ? 1 : 2);
      pos++;
    }
  }
  return val;
}

inline bool in_cone_interior_exact(MyMatrix<T> const &A, std::vector<MyVector<T>> const &l_fct) {
  for (auto &L : l_fct) {
    if (eval_functional<T>(L, A) <= 0) {
      return false;
    }
  }
  return true;
}

// Whether the symmetric matrix W is positive definite (leading minors).
inline bool is_positive_definite_exact(MyMatrix<T> const &W) {
  for (int k = 1; k <= W.rows(); k++) {
    MyMatrix<T> Wk = W.block(0, 0, k, k);
    if (DeterminantMat(Wk) <= 0) {
      return false;
    }
  }
  return true;
}

// The parameters (X, W) of a point of the Siegel space.
template <typename Tf> struct SiegelParam {
  MyMatrix<Tf> X;
  MyMatrix<Tf> W;
};

// The rational parametrization A(X, W) of the Siegel space of J.
class SiegelParametrization {
private:
  int g;
  std::vector<T> l_d;
  std::vector<double> l_d_dbl;

public:
  explicit SiegelParametrization(MyMatrix<Tint> const &J) : g(J.rows() / 2) {
    for (int i = 0; i < g; i++) {
      l_d.push_back(T(J(i, g + i)));
      l_d_dbl.push_back(J(i, g + i).get_d());
    }
  }
  int n_param() const { return g * (g + 1); }
  // A(X, W) times det(W) when scale_det, as a polynomial in the entries.
  MyMatrix<T> point(SiegelParam<T> const &P, bool const &scale_det) const {
    MyMatrix<T> Winv = Inverse(P.W);
    MyMatrix<T> WX = P.W * P.X;
    MyMatrix<T> R = P.X * P.W * P.X + Winv;
    MyMatrix<T> A(2 * g, 2 * g);
    for (int i = 0; i < g; i++) {
      for (int j = 0; j < g; j++) {
        A(i, j) = P.W(i, j);
        A(i, g + j) = WX(i, j) / l_d[j];
        A(g + j, i) = WX(i, j) / l_d[j];
        A(g + i, g + j) = R(i, j) / (l_d[i] * l_d[j]);
      }
    }
    if (scale_det) {
      A *= DeterminantMat(P.W);
    }
    return A;
  }
  // A(X, W) in double precision.
  MyMatrix<double> point_dbl(SiegelParam<double> const &P) const {
    MyMatrix<double> Winv = P.W.inverse();
    MyMatrix<double> WX = P.W * P.X;
    MyMatrix<double> R = P.X * P.W * P.X + Winv;
    MyMatrix<double> A(2 * g, 2 * g);
    for (int i = 0; i < g; i++) {
      for (int j = 0; j < g; j++) {
        A(i, j) = P.W(i, j);
        A(i, g + j) = WX(i, j) / l_d_dbl[j];
        A(g + j, i) = WX(i, j) / l_d_dbl[j];
        A(g + i, g + j) = R(i, j) / (l_d_dbl[i] * l_d_dbl[j]);
      }
    }
    return A;
  }
  // The parameters as a vector (X then W, upper triangles) and back.
  Eigen::VectorXd to_vector(SiegelParam<double> const &Q) const {
    Eigen::VectorXd v(g * (g + 1));
    int pos = 0;
    for (auto *M : {&Q.X, &Q.W}) {
      for (int i = 0; i < g; i++) {
        for (int j = i; j < g; j++) {
          v(pos++) = (*M)(i, j);
        }
      }
    }
    return v;
  }
  SiegelParam<double> from_vector(Eigen::VectorXd const &v) const {
    SiegelParam<double> Q{MyMatrix<double>(g, g), MyMatrix<double>(g, g)};
    int pos = 0;
    for (auto *M : {&Q.X, &Q.W}) {
      for (int i = 0; i < g; i++) {
        for (int j = i; j < g; j++) {
          (*M)(i, j) = v(pos);
          (*M)(j, i) = v(pos);
          pos++;
        }
      }
    }
    return Q;
  }
  // The parameters of a numerical point A of the Siegel space (any scale).
  SiegelParam<double> param_of(MyMatrix<double> const &A) const {
    int n = 2 * g;
    MyMatrix<double> Ap(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double si = (i < g) ? 1.0 : l_d_dbl[i - g];
        double sj = (j < g) ? 1.0 : l_d_dbl[j - g];
        Ap(i, j) = si * A(i, j) * sj;
      }
    }
    // The scale: det of the middle matrix is 1.
    double det = Ap.determinant();
    Ap /= std::pow(det, 1.0 / n);
    MyMatrix<double> W = Ap.block(0, 0, g, g);
    MyMatrix<double> X = W.inverse() * Ap.block(0, g, g, g);
    X = (X + X.transpose()) / 2;
    return {X, W};
  }
  // The rational parameters closest to P with the denominator den.
  static SiegelParam<T> rationalize(SiegelParam<double> const &P, long const &den) {
    auto f_round = [&](MyMatrix<double> const &M) -> MyMatrix<T> {
      int k = M.rows();
      MyMatrix<T> R(k, k);
      for (int i = 0; i < k; i++) {
        for (int j = i; j < k; j++) {
          double x = (M(i, j) + M(j, i)) / 2;
          T val = T(static_cast<long>(std::llround(x * den))) / T(den);
          R(i, j) = val;
          R(j, i) = val;
        }
      }
      return R;
    };
    return {f_round(P.X), f_round(P.W)};
  }
};

/*
  Univariate polynomials with rational coefficients (from the constant
  term), for the Sturm counts.
 */
using Poly = std::vector<T>;

inline void poly_trim(Poly &p) {
  while (!p.empty() && p.back() == 0) {
    p.pop_back();
  }
}

inline T poly_eval(Poly const &p, T const &x) {
  T val(0);
  for (int i = static_cast<int>(p.size()) - 1; i >= 0; i--) {
    val = val * x + p[i];
  }
  return val;
}

inline Poly poly_derivative(Poly const &p) {
  Poly q;
  for (size_t i = 1; i < p.size(); i++) {
    q.push_back(p[i] * T(static_cast<long>(i)));
  }
  poly_trim(q);
  return q;
}

// The remainder of the division of a by b (b nonzero).
inline Poly poly_rem(Poly a, Poly const &b) {
  poly_trim(a);
  int db = b.size() - 1;
  while (static_cast<int>(a.size()) - 1 >= db && !a.empty()) {
    T c = a.back() / b.back();
    int shift = a.size() - 1 - db;
    for (int i = 0; i <= db; i++) {
      a[shift + i] -= c * b[i];
    }
    a.pop_back();
    poly_trim(a);
  }
  return a;
}

inline Poly poly_gcd(Poly a, Poly b) {
  poly_trim(a);
  poly_trim(b);
  while (!b.empty()) {
    Poly r = poly_rem(a, b);
    a = b;
    b = r;
  }
  return a;
}

// The number of distinct roots of p (nonzero, nonzero at x0 and x1) in
// (x0, x1), by the Sturm sequence.
inline int sturm_count(Poly const &p_in, T const &x0, T const &x1) {
  Poly p = p_in;
  poly_trim(p);
  std::vector<Poly> seq{p, poly_derivative(p)};
  while (!seq.back().empty()) {
    Poly r = poly_rem(seq[seq.size() - 2], seq.back());
    for (auto &c : r) {
      c = -c;
    }
    if (r.empty()) {
      break;
    }
    seq.push_back(r);
  }
  auto f_var = [&](T const &x) -> int {
    int n_var = 0;
    int sign_prev = 0;
    for (auto &q : seq) {
      if (q.empty()) {
        continue;
      }
      T val = poly_eval(q, x);
      int sign = (val > 0) ? 1 : ((val < 0) ? -1 : 0);
      if (sign != 0) {
        if (sign_prev != 0 && sign != sign_prev) {
          n_var++;
        }
        sign_prev = sign;
      }
    }
    return n_var;
  };
  return f_var(x0) - f_var(x1);
}

// The polynomial of the values f(t_i) at m points (Newton interpolation).
inline Poly poly_interpolate(std::vector<T> const &l_t, std::vector<T> const &l_val) {
  int m = l_t.size();
  std::vector<T> coef = l_val;
  for (int j = 1; j < m; j++) {
    for (int i = m - 1; i >= j; i--) {
      coef[i] = (coef[i] - coef[i - 1]) / (l_t[i] - l_t[i - j]);
    }
  }
  Poly p{coef[m - 1]};
  for (int i = m - 2; i >= 0; i--) {
    // p = p (t - t_i) + coef[i]
    Poly q(p.size() + 1, T(0));
    for (size_t k = 0; k < p.size(); k++) {
      q[k + 1] += p[k];
      q[k] -= p[k] * l_t[i];
    }
    q[0] += coef[i];
    p = q;
  }
  poly_trim(p);
  return p;
}

/*
  The certificate of a top cell: a rational point of the Siegel space in
  the open cone of facet functionals l_fct, from a numerical point A_num.
 */
inline std::optional<MyMatrix<T>> certify_top_cell(SiegelParametrization const &sp, MyMatrix<double> const &A_num,
                                                   std::vector<MyVector<T>> const &l_fct) {
  SiegelParam<double> P = sp.param_of(A_num);
  for (long den : {100L, 10000L, 1000000L, 100000000L}) {
    SiegelParam<T> Pr = SiegelParametrization::rationalize(P, den);
    if (!is_positive_definite_exact(Pr.W)) {
      continue;
    }
    MyMatrix<T> A = sp.point(Pr, false);
    if (in_cone_interior_exact(A, l_fct)) {
      return A;
    }
  }
  return {};
}

/*
  The certificate of a ridge: the facet i_F of the cone of functionals
  l_fct0 (sigma), and the cone of functionals l_fct1 (sigma') whose facet
  i_F1 is the same facet, A_num a numerical point of the Siegel space on
  the facet. Returns the end A(1) of the curve, an exact point of the
  Siegel space in the open cone sigma', which certifies that top cell.
 */
inline std::optional<MyMatrix<T>> certify_ridge(SiegelParametrization const &sp, MyMatrix<double> const &A_num,
                          std::vector<MyVector<T>> const &l_fct0, int const &i_F,
                          std::vector<MyVector<T>> const &l_fct1, int const &i_F1) {
  SiegelParam<double> P = sp.param_of(A_num);
  int g = P.W.rows();
  // The parameters as a vector (X then W, upper triangles) and back.
  auto f_vec = [&](SiegelParam<double> const &Q) -> Eigen::VectorXd {
    Eigen::VectorXd v(g * (g + 1));
    int pos = 0;
    for (auto *M : {&Q.X, &Q.W}) {
      for (int i = 0; i < g; i++) {
        for (int j = i; j < g; j++) {
          v(pos++) = (*M)(i, j);
        }
      }
    }
    return v;
  };
  auto f_param = [&](Eigen::VectorXd const &v) -> SiegelParam<double> {
    SiegelParam<double> Q{MyMatrix<double>(g, g), MyMatrix<double>(g, g)};
    int pos = 0;
    for (auto *M : {&Q.X, &Q.W}) {
      for (int i = 0; i < g; i++) {
        for (int j = i; j < g; j++) {
          (*M)(i, j) = v(pos);
          (*M)(j, i) = v(pos);
          pos++;
        }
      }
    }
    return Q;
  };
  // The direction across F: the gradient of the functional of F on
  // A(X, W), by central differences (the parameters rounded finely).
  Eigen::VectorXd v0 = f_vec(P);
  int n_par = v0.size();
  auto f_phi_exactpar = [&](Eigen::VectorXd const &v) -> double {
    SiegelParam<double> Q = f_param(v);
    SiegelParam<T> Qr = SiegelParametrization::rationalize(Q, 1000000000L);
    MyMatrix<T> A = sp.point(Qr, false);
    MyMatrix<double> A_dbl = UniversalMatrixConversion<double, T>(A);
    return eval_functional<double>(l_fct0[i_F], A_dbl);
  };
  Eigen::VectorXd grad(n_par);
  double h = 1e-5;
  for (int k = 0; k < n_par; k++) {
    Eigen::VectorXd vp = v0, vm = v0;
    vp(k) += h;
    vm(k) -= h;
    grad(k) = (f_phi_exactpar(vp) - f_phi_exactpar(vm)) / (2 * h);
  }
  if (grad.norm() == 0) {
    return {};
  }
  Eigen::VectorXd dir = grad / grad.norm();
  int n_pt = 2 * g + 8;
  std::vector<T> l_t;
  for (int i = 0; i < n_pt; i++) {
    l_t.push_back(T(i) / T(n_pt - 1));
  }
  for (double eps : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7}) {
    long den = static_cast<long>(1000.0 / eps);
    // sigma is the side where the functional of F is positive.
    SiegelParam<T> P0 = SiegelParametrization::rationalize(f_param(v0 + eps * dir), den);
    SiegelParam<T> P1 = SiegelParametrization::rationalize(f_param(v0 - eps * dir), den);
    if (!is_positive_definite_exact(P0.W) || !is_positive_definite_exact(P1.W)) {
      continue;
    }
    // det(W(t)) A(t) at the interpolation points.
    std::vector<MyMatrix<T>> l_A;
    for (auto &t : l_t) {
      SiegelParam<T> Pt{(1 - t) * P0.X + t * P1.X, (1 - t) * P0.W + t * P1.W};
      l_A.push_back(sp.point(Pt, true));
    }
    auto f_poly = [&](MyVector<T> const &L) -> Poly {
      std::vector<T> l_val;
      for (auto &A : l_A) {
        l_val.push_back(eval_functional<T>(L, A));
      }
      return poly_interpolate(l_t, l_val);
    };
    /*
      A bound on the variation of the curve: the entries of
      B(t) = det(W(t)) A(t) are polynomials, and |B_ab(t) - B_ab(0)| <= E_ab,
      the sum of the absolute values of the coefficients of degree >= 1, on
      [0, 1]. A functional L with L(B(0)) > sum |L_ab| w_ab E_ab (w_ab the
      weights of eval_functional) is then positive on [0, 1], without the
      Sturm sequence.
     */
    int n = l_A[0].rows();
    MyMatrix<T> E(n, n);
    for (int a = 0; a < n; a++) {
      for (int b = a; b < n; b++) {
        std::vector<T> l_val;
        for (auto &A : l_A) {
          l_val.push_back(A(a, b));
        }
        Poly p = poly_interpolate(l_t, l_val);
        T bound(0);
        for (size_t k = 1; k < p.size(); k++) {
          bound += T_abs(p[k]);
        }
        E(a, b) = bound;
        E(b, a) = bound;
      }
    }
    MyMatrix<T> const &B_start = l_A[0];
    auto f_far = [&](MyVector<T> const &L) -> bool {
      T val0 = eval_functional<T>(L, B_start);
      if (val0 <= 0) {
        return false;
      }
      T bound(0);
      int pos = 0;
      for (int a = 0; a < n; a++) {
        for (int b = a; b < n; b++) {
          bound += T_abs(L(pos)) * E(a, b) * (a == b ? 1 : 2);
          pos++;
        }
      }
      return val0 > bound;
    };
    // No root in [0, 1] for the other facets of sigma and sigma'.
    auto f_positive = [&](std::vector<MyVector<T>> const &l_fct, int const &i_skip) -> bool {
      for (size_t i = 0; i < l_fct.size(); i++) {
        if (static_cast<int>(i) == i_skip) {
          continue;
        }
        if (f_far(l_fct[i])) {
          continue;
        }
        Poly p = f_poly(l_fct[i]);
        if (p.empty() || poly_eval(p, T(0)) <= 0 || poly_eval(p, T(1)) <= 0) {
          return false;
        }
        if (sturm_count(p, T(0), T(1)) != 0) {
          return false;
        }
      }
      return true;
    };
    Poly pF = f_poly(l_fct0[i_F]);
    if (pF.empty() || poly_eval(pF, T(0)) <= 0 || poly_eval(pF, T(1)) >= 0) {
      continue;
    }
    if (sturm_count(pF, T(0), T(1)) != 1) {
      continue;
    }
    // The root is simple: gcd(pF, pF') has no root in (0, 1).
    Poly q = poly_gcd(pF, poly_derivative(pF));
    if (q.size() > 1 && sturm_count(q, T(0), T(1)) != 0) {
      continue;
    }
    if (f_positive(l_fct0, i_F) && f_positive(l_fct1, i_F1)) {
      // The end A(1) is in the open cone sigma': its facets other than F
      // are positive on [0, 1], and the functional of F is negative at 1.
      return sp.point(P1, false);
    }
  }
  return {};
}


/*
  The ridges of a top cell found by shooting rays: from the parameters v0
  of an interior point, along random directions of the parameter space,
  the first facet whose functional vanishes (W staying positive definite).
  Returns, for each facet index met, a numerical point of the Siegel space
  on it (the first one found). A ray leaving towards the boundary of the
  parameter space (W degenerating, a cusp) gives nothing.
 */
inline std::map<int, MyMatrix<double>> shoot_rays(SiegelParametrization const &sp, Eigen::VectorXd const &v0,
                                                  std::vector<MyVector<T>> const &l_fct, int const &n_ray,
                                                  std::mt19937 &gen) {
  int n_fct = l_fct.size();
  int dim_sym = l_fct[0].size();
  Eigen::MatrixXd Lmat(n_fct, dim_sym);
  for (int k = 0; k < n_fct; k++) {
    for (int j = 0; j < dim_sym; j++) {
      Lmat(k, j) = l_fct[k](j).get_d();
    }
  }
  auto f_values = [&](Eigen::VectorXd const &v, bool &valid) -> Eigen::VectorXd {
    SiegelParam<double> Q = sp.from_vector(v);
    Eigen::LLT<Eigen::MatrixXd> llt(Q.W);
    valid = (llt.info() == Eigen::Success);
    Eigen::VectorXd vals(n_fct);
    if (!valid) {
      return vals;
    }
    MyMatrix<double> A = sp.point_dbl(Q);
    A /= A.norm();
    int n = A.rows();
    Eigen::VectorXd a(dim_sym);
    int pos = 0;
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        a(pos++) = A(i, j) * (i == j ? 1 : 2);
      }
    }
    vals = Lmat * a;
    return vals;
  };
  std::normal_distribution<double> dist(0.0, 1.0);
  std::map<int, MyMatrix<double>> map_ridge;
  int n_par = v0.size();
  for (int i_ray = 0; i_ray < n_ray; i_ray++) {
    Eigen::VectorXd dir(n_par);
    for (int k = 0; k < n_par; k++) {
      dir(k) = dist(gen);
    }
    dir /= dir.norm();
    double t = 0;
    double h = 1e-3;
    bool valid = true;
    Eigen::VectorXd vals = f_values(v0, valid);
    if (!valid || vals.minCoeff() <= 0) {
      continue;
    }
    while (t < 1e4) {
      double t_new = t + h;
      Eigen::VectorXd vals_new = f_values(v0 + t_new * dir, valid);
      if (!valid) {
        break;
      }
      if (vals_new.minCoeff() > 0) {
        t = t_new;
        h *= 1.5;
        continue;
      }
      // Bisection on the smallest value.
      double t_lo = t, t_hi = t_new;
      for (int it = 0; it < 60; it++) {
        double t_mid = (t_lo + t_hi) / 2;
        Eigen::VectorXd vals_mid = f_values(v0 + t_mid * dir, valid);
        if (valid && vals_mid.minCoeff() > 0) {
          t_lo = t_mid;
        } else {
          t_hi = t_mid;
        }
      }
      Eigen::VectorXd vals_hi = f_values(v0 + t_hi * dir, valid);
      if (!valid) {
        break;
      }
      // The facet crossed, if only one functional is near zero.
      int i_min = -1;
      int n_small = 0;
      for (int k = 0; k < n_fct; k++) {
        if (vals_hi(k) < 1e-9) {
          n_small++;
          i_min = k;
        }
      }
      if (n_small == 1 && map_ridge.count(i_min) == 0) {
        map_ridge[i_min] = sp.point_dbl(sp.from_vector(v0 + t_hi * dir));
      }
      break;
    }
  }
  return map_ridge;
}

}  // namespace psys

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_SIEGEL_VORONOI_H_
// clang-format on
