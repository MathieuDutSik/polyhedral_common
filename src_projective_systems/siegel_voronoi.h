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
  The cells of the Siegel space of a symplectic form J on Z^4 cut out by
  the Voronoi decomposition of the positive definite quaternary forms, as
  MacPherson and McConnell did for Sp_4(Z) (J standard), for any J such as
  the form J_N of the paramodular group K(N) = Sp(J_N, Z).

  The matrices act on the row vectors: a configuration X of vectors of Z^4
  gives the cone C(X) of the forms sum_{v in X} p_v v^T v (p_v > 0), and
  g in Sp(J, Z) maps it to C(X g) by A -> g^T A g. The Siegel space of J
  is the cone of the positive forms A with (A J)^2 = -lambda I, of
  dimension 7 (the symmetric space of dimension 6 times R_{>0}). Each cell
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
};

class SiegelOracle {
private:
  using Mat4 = Eigen::Matrix<double, 4, 4>;
  Mat4 J;
  std::mt19937 gen;
  int n_start;
  // The residual (A J)^2 - tr((A J)^2) / 4 I, as a vector of 16 entries.
  static Eigen::VectorXd residual(Mat4 const &A, Mat4 const &J) {
    Mat4 AJ = A * J;
    Mat4 C = AJ * AJ;
    C -= (C.trace() / 4.0) * Mat4::Identity();
    return Eigen::Map<Eigen::VectorXd>(C.data(), 16);
  }
  // The differential of the residual at A in the direction H.
  static Eigen::VectorXd differential(Mat4 const &A, Mat4 const &H, Mat4 const &J) {
    Mat4 C = H * J * A * J + A * J * H * J;
    C -= (C.trace() / 4.0) * Mat4::Identity();
    return Eigen::Map<Eigen::VectorXd>(C.data(), 16);
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
    starts. Returns the weights.
   */
  std::optional<std::vector<double>> search_point(std::vector<Mat4> const &l_ray,
                                                  std::vector<WeightGroup> const &l_group) {
    int k = l_ray.size();
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    auto f_form = [&](std::vector<double> const &p) -> Mat4 {
      Mat4 A = Mat4::Zero();
      for (int i = 0; i < k; i++) {
        A += p[i] * l_ray[i];
      }
      return A;
    };
    for (int i_start = 0; i_start < n_start; i_start++) {
      Eigen::VectorXd q(k);
      for (int i = 0; i < k; i++) {
        q(i) = dist(gen);
      }
      double mu = 1e-3;
      std::pair<std::vector<double>, Eigen::MatrixXd> pw = weights(q, l_group);
      Mat4 A = f_form(pw.first);
      Eigen::VectorXd r = residual(A, J);
      double cost = r.squaredNorm();
      for (int iter = 0; iter < 400; iter++) {
        if (cost < 1e-28 || mu > 1e14) {
          break;
        }
        Eigen::MatrixXd Jp(16, k);
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
          Mat4 A_new = f_form(pw_new.first);
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
  static std::vector<Mat4> get_rays(MyMatrix<Tint> const &EXT) {
    std::vector<Mat4> l_ray;
    for (int i = 0; i < EXT.rows(); i++) {
      Eigen::Matrix<double, 4, 1> v;
      for (int j = 0; j < 4; j++) {
        v(j) = EXT(i, j).get_d();
      }
      l_ray.push_back(v * v.transpose());
    }
    return l_ray;
  }

public:
  SiegelOracle(MyMatrix<Tint> const &J_in, int const &_n_start, unsigned const &seed)
      : gen(seed), n_start(_n_start) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        J(i, j) = J_in(i, j).get_d();
      }
    }
  }
  /*
    Whether the relative interior of C(X) (X the rows of EXT, up to sign)
    meets the Siegel space, and the dimension of the intersection. The
    margins are tried from the largest: a point found with weights at
    least t is in the relative interior.
   */
  SiegelIntersection intersect(MyMatrix<Tint> const &EXT) {
    int k = EXT.rows();
    std::vector<Mat4> l_ray = get_rays(EXT);
    MyMatrix<T> ScalMat(k, 10);
    for (int i = 0; i < k; i++) {
      int pos = 0;
      for (int a = 0; a < 4; a++) {
        for (int b = a; b < 4; b++) {
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
    for (double t_base : {0.2, 0.05, 0.01, 1e-3, 1e-4}) {
      double t = t_base / k;
      std::vector<int> l_all(k);
      for (int i = 0; i < k; i++) {
        l_all[i] = i;
      }
      std::optional<std::vector<double>> opt = search_point(l_ray, {{l_all, t, 1 - k * t}});
      if (!opt) {
        continue;
      }
      if (t_base < 0.01) {
        return {false, -1, t, true};
      }
      std::vector<double> const &p = *opt;
      Mat4 A = Mat4::Zero();
      for (int i = 0; i < k; i++) {
        A += p[i] * l_ray[i];
      }
      Eigen::MatrixXd Jp(16, k);
      double scale = 0;
      for (int i = 0; i < k; i++) {
        Jp.col(i) = differential(A, l_ray[i], J);
        scale = std::max(scale, A.norm() * l_ray[i].norm() * J.squaredNorm());
      }
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jp);
      Eigen::VectorXd sv = svd.singularValues();
      int rnk = 0;
      for (int i = 0; i < sv.size(); i++) {
        if (sv(i) > 1e-7 * scale) {
          rnk++;
        }
      }
      double margin = *std::min_element(p.begin(), p.end());
#ifdef DEBUG_SIEGEL_VORONOI
      std::cerr << "SIEGEL: k=" << k << " dim_cone=" << dim_cone << " rnk=" << rnk << " t=" << t
                << " margin=" << margin << "\n";
#endif
      return {true, dim_cone - rnk - 1, margin, false};
    }
    return {false, -1, 0.0, false};
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
    std::vector<Mat4> l_ray = get_rays(EXT);
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

}  // namespace psys

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_SIEGEL_VORONOI_H_
// clang-format on
