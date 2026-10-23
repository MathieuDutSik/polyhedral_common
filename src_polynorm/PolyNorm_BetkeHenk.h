// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_BETKEHENK_H_
#define SRC_POLYNORM_POLYNORM_BETKEHENK_H_

// clang-format off
#include "POLY_DirectDualDesc.h"
#include "POLY_DirectTriangulation.h"
#include "POLY_LinearProgramming.h"
#include "PolytopeEquiStab.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <thread>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_BETKE_HENK
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_BETKE_HENK
#endif

#ifdef TIMINGS
#define TIMINGS_POLYNORM_BETKE_HENK
#endif

/*
  The densest lattice packing of a 3-polytope by the algorithm of Betke and
  Henk, Densest lattice packings of 3-polytopes, Comput. Geom. 16 (2000),
  arXiv:math/9909172.

  The lattice packings of P are the admissible lattices of D = P - P
  (no nonzero lattice point in the interior of D), the density is
  vol(P) / det L, and the densest one has det L = Delta(D), the critical
  determinant. Lattice vectors are written u_W = u W with u in Z^3 a row
  and the basis vectors the rows of W.

  Minkowski's theorem (Theorem 2.1 of the paper): some critical lattice has
  a basis W such that a test set U^j_W lies on the boundary of D,
     U^1 = (1,0,0),(0,1,0),(0,0,1),(0,1,-1),(-1,0,1),(1,-1,0)
     U^2 = (1,0,0),(0,1,0),(0,0,1),(0,1,1),(1,0,1),(1,1,0)
     U^3 = U^2 + (1,1,1),
  and, for any choice of supporting hyperplanes H_i of D at the points,
  W is a local minimum of |det| on S_H = {W : u^i_W in H_i}, in one of
  the four cases
     I.   j = 1;  II. j = 2;  III. j = 3;
     IV.  j = 3, with H_6 replaced by the supporting hyperplane of normal
          lambda_1 v_1 + lambda_2 v_2, where lambda_1 v_1 + lambda_2 v_2 +
          lambda_3 v_3 = lambda_6 v_6, lambda_1, lambda_2, lambda_6 > 0
          (v_i the normal of H_i).
  For a polytope the H_i are taken among the facet hyperplanes, so the
  algorithm enumerates the choices of a facet for each test vector, keeps
  those for which the points can lie in their facets (a linear program;
  the additive relations among the test vectors, u_6 = u_1 + u_2 and so
  on, prune most choices beforehand: (F_a + F_b) cap F_c must be nonempty),
  and computes the local minima of det on the affine space S_H.

  On S_H = {C + sum_i lambda_i M_i}, of dimension r = 9 - rank (2 or 3),
  det is a cubic polynomial p(lambda) and its critical points solve r
  quadratic equations. Betke and Henk find them with resultants; here they
  are found numerically by a total degree homotopy (at most 2^3 = 8
  paths), which finds every isolated solution. Every real critical point
  whose test points lie on the boundary of D is a candidate, and its
  admissibility is checked by enumerating the short lattice vectors. The
  smallest admissible |det| is Delta(D): every candidate accepted is an
  admissible lattice, and the critical lattice of the theorem is among
  the candidates.

  The exact field T (rationals, Q(sqrt 2), Q(sqrt 5), a real algebraic
  field) is used for the polytope, the pruning linear programs and the
  affine spaces; the homotopy and the admissibility checks are in double
  precision. Positive dimensional components of critical points (on which
  det is constant, Lemma 4.2 of the paper) are reported when they could
  beat the best value; they occur for tilings.
 */

// The difference body D = P - P with what the enumeration needs.
template <typename T> struct BetkeHenkBody {
  // D = {x : Lmat x <= 1}.
  int m;
  MyMatrix<T> Lmat;
  MyMatrix<double> Lmat_d;
  // The opposite facet.
  std::vector<int> opp;
  // The vertices of D (possibly with some extra points of D).
  std::vector<MyVector<T>> ListVert;
  // The edge-neighbors of each facet: a point of the hyperplane of F lies
  // in F if and only if it satisfies the inequalities of the neighbors.
  std::vector<std::vector<int>> neighbors;
  // The vertices of D on each facet (indices into ListVert).
  std::vector<std::vector<int>> facet_vertices;
  // The bounding boxes of the facets: lower and upper corner.
  std::vector<std::array<double, 6>> box;
  // The largest Euclidean norm of a vertex of D.
  double r_max;
  T vol_P;
  // The orbits of the facets under the linear symmetry group of D (which
  // contains -1): the orbit index of each facet, and whether it is the
  // representative (smallest index) of its orbit.
  std::vector<int> orbit;
  std::vector<uint8_t> is_rep;
  size_t n_orbit;
};

// The linear symmetries of D act on the lattices by W -> W g^T and preserve
// both |det| and admissibility, so the facet of the first test vector can
// be taken up to symmetry.
template <typename T, typename Tgroup>
void BetkeHenk_ComputeOrbits(BetkeHenkBody<T> &body, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int m = body.m;
  std::vector<Telt> LGen =
      LinPolytope_Automorphism_LGen<T, Tgroup>(body.Lmat, os);
  body.orbit.assign(m, -1);
  body.is_rep.assign(m, 0);
  int n_orbit = 0;
  for (int j = 0; j < m; j++) {
    if (body.orbit[j] != -1) {
      continue;
    }
    body.is_rep[j] = 1;
    std::vector<int> stack{j};
    body.orbit[j] = n_orbit;
    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      for (auto &g : LGen) {
        int v = g.at(u);
        if (body.orbit[v] == -1) {
          body.orbit[v] = n_orbit;
          stack.push_back(v);
        }
      }
    }
    n_orbit++;
  }
  body.n_orbit = n_orbit;
}

template <typename T>
BetkeHenkBody<T> BetkeHenk_BuildBody(MyMatrix<T> const &EXT_P,
                                     std::ostream &os) {
  int n = EXT_P.cols() - 1;
  if (n != 3) {
    std::cerr << "BETKE_HENK: the algorithm is for 3-dimensional polytopes\n";
    throw TerminalException{1};
  }
  BetkeHenkBody<T> body;
  body.vol_P = DirectVolumePolytope(EXT_P, TriangulationProgram::lrs, os);
  // The differences of the vertices.
  int n_vert = EXT_P.rows();
  std::set<std::vector<T>> set_diff;
  for (int i = 0; i < n_vert; i++) {
    for (int j = 0; j < n_vert; j++) {
      if (i == j) {
        continue;
      }
      std::vector<T> v(3);
      for (int k = 0; k < 3; k++) {
        v[k] = EXT_P(i, k + 1) - EXT_P(j, k + 1);
      }
      set_diff.insert(v);
    }
  }
  int n_diff = set_diff.size();
  MyMatrix<T> EXTdiff(n_diff, 4);
  int pos = 0;
  for (auto &v : set_diff) {
    EXTdiff(pos, 0) = 1;
    for (int k = 0; k < 3; k++) {
      EXTdiff(pos, k + 1) = v[k];
    }
    pos++;
  }
  // The symmetric polytope D has 0 as interior point, so the facets are
  // ell(x) <= 1.
  // The differences are mostly not vertices of D.
  vectface vf = DirectFacetComputationIncidence_General(EXTdiff, "cdd", os);
  int m = vf.size();
  body.m = m;
  body.Lmat = MyMatrix<T>(m, 3);
  for (int j = 0; j < m; j++) {
    MyVector<T> w = FindFacetInequality(EXTdiff, vf[j]);
    for (int k = 0; k < 3; k++) {
      body.Lmat(j, k) = -w(k + 1) / w(0);
    }
  }
  body.Lmat_d = UniversalMatrixConversion<double, T>(body.Lmat);
  body.opp.resize(m, -1);
  for (int j = 0; j < m; j++) {
    for (int k = 0; k < m; k++) {
      bool is_opp = true;
      for (int u = 0; u < 3; u++) {
        if (body.Lmat(j, u) != -body.Lmat(k, u)) {
          is_opp = false;
        }
      }
      if (is_opp) {
        body.opp[j] = k;
      }
    }
    if (body.opp[j] == -1) {
      std::cerr << "BETKE_HENK: P - P should be centrally symmetric\n";
      throw TerminalException{1};
    }
  }
  // The vertices: the points of EXTdiff on at least 3 facets with rank 3.
  body.r_max = 0;
  for (int i = 0; i < n_diff; i++) {
    MyVector<T> v(3);
    for (int k = 0; k < 3; k++) {
      v(k) = EXTdiff(i, k + 1);
    }
    body.ListVert.push_back(v);
    double nrm = 0;
    for (int k = 0; k < 3; k++) {
      double x = UniversalScalarConversion<double, T>(v(k));
      nrm += x * x;
    }
    body.r_max = std::max(body.r_max, std::sqrt(nrm));
  }
  // The bounding boxes of the facets from the points of EXTdiff on them.
  body.box.resize(m);
  for (int j = 0; j < m; j++) {
    std::array<double, 6> b;
    for (int k = 0; k < 3; k++) {
      b[k] = std::numeric_limits<double>::max();
      b[k + 3] = -std::numeric_limits<double>::max();
    }
    for (auto &v : body.ListVert) {
      T scal(0);
      for (int k = 0; k < 3; k++) {
        scal += body.Lmat(j, k) * v(k);
      }
      if (scal == 1) {
        for (int k = 0; k < 3; k++) {
          double x = UniversalScalarConversion<double, T>(v(k));
          b[k] = std::min(b[k], x);
          b[k + 3] = std::max(b[k + 3], x);
        }
      }
    }
    body.box[j] = b;
  }
  // Two facets are neighbors when they share at least two distinct points
  // of EXTdiff, which then span their common edge.
  std::vector<std::vector<int>> facet_points(m);
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n_diff; i++) {
      T scal(0);
      for (int k = 0; k < 3; k++) {
        scal += body.Lmat(j, k) * body.ListVert[i](k);
      }
      if (scal == 1) {
        facet_points[j].push_back(i);
      }
    }
  }
  // The vertices: the points on facets whose normals have rank 3.
  body.facet_vertices.resize(m);
  for (int i = 0; i < n_diff; i++) {
    std::vector<int> inc;
    for (int j = 0; j < m; j++) {
      if (std::find(facet_points[j].begin(), facet_points[j].end(), i) !=
          facet_points[j].end()) {
        inc.push_back(j);
      }
    }
    if (inc.size() >= 3 && RankMat(SelectRow(body.Lmat, inc)) == 3) {
      for (auto &j : inc) {
        body.facet_vertices[j].push_back(i);
      }
    }
  }
  body.neighbors.resize(m);
  for (int j = 0; j < m; j++) {
    for (int k = 0; k < m; k++) {
      if (k == j) {
        continue;
      }
      int n_common = 0;
      for (auto &i : facet_points[j]) {
        if (std::find(facet_points[k].begin(), facet_points[k].end(), i) !=
            facet_points[k].end()) {
          n_common++;
        }
      }
      if (n_common >= 2) {
        body.neighbors[j].push_back(k);
      }
    }
  }
  return body;
}

/*
  The pruning: (s_a F_a + s_b F_b) cap F_c nonempty, with s F_a = F_a or
  F_opp(a). In terms of facets: (F_a' + F_b') cap F_c, with a' = a or
  opp(a). The test is a box filter followed by an exact linear program in
  (x, y) in R^6: x in F_a', y in F_b', x + y in F_c.
 */
template <typename T> struct BetkeHenkTriples {
  BetkeHenkBody<T> const &body;
  std::ostream &os;
  std::unordered_map<uint64_t, bool> memo;
  size_t n_lp = 0;
  size_t n_box_reject = 0;

  BetkeHenkTriples(BetkeHenkBody<T> const &_body, std::ostream &_os)
      : body(_body), os(_os) {}

  bool BoxTest(int a, int b, int c) const {
    double tol = 1.0e-9;
    for (int k = 0; k < 3; k++) {
      double lo = body.box[a][k] + body.box[b][k];
      double hi = body.box[a][k + 3] + body.box[b][k + 3];
      if (hi < body.box[c][k] - tol || body.box[c][k + 3] < lo - tol) {
        return false;
      }
    }
    return true;
  }

  bool Feasible(int a, int b, int c) {
    uint64_t key = (uint64_t(a) * 4096 + uint64_t(b)) * 4096 + uint64_t(c);
    auto iter = memo.find(key);
    if (iter != memo.end()) {
      return iter->second;
    }
    bool result;
    if (!BoxTest(a, b, c)) {
      n_box_reject++;
      result = false;
    } else {
      n_lp++;
      result = LinearProgram(a, b, c);
    }
    memo[key] = result;
    return result;
  }

  bool LinearProgram(int a, int b, int c) const {
    // x in F_a, y in F_b, x + y in F_c, each by its hyperplane and the
    // inequalities of the neighbors of the facet.
    std::vector<std::array<int, 3>> rows; // (facet, which, is_eq)
    // which: 0 for x, 1 for y, 2 for x + y.
    auto add_facet = [&](int f, int which) {
      rows.push_back({f, which, 1});
      for (auto &g : body.neighbors[f]) {
        rows.push_back({g, which, 0});
      }
    };
    add_facet(a, 0);
    add_facet(b, 1);
    add_facet(c, 2);
    MyMatrix<T> ListIneq = ZeroMatrix<T>(rows.size(), 7);
    for (size_t pos = 0; pos < rows.size(); pos++) {
      int f = rows[pos][0];
      int which = rows[pos][1];
      bool is_eq = rows[pos][2] == 1;
      // Inequality 1 - d_f z >= 0, or for the facet itself d_f z - 1 >= 0
      // (with the inequality d_f z <= 1 implied by the neighbors' ones
      // only together with it, it is imposed as an equality through the
      // two rows below).
      T sign = is_eq ? T(1) : T(-1);
      ListIneq(pos, 0) = is_eq ? T(-1) : T(1);
      for (int u = 0; u < 3; u++) {
        if (which == 0 || which == 2) {
          ListIneq(pos, 1 + u) = sign * body.Lmat(f, u);
        }
        if (which == 1 || which == 2) {
          ListIneq(pos, 4 + u) = sign * body.Lmat(f, u);
        }
      }
    }
    // The equalities also need d_f z <= 1.
    int n_eq = 3;
    MyMatrix<T> Full = ZeroMatrix<T>(rows.size() + n_eq, 7);
    for (size_t pos = 0; pos < rows.size(); pos++) {
      for (int u = 0; u < 7; u++) {
        Full(pos, u) = ListIneq(pos, u);
      }
    }
    int pos2 = rows.size();
    for (size_t pos = 0; pos < rows.size(); pos++) {
      if (rows[pos][2] == 1) {
        for (int u = 0; u < 7; u++) {
          Full(pos2, u) = -ListIneq(pos, u);
        }
        pos2++;
      }
    }
    MyVector<T> ToBeMinimized = ZeroVector<T>(7);
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(Full, ToBeMinimized, os);
    return eSol.DirectSolution.has_value() && eSol.DualSolution.has_value();
  }
};

/*
  The polynomials of degree at most 3 in r <= 3 variables, for det on the
  affine space.
 */
template <typename T> using BetkeHenkPoly = std::map<std::array<int, 3>, T>;

template <typename T>
BetkeHenkPoly<T> BetkeHenk_Mult(BetkeHenkPoly<T> const &P,
                                BetkeHenkPoly<T> const &Q) {
  BetkeHenkPoly<T> R;
  for (auto &kv1 : P) {
    for (auto &kv2 : Q) {
      std::array<int, 3> e;
      for (int i = 0; i < 3; i++) {
        e[i] = kv1.first[i] + kv2.first[i];
      }
      R[e] += kv1.second * kv2.second;
    }
  }
  return R;
}

// A quadratic system F_i(x) = c_i + a_i . x + x^T B_i x, i < r, in double.
struct BetkeHenkQuadSystem {
  int r;
  std::vector<double> c;
  std::vector<std::vector<double>> a;
  std::vector<std::vector<std::vector<double>>> B;
};

// The gradient of p as a quadratic system, in double.
template <typename T>
BetkeHenkQuadSystem BetkeHenk_Gradient(BetkeHenkPoly<T> const &p, int r) {
  BetkeHenkQuadSystem sys;
  sys.r = r;
  sys.c.assign(r, 0);
  sys.a.assign(r, std::vector<double>(r, 0));
  sys.B.assign(r,
               std::vector<std::vector<double>>(r, std::vector<double>(r, 0)));
  for (auto &kv : p) {
    double coef = UniversalScalarConversion<double, T>(kv.second);
    std::array<int, 3> const &e = kv.first;
    for (int i = 0; i < r; i++) {
      if (e[i] == 0) {
        continue;
      }
      // d/dx_i of coef x^e = coef e_i x^{e - 1_i}.
      std::array<int, 3> f = e;
      f[i]--;
      double c2 = coef * e[i];
      int deg = f[0] + f[1] + f[2];
      if (deg == 0) {
        sys.c[i] += c2;
      } else if (deg == 1) {
        for (int k = 0; k < r; k++) {
          if (f[k] == 1) {
            sys.a[i][k] += c2;
          }
        }
      } else {
        // deg 2: either x_k^2 or x_k x_l.
        for (int k = 0; k < r; k++) {
          if (f[k] == 2) {
            sys.B[i][k][k] += c2;
          }
          for (int l = k + 1; l < r; l++) {
            if (f[k] == 1 && f[l] == 1) {
              sys.B[i][k][l] += c2 / 2;
              sys.B[i][l][k] += c2 / 2;
            }
          }
        }
      }
    }
  }
  // Each equation is divided by its largest coefficient: the solutions do
  // not change and the homotopy is better conditioned.
  for (int i = 0; i < r; i++) {
    double mx = std::fabs(sys.c[i]);
    for (int k = 0; k < r; k++) {
      mx = std::max(mx, std::fabs(sys.a[i][k]));
      for (int l = 0; l < r; l++) {
        mx = std::max(mx, std::fabs(sys.B[i][k][l]));
      }
    }
    if (mx > 0) {
      sys.c[i] /= mx;
      for (int k = 0; k < r; k++) {
        sys.a[i][k] /= mx;
        for (int l = 0; l < r; l++) {
          sys.B[i][k][l] /= mx;
        }
      }
    }
  }
  return sys;
}

using BetkeHenkCVec = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using BetkeHenkCMat =
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

inline BetkeHenkCVec BetkeHenk_EvalF(BetkeHenkQuadSystem const &sys,
                                     BetkeHenkCVec const &x) {
  int r = sys.r;
  BetkeHenkCVec F(r);
  for (int i = 0; i < r; i++) {
    std::complex<double> v = sys.c[i];
    for (int k = 0; k < r; k++) {
      v += sys.a[i][k] * x(k);
      for (int l = 0; l < r; l++) {
        v += sys.B[i][k][l] * x(k) * x(l);
      }
    }
    F(i) = v;
  }
  return F;
}

inline BetkeHenkCMat BetkeHenk_EvalJ(BetkeHenkQuadSystem const &sys,
                                     BetkeHenkCVec const &x) {
  int r = sys.r;
  BetkeHenkCMat J(r, r);
  for (int i = 0; i < r; i++) {
    for (int k = 0; k < r; k++) {
      std::complex<double> v = sys.a[i][k];
      for (int l = 0; l < r; l++) {
        v += 2.0 * sys.B[i][k][l] * x(l);
      }
      J(i, k) = v;
    }
  }
  return J;
}

// The endpoints of the total degree homotopy H = (1 - t) gamma G + t F with
// G_i = x_i^2 - 1, tracked in projective space: X = (X_0, X_1, ..., X_r)
// with x = X / X_0, the homogenized equations
//    F^h_i(X) = c_i X_0^2 + sum_k a_ik X_k X_0 + sum_kl B_ikl X_k X_l,
//    G^h_i(X) = X_i^2 - X_0^2,
// and a random affine chart v . X = 1. The paths whose affine solution goes
// to infinity stay bounded and end with X_0 -> 0, so no path is lost to
// the ill-conditioning of large affine coordinates. Each endpoint carries
// whether the Jacobian of F is (numerically) singular there.
struct BetkeHenkRoot {
  std::vector<double> x;
  bool singular;
};

inline std::vector<BetkeHenkRoot>
BetkeHenk_SolveQuadraticGamma(BetkeHenkQuadSystem const &sys,
                              std::complex<double> const &gamma,
                              size_t &n_stalled_mid, size_t &n_stalled_end) {
  int r = sys.r;
  int N = r + 1;
  using C = std::complex<double>;
  BetkeHenkCVec chart(N);
  double chart_re[4] = {0.8147236, -0.1270700, 0.6323592, 0.2784982};
  double chart_im[4] = {0.9057919, 0.9133759, -0.0975404, 0.5468815};
  for (int i = 0; i < N; i++) {
    chart(i) = C(chart_re[i], chart_im[i]);
  }
  auto evalFh = [&](BetkeHenkCVec const &X) -> BetkeHenkCVec {
    BetkeHenkCVec F(r);
    for (int i = 0; i < r; i++) {
      C v = sys.c[i] * X(0) * X(0);
      for (int k = 0; k < r; k++) {
        v += sys.a[i][k] * X(k + 1) * X(0);
        for (int l = 0; l < r; l++) {
          v += sys.B[i][k][l] * X(k + 1) * X(l + 1);
        }
      }
      F(i) = v;
    }
    return F;
  };
  auto evalJFh = [&](BetkeHenkCVec const &X) -> BetkeHenkCMat {
    BetkeHenkCMat J = BetkeHenkCMat::Zero(r, N);
    for (int i = 0; i < r; i++) {
      C d0 = 2.0 * sys.c[i] * X(0);
      for (int k = 0; k < r; k++) {
        d0 += sys.a[i][k] * X(k + 1);
        C dk = sys.a[i][k] * X(0);
        for (int l = 0; l < r; l++) {
          dk += 2.0 * sys.B[i][k][l] * X(l + 1);
        }
        J(i, k + 1) = dk;
      }
      J(i, 0) = d0;
    }
    return J;
  };
  auto evalGh = [&](BetkeHenkCVec const &X) -> BetkeHenkCVec {
    BetkeHenkCVec G(r);
    for (int i = 0; i < r; i++) {
      G(i) = X(i + 1) * X(i + 1) - X(0) * X(0);
    }
    return G;
  };
  auto evalJGh = [&](BetkeHenkCVec const &X) -> BetkeHenkCMat {
    BetkeHenkCMat J = BetkeHenkCMat::Zero(r, N);
    for (int i = 0; i < r; i++) {
      J(i, 0) = -2.0 * X(0);
      J(i, i + 1) = 2.0 * X(i + 1);
    }
    return J;
  };
  // The full system (homotopy and chart) and its Jacobian.
  auto evalH = [&](BetkeHenkCVec const &X, double t) -> BetkeHenkCVec {
    BetkeHenkCVec H(N);
    H.head(r) = (1 - t) * gamma * evalGh(X) + t * evalFh(X);
    H(r) = chart.dot(X) - 1.0;
    return H;
  };
  auto evalHx = [&](BetkeHenkCVec const &X, double t) -> BetkeHenkCMat {
    BetkeHenkCMat J(N, N);
    J.topRows(r) = (1 - t) * gamma * evalJGh(X) + t * evalJFh(X);
    for (int i = 0; i < N; i++) {
      // chart . X with dot conjugating the first argument: use transpose.
      J(r, i) = std::conj(chart(i));
    }
    return J;
  };
  auto evalHt = [&](BetkeHenkCVec const &X) -> BetkeHenkCVec {
    BetkeHenkCVec Ht = BetkeHenkCVec::Zero(N);
    Ht.head(r) = evalFh(X) - gamma * evalGh(X);
    return Ht;
  };
  std::vector<BetkeHenkRoot> ListRoot;
  int n_start = 1 << r;
  for (int s = 0; s < n_start; s++) {
    BetkeHenkCVec X(N);
    X(0) = 1.0;
    for (int i = 0; i < r; i++) {
      X(i + 1) = (s & (1 << i)) ? -1.0 : 1.0;
    }
    X /= chart.dot(X);
    double t = 0;
    double dt = 0.02;
    int n_step = 0;
    while (t < 1 && n_step < 100000) {
      n_step++;
      double dt_eff = std::min(dt, 1 - t);
      auto deriv = [&](BetkeHenkCVec const &Y, double tt) -> BetkeHenkCVec {
        return -evalHx(Y, tt).fullPivLu().solve(evalHt(Y));
      };
      BetkeHenkCVec k1 = deriv(X, t);
      BetkeHenkCVec k2 = deriv(X + 0.5 * dt_eff * k1, t + 0.5 * dt_eff);
      BetkeHenkCVec k3 = deriv(X + 0.5 * dt_eff * k2, t + 0.5 * dt_eff);
      BetkeHenkCVec k4 = deriv(X + dt_eff * k3, t + dt_eff);
      BetkeHenkCVec Y = X + dt_eff / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
      double t_new = t + dt_eff;
      // The corrector. In the second half of the path the Jacobian can be
      // ill-conditioned when the endpoint is singular (the convergence of
      // Newton's method is then only linear), and a looser tolerance is
      // accepted; the step size control still prevents path jumping.
      double tol_corr = t_new > 0.5 ? 1.0e-8 : 1.0e-11;
      bool ok = false;
      for (int it = 0; it < 8; it++) {
        BetkeHenkCVec delta =
            evalHx(Y, t_new).fullPivLu().solve(-evalH(Y, t_new));
        Y += delta;
        if (delta.norm() < tol_corr * (1 + Y.norm())) {
          ok = true;
          break;
        }
      }
      if (!ok || !(Y.norm() < 1.0e12)) {
        dt /= 2;
        if (dt < 1.0e-13) {
          break;
        }
        continue;
      }
      X = Y;
      t = t_new;
      dt = std::min(2 * dt, 0.05);
    }
    // A path that stops before t = 1 (the step size fell to 1e-13). Close
    // to t = 1 this is the endgame of a singular endpoint (a multiple root,
    // a point of a positive dimensional component or a point at infinity,
    // all frequent here), which Newton's method below refines; with the
    // random gamma a stop in the middle of a path has probability zero and
    // could lose a root, so those are counted separately.
    if (t <= 0.9) {
      n_stalled_mid++;
#ifdef DEBUG_POLYNORM_BETKE_HENK
      std::cerr << "BETKE_HENK: path stalled at t=" << t << " |X|=" << X.norm()
                << " |X_0|=" << std::abs(X(0)) << "\n";
#endif
    } else if (t < 1 - 1.0e-6) {
      n_stalled_end++;
    }
    // At infinity: X_0 = 0 relative to X.
    if (std::abs(X(0)) < 1.0e-8 * X.norm()) {
      continue;
    }
    BetkeHenkCVec x = X.tail(r) / X(0);
    // Final Newton steps on F in the affine coordinates.
    for (int it = 0; it < 20; it++) {
      BetkeHenkCMat J = BetkeHenk_EvalJ(sys, x);
      BetkeHenkCVec delta =
          J.completeOrthogonalDecomposition().solve(-BetkeHenk_EvalF(sys, x));
      x += delta;
      if (delta.norm() < 1.0e-14 * (1 + x.norm())) {
        break;
      }
    }
    if (x.norm() > 1.0e8) {
      continue;
    }
    double max_imag = 0;
    for (int i = 0; i < r; i++) {
      max_imag = std::max(max_imag, std::fabs(x(i).imag()));
    }
    if (max_imag > 1.0e-6 * (1 + x.norm())) {
      // A non-real endpoint. It can be a point of a positive dimensional
      // component; the components that matter (those containing a local
      // extremum) are real affine subspaces (Lemma 4.2 of Betke and Henk),
      // so the real part of the point is on them as well. It is kept when
      // it is a root.
      BetkeHenkCVec x_re(r);
      for (int i = 0; i < r; i++) {
        x_re(i) = x(i).real();
      }
      if (BetkeHenk_EvalF(sys, x_re).norm() >
          1.0e-7 * (1 + x_re.squaredNorm())) {
        continue;
      }
      x = x_re;
    }
    if (BetkeHenk_EvalF(sys, x).norm() > 1.0e-7 * (1 + x.squaredNorm())) {
      continue;
    }
    BetkeHenkRoot root;
    for (int i = 0; i < r; i++) {
      root.x.push_back(x(i).real());
    }
    Eigen::JacobiSVD<BetkeHenkCMat> svd(BetkeHenk_EvalJ(sys, x));
    double smax = svd.singularValues()(0);
    double smin = svd.singularValues()(r - 1);
    root.singular = smin < 1.0e-7 * (1 + smax);
    ListRoot.push_back(root);
  }
  return ListRoot;
}

// The homotopy with the gamma trick: a stop in the middle of a path is a
// probability zero event for a random gamma, and when it happens the system
// is solved again with another gamma (the endpoints do not depend on it).
// The attempt with the fewest such stops is kept.
inline std::vector<BetkeHenkRoot>
BetkeHenk_SolveQuadratic(BetkeHenkQuadSystem const &sys, size_t &n_stalled_mid,
                         size_t &n_stalled_end) {
  std::complex<double> ListGamma[4] = {{0.6129831, 0.7900987},
                                       {-0.4161468, 0.9092974},
                                       {0.9800666, -0.1986693},
                                       {-0.8775826, -0.4794255}};
  std::vector<BetkeHenkRoot> best;
  size_t best_mid = std::numeric_limits<size_t>::max();
  size_t best_end = 0;
  for (auto &gamma : ListGamma) {
    size_t mid = 0, end = 0;
    std::vector<BetkeHenkRoot> L =
        BetkeHenk_SolveQuadraticGamma(sys, gamma, mid, end);
    if (mid < best_mid) {
      best = std::move(L);
      best_mid = mid;
      best_end = end;
    }
    if (best_mid == 0) {
      break;
    }
  }
  n_stalled_mid += best_mid;
  n_stalled_end += best_end;
  return best;
}

/*
  The admissibility in double precision: the minimum of ||z W||_D over the
  nonzero z. Since ||x||_D >= |x|_2 / r_max, only the z with
  |z W|_2 < r_max can have ||zW||_D < 1, and they lie in a box given by
  the norms of the columns of W^{-1}.
 */
template <typename T>
double BetkeHenk_MinGauge(BetkeHenkBody<T> const &body,
                          MyMatrix<double> const &W) {
  MyMatrix<double> Winv = W.inverse();
  double R = body.r_max * 1.0001;
  int bound[3];
  for (int i = 0; i < 3; i++) {
    double nrm = 0;
    for (int k = 0; k < 3; k++) {
      nrm += Winv(k, i) * Winv(k, i);
    }
    bound[i] = static_cast<int>(std::floor(R * std::sqrt(nrm))) + 1;
    if (bound[i] > 200) {
      return -1;
    }
  }
  double best = std::numeric_limits<double>::max();
  for (int z0 = -bound[0]; z0 <= bound[0]; z0++) {
    for (int z1 = -bound[1]; z1 <= bound[1]; z1++) {
      for (int z2 = -bound[2]; z2 <= bound[2]; z2++) {
        if (z0 == 0 && z1 == 0 && z2 == 0) {
          continue;
        }
        double x[3];
        for (int k = 0; k < 3; k++) {
          x[k] = z0 * W(0, k) + z1 * W(1, k) + z2 * W(2, k);
        }
        double g = -std::numeric_limits<double>::max();
        for (int j = 0; j < body.m; j++) {
          double s = body.Lmat_d(j, 0) * x[0] + body.Lmat_d(j, 1) * x[1] +
                     body.Lmat_d(j, 2) * x[2];
          g = std::max(g, s);
        }
        best = std::min(best, g);
      }
    }
  }
  return best;
}

/*
  The enumeration.
 */

struct BetkeHenkCase {
  // 1, 2, 3 or 4.
  int case_id;
  std::vector<std::array<int, 3>> U;
};

inline std::vector<BetkeHenkCase> BetkeHenk_Cases() {
  std::vector<std::array<int, 3>> U1{{1, 0, 0},  {0, 1, 0},  {0, 0, 1},
                                     {0, 1, -1}, {-1, 0, 1}, {1, -1, 0}};
  std::vector<std::array<int, 3>> U2{{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                     {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
  std::vector<std::array<int, 3>> U3 = U2;
  U3.push_back({1, 1, 1});
  return {{1, U1}, {2, U2}, {3, U3}, {4, U3}};
}

// A relation u_c = s_a u_a + s_b u_b with a < b < c.
struct BetkeHenkRelation {
  int a, b, c, s_a, s_b;
};

inline std::vector<BetkeHenkRelation>
BetkeHenk_Relations(std::vector<std::array<int, 3>> const &U) {
  std::vector<BetkeHenkRelation> ListRel;
  int k = U.size();
  for (int c = 0; c < k; c++) {
    for (int a = 0; a < c; a++) {
      for (int b = a + 1; b < c; b++) {
        for (int s_a : {1, -1}) {
          for (int s_b : {1, -1}) {
            bool match = true;
            for (int u = 0; u < 3; u++) {
              if (U[c][u] != s_a * U[a][u] + s_b * U[b][u]) {
                match = false;
              }
            }
            if (match) {
              ListRel.push_back({a, b, c, s_a, s_b});
            }
          }
        }
      }
    }
  }
  return ListRel;
}

// A real critical point of det on some S_H, with its test points on the
// boundary of D.
struct BetkeHenkCandidate {
  double det;
  int case_id;
  std::vector<int> facets;
  MyMatrix<double> W;
  bool singular;
  // The minimum of ||z W||_D over z != 0 (1 for an admissible lattice
  // touching D), -1 if not computed.
  double min_gauge;
};

struct BetkeHenkStat {
  size_t n_assignment = 0;
  size_t n_lp_feasible = 0;
  size_t n_rank_ok = 0;
  size_t n_root = 0;
  size_t n_candidate = 0;
  size_t n_case4_degenerate = 0;
  size_t n_partial_lp = 0;
  size_t n_singular = 0;
  size_t n_singular_unresolved = 0;
  size_t n_path_stalled_mid = 0;
  size_t n_path_stalled_end = 0;
  // Time spent (microseconds) in the linear programs, the homotopies and
  // the admissibility checks.
  double t_lp = 0;
  double t_homotopy = 0;
  double t_admissible = 0;
};

template <typename T> struct BetkeHenkSolver {
  BetkeHenkBody<T> const &body;
  std::ostream &os;
  BetkeHenkTriples<T> triples;
  BetkeHenkStat stat;
  std::vector<BetkeHenkCandidate> ListCand;
  // Candidates with |det| above this bound are discarded (the best
  // admissible lattice known, times 1 + tolerance).
  double det_bound;

  BetkeHenkSolver(BetkeHenkBody<T> const &_body, double _det_bound,
                  std::ostream &_os)
      : body(_body), os(_os), triples(_body, _os), det_bound(_det_bound) {}

  int Signed(int j, int s) const { return s == 1 ? j : body.opp[j]; }

  // The rows (c_0, c) of c_0 + c . vec(W) >= 0 expressing u W in the facet
  // F: d_F (u W) = 1 and d_G (u W) <= 1 for the neighbors G of F. With
  // a normal nt (level 1), the additional equality nt (u W) = 1.
  // The cache of the rows d (u W) <= 1 and d (u W) >= 1 for the integer
  // vectors u met so far and the facets d: the arithmetic of the number
  // fields makes recomputing them expensive.
  std::map<std::array<int, 3>, std::vector<std::array<MyVector<T>, 2>>>
      row_cache;

  std::vector<std::array<MyVector<T>, 2>> const &
  GetRows(std::array<int, 3> const &u) {
    auto iter = row_cache.find(u);
    if (iter != row_cache.end()) {
      return iter->second;
    }
    std::vector<std::array<MyVector<T>, 2>> L(body.m);
    for (int j = 0; j < body.m; j++) {
      for (int sign : {1, -1}) {
        MyVector<T> R(10);
        R(0) = T(-sign);
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            R(1 + 3 * a + b) = T(sign * u[a]) * body.Lmat(j, b);
          }
        }
        L[j][sign == 1 ? 0 : 1] = R;
      }
    }
    return row_cache.emplace(u, std::move(L)).first->second;
  }

  // The rows (c_0, c) of c_0 + c . vec(W) >= 0 expressing u W in the facet
  // F: d_F (u W) = 1 and d_G (u W) <= 1 for the neighbors G of F, and the
  // additional equalities extra_eq (normal (u W) = level).
  void AddPointInFacet(std::vector<MyVector<T>> &rows,
                       std::array<int, 3> const &u, int f,
                       std::vector<std::pair<MyVector<T>, T>> const &extra_eq) {
    std::vector<std::array<MyVector<T>, 2>> const &L = GetRows(u);
    rows.push_back(L[f][0]);
    rows.push_back(L[f][1]);
    for (auto &g : body.neighbors[f]) {
      rows.push_back(L[g][1]);
    }
    for (auto &eq : extra_eq) {
      for (int sign : {1, -1}) {
        MyVector<T> R(10);
        R(0) = T(-sign) * eq.second;
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            R(1 + 3 * a + b) = T(sign * u[a]) * eq.first(b);
          }
        }
        rows.push_back(R);
      }
    }
  }

  bool FeasibleRows(std::vector<MyVector<T>> const &rows) {
    MicrosecondTime time;
    bool result = FeasibleRowsKernel(rows);
    stat.t_lp += time.const_eval_int64();
    return result;
  }

  bool FeasibleRowsKernel(std::vector<MyVector<T>> const &rows) {
    MyMatrix<T> ListIneq(rows.size(), 10);
    for (size_t i = 0; i < rows.size(); i++) {
      for (int u = 0; u < 10; u++) {
        ListIneq(i, u) = rows[i](u);
      }
    }
    MyVector<T> ToBeMinimized = ZeroVector<T>(10);
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
    return eSol.DirectSolution.has_value() && eSol.DualSolution.has_value();
  }

  // The rows of the equations u W d^T = h (for the normal d and the level
  // h) in the 9 entries of W, and the inequalities u W in D.
  // The equalities of the affine space: for each test vector u_i a list of
  // (normal, level) with normal (u_i W) = level.
  using EqList = std::vector<std::vector<std::pair<MyVector<T>, T>>>;

  void Process(BetkeHenkCase const &eCase, std::vector<int> const &facets) {
    stat.n_assignment++;
    std::vector<std::array<int, 3>> const &U = eCase.U;
    int k = U.size();
    EqList eqs(k);
    for (int i = 0; i < k; i++) {
      eqs[i].push_back({GetMatrixRow(body.Lmat, facets[i]), T(1)});
    }
    if (eCase.case_id != 4) {
      ProcessAffine(eCase, facets, eqs);
      return;
    }
    // Case IV: lambda_1 v_1 + lambda_2 v_2 + lambda_3 v_3 = lambda_6 v_6.
    std::vector<MyVector<T>> v(k);
    for (int i = 0; i < k; i++) {
      v[i] = GetMatrixRow(body.Lmat, facets[i]);
    }
    MyMatrix<T> M(4, 3);
    for (int u = 0; u < 3; u++) {
      M(0, u) = v[0](u);
      M(1, u) = v[1](u);
      M(2, u) = v[2](u);
      M(3, u) = -v[5](u);
    }
    MyMatrix<T> NSP = NullspaceTrMat(TransposedMat(M));
    if (NSP.rows() == 1) {
      MyVector<T> lambda = GetMatrixRow(NSP, 0);
      if (lambda(3) < 0) {
        lambda = -lambda;
      }
      if (!(lambda(0) > 0 && lambda(1) > 0 && lambda(3) > 0)) {
        return;
      }
      MyVector<T> nt = lambda(0) * v[0] + lambda(1) * v[1];
      // The supporting hyperplane of normal nt.
      std::optional<T> h;
      for (auto &x : body.ListVert) {
        T sc = nt.dot(x);
        if (!h || sc > *h) {
          h = sc;
        }
      }
      if (!h || *h <= 0) {
        return;
      }
      MyVector<T> ntn = nt / *h;
      if (ntn == v[5]) {
        // The facet hyperplane itself: this is case III.
        return;
      }
      eqs[5] = {{v[5], T(1)}, {ntn, T(1)}};
      // Only H_6 tilde defines the affine space, the facet F_6 only
      // constrains the feasibility: the first equation is dropped from
      // the affine space but kept for the feasibility.
      ProcessAffine(eCase, facets, eqs, 1);
      return;
    }
    // Degenerate case IV: v_1, v_2, v_3, v_6 in a plane, the normal of
    // H_6 tilde is not unique. The critical lattice is a local minimum on
    // S_H tilde for some of them, hence on the smaller space where u_6 W
    // lies in the affine hull of the face G of D exposed by that normal,
    // and G meets F_6. So it is enough to take for G every edge and every
    // vertex of F_6: a superset of the candidates.
    stat.n_case4_degenerate++;
    int f6 = facets[5];
    for (auto &g : body.neighbors[f6]) {
      EqList eqs2 = eqs;
      eqs2[5] = {{v[5], T(1)}, {GetMatrixRow(body.Lmat, g), T(1)}};
      ProcessAffine(eCase, facets, eqs2);
    }
    for (auto &iv : body.facet_vertices[f6]) {
      EqList eqs2 = eqs;
      // The facet equation first (as AddPointInFacet expects), then u_6 W
      // equal to the vertex.
      eqs2[5] = {{v[5], T(1)}};
      for (int b = 0; b < 3; b++) {
        MyVector<T> e = ZeroVector<T>(3);
        e(b) = 1;
        eqs2[5].push_back({e, body.ListVert[iv](b)});
      }
      ProcessAffine(eCase, facets, eqs2);
    }
  }

  // The test points in D for the lattice W (they are on their hyperplanes).
  bool TestPointsInD(std::vector<std::array<int, 3>> const &U,
                     MyMatrix<double> const &W) const {
    for (auto &u : U) {
      double x[3];
      for (int b = 0; b < 3; b++) {
        x[b] = u[0] * W(0, b) + u[1] * W(1, b) + u[2] * W(2, b);
      }
      for (int j = 0; j < body.m; j++) {
        double sc = body.Lmat_d(j, 0) * x[0] + body.Lmat_d(j, 1) * x[1] +
                    body.Lmat_d(j, 2) * x[2];
        if (sc > 1 + 1.0e-7) {
          return false;
        }
      }
    }
    return true;
  }

  void AddCandidate(BetkeHenkCase const &eCase, std::vector<int> const &facets,
                    MyMatrix<double> const &W, bool singular) {
    double dt = std::fabs(W.determinant());
    stat.n_candidate++;
    MicrosecondTime time_a;
    double mg = BetkeHenk_MinGauge(body, W);
    stat.t_admissible += time_a.const_eval_int64();
    ListCand.push_back({dt, eCase.case_id, facets, W, singular, mg});
  }

  // The affine space of the equalities eqs (skipping the first n_skip_6
  // equalities of the test vector 6 for the affine space, which are only
  // used for the feasibility), its feasibility, and the critical points of
  // det on it.
  void ProcessAffine(BetkeHenkCase const &eCase, std::vector<int> const &facets,
                     EqList const &eqs, int n_skip_6 = 0) {
    std::vector<std::array<int, 3>> const &U = eCase.U;
    int k = U.size();
    // Feasibility: u_i W in the facet F_i and on its other hyperplanes.
    std::vector<MyVector<T>> rows;
    for (int i = 0; i < k; i++) {
      std::vector<std::pair<MyVector<T>, T>> extra(eqs[i].begin() + 1,
                                                   eqs[i].end());
      AddPointInFacet(rows, U[i], facets[i], extra);
    }
    if (!FeasibleRows(rows)) {
      return;
    }
    stat.n_lp_feasible++;
    // The affine space.
    std::vector<MyVector<T>> ListRow;
    std::vector<T> ListLevel;
    for (int i = 0; i < k; i++) {
      for (size_t e = 0; e < eqs[i].size(); e++) {
        if (i == 5 && static_cast<int>(e) < n_skip_6) {
          continue;
        }
        MyVector<T> R(9);
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            R(3 * a + b) = U[i][a] * eqs[i][e].first(b);
          }
        }
        ListRow.push_back(R);
        ListLevel.push_back(eqs[i][e].second);
      }
    }
    int n_eq = ListRow.size();
    MyMatrix<T> A(n_eq, 9);
    MyVector<T> bvec(n_eq);
    for (int i = 0; i < n_eq; i++) {
      for (int u = 0; u < 9; u++) {
        A(i, u) = ListRow[i](u);
      }
      bvec(i) = ListLevel[i];
    }
    int rnk = RankMat(A);
    if (rnk < 6) {
      // Lemma 4.1: no local minimum with det != 0.
      return;
    }
    stat.n_rank_ok++;
    MyMatrix<T> AT = TransposedMat(A);
    std::optional<MyVector<T>> opt = SolutionMat(AT, bvec);
    if (!opt) {
      return;
    }
    MyMatrix<T> NSP = NullspaceTrMat(A);
    int r = NSP.rows();
    // The directions are scaled by powers of 2 (exact in every field) to
    // have a norm close to that of the particular solution: the integer
    // scaled kernel vectors can be huge, which makes the homotopy badly
    // conditioned.
    {
      double nrm0 = 1;
      for (int u = 0; u < 9; u++) {
        double x = UniversalScalarConversion<double, T>((*opt)(u));
        nrm0 = std::max(nrm0, std::fabs(x));
      }
      for (int i = 0; i < r; i++) {
        double nrm = 0;
        for (int u = 0; u < 9; u++) {
          double x = UniversalScalarConversion<double, T>(NSP(i, u));
          nrm = std::max(nrm, std::fabs(x));
        }
        int e = static_cast<int>(std::lround(std::log2(nrm / nrm0)));
        T fact(1);
        for (int j = 0; j < std::abs(e); j++) {
          fact *= 2;
        }
        for (int u = 0; u < 9; u++) {
          if (e > 0) {
            NSP(i, u) /= fact;
          } else {
            NSP(i, u) *= fact;
          }
        }
      }
    }
    MyVector<double> w0 = UniversalVectorConversion<double, T>(*opt);
    MyMatrix<double> NSP_d = UniversalMatrixConversion<double, T>(NSP);
    auto get_W = [&](std::vector<double> const &x) -> MyMatrix<double> {
      MyMatrix<double> W(3, 3);
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          double val = w0(3 * a + b);
          for (int i = 0; i < r; i++) {
            val += x[i] * NSP_d(i, 3 * a + b);
          }
          W(a, b) = val;
        }
      }
      return W;
    };
    if (r == 0) {
      MyMatrix<double> W = get_W({});
      double dt = std::fabs(W.determinant());
      if (dt > 1.0e-9 && dt <= det_bound && TestPointsInD(U, W)) {
        AddCandidate(eCase, facets, W, false);
      }
      return;
    }
    // W(lambda) entries as affine forms.
    std::vector<std::vector<BetkeHenkPoly<T>>> Wpoly(
        3, std::vector<BetkeHenkPoly<T>>(3));
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        int idx = 3 * a + b;
        BetkeHenkPoly<T> P;
        if ((*opt)(idx) != 0) {
          P[{0, 0, 0}] = (*opt)(idx);
        }
        for (int i = 0; i < r; i++) {
          if (NSP(i, idx) != 0) {
            std::array<int, 3> e{0, 0, 0};
            e[i] = 1;
            P[e] = NSP(i, idx);
          }
        }
        Wpoly[a][b] = P;
      }
    }
    // det by the Leibniz formula.
    BetkeHenkPoly<T> det;
    int perms[6][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1},
                       {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
    int signs[6] = {1, 1, 1, -1, -1, -1};
    for (int sg = 0; sg < 6; sg++) {
      BetkeHenkPoly<T> P = BetkeHenk_Mult(
          BetkeHenk_Mult(Wpoly[0][perms[sg][0]], Wpoly[1][perms[sg][1]]),
          Wpoly[2][perms[sg][2]]);
      for (auto &kv : P) {
        det[kv.first] += T(signs[sg]) * kv.second;
      }
    }
    BetkeHenkQuadSystem sys = BetkeHenk_Gradient(det, r);
    MicrosecondTime time_h;
    std::vector<BetkeHenkRoot> ListRoot = BetkeHenk_SolveQuadratic(
        sys, stat.n_path_stalled_mid, stat.n_path_stalled_end);
    stat.t_homotopy += time_h.const_eval_int64();
    stat.n_root += ListRoot.size();
    for (auto &root : ListRoot) {
      MyMatrix<double> W = get_W(root.x);
      double dt = std::fabs(W.determinant());
      if (dt < 1.0e-9 || dt > det_bound) {
        continue;
      }
      bool in_D = TestPointsInD(U, W);
      if (root.singular) {
        // A point of a positive dimensional component of critical points,
        // on which det is constant (Lemma 4.2). Another point of the
        // component may have its test points in D and be admissible.
        std::optional<MyMatrix<double>> optW =
            ResolveComponent(U, sys, root.x, get_W);
        if (optW) {
          AddCandidate(eCase, facets, *optW, false);
          continue;
        }
      }
      if (!in_D) {
        continue;
      }
      AddCandidate(eCase, facets, W, root.singular);
    }
  }

  // For a singular root x0 of the gradient system: the affine components
  // x0 + span(V) of the critical set through x0 (V in the kernel of the
  // Jacobian, with the quadratic parts vanishing on V), and on each a point
  // whose test points lie in D, which is returned if its lattice is
  // admissible.
  template <typename Fget>
  std::optional<MyMatrix<double>>
  ResolveComponent(std::vector<std::array<int, 3>> const &U,
                   BetkeHenkQuadSystem const &sys,
                   std::vector<double> const &x0, Fget get_W) {
    stat.n_singular++;
    int r = sys.r;
    MyMatrix<double> J(r, r);
    for (int i = 0; i < r; i++) {
      for (int kk = 0; kk < r; kk++) {
        double val = sys.a[i][kk];
        for (int l = 0; l < r; l++) {
          val += 2 * sys.B[i][kk][l] * x0[l];
        }
        J(i, kk) = val;
      }
    }
    Eigen::JacobiSVD<MyMatrix<double>> svd(J, Eigen::ComputeFullV);
    MyVector<double> sv = svd.singularValues();
    double smax = std::max(sv(0), 1.0);
    std::vector<MyVector<double>> ker;
    for (int i = 0; i < r; i++) {
      if (sv(i) < 1.0e-7 * smax) {
        ker.push_back(svd.matrixV().col(i));
      }
    }
    auto quad = [&](int i, MyVector<double> const &p,
                    MyVector<double> const &q) -> double {
      double val = 0;
      for (int kk = 0; kk < r; kk++) {
        for (int l = 0; l < r; l++) {
          val += sys.B[i][kk][l] * p(kk) * q(l);
        }
      }
      return val;
    };
    double tolq = 1.0e-7;
    // The candidate subspaces V.
    std::vector<std::vector<MyVector<double>>> ListV;
    int dk = ker.size();
    if (dk == 1) {
      bool ok = true;
      for (int i = 0; i < r; i++) {
        if (std::fabs(quad(i, ker[0], ker[0])) > tolq) {
          ok = false;
        }
      }
      if (ok) {
        ListV.push_back({ker[0]});
      }
    } else if (dk == 2) {
      bool all_zero = true;
      for (int i = 0; i < r; i++) {
        if (std::fabs(quad(i, ker[0], ker[0])) > tolq ||
            std::fabs(quad(i, ker[0], ker[1])) > tolq ||
            std::fabs(quad(i, ker[1], ker[1])) > tolq) {
          all_zero = false;
        }
      }
      if (all_zero) {
        ListV.push_back({ker[0], ker[1]});
      } else {
        // The directions cos(t) k0 + sin(t) k1 isotropic for every form:
        // the roots of the first nonzero binary form, checked on the
        // others.
        for (int i = 0; i < r; i++) {
          double qa = quad(i, ker[0], ker[0]);
          double qb = 2 * quad(i, ker[0], ker[1]);
          double qc = quad(i, ker[1], ker[1]);
          if (std::fabs(qa) + std::fabs(qb) + std::fabs(qc) <= 3 * tolq) {
            continue;
          }
          std::vector<MyVector<double>> ListDir;
          if (std::fabs(qa) <= tolq) {
            ListDir.push_back(ker[0]);
            if (std::fabs(qb) > tolq) {
              ListDir.push_back(-qc * ker[0] + qb * ker[1]);
            }
          } else {
            double disc = qb * qb - 4 * qa * qc;
            if (disc >= -tolq) {
              double sq = std::sqrt(std::max(disc, 0.0));
              for (double sgn : {1.0, -1.0}) {
                double tt = (-qb + sgn * sq) / (2 * qa);
                ListDir.push_back(tt * ker[0] + ker[1]);
              }
            }
          }
          for (auto &dir : ListDir) {
            MyVector<double> dn = dir / dir.norm();
            bool ok = true;
            for (int i2 = 0; i2 < r; i2++) {
              if (std::fabs(quad(i2, dn, dn)) > 1.0e-6) {
                ok = false;
              }
            }
            if (ok) {
              ListV.push_back({dn});
            }
          }
          break;
        }
      }
    } else if (dk >= 3) {
      stat.n_singular_unresolved++;
      return {};
    }
    for (auto &V : ListV) {
      std::optional<std::vector<double>> opt_mu =
          FeasibleOnComponent(U, x0, V, get_W);
      if (!opt_mu) {
        continue;
      }
      std::vector<double> x = x0;
      for (size_t c = 0; c < V.size(); c++) {
        for (int i = 0; i < r; i++) {
          x[i] += (*opt_mu)[c] * V[c](i);
        }
      }
      MyMatrix<double> W = get_W(x);
      double mg = BetkeHenk_MinGauge(body, W);
      if (mg >= 1 - 1.0e-7) {
        return W;
      }
    }
    if (ListV.empty()) {
      // An isolated singular root (a multiple root).
      return {};
    }
    return {};
  }

  // A point mu (of dimension 1 or 2) with the test points of
  // W(x0 + V mu) in D, as central as possible: the midpoint of the
  // interval, or the average of the vertices of the polygon.
  template <typename Fget>
  std::optional<std::vector<double>>
  FeasibleOnComponent(std::vector<std::array<int, 3>> const &U,
                      std::vector<double> const &x0,
                      std::vector<MyVector<double>> const &V, Fget get_W) {
    int r = x0.size();
    int dv = V.size();
    MyMatrix<double> W0 = get_W(x0);
    std::vector<MyMatrix<double>> ListDW;
    for (int c = 0; c < dv; c++) {
      std::vector<double> x(r);
      for (int i = 0; i < r; i++) {
        x[i] = x0[i] + V[c](i);
      }
      ListDW.push_back(get_W(x) - W0);
    }
    // Constraints c0 + c . mu <= 1 for each test point and facet.
    std::vector<std::pair<double, std::vector<double>>> cons;
    for (auto &u : U) {
      for (int j = 0; j < body.m; j++) {
        auto lin = [&](MyMatrix<double> const &M) -> double {
          double val = 0;
          for (int b = 0; b < 3; b++) {
            double xb = u[0] * M(0, b) + u[1] * M(1, b) + u[2] * M(2, b);
            val += body.Lmat_d(j, b) * xb;
          }
          return val;
        };
        std::vector<double> cc(dv);
        for (int c = 0; c < dv; c++) {
          cc[c] = lin(ListDW[c]);
        }
        cons.push_back({lin(W0), cc});
      }
    }
    double tol = 1.0e-9;
    if (dv == 1) {
      double lo = -1.0e30, hi = 1.0e30;
      for (auto &eC : cons) {
        double c0 = eC.first, c1 = eC.second[0];
        if (std::fabs(c1) < 1.0e-14) {
          if (c0 > 1 + tol) {
            return {};
          }
        } else if (c1 > 0) {
          hi = std::min(hi, (1 - c0) / c1);
        } else {
          lo = std::max(lo, (1 - c0) / c1);
        }
      }
      if (lo > hi + tol || lo < -1.0e29 || hi > 1.0e29) {
        return {};
      }
      return std::vector<double>{(lo + hi) / 2};
    }
    // dv == 2: the vertices of the polygon.
    std::vector<std::array<double, 2>> ListVertex;
    int nc = cons.size();
    for (int i1 = 0; i1 < nc; i1++) {
      for (int i2 = i1 + 1; i2 < nc; i2++) {
        double a11 = cons[i1].second[0], a12 = cons[i1].second[1];
        double a21 = cons[i2].second[0], a22 = cons[i2].second[1];
        double dd = a11 * a22 - a12 * a21;
        if (std::fabs(dd) < 1.0e-12) {
          continue;
        }
        double r1 = 1 - cons[i1].first, r2 = 1 - cons[i2].first;
        double m0 = (r1 * a22 - r2 * a12) / dd;
        double m1 = (a11 * r2 - a21 * r1) / dd;
        bool ok = true;
        for (auto &eC : cons) {
          if (eC.first + eC.second[0] * m0 + eC.second[1] * m1 > 1 + 1.0e-7) {
            ok = false;
            break;
          }
        }
        if (ok) {
          ListVertex.push_back({m0, m1});
        }
      }
    }
    if (ListVertex.empty()) {
      return {};
    }
    double s0 = 0, s1 = 0;
    for (auto &v : ListVertex) {
      s0 += v[0];
      s1 += v[1];
    }
    return std::vector<double>{s0 / ListVertex.size(), s1 / ListVertex.size()};
  }

  // Whether some W puts u_0, ..., u_{len-1} in their facets.
  bool PartialFeasible(std::vector<std::array<int, 3>> const &U,
                       std::vector<int> const &facets, int len) {
    stat.n_partial_lp++;
    std::vector<MyVector<T>> rows;
    for (int i = 0; i < len; i++) {
      AddPointInFacet(rows, U[i], facets[i], {});
    }
    return FeasibleRows(rows);
  }

  // The backtracking over the facets of the test vectors.
  // With several threads, the depth 2 prefixes (facets of the first two
  // test vectors) are claimed dynamically: every thread enumerates the same
  // sequence of prefixes and processes the prefix of index my_next, then
  // claims the next free index from the shared counter. So each prefix is
  // processed by exactly one thread, whichever is free.
  std::atomic<int> *shared_counter = nullptr;
  int my_next = -1;

  void Enumerate(BetkeHenkCase const &eCase) {
    int prefix_index = 0;
    if (shared_counter) {
      my_next = shared_counter->fetch_add(1);
    }
    std::vector<std::array<int, 3>> const &U = eCase.U;
    int k = U.size();
    std::vector<BetkeHenkRelation> ListRel = BetkeHenk_Relations(U);
    std::vector<int> facets(k, -1);
    int m = body.m;
    int n_sym = eCase.case_id == 4 ? 2 : 3;
    // The relations whose last vector is c.
    std::vector<std::vector<BetkeHenkRelation>> rel_by_c(k);
    for (auto &rel : ListRel) {
      rel_by_c[rel.c].push_back(rel);
    }
    // Pair compatibility for the first three vectors: for a relation
    // u_c = s_a u_a + s_b u_b with c >= 3 and a, b < 3, some facet F_c
    // must exist. Computed lazily.
    std::map<std::array<int, 4>, bool> pair_ok;
    auto pair_feasible = [&](int ja, int jb) -> bool {
      // ja, jb already signed.
      std::array<int, 4> key{ja, jb, 0, 0};
      auto iter = pair_ok.find(key);
      if (iter != pair_ok.end()) {
        return iter->second;
      }
      bool ok = false;
      for (int c = 0; c < m && !ok; c++) {
        if (triples.Feasible(ja, jb, c)) {
          ok = true;
        }
      }
      pair_ok[key] = ok;
      return ok;
    };
    std::function<void(int)> rec = [&](int i) {
      if (i == k) {
        Process(eCase, facets);
        return;
      }
      for (int j = 0; j < m; j++) {
        facets[i] = j;
        // Symmetry. The permutations of the basis vectors that preserve
        // the case (all of them for the cases I to III, the exchange of
        // the first two for the case IV, where u_6 = u_1 + u_2 plays a
        // special role) permute the first n_sym test vectors, and the
        // symmetries of D act on the facets. So the first vector can be
        // taken with the smallest orbit among the first n_sym, and its
        // facet the representative of that orbit.
        if (i == 0 && !body.is_rep[j]) {
          continue;
        }
        if (i > 0 && i < n_sym && body.orbit[j] < body.orbit[facets[0]]) {
          continue;
        }
        bool ok = true;
        if (i < 3) {
          // The relations with a, b among the first three vectors and
          // a later c give pair conditions.
          for (auto &rel : ListRel) {
            if (rel.b == i && rel.c >= 3) {
              int ja = Signed(facets[rel.a], rel.s_a);
              int jb = Signed(facets[rel.b], rel.s_b);
              if (!pair_feasible(ja, jb)) {
                ok = false;
                break;
              }
            }
          }
        }
        if (ok) {
          for (auto &rel : rel_by_c[i]) {
            int ja = Signed(facets[rel.a], rel.s_a);
            int jb = Signed(facets[rel.b], rel.s_b);
            if (!triples.Feasible(ja, jb, j)) {
              ok = false;
              break;
            }
          }
        }
        bool claimed = false;
        if (ok && i == 1 && shared_counter) {
          int idx = prefix_index++;
          if (idx != my_next) {
            ok = false;
          } else {
            claimed = true;
          }
        }
        if (ok && (i == 4 || i == 5) && i < k - 1) {
          // Partial feasibility: the first i + 1 points in their facets.
          ok = PartialFeasible(U, facets, i + 1);
        }
        if (ok) {
          rec(i + 1);
        }
        if (claimed) {
          my_next = shared_counter->fetch_add(1);
        }
      }
      facets[i] = -1;
    };
    rec(0);
  }
};

// The enumeration of a case by n_thread threads, each with its own solver
// (memoized triples, cached rows, candidates), merged at the end. The
// program must use the thread-local GMP pool (maybe_install_tls_gmp_pool).
template <typename T>
void BetkeHenk_EnumerateParallel(BetkeHenkSolver<T> &solver,
                                 BetkeHenkCase const &eCase, int n_thread) {
  if (n_thread <= 1) {
    solver.Enumerate(eCase);
    return;
  }
  std::atomic<int> counter(0);
  std::vector<BetkeHenkSolver<T>> ListSolver;
  ListSolver.reserve(n_thread);
  for (int t = 0; t < n_thread; t++) {
    ListSolver.emplace_back(solver.body, solver.det_bound, solver.os);
    ListSolver.back().shared_counter = &counter;
  }
  std::vector<std::thread> ListThread;
  for (int t = 0; t < n_thread; t++) {
    ListThread.emplace_back(
        [&ListSolver, &eCase, t]() { ListSolver[t].Enumerate(eCase); });
  }
  for (auto &th : ListThread) {
    th.join();
  }
  for (auto &sub : ListSolver) {
    for (auto &c : sub.ListCand) {
      solver.ListCand.push_back(c);
    }
    BetkeHenkStat &a = solver.stat;
    BetkeHenkStat const &b = sub.stat;
    a.n_assignment += b.n_assignment;
    a.n_lp_feasible += b.n_lp_feasible;
    a.n_rank_ok += b.n_rank_ok;
    a.n_root += b.n_root;
    a.n_candidate += b.n_candidate;
    a.n_case4_degenerate += b.n_case4_degenerate;
    a.n_partial_lp += b.n_partial_lp;
    a.n_singular += b.n_singular;
    a.n_singular_unresolved += b.n_singular_unresolved;
    a.n_path_stalled_mid += b.n_path_stalled_mid;
    a.n_path_stalled_end += b.n_path_stalled_end;
    a.t_lp += b.t_lp;
    a.t_homotopy += b.t_homotopy;
    a.t_admissible += b.t_admissible;
    solver.triples.n_lp += sub.triples.n_lp;
    solver.triples.n_box_reject += sub.triples.n_box_reject;
  }
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_BETKEHENK_H_
// clang-format on
