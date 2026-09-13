// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_JOINTCOVERINGDOUBLE_H_
#define SRC_DELAUNAY_JOINTCOVERINGDOUBLE_H_

// clang-format off
#include <Eigen/Dense>
#include <libqhull_r/qhull_ra.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  Joint optimization of the covering density of a periodic point set
  Z^n + {C_0 = 0, C_1, ..., C_{m-1}} over the Gram matrix Q and the
  continuous coset positions C, entirely in double precision.

  This is the search tool: fast and numerical, complementing the exact
  machinery of PeriodicDelaunay.h. It is meant for configurations too rich
  for the exact side -- ten cosets and thousands of Delaunay simplices --
  where an interesting candidate, once found, is to be re-verified exactly.

  The tessellation is the Delaunay triangulation, in the metric Q, of the
  points of the set of Q-norm at most R, computed by qhull. A retained cell
  is provably a Delaunay cell of the full periodic set once its whole
  circumsphere lies inside the enumerated ball (|center|_Q + radius <= R);
  when a cell incident to the central copy fails that test, R grows and
  everything is redone. The enumeration runs in an LLL-reduced basis so that
  the ball needs few lattice shells.

  The optimization minimizes the smooth surrogate

     F(Q, C) = (n/2) log smax_beta { R^2_S(Q, C) } - (1/2) log det Q,

  smax_beta the soft maximum over the translation classes S of Delaunay
  cells, by limited-memory BFGS with analytic gradients; beta is annealed
  and the cell list refreshed between stages, each stage's result being
  verified against a fresh tessellation. As beta -> infinity, F tends to
  log Theta up to an additive constant.
 */

#ifdef DEBUG
#define DEBUG_JOINT_COVERING_DOUBLE
#endif

#ifdef TIMINGS
#define TIMINGS_JOINT_COVERING_DOUBLE
#endif

namespace joint_covering_double {

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

// ---------------------------------------------------------------------------
// Configurations
// ---------------------------------------------------------------------------

struct PeriodicConfig {
  int n;
  int m;
  MatrixXd Q;   // n x n positive definite
  MatrixXd C;   // m x n, first row zero
};

// One translation class of Delaunay cells: the lattice parts of its n+1
// vertices and their coset indices.
struct CellClass {
  MatrixXi lat;   // (n+1) x n
  VectorXi cos;   // n+1
};

// ---------------------------------------------------------------------------
// LLL reduction (double precision, on the Gram matrix)
// ---------------------------------------------------------------------------

inline MatrixXi LLLReduceDouble(MatrixXd const &Q, double delta = 0.75) {
  int n = Q.rows();
  Eigen::LLT<MatrixXd> llt(Q);
  MatrixXd B = llt.matrixL();          // rows of B: basis vectors
  B.transposeInPlace();                // now B(i,:) . B(j,:) = Q(i,j)? no:
  B = MatrixXd(llt.matrixL());         // L with Q = L L^T: rows of L work:
  // (L L^T)(i,j) = row_i(L) . row_j(L) = Q(i,j), so rows of L are the basis.
  MatrixXi U = MatrixXi::Identity(n, n);
  auto gso = [&](MatrixXd const &Bc, MatrixXd &Bs, MatrixXd &mu) -> void {
    Bs = Bc;
    mu = MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        mu(i, j) = Bc.row(i).dot(Bs.row(j)) / Bs.row(j).squaredNorm();
        Bs.row(i) -= mu(i, j) * Bs.row(j);
      }
    }
  };
  MatrixXd Bs, mu;
  int k = 1;
  int guard = 0;
  while (k < n && guard++ < 10000) {
    gso(B, Bs, mu);
    for (int j = k - 1; j >= 0; j--) {
      long r = std::lround(mu(k, j));
      if (r != 0) {
        B.row(k) -= double(r) * B.row(j);
        U.row(k) -= int(r) * U.row(j);
        gso(B, Bs, mu);
      }
    }
    if (Bs.row(k).squaredNorm() >=
        (delta - mu(k, k - 1) * mu(k, k - 1)) * Bs.row(k - 1).squaredNorm()) {
      k++;
    } else {
      B.row(k).swap(B.row(k - 1));
      U.row(k).swap(U.row(k - 1));
      k = std::max(k - 1, 1);
    }
  }
  return U;
}

// ---------------------------------------------------------------------------
// The tessellation
// ---------------------------------------------------------------------------

/*
  Circumcenter (in the coordinates of the point set) and squared
  circumradius, in the metric Q, of the simplex with vertex rows P.
 */
inline void CircumcenterRadius2(MatrixXd const &P, MatrixXd const &Q,
                                VectorXd &center, double &R2) {
  int n = Q.rows();
  MatrixXd V(n, n);
  for (int k = 0; k < n; k++) {
    V.row(k) = P.row(k + 1) - P.row(0);
  }
  MatrixXd G = V * Q * V.transpose();
  VectorXd q = G.diagonal();
  VectorXd a = G.ldlt().solve(q);
  R2 = 0.25 * q.dot(a);
  center = P.row(0).transpose() + 0.5 * V.transpose() * a;
}

/*
  All translation classes of Delaunay cells of the periodic set at (Q, C).
  Throws std::runtime_error when the ball never stabilizes.
 */
inline std::vector<CellClass> DelaunayCellClasses(PeriodicConfig const &conf,
                                                  int max_iter = 6,
                                                  double growth = 1.3) {
  int n = conf.n;
  int m = conf.m;
  MatrixXd const &Q = conf.Q;
  // A non positive definite form has no covering geometry and would send the
  // ball-growth loop below into an unbounded spin; reject it at once. Eigen's
  // LLT can accept a mildly indefinite matrix, so the smallest eigenvalue is
  // tested directly.
  {
    Eigen::SelfAdjointEigenSolver<MatrixXd> esg(Q, Eigen::EigenvaluesOnly);
    if (esg.eigenvalues()(0) <= 1e-12 * esg.eigenvalues()(n - 1)) {
      throw std::runtime_error(
          "DelaunayCellClasses: Q is not positive definite");
    }
  }
  MatrixXi U = LLLReduceDouble(Q);
  MatrixXd Ud = U.cast<double>();
  MatrixXd Qr = Ud * Q * Ud.transpose();
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(Qr);
  double sig_min = std::sqrt(std::max(es.eigenvalues()(0), 1e-14));
  // covering radius of the sublattice Z^n bounds the one of the union
  double mu_bound = 0.5 * std::sqrt(Qr.trace());
  double R = 2.5 * mu_bound;
  Eigen::LLT<MatrixXd> llt(Q);
  MatrixXd L = llt.matrixL();
  // wall-clock deadline: a pathological (Q, c) -- near-cocircular cosets --
  // makes the joggled qhull and the ball-emptiness certificate thrash. Rather
  // than stall the whole search, such a configuration is abandoned.
  auto t_start = std::chrono::steady_clock::now();
  for (int iter = 0; iter < max_iter; iter++) {
    if (std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - t_start).count() > 20) {
      throw std::runtime_error("DelaunayCellClasses: tessellation deadline");
    }
    // enumerate the points of Q-norm <= R through the reduced basis
    int B0 = int(std::ceil(R / sig_min)) + 1;
    std::vector<VectorXi> lats;
    std::vector<int> coss;
    std::vector<double> coords;   // qhull input, transformed x -> L^T x
    VectorXi b = VectorXi::Constant(n, -B0);
    while (true) {
      VectorXi a = U.transpose() * b;   // a = b U in row convention
      VectorXd af = a.cast<double>();
      for (int t = 0; t < m; t++) {
        VectorXd x = af + conf.C.row(t).transpose();
        double nrm2 = x.dot(Q * x);
        if (nrm2 <= R * R) {
          lats.push_back(a);
          coss.push_back(t);
          VectorXd y = L.transpose() * x;
          for (int i = 0; i < n; i++) {
            coords.push_back(y(i));
          }
        }
      }
      int pos = 0;
      while (pos < n && b(pos) == B0) {
        b(pos) = -B0;
        pos++;
      }
      if (pos == n) {
        break;
      }
      b(pos)++;
    }
    int n_pts = lats.size();
    // qhull Delaunay of the transformed cloud
    qhT qh_qh;
    qhT *qh = &qh_qh;
    qh_zero(qh, stderr);
    FILE *devnull = fopen("/dev/null", "w");
    char options[] = "qhull d QJ Qbb";
    int exitcode = qh_new_qhull(qh, n, n_pts, coords.data(), False, options,
                                nullptr, devnull);
    std::map<std::vector<int>, CellClass> classes;
    bool ok = (exitcode == 0);
    if (ok) {
      facetT *facet;
      vertexT *vertex, **vertexp;
      MatrixXd P(n + 1, n);
      MatrixXi lat(n + 1, n);
      VectorXi cos(n + 1);
      FORALLfacets {
        if (facet->upperdelaunay) {
          continue;
        }
        int idx = 0;
        bool has_central = false;
        FOREACHvertex_(facet->vertices) {
          int id = qh_pointid(qh, vertex->point);
          if (idx > n) {
            break;
          }
          lat.row(idx) = lats[id].transpose();
          cos(idx) = coss[id];
          P.row(idx) = lats[id].cast<double>().transpose() +
                       conf.C.row(coss[id]);
          if (lats[id].isZero()) {
            has_central = true;
          }
          idx++;
        }
        if (idx != n + 1 || !has_central) {
          continue;
        }
        // joggle slivers: flat pieces of a split cocircular cell
        MatrixXd V(n, n);
        for (int k = 0; k < n; k++) {
          V.row(k) = P.row(k + 1) - P.row(0);
        }
        if (std::abs(V.determinant()) < 1e-9) {
          continue;
        }
        VectorXd center;
        double R2;
        CircumcenterRadius2(P, Q, center, R2);
        // the certificate making the cell a Delaunay cell of the full set
        if (std::sqrt(center.dot(Q * center)) + std::sqrt(std::max(R2, 0.0)) >
            R) {
          ok = false;
          break;
        }
        // canonical key: minimum over all vertex anchors of the sorted rows
        std::vector<int> best_key;
        CellClass best_val;
        for (int i0 = 0; i0 <= n; i0++) {
          std::vector<std::vector<int>> rows(n + 1,
                                             std::vector<int>(n + 1));
          for (int k = 0; k <= n; k++) {
            for (int j = 0; j < n; j++) {
              rows[k][j] = lat(k, j) - lat(i0, j);
            }
            rows[k][n] = cos(k);
          }
          std::sort(rows.begin(), rows.end());
          std::vector<int> key;
          key.reserve((n + 1) * (n + 1));
          for (auto &r : rows) {
            key.insert(key.end(), r.begin(), r.end());
          }
          if (best_key.empty() || key < best_key) {
            best_key = key;
            CellClass cc;
            cc.lat.resize(n + 1, n);
            cc.cos.resize(n + 1);
            for (int k = 0; k <= n; k++) {
              for (int j = 0; j < n; j++) {
                cc.lat(k, j) = rows[k][j];
              }
              cc.cos(k) = rows[k][n];
            }
            best_val = cc;
          }
        }
        classes[best_key] = best_val;
      }
    }
    qh_freeqhull(qh, !qh_ALL);
    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);
    fclose(devnull);
    if (ok && !classes.empty()) {
      std::vector<CellClass> ret;
      ret.reserve(classes.size());
      for (auto &kv : classes) {
        ret.push_back(kv.second);
      }
      return ret;
    }
    R *= growth;
  }
  throw std::runtime_error("DelaunayCellClasses: ball never stabilized");
}

// ---------------------------------------------------------------------------
// Covering density
// ---------------------------------------------------------------------------

inline double VolumeUnitBall(int n) {
  return std::pow(M_PI, 0.5 * n) / std::tgamma(0.5 * n + 1.0);
}

inline MatrixXd CellPositions(CellClass const &cl, MatrixXd const &C) {
  int n = C.cols();
  MatrixXd P(n + 1, n);
  for (int k = 0; k <= n; k++) {
    P.row(k) = cl.lat.row(k).cast<double>() + C.row(cl.cos(k));
  }
  return P;
}

struct DensityResult {
  double theta;
  double mu2;
  std::vector<CellClass> cells;
};

inline DensityResult CoveringDensity(PeriodicConfig const &conf) {
  std::vector<CellClass> cells = DelaunayCellClasses(conf);
  double mu2 = 0;
  for (auto &cl : cells) {
    MatrixXd P = CellPositions(cl, conf.C);
    VectorXd center;
    double R2;
    CircumcenterRadius2(P, conf.Q, center, R2);
    mu2 = std::max(mu2, R2);
  }
  double det = conf.Q.determinant();
  double theta = conf.m * VolumeUnitBall(conf.n) *
                 std::pow(mu2, 0.5 * conf.n) / std::sqrt(det);
  return {theta, mu2, std::move(cells)};
}

// ---------------------------------------------------------------------------
// The Q-step: exact convex SDP at fixed cosets
// ---------------------------------------------------------------------------
//
// With the cosets and the cell list fixed, the vertices are fixed, so
// r(Delta)^2 <= 1 is an LMI linear in Q (Delone-Dolbilin-Ryshkov-Stogrin;
// see CoveringMaxdet.h): the block
//
//   BR_i(Q) = [ 4                diag(G_i)^T ]  >= 0,   G_i = V_i Q V_i^T,
//             [ diag(G_i)        G_i         ]
//
// with V_i the n x n matrix of edge vectors of simplex i. Minimizing
// -log det Q under these blocks is a determinant maximization problem,
// convex in Q; its optimum is the least dense form whose covering radius is
// at most 1 for this coset configuration. We solve it by a path-following
// barrier, mirroring the exact solver of CoveringMaxdet.h in double.
//
// Q is parametrized by its n(n+1)/2 independent entries in the basis A_u
// (E_kk on the diagonal, E_kl + E_lk off it), Q = sum_u x_u A_u.

// The symmetric basis: index pairs (k <= l).
inline std::vector<std::pair<int,int>> SymBasis(int n) {
  std::vector<std::pair<int,int>> b;
  for (int k = 0; k < n; k++) {
    for (int l = k; l < n; l++) {
      b.emplace_back(k, l);
    }
  }
  return b;
}

inline MatrixXd QFromX(int n, std::vector<std::pair<int,int>> const &basis,
                       VectorXd const &x) {
  MatrixXd Q = MatrixXd::Zero(n, n);
  for (size_t u = 0; u < basis.size(); u++) {
    int k = basis[u].first, l = basis[u].second;
    Q(k, l) += x(u);
    if (k != l) {
      Q(l, k) += x(u);
    }
  }
  return Q;
}

// The n x n edge matrices V_i (rows v_1 - v_0, ..., v_n - v_0) of the cells.
inline std::vector<MatrixXd> EdgeMatrices(std::vector<CellClass> const &cells,
                                          MatrixXd const &C) {
  int n = C.cols();
  std::vector<MatrixXd> Vs;
  Vs.reserve(cells.size());
  for (auto &cl : cells) {
    MatrixXd P = CellPositions(cl, C);
    MatrixXd V(n, n);
    for (int k = 0; k < n; k++) {
      V.row(k) = P.row(k + 1) - P.row(0);
    }
    Vs.push_back(std::move(V));
  }
  return Vs;
}

// Embed an n x n symmetric block M as the (n+1) x (n+1) circumradius block
// with border diag(M) and the given corner value.
inline MatrixXd EmbedBR(MatrixXd const &M, double corner) {
  int n = M.rows();
  MatrixXd F = MatrixXd::Zero(n + 1, n + 1);
  F(0, 0) = corner;
  for (int i = 0; i < n; i++) {
    F(0, i + 1) = M(i, i);
    F(i + 1, 0) = M(i, i);
    for (int j = 0; j < n; j++) {
      F(i + 1, j + 1) = M(i, j);
    }
  }
  return F;
}

struct QStepData {
  int n;
  int dim;
  std::vector<std::pair<int,int>> basis;
  std::vector<MatrixXd> Vs;           // edge matrices, per cell
  // per cell, per basis index u: the n x n matrix V_i A_u V_i^T
  std::vector<std::vector<MatrixXd>> Muc;
  double nu;                          // barrier parameter
};

inline QStepData BuildQStepData(std::vector<CellClass> const &cells,
                                MatrixXd const &C) {
  int n = C.cols();
  QStepData qd;
  qd.n = n;
  qd.basis = SymBasis(n);
  qd.dim = qd.basis.size();
  qd.Vs = EdgeMatrices(cells, C);
  int n_cell = qd.Vs.size();
  qd.Muc.resize(n_cell);
  for (int i = 0; i < n_cell; i++) {
    MatrixXd const &V = qd.Vs[i];
    qd.Muc[i].resize(qd.dim);
    for (int u = 0; u < qd.dim; u++) {
      int k = qd.basis[u].first, l = qd.basis[u].second;
      // V A_u V^T = V(:,k) V(:,l)^T + (k!=l) V(:,l) V(:,k)^T
      MatrixXd Vk = V.col(k);
      MatrixXd Vl = V.col(l);
      MatrixXd M = Vk * Vl.transpose();
      if (k != l) {
        M += Vl * Vk.transpose();
      }
      qd.Muc[i][u] = std::move(M);
    }
  }
  qd.nu = double(n) + double(n_cell) * double(n + 1);
  return qd;
}

// value / gradient / Hessian of the weighted barrier
//   phi_t(x) = -(1 + 1/t) log det Q - (1/t) sum_i log det BR_i(Q)
// returns false when x is not strictly feasible.
inline bool QStepBarrier(QStepData const &qd, VectorXd const &x, double t,
                         double &val, VectorXd *grad, MatrixXd *hess) {
  int n = qd.n, d = qd.dim;
  double wbar = 1.0 / t;
  MatrixXd Q = QFromX(n, qd.basis, x);
  Eigen::LLT<MatrixXd> lltQ(Q);
  if (lltQ.info() != Eigen::Success) {
    return false;
  }
  double logdetQ = 2.0 * lltQ.matrixLLT().diagonal().array().abs().log().sum();
  MatrixXd Qinv = lltQ.solve(MatrixXd::Identity(n, n));
  val = -(1.0 + wbar) * logdetQ;
  if (grad) {
    grad->setZero(d);
  }
  if (hess) {
    hess->setZero(d, d);
  }
  // objective term
  if (grad || hess) {
    std::vector<MatrixXd> WQ(d);
    for (int u = 0; u < d; u++) {
      int k = qd.basis[u].first, l = qd.basis[u].second;
      MatrixXd Au = MatrixXd::Zero(n, n);
      Au(k, l) = 1;
      if (k != l) {
        Au(l, k) = 1;
      }
      WQ[u] = Qinv * Au;
    }
    double w = 1.0 + wbar;
    for (int u = 0; u < d; u++) {
      if (grad) {
        (*grad)(u) -= w * WQ[u].trace();
      }
      if (hess) {
        for (int v = u; v < d; v++) {
          double val2 = w * (WQ[u].array() * WQ[v].transpose().array()).sum();
          (*hess)(u, v) += val2;
          if (v > u) {
            (*hess)(v, u) += val2;
          }
        }
      }
    }
  }
  // barrier blocks
  int n_cell = qd.Vs.size();
  for (int i = 0; i < n_cell; i++) {
    MatrixXd G = qd.Vs[i] * Q * qd.Vs[i].transpose();
    MatrixXd F = EmbedBR(G, 4.0);
    Eigen::LLT<MatrixXd> lltF(F);
    if (lltF.info() != Eigen::Success) {
      return false;
    }
    double logdetF = 2.0 * lltF.matrixLLT().diagonal().array().abs().log().sum();
    val -= wbar * logdetF;
    if (grad || hess) {
      MatrixXd Finv = lltF.solve(MatrixXd::Identity(n + 1, n + 1));
      std::vector<MatrixXd> WF(d);
      for (int u = 0; u < d; u++) {
        WF[u] = Finv * EmbedBR(qd.Muc[i][u], 0.0);
      }
      for (int u = 0; u < d; u++) {
        if (grad) {
          (*grad)(u) -= wbar * WF[u].trace();
        }
        if (hess) {
          for (int v = u; v < d; v++) {
            double val2 =
                wbar * (WF[u].array() * WF[v].transpose().array()).sum();
            (*hess)(u, v) += val2;
            if (v > u) {
              (*hess)(v, u) += val2;
            }
          }
        }
      }
    }
  }
  return true;
}

// Solve the Q-step SDP: returns the optimizer Q (normalized so that the
// covering radius squared is 1), from a strictly feasible start Q_init.
inline MatrixXd QStep(std::vector<CellClass> const &cells, MatrixXd const &C,
                      MatrixXd const &Q_init, int max_outer = 60,
                      double tol_gap = 1e-9, double mu = 15.0) {
  QStepData qd = BuildQStepData(cells, C);
  int n = qd.n, d = qd.dim;
  // strictly feasible start: scale Q so max circumradius^2 < 1
  double maxR2 = 0;
  for (auto &V : qd.Vs) {
    MatrixXd G = V * Q_init * V.transpose();
    VectorXd q = G.diagonal();
    double R2 = 0.25 * q.dot(G.ldlt().solve(q));
    maxR2 = std::max(maxR2, R2);
  }
  MatrixXd Q = Q_init / (2.0 * maxR2);
  VectorXd x(d);
  for (int u = 0; u < d; u++) {
    x(u) = Q(qd.basis[u].first, qd.basis[u].second);
  }
  double t = 1.0;
  VectorXd grad(d);
  MatrixXd hess(d, d);
  double val;
  for (int outer = 0; outer < max_outer; outer++) {
    if (qd.nu / t < tol_gap) {
      break;
    }
    for (int nit = 0; nit < 100; nit++) {
      if (!QStepBarrier(qd, x, t, val, &grad, &hess)) {
        break;
      }
      Eigen::LDLT<MatrixXd> ldlt(hess);
      VectorXd step = ldlt.solve(-grad);
      double dec2 = -grad.dot(step);
      if (!(dec2 > 0) || dec2 / 2 < 1e-13) {
        break;
      }
      double s = 1.0;
      bool ok = false;
      for (int bt = 0; bt < 60; bt++) {
        double vnew;
        if (QStepBarrier(qd, x + s * step, t, vnew, nullptr, nullptr) &&
            vnew <= val + 0.01 * s * grad.dot(step)) {
          x += s * step;
          ok = true;
          break;
        }
        s *= 0.5;
      }
      if (!ok) {
        break;
      }
    }
    t *= mu;
  }
  Q = QFromX(n, qd.basis, x);
  // renormalize so max circumradius^2 = 1
  double mR2 = 0;
  for (auto &V : qd.Vs) {
    MatrixXd G = V * Q * V.transpose();
    VectorXd q = G.diagonal();
    mR2 = std::max(mR2, 0.25 * q.dot(G.ldlt().solve(q)));
  }
  return Q / mR2;
}

// ---------------------------------------------------------------------------
// The smooth objective and its analytic gradient
// ---------------------------------------------------------------------------

/*
  Variables: the lower triangle of the Cholesky factor L of Q (n(n+1)/2
  entries, diagonal kept positive implicitly by the objective), then the
  cosets C_1, ..., C_{m-1} row by row.
 */
struct Packing {
  int n;
  int m;
  int dimL;
  int dim() const { return dimL + n * (m - 1); }
};

inline void Unpack(Packing const &pk, VectorXd const &x, MatrixXd &Lo,
                   MatrixXd &C) {
  int n = pk.n;
  Lo = MatrixXd::Zero(n, n);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      Lo(i, j) = x(pos++);
    }
  }
  C = MatrixXd::Zero(pk.m, n);
  for (int t = 1; t < pk.m; t++) {
    for (int j = 0; j < n; j++) {
      C(t, j) = x(pos++);
    }
  }
}

inline VectorXd Pack(Packing const &pk, MatrixXd const &Lo, MatrixXd const &C) {
  VectorXd x(pk.dim());
  int pos = 0;
  for (int i = 0; i < pk.n; i++) {
    for (int j = 0; j <= i; j++) {
      x(pos++) = Lo(i, j);
    }
  }
  for (int t = 1; t < pk.m; t++) {
    for (int j = 0; j < pk.n; j++) {
      x(pos++) = C(t, j);
    }
  }
  return x;
}

/*
  F(x) = (n/2) log smax_beta(R^2_S) - (1/2) log det Q and its gradient.
  The R^2 derivative in the Gram matrix G of the edge vectors is
     dR^2 = tr(W dG),   W = (2 diag(a) - a a^T)/4,   a = G^{-1} q,
  which chains to Q through G = V Q V^T, to the cosets through the vertex
  positions, and to L through Q = L L^T.
 */
inline double ObjectiveGradient(Packing const &pk,
                                std::vector<CellClass> const &cells,
                                double beta, VectorXd const &x,
                                VectorXd &grad) {
  int n = pk.n;
  MatrixXd Lo, C;
  Unpack(pk, x, Lo, C);
  MatrixXd Q = Lo * Lo.transpose();
  int n_cell = cells.size();
  std::vector<double> R2(n_cell);
  std::vector<MatrixXd> Wc(n_cell), Vc(n_cell);
  double maxR2 = -1;
  for (int s = 0; s < n_cell; s++) {
    MatrixXd P = CellPositions(cells[s], C);
    MatrixXd V(n, n);
    for (int k = 0; k < n; k++) {
      V.row(k) = P.row(k + 1) - P.row(0);
    }
    MatrixXd G = V * Q * V.transpose();
    VectorXd q = G.diagonal();
    VectorXd a = G.ldlt().solve(q);
    R2[s] = 0.25 * q.dot(a);
    MatrixXd W = -0.25 * (a * a.transpose());
    W.diagonal() += 0.5 * a;
    Wc[s] = std::move(W);
    Vc[s] = std::move(V);
    maxR2 = std::max(maxR2, R2[s]);
  }
  // soft maximum and its weights
  double Z = 0;
  std::vector<double> w(n_cell);
  for (int s = 0; s < n_cell; s++) {
    w[s] = std::exp(beta * (R2[s] - maxR2));
    Z += w[s];
  }
  double smax = maxR2 + std::log(Z) / beta;
  for (int s = 0; s < n_cell; s++) {
    w[s] /= Z;
  }
  double logdet = 0;
  for (int i = 0; i < n; i++) {
    logdet += 2 * std::log(std::abs(Lo(i, i)) + 1e-300);
  }
  double F = 0.5 * n * std::log(smax) - 0.5 * logdet;
  // gradient assembly
  double c_out = 0.5 * n / smax;
  MatrixXd dFdQ = MatrixXd::Zero(n, n);
  MatrixXd dFdC = MatrixXd::Zero(pk.m, n);
  for (int s = 0; s < n_cell; s++) {
    double cw = c_out * w[s];
    if (cw < 1e-16) {
      continue;
    }
    MatrixXd const &V = Vc[s];
    MatrixXd const &W = Wc[s];
    dFdQ += cw * (V.transpose() * W * V);
    // dR^2/dV = 2 W V Q ; vertex j>=1 gets row j-1, vertex 0 minus the sum
    MatrixXd gV = 2.0 * (W * V * Q);
    VectorXd rowsum = gV.colwise().sum();
    for (int k = 0; k <= n; k++) {
      int t = cells[s].cos(k);
      if (t == 0) {
        continue;
      }
      if (k == 0) {
        dFdC.row(t) -= cw * rowsum.transpose();
      } else {
        dFdC.row(t) += cw * gV.row(k - 1);
      }
    }
  }
  dFdQ = 0.5 * (dFdQ + dFdQ.transpose());
  // chain to L: dQ = dL L^T + L dL^T  ->  grad_L = 2 dFdQ L (lower triangle)
  MatrixXd gL = 2.0 * dFdQ * Lo;
  for (int i = 0; i < n; i++) {
    gL(i, i) -= 1.0 / Lo(i, i);
  }
  grad.resize(pk.dim());
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      grad(pos++) = gL(i, j);
    }
  }
  for (int t = 1; t < pk.m; t++) {
    for (int j = 0; j < n; j++) {
      grad(pos++) = dFdC(t, j);
    }
  }
  return F;
}

// ---------------------------------------------------------------------------
// L-BFGS with Armijo backtracking
// ---------------------------------------------------------------------------

template <typename Fobj>
void LBFGS(Fobj f, VectorXd &x, int max_iter, double gtol,
           [[maybe_unused]] std::ostream &os) {
  int dim = x.size();
  int mem = 10;
  std::vector<VectorXd> s_hist, y_hist;
  std::vector<double> rho;
  VectorXd g(dim), g_new(dim);
  double F = f(x, g);
  for (int it = 0; it < max_iter; it++) {
    if (g.norm() < gtol) {
      break;
    }
    // two-loop recursion
    VectorXd d = -g;
    int k = s_hist.size();
    std::vector<double> alpha(k);
    for (int i = k - 1; i >= 0; i--) {
      alpha[i] = rho[i] * s_hist[i].dot(d);
      d -= alpha[i] * y_hist[i];
    }
    if (k > 0) {
      double gamma = s_hist[k - 1].dot(y_hist[k - 1]) /
                     y_hist[k - 1].squaredNorm();
      d *= gamma;
    }
    for (int i = 0; i < k; i++) {
      double bcoef = rho[i] * y_hist[i].dot(d);
      d += (alpha[i] - bcoef) * s_hist[i];
    }
    double gd = g.dot(d);
    if (gd > 0) {
      d = -g;
      gd = -g.squaredNorm();
    }
    double step = 1.0;
    double F_new = F;
    bool moved = false;
    for (int bt = 0; bt < 50; bt++) {
      VectorXd x_new = x + step * d;
      F_new = f(x_new, g_new);
      if (std::isfinite(F_new) && F_new <= F + 1e-4 * step * gd) {
        VectorXd s_vec = x_new - x;
        VectorXd y_vec = g_new - g;
        double sy = s_vec.dot(y_vec);
        if (sy > 1e-12) {
          s_hist.push_back(s_vec);
          y_hist.push_back(y_vec);
          rho.push_back(1.0 / sy);
          if (int(s_hist.size()) > mem) {
            s_hist.erase(s_hist.begin());
            y_hist.erase(y_hist.begin());
            rho.erase(rho.begin());
          }
        }
        x = x_new;
        g = g_new;
        F = F_new;
        moved = true;
        break;
      }
      step *= 0.5;
    }
    if (!moved) {
      break;
    }
  }
}

// ---------------------------------------------------------------------------
// The c-step: trust-region SLP on the coset minimax
// ---------------------------------------------------------------------------
//
// With Q and the cell list fixed, minimize max_i r(Delta_i)^2 over the
// cosets C by sequential linear programming. At the current C each active
// simplex contributes r_i^2 and its gradient dr_i^2/dC (analytic, the coset
// chain of ObjectiveGradient); the step solves the tiny LP
//
//   min_{dc, u}  u   s.t.  r_i^2 + g_i . dc <= u,   |dc|_inf <= rho,
//
// over the simplices within a band of the current maximum, then C += dc,
// with a trust radius rho grown on success and shrunk on failure. This is
// the "inequality matching" of the minimax: the binding simplices are the
// active LP rows. The LP is solved by its dual with a projected subgradient,
// which suffices at this size (5(m-1)+1 variables) and needs no external
// solver.

// r_i^2 and its gradient in the free cosets C_1..C_{m-1}, packed row-major.
inline void CellR2Grad(CellClass const &cl, MatrixXd const &Q, MatrixXd const &C,
                       double &R2, VectorXd &g) {
  int n = Q.rows();
  int m = C.rows();
  MatrixXd P = CellPositions(cl, C);
  MatrixXd V(n, n);
  for (int k = 0; k < n; k++) {
    V.row(k) = P.row(k + 1) - P.row(0);
  }
  MatrixXd G = V * Q * V.transpose();
  VectorXd q = G.diagonal();
  Eigen::LDLT<MatrixXd> ldlt(G);
  VectorXd a = ldlt.solve(q);
  R2 = 0.25 * q.dot(a);
  // dR2/dV = 2 W V Q, W = (2 diag(a) - a a^T)/4
  MatrixXd W = -0.25 * (a * a.transpose());
  W.diagonal() += 0.5 * a;
  MatrixXd gV = 2.0 * (W * V * Q);      // rows = d/d(edge_k)
  VectorXd rowsum = gV.colwise().sum();
  g = VectorXd::Zero(n * (m - 1));
  for (int k = 0; k <= n; k++) {
    int t = cl.cos(k);
    if (t == 0) {
      continue;
    }
    int base = (t - 1) * n;
    if (k == 0) {
      for (int j = 0; j < n; j++) {
        g(base + j) -= rowsum(j);
      }
    } else {
      for (int j = 0; j < n; j++) {
        g(base + j) += gV(k - 1, j);
      }
    }
  }
}

// One c-step: minimize max_i r_i^2 over the cosets at fixed Q, by steepest
// descent for the maximum. The descent direction is the negative of the
// min-norm element of the convex hull of the active gradients (the correct
// steepest-descent direction for a minimax, found by Frank-Wolfe), and the
// step is chosen by Armijo backtracking on the true maximum. Returns the
// new maximum r^2.
inline double CStepMinimax(std::vector<CellClass> const &cells,
                           MatrixXd const &Q, MatrixXd &C, int iters = 20) {
  int n = Q.cols();
  int m = C.rows();
  int dc_dim = n * (m - 1);
  int n_cell = cells.size();
  double cur_max = -1;
  for (int it = 0; it < iters; it++) {
    std::vector<double> R2(n_cell);
    std::vector<VectorXd> g(n_cell);
    double mx = -1;
    for (int i = 0; i < n_cell; i++) {
      CellR2Grad(cells[i], Q, C, R2[i], g[i]);
      mx = std::max(mx, R2[i]);
    }
    cur_max = mx;
    std::vector<int> act;
    for (int i = 0; i < n_cell; i++) {
      if (R2[i] > mx - 0.02 * mx - 1e-12) {
        act.push_back(i);
      }
    }
    // min-norm element of conv{ g_i : i active } by Frank-Wolfe
    VectorXd d = g[act[0]];
    for (int k = 0; k < 80; k++) {
      int jmin = act[0];
      double best = d.dot(g[act[0]]);
      for (int idx : act) {
        double val = d.dot(g[idx]);
        if (val < best) {
          best = val;
          jmin = idx;
        }
      }
      double gamma = 2.0 / (k + 2.0);
      d = (1.0 - gamma) * d + gamma * g[jmin];
    }
    double dn = d.norm();
    if (dn < 1e-10) {
      break;   // stationary: 0 is in the hull of active gradients
    }
    VectorXd dir = -d / dn;
    // Armijo line search on the true maximum
    double step = 0.25;
    bool moved = false;
    for (int bt = 0; bt < 40; bt++) {
      MatrixXd Ctrial = C;
      for (int t = 1; t < m; t++) {
        for (int j = 0; j < n; j++) {
          Ctrial(t, j) += step * dir((t - 1) * n + j);
        }
      }
      double mx_new = -1;
      for (int i = 0; i < n_cell; i++) {
        MatrixXd P = CellPositions(cells[i], Ctrial);
        MatrixXd V(n, n);
        for (int k = 0; k < n; k++) {
          V.row(k) = P.row(k + 1) - P.row(0);
        }
        MatrixXd G = V * Q * V.transpose();
        VectorXd q = G.diagonal();
        mx_new = std::max(mx_new, 0.25 * q.dot(G.ldlt().solve(q)));
      }
      if (mx_new < mx - 1e-4 * step * dn) {
        C = Ctrial;
        cur_max = mx_new;
        moved = true;
        break;
      }
      step *= 0.5;
    }
    if (!moved) {
      break;
    }
    (void)dc_dim;
  }
  return cur_max;
}

// ---------------------------------------------------------------------------
// The descent
// ---------------------------------------------------------------------------

struct DescendResult {
  bool success = false;
  double theta = 0;
  PeriodicConfig conf;
  int n_retessellations = 0;
};

// Alternating block descent: exact SDP Q-step, trust-region SLP c-step,
// re-tessellating between outer rounds. This is the fast successor of the
// soft-max L-BFGS Descend below; it exploits the convexity of the problem
// in Q and treats the coset minimax by its active set.
inline DescendResult DescendAlt(PeriodicConfig conf, int rounds,
                                std::ostream &os, bool verbose = false,
                                int deadline_sec = 120) {
  DescendResult best;
  int stall = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int it = 0; it < rounds; it++) {
    if (std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - t0).count() > deadline_sec) {
      return best;
    }
    // one tessellation per round: it both scores the current config and
    // supplies the cell list for this round's Q-step and c-step. The next
    // round's tessellation verifies the move just made.
    DensityResult dr;
    try {
      dr = CoveringDensity(conf);
    } catch (std::runtime_error &e) {
      return best;
    }
    if (verbose) {
      os << "  round " << it << ": theta = " << dr.theta << "\n";
    }
    if (!best.success || dr.theta < best.theta - 1e-11) {
      best.success = true;
      best.theta = dr.theta;
      best.conf = conf;
      stall = 0;
    } else {
      stall++;
      if (stall >= 3) {
        return best;
      }
    }
    best.n_retessellations = it + 1;
    // Q-step: exact SDP optimum for the current cosets and cell list
    MatrixXd Qnew;
    try {
      Qnew = QStep(dr.cells, conf.C, conf.Q);
    } catch (std::runtime_error &e) {
      return best;
    }
    // guard: a bad cell list can give an anisotropic SDP optimum whose
    // tessellation would explode. If the new form is too skewed, keep the
    // best so far and stop rather than crawl.
    {
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(Qnew, Eigen::EigenvaluesOnly);
      double cond = es.eigenvalues()(conf.n - 1) / es.eigenvalues()(0);
      if (!(cond < 300.0)) {
        return best;
      }
    }
    conf.Q = Qnew;
    // c-step: minimize the covering radius over the cosets at this Q
    CStepMinimax(dr.cells, conf.Q, conf.C);
    conf.C.row(0).setZero();
  }
  // final scoring after the last move
  try {
    DensityResult dr = CoveringDensity(conf);
    if (dr.theta < best.theta - 1e-11) {
      best.theta = dr.theta;
      best.conf = conf;
    }
  } catch (std::runtime_error &e) {
  }
  return best;
}

inline DescendResult Descend(PeriodicConfig conf, int rounds, std::ostream &os,
                             bool verbose = false) {
  Packing pk{conf.n, conf.m, conf.n * (conf.n + 1) / 2};
  DescendResult best;
  int stall = 0;
  for (int it = 0; it < rounds; it++) {
    DensityResult dr;
    try {
      dr = CoveringDensity(conf);
    } catch (std::runtime_error &e) {
      return best;
    }
    // the input itself is the first incumbent: a descent must never
    // return something worse than what it was given
    if (!best.success || dr.theta < best.theta - 1e-11) {
      best.success = true;
      best.theta = dr.theta;
      best.conf = conf;
    }
    // normalize the scale so that mu^2 = 1
    conf.Q /= dr.mu2;
    Eigen::LLT<MatrixXd> llt(conf.Q);
    if (llt.info() != Eigen::Success) {
      return best;
    }
    MatrixXd Lo = llt.matrixL();
    VectorXd x = Pack(pk, Lo, conf.C);
    double beta = 200.0;
    for (int stage = 0; stage < 6; stage++) {
      auto f = [&](VectorXd const &xv, VectorXd &g) -> double {
        return ObjectiveGradient(pk, dr.cells, beta, xv, g);
      };
      LBFGS(f, x, 150, 1e-10, os);
      beta *= 5.0;
    }
    MatrixXd C_new;
    Unpack(pk, x, Lo, C_new);
    conf.Q = Lo * Lo.transpose();
    conf.C = C_new;
    conf.C.row(0).setZero();
    DensityResult ver;
    try {
      ver = CoveringDensity(conf);
    } catch (std::runtime_error &e) {
      return best;
    }
    best.n_retessellations = it + 1;
    if (verbose) {
      os << "  round " << it << ": verified theta = " << ver.theta << "\n";
    }
    if (!best.success || ver.theta < best.theta - 1e-11) {
      best.success = true;
      best.theta = ver.theta;
      best.conf = conf;
      stall = 0;
    } else {
      stall++;
      if (stall >= 2) {
        return best;
      }
    }
  }
  return best;
}

// ---------------------------------------------------------------------------
// Configuration file I/O
// ---------------------------------------------------------------------------

inline PeriodicConfig ReadConfigFile(std::string const &FileName) {
  FILE *f = fopen(FileName.c_str(), "r");
  if (!f) {
    throw std::runtime_error("cannot open " + FileName);
  }
  PeriodicConfig conf;
  if (fscanf(f, "%d %d", &conf.n, &conf.m) != 2) {
    throw std::runtime_error("bad header in " + FileName);
  }
  conf.Q.resize(conf.n, conf.n);
  for (int i = 0; i < conf.n; i++) {
    for (int j = 0; j < conf.n; j++) {
      if (fscanf(f, "%lf", &conf.Q(i, j)) != 1) {
        throw std::runtime_error("bad Q in " + FileName);
      }
    }
  }
  conf.C.resize(conf.m, conf.n);
  for (int t = 0; t < conf.m; t++) {
    for (int j = 0; j < conf.n; j++) {
      if (fscanf(f, "%lf", &conf.C(t, j)) != 1) {
        throw std::runtime_error("bad C in " + FileName);
      }
    }
  }
  fclose(f);
  return conf;
}

inline void WriteConfigFile(std::string const &FileName,
                            PeriodicConfig const &conf) {
  FILE *f = fopen(FileName.c_str(), "w");
  fprintf(f, "%d %d\n", conf.n, conf.m);
  for (int i = 0; i < conf.n; i++) {
    for (int j = 0; j < conf.n; j++) {
      fprintf(f, " %.17g", conf.Q(i, j));
    }
    fprintf(f, "\n");
  }
  for (int t = 0; t < conf.m; t++) {
    for (int j = 0; j < conf.n; j++) {
      fprintf(f, " %.17g", conf.C(t, j));
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

// clang-format off
}  // namespace joint_covering_double
// clang-format on

// clang-format off
#endif  // SRC_DELAUNAY_JOINTCOVERINGDOUBLE_H_
// clang-format on
