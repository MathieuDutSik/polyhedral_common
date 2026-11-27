// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_LOCALDELAUNAYCONFIG_H_
#define SRC_DELAUNAY_LOCALDELAUNAYCONFIG_H_

// clang-format off
#include "CoveringMaxdet.h"
#include "SignatureSymmetric.h"
#include <cmath>
#include <map>
#include <functional>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  Local configurations of Delaunay simplices for the lattice covering problem,
  after A. Garber, "The Schuermann-Vallentin lattice L^c_6 gives the thinnest
  lattice sphere covering in dimension six", arXiv:2611.28562.

  Let Q be a form with covering radius R(Q) = 1 and covering density at most
  Theta_rec. Then det Q >= dstar = kappa_n^2 / Theta_rec^2, every Delaunay
  simplex has circumradius at most 1, and a configuration L of simplices of
  the Delaunay triangulation of Q satisfies

      det Q <= max { det Q' : BR_S(Q') >= 0 for S in L },

  the maximum being a determinant maximization (CoveringMaxdet.h). A
  configuration whose maximum is below dstar cannot occur.

  --- The neighbours of a simplex ---

  Let S = conv(v_0, ..., v_n) be a Delaunay simplex of normalized volume V.
  Its inertia matrix has smallest eigenvalue lambda_S with

      lambda (n + 1 - lambda)^(n-1) / (n-1)^(n-1) >= V^2 dstar / (n + 1)

  (Garber, Lemma 5.6), and a lattice vector x = sum_i alpha_i v_i with
  sum_i alpha_i = 0 has Q[x] >= lambda_S sum_i alpha_i^2 (Lemma 5.7). An edge
  has Q[x] <= 4, so its coefficients satisfy sum_i alpha_i^2 <= 4 / lambda_S.

  The neighbour of S across the facet F_k opposite to v_k is F_k + w for a
  lattice point w whose barycentric coordinates beta satisfy beta_k < 0. The
  edges w - v_j, j != k, have coefficients beta - e_j, which bounds
  sum_i beta_i^2 by (1 + sqrt(4 / lambda_S))^2, a finite enumeration. The
  candidates are then filtered by
  * the volume V |beta_k| of the neighbour, at most a given Vmax;
  * the edge bound for each w - v_j;
  * no edge in 2 Z^n (its midpoint would be a lattice point);
  * the parity frame: two edges in the same class of Z^n / 2 Z^n are equal up
    to sign (Garber, Lemma 5.2).

  --- Emptiness ---

  Optionally, the linear constraint that w is outside the circumsphere of S,
  f_{S,w}(Q) = Q[w - v_0] - sum_i beta_i Q[v_i - v_0] >= 0, is added, and at
  the optimum any lattice point found strictly inside the circumsphere of S
  or of the neighbour gives a further valid linear constraint. Every
  constraint added is valid, so the bound stays an upper bound whether or not
  the cuts are found. When no strictly feasible starting point is found for
  the extra constraints, they are dropped, which also keeps the bound valid.

  The determinant maximization is numerical. The bound reported is
  det * exp(gap_bound), the bound of the solver on the optimum; it is not a
  rational certificate.
 */

#ifdef DEBUG
#define DEBUG_LOCAL_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LOCAL_DELAUNAY
#endif

#ifdef TIMINGS
#define TIMINGS_LOCAL_DELAUNAY
#endif

namespace local_delaunay {

inline double VolumeUnitBall(int n) {
  return std::pow(M_PI, double(n) / 2.0) / std::tgamma(double(n) / 2.0 + 1.0);
}

// The threshold kappa_n^2 / Theta_rec^2 on det Q.
inline double DeterminantThreshold(int n, double ThetaRec) {
  double kappa = VolumeUnitBall(n);
  return (kappa / ThetaRec) * (kappa / ThetaRec);
}

/*
  The lower bound on the smallest eigenvalue of the inertia matrix of a
  simplex of normalized volume V, or nothing if no such simplex is possible.
  The function lambda -> lambda ((n + 1 - lambda) / (n - 1))^(n - 1) is
  increasing on [0, (n + 1) / n].
 */
inline std::optional<double> InertiaLowerBound(int n, int V, double dstar) {
  double rhs = double(V) * double(V) * dstar / double(n + 1);
  auto f = [&](double lambda) -> double {
    return lambda * std::pow((double(n + 1) - lambda) / double(n - 1), n - 1);
  };
  double hi = double(n + 1) / double(n);
  if (f(hi) < rhs) {
    return {};
  }
  double lo = 0;
  for (int iter = 0; iter < 200; iter++) {
    double mid = (lo + hi) / 2;
    if (f(mid) < rhs) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  // lo is below the root: a smaller lambda gives a larger, still valid,
  // edge bound.
  return lo;
}

template <typename T> struct SimplexData {
  int n;
  // The vertices, one per row.
  MyMatrix<T> Vert;
  // The rows v_i - v_0, i = 1, ..., n, and the inverse of that matrix.
  MyMatrix<T> M;
  MyMatrix<T> Minv;
  int V;
};

template <typename T> SimplexData<T> GetSimplexData(MyMatrix<T> const &Vert) {
  int n = Vert.cols();
  MyMatrix<T> M(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M(i, j) = Vert(i + 1, j) - Vert(0, j);
    }
  }
  T det = DeterminantMat(M);
  T absdet = T_abs(det);
  int V = UniversalScalarConversion<int, T>(absdet);
  return {n, Vert, M, Inverse(M), V};
}

// The barycentric coordinates (beta_0, ..., beta_n) of the point w.
template <typename T>
MyVector<T> Barycentric(SimplexData<T> const &sd, MyVector<T> const &w) {
  int n = sd.n;
  MyVector<T> z(n);
  for (int i = 0; i < n; i++) {
    z(i) = w(i) - sd.Vert(0, i);
  }
  MyVector<T> betap = sd.Minv.transpose() * z;
  MyVector<T> beta(n + 1);
  T sum(0);
  for (int i = 0; i < n; i++) {
    beta(i + 1) = betap(i);
    sum += betap(i);
  }
  beta(0) = T(1) - sum;
  return beta;
}

template <typename T> MyVector<T> GetRow(MyMatrix<T> const &M, int i) {
  return GetMatrixRow(M, i);
}

/*
  The parity frame of a family of edges: the class of each edge in
  Z^n / 2 Z^n and the edge itself, up to sign.
 */
template <typename T> struct ParityFrame {
  std::map<std::vector<int>, MyVector<T>> map;

  static std::vector<int> Class(MyVector<T> const &e) {
    std::vector<int> c(e.size());
    for (int i = 0; i < e.size(); i++) {
      // e is integral with small coordinates: its class mod 2.
      int v = UniversalScalarConversion<int, T>(e(i));
      c[i] = ((v % 2) + 2) % 2;
    }
    return c;
  }
  static bool IsZeroClass(std::vector<int> const &c) {
    for (int x : c) {
      if (x != 0) {
        return false;
      }
    }
    return true;
  }
  // Inserts e; returns false if e is even or conflicts with the frame.
  bool Insert(MyVector<T> const &e) {
    std::vector<int> c = Class(e);
    if (IsZeroClass(c)) {
      return false;
    }
    auto iter = map.find(c);
    if (iter == map.end()) {
      map[c] = e;
      return true;
    }
    MyVector<T> const &f = iter->second;
    return f == e || f == MyVector<T>(-e);
  }
};

template <typename T> ParityFrame<T> GetSimplexFrame(MyMatrix<T> const &Vert) {
  ParityFrame<T> pf;
  int nbVert = Vert.rows();
  for (int i = 0; i < nbVert; i++) {
    for (int j = i + 1; j < nbVert; j++) {
      MyVector<T> e = GetRow(Vert, j) - GetRow(Vert, i);
      (void)pf.Insert(e);
    }
  }
  return pf;
}

template <typename T> struct NeighborCandidate {
  int k;
  MyVector<T> w;
  MyVector<T> beta;
  // The normalized volume of the neighbour.
  int V;
};

/*
  The classes of Z^n modulo the lattice spanned by the edges v_i - v_0 of S,
  as the vectors a in [0, V)^n with V beta' = a mod V, beta' the barycentric
  coordinates beta_1, ..., beta_n. V Minv is integral, so a breadth first
  search by the unit vectors reaches the V classes.
 */
template <typename T>
std::vector<std::vector<int64_t>>
CosetRepresentatives(SimplexData<T> const &sd) {
  int n = sd.n;
  int64_t V = sd.V;
  std::vector<std::vector<int64_t>> ListStep(n, std::vector<int64_t>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T val = T(sd.V) * sd.Minv(i, j);
      int64_t v = UniversalScalarConversion<int64_t, T>(val);
      ListStep[i][j] = ((v % V) + V) % V;
    }
  }
  std::vector<std::vector<int64_t>> ListRep{std::vector<int64_t>(n, 0)};
  std::set<std::vector<int64_t>> SetRep(ListRep.begin(), ListRep.end());
  size_t pos = 0;
  while (pos < ListRep.size()) {
    std::vector<int64_t> a = ListRep[pos];
    pos++;
    for (int i = 0; i < n; i++) {
      std::vector<int64_t> b(n);
      for (int j = 0; j < n; j++) {
        b[j] = (a[j] + ListStep[i][j]) % V;
      }
      if (SetRep.insert(b).second) {
        ListRep.push_back(b);
      }
    }
  }
#ifdef SANITY_CHECK_LOCAL_DELAUNAY
  if (static_cast<int64_t>(ListRep.size()) != V) {
    std::cerr << "LOCAL_DELAUNAY: found " << ListRep.size()
              << " classes instead of V=" << V << "\n";
    throw TerminalException{1};
  }
#endif
  return ListRep;
}

/*
  The candidate apexes of the neighbours of S across the facet opposite to
  v_k, for each k in ListFacet.

  With b = V beta (integers, sum V), the edge w - v_j has coefficients
  beta - e_j and the edge bound reads, for every j != k,
      sum_i b_i^2 <= V^2 (EdgeBound - 1) + 2 V b_j,
  and the volume of the neighbour is -b_k, in [1, Vmax]. The recursion
  fixes b_k, then assigns b_1, ..., b_n (b_0 being determined by the sum)
  in the class of a coset representative, pruning with the partial sums of
  squares against the bound for the smallest b_j assigned so far.
 */
template <typename T, typename Tint>
std::vector<NeighborCandidate<T>>
EnumerateNeighborCandidates(SimplexData<T> const &sd,
                            std::vector<int> const &ListFacet, double EdgeBound,
                            int Vmax, [[maybe_unused]] std::ostream &os) {
  int n = sd.n;
  int64_t V = sd.V;
  double VV = double(V);
  double base = VV * VV * (EdgeBound - 1) + 1e-6;
  std::vector<std::vector<int64_t>> ListRep = CosetRepresentatives(sd);
  ParityFrame<T> frame = GetSimplexFrame(sd.Vert);
  std::vector<NeighborCandidate<T>> ListCand;
  // The smallest b_j, j != k, is at most the mean (V - b_k) / n.
  double mean_bound = (VV + double(Vmax)) / double(n);
  double cap = base + 2 * VV * mean_bound;
  int64_t bmax = static_cast<int64_t>(std::floor(std::sqrt(cap))) + 1;
  std::vector<int64_t> b(n + 1);
#ifdef DEBUG_LOCAL_DELAUNAY
  size_t n_leaf = 0;
#endif
  auto f_leaf = [&](int k) -> void {
#ifdef DEBUG_LOCAL_DELAUNAY
    n_leaf++;
#endif
    int64_t sum = 0;
    for (int i = 1; i <= n; i++) {
      sum += b[i];
    }
    b[0] = V - sum;
    if (k == 0 && (b[0] >= 0 || -b[0] > Vmax)) {
      return;
    }
    double sumsq = 0;
    for (int i = 0; i <= n; i++) {
      sumsq += double(b[i]) * double(b[i]);
    }
    for (int j = 0; j <= n; j++) {
      if (j != k && sumsq > base + 2 * VV * double(b[j])) {
        return;
      }
    }
    MyVector<T> beta(n + 1);
    for (int i = 0; i <= n; i++) {
      beta(i) = T(static_cast<long>(b[i])) / T(static_cast<long>(V));
    }
    MyVector<T> w(n);
    for (int j = 0; j < n; j++) {
      T val = sd.Vert(0, j);
      for (int i = 1; i <= n; i++) {
        val += beta(i) * sd.M(i - 1, j);
      }
      w(j) = val;
    }
    ParityFrame<T> frame_k = frame;
    for (int j = 0; j <= n; j++) {
      if (j == k) {
        continue;
      }
      MyVector<T> e = w - GetRow(sd.Vert, j);
      if (!frame_k.Insert(e)) {
        return;
      }
    }
    int vol = static_cast<int>(-b[k]);
    ListCand.push_back({k, w, beta, vol});
  };
  for (int k : ListFacet) {
    for (auto &a : ListRep) {
      // The free indices: 1..n except k (b_k is fixed when k >= 1).
      std::vector<int> ListFree;
      for (int i = 1; i <= n; i++) {
        if (i != k) {
          ListFree.push_back(i);
        }
      }
      std::function<void(size_t, double, double)> f_rec =
          [&](size_t pos, double partial, double min_b) -> void {
        double cap_now = base + 2 * VV * std::min(min_b, mean_bound);
        if (partial > cap_now) {
          return;
        }
        if (pos == ListFree.size()) {
          f_leaf(k);
          return;
        }
        int i = ListFree[pos];
        // b_i = a_{i-1} mod V, |b_i| <= bmax.
        int64_t r = a[i - 1];
        // The smallest v >= -bmax with v = r mod V.
        int64_t v_first = -bmax + (((r + bmax) % V) + V) % V;
        for (int64_t v = v_first; v <= bmax; v += V) {
          double nv = partial + double(v) * double(v);
          if (nv > cap_now) {
            continue;
          }
          b[i] = v;
          f_rec(pos + 1, nv, std::min(min_b, double(v)));
        }
      };
      if (k == 0) {
        f_rec(0, 0.0, 1e100);
      } else {
        int64_t r = a[k - 1];
        for (int64_t vk = -1; vk >= -Vmax; vk--) {
          if ((((vk - r) % V) + V) % V != 0) {
            continue;
          }
          b[k] = vk;
          f_rec(0, double(vk) * double(vk), 1e100);
        }
      }
    }
  }
#ifdef DEBUG_LOCAL_DELAUNAY
  os << "LOCAL_DELAUNAY: n_leaf=" << n_leaf << " |ListCand|=" << ListCand.size()
     << "\n";
#endif
  return ListCand;
}

// The basis E_ii, E_ij + E_ji of the symmetric matrices.
template <typename T> std::vector<MyMatrix<T>> GetSymmetricBasis(int n) {
  std::vector<MyMatrix<T>> ListA;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      MyMatrix<T> A = ZeroMatrix<T>(n, n);
      A(i, j) = T(1);
      A(j, i) = T(1);
      ListA.push_back(A);
    }
  }
  return ListA;
}

template <typename T>
MyVector<T> GramToSymmetricCoordinates(MyMatrix<T> const &Q) {
  int n = Q.rows();
  MyVector<T> x(n * (n + 1) / 2);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      x(pos) = (i == j) ? Q(i, i) : Q(i, j);
      pos++;
    }
  }
  return x;
}

// The rows v_i - v_0 of a simplex given by its vertices.
template <typename T> MyMatrix<T> DifferenceMatrix(MyMatrix<T> const &Vert) {
  int n = Vert.cols();
  MyMatrix<T> M(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M(i, j) = Vert(i + 1, j) - Vert(0, j);
    }
  }
  return M;
}

/*
  The linear form f_{S,w}(Q) >= 0 in the symmetric coordinates: w outside
  the circumsphere of S.
 */
template <typename T>
MyVector<T> EmptinessInequality(MyMatrix<T> const &Vert, MyVector<T> const &w) {
  SimplexData<T> sd = GetSimplexData(Vert);
  int n = sd.n;
  MyVector<T> beta = Barycentric(sd, w);
  MyMatrix<T> F = ZeroMatrix<T>(n, n);
  MyVector<T> y(n);
  for (int i = 0; i < n; i++) {
    y(i) = w(i) - Vert(0, i);
  }
  F += y * y.transpose();
  for (int i = 1; i <= n; i++) {
    MyVector<T> d = GetRow(sd.M, i - 1);
    F -= beta(i) * (d * d.transpose());
  }
  // f(Q) = trace(F Q) = sum_{i<=j} x_ij coef_ij with x_ii = Q_ii, x_ij = Q_ij.
  MyVector<T> a(n * (n + 1) / 2);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      a(pos) = (i == j) ? F(i, i) : T(2) * F(i, j);
      pos++;
    }
  }
  return a;
}

/*
  A rank one positive semidefinite form D with f_{S,w}(D) > 0, for w outside
  S: with beta_k < 0, the dual vector u of the vertex v_k with respect to
  another vertex gives f_{S,w}(u u^T) = beta_k^2 - beta_k > 0.
 */
template <typename T>
std::optional<MyMatrix<T>> PositiveDirection(MyMatrix<T> const &Vert,
                                             MyVector<T> const &w) {
  int n = Vert.cols();
  SimplexData<T> sd0 = GetSimplexData(Vert);
  MyVector<T> beta = Barycentric(sd0, w);
  int k = -1;
  for (int i = 0; i <= n; i++) {
    if (beta(i) < 0) {
      k = i;
      break;
    }
  }
  if (k == -1) {
    return {};
  }
  // Base vertex b != k, the vertices reordered so that b comes first.
  int b = (k == 0) ? 1 : 0;
  MyMatrix<T> Vert2 = Vert;
  for (int j = 0; j < n; j++) {
    Vert2(0, j) = Vert(b, j);
    Vert2(b, j) = Vert(0, j);
  }
  // In the reordered simplex the index of v_k is k (k != 0, b) or b (k == 0).
  int k2 = (k == 0) ? b : k;
  SimplexData<T> sd = GetSimplexData(Vert2);
  // u is column k2 - 1 of Minv^T, i.e. (v_i - v_b) . u = delta_{i,k2}.
  MyVector<T> u(n);
  for (int j = 0; j < n; j++) {
    u(j) = sd.Minv(j, k2 - 1);
  }
  return MyMatrix<T>(u * u.transpose());
}

struct PairBound {
  double det;
  double bound;
  bool success;
  int n_cut;
};

/*
  The determinant maximization over the simplices ListVert, with the linear
  constraints ListIneq. Returns the bound det * exp(gap) of the solver.
 */
template <typename T>
std::optional<std::pair<double, double>>
SolveConfiguration(std::vector<MyMatrix<T>> const &ListVert,
                   std::vector<MyVector<T>> const &ListIneq,
                   std::vector<MyMatrix<T>> const &ListDir,
                   MyVector<double> &x_opt, std::ostream &os) {
  using namespace covering_maxdet;
  int n = ListVert[0].cols();
  int dim = n * (n + 1) / 2;
  CoveringData<T> cd;
  cd.n = n;
  cd.dim = dim;
  cd.ListA = GetSymmetricBasis<T>(n);
  cd.Alin = MyMatrix<T>(ListIneq.size(), dim);
  for (size_t i = 0; i < ListIneq.size(); i++) {
    for (int u = 0; u < dim; u++) {
      cd.Alin(i, u) = ListIneq[i](u);
    }
  }
  for (auto &Vert : ListVert) {
    cd.ListSimplex.push_back(DifferenceMatrix(Vert));
  }
  cd.point_density = T(1);
  // The starting point: the identity plus the positive directions, until
  // every linear constraint is strict, then scaled to circumradii <= 1/2.
  MyMatrix<T> Q0 = IdentityMat<T>(n);
  for (int iter = 0; iter < 30; iter++) {
    bool all_strict = true;
    for (size_t i = 0; i < ListIneq.size(); i++) {
      T val = ListIneq[i].dot(GramToSymmetricCoordinates(Q0));
      if (val <= 0) {
        all_strict = false;
        Q0 += (T(1) - val) * ListDir[i];
      }
    }
    if (all_strict) {
      break;
    }
    if (iter == 29) {
      return {};
    }
  }
  T max_sq = MaxSquaredCircumRadius(cd, Q0);
  MyMatrix<T> Q_start = Q0 / (T(2) * max_sq);
  MyVector<T> x_T = GramToSymmetricCoordinates(Q_start);
  for (size_t i = 0; i < ListIneq.size(); i++) {
    if (ListIneq[i].dot(x_T) <= 0) {
      return {};
    }
  }
  CoveringData<double> cd_f = ConvertCoveringData<double, T>(cd);
  MyVector<double> x_start(dim);
  for (int u = 0; u < dim; u++) {
    x_start(u) = UniversalScalarConversion<double, T>(x_T(u));
  }
  MaxdetOptions<double> opts = GetDefaultMaxdetOptions<double>();
  MaxdetResult<double> res = SolveCoveringMaxdet(cd_f, x_start, opts, os);
  // Only a converged solve gives a bound: after a stall the gap bound
  // nu / t is not valid, the point not being centered.
  if (!res.success || res.gap_bound < 0) {
    return {};
  }
  x_opt = res.x;
  return std::make_pair(res.det, res.det * std::exp(res.gap_bound));
}

/*
  The lattice points strictly inside the circumsphere of a simplex for the
  Gram matrix Q (a rational approximation of the optimum).
 */
template <typename T, typename Tint>
std::optional<MyVector<T>> FindInteriorPoint(MyMatrix<T> const &Q,
                                             MyMatrix<T> const &Vert,
                                             std::ostream &os) {
  int n = Vert.cols();
  MyMatrix<T> EXT(n + 1, n + 1);
  for (int i = 0; i <= n; i++) {
    EXT(i, 0) = T(1);
    for (int j = 0; j < n; j++) {
      EXT(i, j + 1) = Vert(i, j);
    }
  }
  CP<T> cp = CenterRadiusDelaunayPolytopeGeneral(Q, EXT);
  MyVector<T> c(n);
  for (int i = 0; i < n; i++) {
    c(i) = cp.eCent(i + 1);
  }
  CVPSolver<T, Tint> solver(Q, os);
  ResultShortest<Tint> res = solver.ShortestAtDistance(c, cp.SquareRadius);
  if (res.better_vector) {
    return UniversalVectorConversion<T, Tint>(*res.better_vector);
  }
  return {};
}

/*
  The bound for the pair of simplices S and S' = F_k + w. With emptiness,
  the local Delaunay constraint and up to max_cut rounds of cuts are added.
 */
template <typename T, typename Tint>
PairBound BoundPair(MyMatrix<T> const &VertS, NeighborCandidate<T> const &cand,
                    bool emptiness, int max_cut, std::ostream &os) {
  int n = VertS.cols();
  MyMatrix<T> VertS2 = VertS;
  for (int j = 0; j < n; j++) {
    VertS2(cand.k, j) = cand.w(j);
  }
  std::vector<MyMatrix<T>> ListVert{VertS, VertS2};
  std::vector<MyVector<T>> ListIneq;
  std::vector<MyMatrix<T>> ListDir;
  if (emptiness) {
    ListIneq.push_back(EmptinessInequality(VertS, cand.w));
    std::optional<MyMatrix<T>> opt = PositiveDirection(VertS, cand.w);
    if (opt) {
      ListDir.push_back(*opt);
    } else {
      ListIneq.clear();
    }
  }
  MyVector<double> x_opt;
  std::optional<std::pair<double, double>> opt =
      SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
  if (!opt && !ListIneq.empty()) {
    ListIneq.clear();
    ListDir.clear();
    opt = SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
  }
  if (!opt) {
    return {0, 0, false, 0};
  }
  int n_cut = 0;
  for (int i_cut = 0; i_cut < max_cut && emptiness; i_cut++) {
    MyMatrix<double> Qd = ZeroMatrix<double>(n, n);
    std::vector<MyMatrix<double>> ListA = GetSymmetricBasis<double>(n);
    for (size_t u = 0; u < ListA.size(); u++) {
      Qd += x_opt(u) * ListA[u];
    }
    MyMatrix<T> Q(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        // Rounded to 1e-9: any form serves to find cuts.
        long v = std::lround(Qd(i, j) * 1e9);
        Q(i, j) = T(v) / T(1000000000);
      }
    }
    bool found = false;
    for (auto &Vert : ListVert) {
      std::optional<MyVector<T>> opt_w =
          FindInteriorPoint<T, Tint>(Q, Vert, os);
      if (opt_w) {
        std::optional<MyMatrix<T>> opt_d = PositiveDirection(Vert, *opt_w);
        if (opt_d) {
          ListIneq.push_back(EmptinessInequality(Vert, *opt_w));
          ListDir.push_back(*opt_d);
          found = true;
        }
      }
    }
    if (!found) {
      break;
    }
    MyVector<double> x_new;
    std::optional<std::pair<double, double>> opt_new =
        SolveConfiguration(ListVert, ListIneq, ListDir, x_new, os);
    if (!opt_new) {
      break;
    }
    n_cut++;
    opt = opt_new;
    x_opt = x_new;
  }
  return {opt->first, opt->second, true, n_cut};
}

// ---------------------------------------------------------------------------
// The packing-covering problem
// ---------------------------------------------------------------------------

/*
  For the packing-covering problem, normalize R(Q) = 1. A configuration L
  of Delaunay simplices and any set W of non-zero lattice vectors give

      lambda(Q) <= max { t : BR_S(Q') >= 0 for S in L, Q'[w] >= t for w in W },

  and gamma(Q) = 2 / sqrt(lambda(Q)). The semidefinite program is solved by
  a primal barrier method on (x, t):

      minimize  -s t - sum_S log det BR_S(x) - log det Q(x)
                - sum_w log(Q(x)[w] - t) - sum_i log(a_i . x),

  s increasing. The variable t being free, (c I, t0) with c small and t0
  below the norms is strictly feasible: there is no phase one. The barrier
  parameter nu bounds the gap on t by nu / s.
 */
struct SdpResult {
  bool success;
  double t;
  double gap;
  MyVector<double> x;
};

// One term -log det(M0 + sum_u y_u M_u) of the barrier.
struct LogDetTerm {
  MyMatrix<double> M0;
  std::vector<MyMatrix<double>> ListM;
};

inline std::optional<double>
BarrierPackCov(std::vector<LogDetTerm> const &ListTerm,
               MyMatrix<double> const &Lin, MyVector<double> const &y, double s,
               MyVector<double> *grad, MyMatrix<double> *hess) {
  int dimy = y.size();
  double val = -s * y(dimy - 1);
  if (grad) {
    *grad = ZeroVector<double>(dimy);
    (*grad)(dimy - 1) = -s;
  }
  if (hess) {
    *hess = ZeroMatrix<double>(dimy, dimy);
  }
  for (auto &term : ListTerm) {
    MyMatrix<double> M = term.M0;
    int nvar = term.ListM.size();
    for (int u = 0; u < nvar; u++) {
      M += y(u) * term.ListM[u];
    }
    Eigen::LLT<MyMatrix<double>> llt(M);
    if (llt.info() != Eigen::Success) {
      return {};
    }
    MyMatrix<double> L = llt.matrixL();
    double logdet = 0;
    for (int i = 0; i < M.rows(); i++) {
      if (L(i, i) <= 0) {
        return {};
      }
      logdet += 2 * std::log(L(i, i));
    }
    val -= logdet;
    if (grad || hess) {
      MyMatrix<double> Minv = llt.solve(IdentityMat<double>(M.rows()));
      std::vector<MyMatrix<double>> ListP(nvar);
      for (int u = 0; u < nvar; u++) {
        ListP[u] = Minv * term.ListM[u];
        if (grad) {
          (*grad)(u) -= ListP[u].trace();
        }
      }
      if (hess) {
        for (int u = 0; u < nvar; u++) {
          for (int v = u; v < nvar; v++) {
            double h = (ListP[u] * ListP[v]).trace();
            (*hess)(u, v) += h;
            if (u != v) {
              (*hess)(v, u) += h;
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < Lin.rows(); i++) {
    double l = 0;
    for (int u = 0; u < dimy; u++) {
      l += Lin(i, u) * y(u);
    }
    if (l <= 0) {
      return {};
    }
    val -= std::log(l);
    if (grad || hess) {
      for (int u = 0; u < dimy; u++) {
        if (grad) {
          (*grad)(u) -= Lin(i, u) / l;
        }
        if (hess) {
          for (int v = 0; v < dimy; v++) {
            (*hess)(u, v) += Lin(i, u) * Lin(i, v) / (l * l);
          }
        }
      }
    }
  }
  return val;
}

/*
  The semidefinite program over a covering data cd: the T-space basis
  ListA, the homogeneous linear constraints Alin (an iso-Delaunay cone, or
  emptiness constraints), the simplices whose circumradius is at most 1, and
  the lattice vectors ListW with Q[w] >= t. x_start has to be strictly
  feasible for Alin and the circumradii; t starts below the norms.
 */
inline SdpResult
SolvePackCovData(covering_maxdet::CoveringData<double> const &cd,
                 std::vector<MyVector<double>> const &ListW,
                 MyVector<double> const &x_start,
                 [[maybe_unused]] std::ostream &os) {
  using namespace covering_maxdet;
  int n = cd.n;
  int dim = cd.dim;
  int dimy = dim + 1;
  std::vector<LogDetTerm> ListTerm;
  for (size_t i = 0; i < cd.ListSimplex.size(); i++) {
    std::vector<MyMatrix<double>> ListB = GetLmiBlock(cd, i);
    LogDetTerm term;
    term.M0 = ListB[0];
    for (int u = 0; u < dim; u++) {
      term.ListM.push_back(ListB[u + 1]);
    }
    ListTerm.push_back(term);
  }
  {
    LogDetTerm term;
    term.M0 = ZeroMatrix<double>(n, n);
    term.ListM = cd.ListA;
    ListTerm.push_back(term);
  }
  int n_lin = ListW.size() + cd.Alin.rows();
  MyMatrix<double> Lin = ZeroMatrix<double>(n_lin, dimy);
  int pos = 0;
  for (auto &w : ListW) {
    for (int u = 0; u < dim; u++) {
      Lin(pos, u) = w.dot(cd.ListA[u] * w);
    }
    Lin(pos, dim) = -1;
    pos++;
  }
  for (int i = 0; i < cd.Alin.rows(); i++) {
    for (int u = 0; u < dim; u++) {
      Lin(pos, u) = cd.Alin(i, u);
    }
    pos++;
  }
  double nu =
      double(n) + double(cd.ListSimplex.size()) * double(n + 1) + double(n_lin);
  MyVector<double> y = ZeroVector<double>(dimy);
  for (int u = 0; u < dim; u++) {
    y(u) = x_start(u);
  }
  double min_norm = 1e100;
  for (int i = 0; i < static_cast<int>(ListW.size()); i++) {
    double l = 0;
    for (int u = 0; u < dim; u++) {
      l += Lin(i, u) * y(u);
    }
    min_norm = std::min(min_norm, l);
  }
  y(dim) = min_norm / 2;
  if (!BarrierPackCov(ListTerm, Lin, y, 1.0, nullptr, nullptr)) {
    // The start is not strictly feasible.
    return {false, 0, 0, y};
  }
  double s = 1;
  // The gap bounds t* - t only near the central path: at a point whose
  // Newton decrement delta = sqrt(g^T H^{-1} g) for the parameter s is below
  // 1, the gap is at most (nu + (delta + sqrt(nu)) delta / (1 - delta)) / s
  // (Nesterov and Nemirovski); we accept delta <= 1/2 and report it. A
  // point further from the path gives no bound and the solve is reported
  // as failed, which keeps the configuration.
  double delta = 1e100;
  for (int outer = 0; outer < 200; outer++) {
    for (int newton = 0; newton < 200; newton++) {
      MyVector<double> g;
      MyMatrix<double> H;
      std::optional<double> val = BarrierPackCov(ListTerm, Lin, y, s, &g, &H);
      MyVector<double> dy = -H.ldlt().solve(g);
      double decrement = -g.dot(dy);
      delta = std::sqrt(std::max(decrement, 0.0));
      if (decrement / 2 < 1e-8) {
        break;
      }
      double step = 1;
      bool moved = false;
      for (int bt = 0; bt < 80; bt++) {
        MyVector<double> ynew = y + step * dy;
        std::optional<double> vnew =
            BarrierPackCov(ListTerm, Lin, ynew, s, nullptr, nullptr);
        if (vnew && *vnew <= *val - 0.01 * step * decrement) {
          y = ynew;
          moved = true;
          break;
        }
        step /= 2;
      }
      if (!moved) {
        break;
      }
    }
    double gap = delta < 1
                     ? (nu + (delta + std::sqrt(nu)) * delta / (1 - delta)) / s
                     : 1e100;
    if (gap < 1e-8 * std::max(1.0, std::abs(y(dim)))) {
      MyVector<double> x(dim);
      for (int u = 0; u < dim; u++) {
        x(u) = y(u);
      }
      return {delta <= 0.5, y(dim), gap, x};
    }
    s *= 10;
  }
#ifdef DEBUG_LOCAL_DELAUNAY
  os << "LOCAL_DELAUNAY: SolvePackCovConfiguration did not converge\n";
#endif
  MyVector<double> x(dim);
  for (int u = 0; u < dim; u++) {
    x(u) = y(u);
  }
  return {false, y(dim), nu / s, x};
}

/*
  The semidefinite program for the simplices ListVert, the lattice vectors
  ListW and the homogeneous linear constraints ListIneq (in the symmetric
  coordinates), started from a multiple of the identity.
 */
template <typename T>
SdpResult SolvePackCovConfiguration(std::vector<MyMatrix<T>> const &ListVert,
                                    std::vector<MyVector<T>> const &ListW,
                                    std::vector<MyVector<T>> const &ListIneq,
                                    std::ostream &os) {
  using namespace covering_maxdet;
  int n = ListVert[0].cols();
  int dim = n * (n + 1) / 2;
  CoveringData<double> cd;
  cd.n = n;
  cd.dim = dim;
  cd.ListA = GetSymmetricBasis<double>(n);
  cd.Alin = MyMatrix<double>(ListIneq.size(), dim);
  for (size_t i = 0; i < ListIneq.size(); i++) {
    for (int u = 0; u < dim; u++) {
      cd.Alin(i, u) = UniversalScalarConversion<double, T>(ListIneq[i](u));
    }
  }
  cd.point_density = 1;
  for (auto &Vert : ListVert) {
    cd.ListSimplex.push_back(
        UniversalMatrixConversion<double, T>(DifferenceMatrix(Vert)));
  }
  std::vector<MyVector<double>> ListW_d;
  for (auto &w : ListW) {
    ListW_d.push_back(UniversalVectorConversion<double, T>(w));
  }
  // Start: c I with the circumradii 1/2.
  MyMatrix<double> Id = IdentityMat<double>(n);
  double max_sq = MaxSquaredCircumRadius(cd, Id);
  MyVector<double> x0 =
      GramToSymmetricCoordinates(MyMatrix<double>(Id / (2 * max_sq)));
  return SolvePackCovData(cd, ListW_d, x0, os);
}

/*
  The packing-covering bound of a configuration of simplices: an upper bound
  on lambda for R = 1, from W = all differences of the vertices of the
  simplices, plus up to max_cut rounds of short vectors of the optimum.
  Returns lambda at the optimum and its upper bound t + gap.
 */
template <typename T, typename Tint>
PairBound BoundConfigPackCov(std::vector<MyMatrix<T>> const &ListVert,
                             std::vector<MyVector<T>> ListIneq, int max_cut,
                             std::ostream &os) {
  int n = ListVert[0].cols();
  std::vector<MyVector<T>> ListW;
  std::set<std::vector<T>> SetW;
  auto f_insert_w = [&](MyVector<T> const &w) -> void {
    MyVector<T> wc = w;
    // One of w, -w: the first non-zero coordinate positive.
    for (int i = 0; i < n; i++) {
      if (wc(i) != 0) {
        if (wc(i) < 0) {
          wc = -wc;
        }
        break;
      }
    }
    std::vector<T> key(wc.data(), wc.data() + n);
    if (SetW.insert(key).second) {
      ListW.push_back(wc);
    }
  };
  std::vector<MyVector<T>> ListPt;
  std::set<std::vector<T>> SetPt;
  for (auto &Vert : ListVert) {
    for (int i = 0; i <= n; i++) {
      MyVector<T> v = GetRow(Vert, i);
      std::vector<T> key(v.data(), v.data() + n);
      if (SetPt.insert(key).second) {
        ListPt.push_back(v);
      }
    }
  }
  for (size_t i = 0; i < ListPt.size(); i++) {
    for (size_t j = i + 1; j < ListPt.size(); j++) {
      f_insert_w(MyVector<T>(ListPt[j] - ListPt[i]));
    }
  }
  // The vectors with coordinates in {-1, 0, 1}: they keep the optimum away
  // from the singular forms, which the vertex differences alone do not.
  {
    std::vector<int> c(n, -1);
    while (true) {
      MyVector<T> v(n);
      bool nonzero = false;
      for (int i = 0; i < n; i++) {
        v(i) = T(c[i]);
        if (c[i] != 0) {
          nonzero = true;
        }
      }
      if (nonzero) {
        f_insert_w(v);
      }
      int pos = 0;
      while (pos < n && c[pos] == 1) {
        c[pos] = -1;
        pos++;
      }
      if (pos == n) {
        break;
      }
      c[pos]++;
    }
  }
  SdpResult res = SolvePackCovConfiguration(ListVert, ListW, ListIneq, os);
  if (!res.success && !ListIneq.empty()) {
    SdpResult res2 = SolvePackCovConfiguration(ListVert, ListW, {}, os);
    if (res2.success) {
      ListIneq.clear();
      res = res2;
    }
  }
  int n_cut = 0;
  // The cuts only use the point found, which need not be centered; the
  // success is required of the last solve.
  for (int i_cut = 0; i_cut < max_cut; i_cut++) {
    MyMatrix<double> Qd = ZeroMatrix<double>(n, n);
    std::vector<MyMatrix<double>> ListA = GetSymmetricBasis<double>(n);
    for (size_t u = 0; u < ListA.size(); u++) {
      Qd += res.x(u) * ListA[u];
    }
    MyMatrix<T> Q(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        long v = std::lround(Qd(i, j) * 1e9);
        Q(i, j) = T(v) / T(1000000000);
      }
    }
    // Short vectors of the optimum: the optimum can be close to singular
    // (lambda is only a minimum over W), so an LLL reduction of Q + eps I
    // is used and its basis vectors and their sums and differences are
    // added when their Q-norm is below t. Any lattice vector gives a valid
    // constraint.
    long tv = std::lround(res.t * 1e9);
    T bound = T(tv) / T(1000000000);
    MyMatrix<T> Qreg = Q;
    T eps = T(1) / T(1000000);
    for (int i = 0; i < n; i++) {
      Qreg(i, i) += eps;
    }
    LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(Qreg, os);
    std::vector<MyVector<T>> ListShort;
    for (int i = 0; i < n; i++) {
      MyVector<T> v_i =
          UniversalVectorConversion<T, Tint>(GetMatrixRow(lll.Pmat, i));
      ListShort.push_back(v_i);
      for (int j = i + 1; j < n; j++) {
        MyVector<T> v_j =
            UniversalVectorConversion<T, Tint>(GetMatrixRow(lll.Pmat, j));
        ListShort.push_back(MyVector<T>(v_i + v_j));
        ListShort.push_back(MyVector<T>(v_i - v_j));
      }
    }
    size_t n_before = ListW.size();
    for (auto &v_T : ListShort) {
      if (EvaluationQuadForm(Q, v_T) < bound) {
        f_insert_w(v_T);
      }
    }
    if (ListW.size() == n_before) {
      break;
    }
    res = SolvePackCovConfiguration(ListVert, ListW, ListIneq, os);
    n_cut++;
  }
  if (!res.success) {
    return {0, 0, false, n_cut};
  }
  return {res.t, res.t + res.gap, true, n_cut};
}

/*
  The packing-covering constant over an iso-Delaunay domain: with R = 1,
  the maximum of lambda over the closed domain, lambda being the minimum of
  Q over the edges of the Delaunay subdivision (a shortest vector is an
  edge, and the edges of a coarser subdivision are edges of the domain's).
  ListW0 holds the edge vectors of the orbit representatives of the
  Delaunay polytopes; the other edges are their images, which are added as
  cuts when shorter than t at the optimum, with all the vectors shorter
  than t found by an exact enumeration on a rational rounding of the
  optimum. The result gives the lower bound 2 / sqrt(lambda + gap) of gamma
  on the domain, and the optimum as a Gram matrix.
 */
struct PackCovDomainResult {
  bool success;
  std::string message;
  double lambda;
  double gap;
  double gamma_lower;
  int n_cut;
  MyMatrix<double> Q;
};

template <typename T, typename Tint>
PackCovDomainResult
SolvePackingCoveringDomain(covering_maxdet::CoveringData<T> const &cd,
                           std::vector<MyVector<T>> const &ListW0, int max_cut,
                           std::ostream &os) {
  using namespace covering_maxdet;
  int n = cd.n;
  PackCovDomainResult res{false, "", 0, 0, 0, 0, MyMatrix<double>()};
  CoveringData<double> cd_f = ConvertCoveringData<double, T>(cd);
  std::optional<MyVector<double>> opt_start =
      GetStartingPoint<T, double>(cd, cd_f, os);
  if (!opt_start) {
    res.message = "the domain is too thin to be optimized in floating point";
    return res;
  }
  std::vector<MyVector<double>> ListW;
  std::set<std::vector<long>> SetW;
  auto f_insert = [&](MyVector<T> const &w_T) -> void {
    std::vector<long> key(n);
    for (int i = 0; i < n; i++) {
      key[i] = UniversalScalarConversion<long, T>(w_T(i));
    }
    for (int i = 0; i < n; i++) {
      if (key[i] != 0) {
        if (key[i] < 0) {
          for (auto &x : key) {
            x = -x;
          }
        }
        break;
      }
    }
    if (SetW.insert(key).second) {
      ListW.push_back(UniversalVectorConversion<double, T>(w_T));
    }
  };
  for (auto &w : ListW0) {
    f_insert(w);
  }
  SdpResult sres = SolvePackCovData(cd_f, ListW, *opt_start, os);
  for (int i_cut = 0; i_cut < max_cut; i_cut++) {
    MyMatrix<double> Qd = ZeroMatrix<double>(n, n);
    for (int u = 0; u < cd.dim; u++) {
      Qd += sres.x(u) * cd_f.ListA[u];
    }
    MyMatrix<T> Q(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        long v = std::lround(Qd(i, j) * 1e9);
        Q(i, j) = T(v) / T(1000000000);
      }
    }
    long tv = std::lround(sres.t * 1e9);
    T bound = T(tv) / T(1000000000);
    size_t n_before = ListW.size();
    if (IsPositiveDefinite(Q, os)) {
      CVPSolver<T, Tint> solver(Q, os);
      std::vector<MyVector<Tint>> ListShort =
          solver.at_most_norm_vectors(bound);
      for (auto &v : ListShort) {
        MyVector<T> v_T = UniversalVectorConversion<T, Tint>(v);
        if (EvaluationQuadForm(Q, v_T) < bound) {
          f_insert(v_T);
        }
      }
    }
    if (ListW.size() == n_before) {
      break;
    }
    sres = SolvePackCovData(cd_f, ListW, *opt_start, os);
    res.n_cut++;
  }
  res.success = sres.success;
  res.message = sres.success ? "converged" : "no centered point reached";
  res.lambda = sres.t;
  res.gap = sres.gap;
  res.gamma_lower = 2.0 / std::sqrt(sres.t + sres.gap);
  res.Q = ZeroMatrix<double>(n, n);
  for (int u = 0; u < cd.dim; u++) {
    res.Q += sres.x(u) * cd_f.ListA[u];
  }
  return res;
}

// The packing-covering bound of the pair S, S' = F_k + w.
template <typename T, typename Tint>
PairBound BoundPairPackCov(MyMatrix<T> const &VertS,
                           NeighborCandidate<T> const &cand, bool emptiness,
                           int max_cut, std::ostream &os) {
  int n = VertS.cols();
  MyMatrix<T> VertS2 = VertS;
  for (int j = 0; j < n; j++) {
    VertS2(cand.k, j) = cand.w(j);
  }
  std::vector<MyMatrix<T>> ListVert{VertS, VertS2};
  std::vector<MyVector<T>> ListIneq;
  if (emptiness) {
    ListIneq.push_back(EmptinessInequality(VertS, cand.w));
  }
  return BoundConfigPackCov<T, Tint>(ListVert, ListIneq, max_cut, os);
}

// ---------------------------------------------------------------------------
// The star of a simplex
// ---------------------------------------------------------------------------

/*
  The affine integral maps permuting the vertices of S: a permutation sigma
  of the vertices defines the affine map v_i -> v_sigma(i), which preserves
  Z^n if and only if its linear part M_sigma M^{-1} is integral (its
  determinant is then +-1, the volume being preserved).
 */
template <typename T> struct AffineMap {
  MyMatrix<T> L;
  MyVector<T> b;
  std::vector<int> perm;
};

template <typename T>
std::vector<AffineMap<T>> SimplexStabilizer(SimplexData<T> const &sd) {
  int n = sd.n;
  std::vector<int> perm(n + 1);
  for (int i = 0; i <= n; i++) {
    perm[i] = i;
  }
  std::vector<AffineMap<T>> ListMap;
  do {
    // The images of v_i - v_0 are v_perm(i) - v_perm(0).
    MyMatrix<T> Mp(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Mp(i, j) = sd.Vert(perm[i + 1], j) - sd.Vert(perm[0], j);
      }
    }
    // Row convention: x -> x L + b with M L = Mp.
    MyMatrix<T> L = sd.Minv * Mp;
    if (!IsIntegralMatrix(L)) {
      continue;
    }
    MyVector<T> b(n);
    MyVector<T> v0 = GetRow(sd.Vert, 0);
    MyVector<T> img = L.transpose() * v0;
    for (int j = 0; j < n; j++) {
      b(j) = sd.Vert(perm[0], j) - img(j);
    }
    ListMap.push_back({L, b, perm});
  } while (std::next_permutation(perm.begin(), perm.end()));
  return ListMap;
}

template <typename T>
MyVector<T> ApplyAffine(AffineMap<T> const &am, MyVector<T> const &w) {
  MyVector<T> img = am.L.transpose() * w;
  return img + am.b;
}

/*
  The orbits of the candidate neighbours under the stabilizer of S: g sends
  the neighbour across the facet opposite to v_k with apex w to the
  neighbour across the facet opposite to v_perm(k) with apex g(w). The set
  of candidates over all facets is invariant (the edge bound, the volume and
  the parity frame are). Returns for each candidate the index of its orbit
  representative.
 */
template <typename T>
std::vector<size_t>
OrbitRepresentatives(std::vector<NeighborCandidate<T>> const &ListCand,
                     std::vector<AffineMap<T>> const &ListMap) {
  size_t n_cand = ListCand.size();
  std::vector<size_t> ListRep(n_cand, n_cand);
  if (n_cand == 0) {
    return ListRep;
  }
  int n = ListCand[0].w.size();
  auto key_of = [&](int k, std::vector<long> const &w) -> std::vector<long> {
    std::vector<long> key(n + 1);
    key[0] = k;
    for (int i = 0; i < n; i++) {
      key[i + 1] = w[i];
    }
    return key;
  };
  std::map<std::vector<long>, size_t> map_idx;
  std::vector<std::vector<long>> ListW(n_cand, std::vector<long>(n));
  for (size_t i = 0; i < n_cand; i++) {
    for (int j = 0; j < n; j++) {
      ListW[i][j] = UniversalScalarConversion<long, T>(ListCand[i].w(j));
    }
    map_idx[key_of(ListCand[i].k, ListW[i])] = i;
  }
  // The maps in integers: x -> x L + b.
  struct IntMap {
    std::vector<std::vector<long>> L;
    std::vector<long> b;
    std::vector<int> perm;
  };
  std::vector<IntMap> ListIntMap;
  for (auto &am : ListMap) {
    IntMap im;
    im.L.assign(n, std::vector<long>(n));
    im.b.assign(n, 0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        im.L[i][j] = UniversalScalarConversion<long, T>(am.L(i, j));
      }
      im.b[i] = UniversalScalarConversion<long, T>(am.b(i));
    }
    im.perm = am.perm;
    ListIntMap.push_back(im);
  }
  std::vector<long> img(n);
  for (size_t i = 0; i < n_cand; i++) {
    if (ListRep[i] != n_cand) {
      continue;
    }
    for (auto &im : ListIntMap) {
      for (int j = 0; j < n; j++) {
        long v = im.b[j];
        for (int l = 0; l < n; l++) {
          v += ListW[i][l] * im.L[l][j];
        }
        img[j] = v;
      }
      // When only some facets are enumerated, the images across the other
      // facets are not candidates; the maps fixing the facets give the
      // orbits.
      auto iter = map_idx.find(key_of(im.perm[ListCand[i].k], img));
      if (iter != map_idx.end()) {
        ListRep[iter->second] = i;
      }
    }
  }
  return ListRep;
}

/*
  Two neighbours F_a + w_a and F_b + w_b of S (a != b) around the ridge
  F_a cap F_b. The affine map x -> (beta_a(x), beta_b(x)) sends the ridge to
  the origin, v_a to (1, 0) and v_b to (0, 1); the three simplices project
  to the cones spanned by (v_a, v_b), (v_b, w_a) and (v_a, w_b). The two
  neighbours overlap if and only if the angle of w_a, in (90, 270) degrees,
  exceeds the angle of w_b taken in (180, 360) degrees.
 */
template <typename T>
bool NeighborsOverlap(MyVector<T> const &beta_wa, int a,
                      MyVector<T> const &beta_wb, int b) {
  double xa = UniversalScalarConversion<double, T>(beta_wa(a));
  double ya = UniversalScalarConversion<double, T>(beta_wa(b));
  double xb = UniversalScalarConversion<double, T>(beta_wb(a));
  double yb = UniversalScalarConversion<double, T>(beta_wb(b));
  double tha = std::atan2(ya, xa);
  if (tha < 0) {
    tha += 2 * M_PI;
  }
  double thb = std::atan2(yb, xb);
  if (thb < M_PI / 2) {
    thb += 2 * M_PI;
  }
  return tha > thb + 1e-12;
}

template <typename T> struct StarCandidate {
  int k;
  MyVector<T> w;
  MyVector<T> beta;
  // The vertices of the neighbour F_k + w.
  MyMatrix<T> Vert;
  double pair_bound;
  // Filled by StarSearch: the edges w - v_j as (parity class, edge id) and
  // beta in floating point.
  std::vector<std::pair<int, int>> edges;
  std::vector<double> beta_d;
};

struct StarSearchResult {
  size_t n_node;
  size_t n_maxdet;
  size_t n_full;
  bool complete;
  // The smallest and largest bound among the complete stars.
  double min_full_bound;
  double max_full_bound;
};

/*
  The edge ids: one id per edge vector up to sign, and its class in
  Z^n / 2 Z^n as an integer in [1, 2^n).
 */
template <typename T> struct EdgeIndex {
  std::map<std::vector<long>, int> map;
  std::pair<int, int> get(MyVector<T> const &e) {
    int n = e.size();
    std::vector<long> key(n);
    int cls = 0;
    for (int i = 0; i < n; i++) {
      long v = UniversalScalarConversion<long, T>(e(i));
      key[i] = v;
      if (((v % 2) + 2) % 2 == 1) {
        cls |= (1 << i);
      }
    }
    for (int i = 0; i < n; i++) {
      if (key[i] != 0) {
        if (key[i] < 0) {
          for (auto &x : key) {
            x = -x;
          }
        }
        break;
      }
    }
    auto iter = map.find(key);
    if (iter != map.end()) {
      return {cls, iter->second};
    }
    int id = map.size();
    map[key] = id;
    return {cls, id};
  }
};

/*
  The depth first search over the stars of S: one neighbour per facet, with
  the parity frame, the non-overlap around the ridges and the determinant
  maximization of S and the chosen neighbours with all the emptiness
  constraints between them. A star surviving all tests is reported; if
  there is none and the search is complete, S does not occur in the
  Delaunay triangulation of a form below the threshold.
 */
template <typename T, typename Tint>
StarSearchResult
StarSearch(SimplexData<T> const &sd,
           std::vector<std::vector<StarCandidate<T>>> &ListCandFacet,
           std::vector<StarCandidate<T>> &RootList, int root_facet,
           double threshold, bool is_pc, size_t node_limit, size_t max_print,
           std::ostream &os_out, std::ostream &os) {
  int n = sd.n;
  EdgeIndex<T> eidx;
  std::vector<int> slot0(1 << n, -1);
  for (int i = 0; i <= n; i++) {
    for (int j = i + 1; j <= n; j++) {
      std::pair<int, int> pr =
          eidx.get(MyVector<T>(GetRow(sd.Vert, j) - GetRow(sd.Vert, i)));
      slot0[pr.first] = pr.second;
    }
  }
  auto prepare = [&](StarCandidate<T> &c) -> void {
    c.edges.clear();
    for (int j = 0; j <= n; j++) {
      if (j != c.k) {
        c.edges.push_back(eidx.get(MyVector<T>(c.w - GetRow(sd.Vert, j))));
      }
    }
    c.beta_d.resize(n + 1);
    for (int i = 0; i <= n; i++) {
      c.beta_d[i] = UniversalScalarConversion<double, T>(c.beta(i));
    }
  };
  for (auto &list : ListCandFacet) {
    for (auto &c : list) {
      prepare(c);
    }
  }
  for (auto &c : RootList) {
    prepare(c);
  }
  StarSearchResult res{0, 0, 0, true, 0, 0};
  std::vector<StarCandidate<T> const *> chosen(n + 1, nullptr);
  // The bound of the configuration of the current node.
  double node_bound = 0;
  // The overlap test of NeighborsOverlap, on the floating point betas.
  auto overlap = [&](StarCandidate<T> const &ca,
                     StarCandidate<T> const &cb) -> bool {
    int a = ca.k, b = cb.k;
    double tha = std::atan2(ca.beta_d[b], ca.beta_d[a]);
    if (tha < 0) {
      tha += 2 * M_PI;
    }
    double thb = std::atan2(cb.beta_d[b], cb.beta_d[a]);
    if (thb < M_PI / 2) {
      thb += 2 * M_PI;
    }
    return tha > thb + 1e-12;
  };
  auto compatible = [&](StarCandidate<T> const &c,
                        std::vector<int> const &slot) -> bool {
    for (auto &pr : c.edges) {
      if (slot[pr.first] != -1 && slot[pr.first] != pr.second) {
        return false;
      }
    }
    for (int a = 0; a <= n; a++) {
      if (chosen[a] != nullptr && overlap(*chosen[a], c)) {
        return false;
      }
    }
    return true;
  };
  auto bound_config = [&]() -> double {
    std::vector<MyMatrix<T>> ListVert{sd.Vert};
    std::vector<MyVector<T>> ListIneq;
    std::vector<MyMatrix<T>> ListDir;
    std::vector<StarCandidate<T> const *> ListC;
    for (int a = 0; a <= n; a++) {
      if (chosen[a] != nullptr) {
        ListC.push_back(chosen[a]);
        ListVert.push_back(chosen[a]->Vert);
      }
    }
    auto add = [&](MyMatrix<T> const &Vert, MyVector<T> const &w) -> void {
      std::optional<MyMatrix<T>> opt = PositiveDirection(Vert, w);
      if (opt) {
        ListIneq.push_back(EmptinessInequality(Vert, w));
        ListDir.push_back(*opt);
      }
    };
    for (auto *c : ListC) {
      add(sd.Vert, c->w);
      for (auto *d : ListC) {
        if (c != d && !(c->w == d->w)) {
          add(c->Vert, d->w);
        }
      }
    }
    res.n_maxdet++;
    if (is_pc) {
      PairBound pb = BoundConfigPackCov<T, Tint>(ListVert, ListIneq, 3, os);
      return pb.success ? pb.bound : 1e100;
    }
    MyVector<double> x_opt;
    std::optional<std::pair<double, double>> opt =
        SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
    if (!opt) {
      ListIneq.clear();
      ListDir.clear();
      opt = SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
    }
    if (!opt) {
      // No bound: keep the configuration.
      return 1e100;
    }
    return opt->second;
  };
  std::function<void(std::vector<int> const &, int)> f_rec =
      [&](std::vector<int> const &slot, int depth) -> void {
    res.n_node++;
    if (res.n_node % 1000 == 0) {
      os_out << "progress: n_node=" << res.n_node
             << " n_maxdet=" << res.n_maxdet << " n_full=" << res.n_full
             << " depth=" << depth << std::endl;
    }
    if (res.n_node > node_limit) {
      res.complete = false;
      return;
    }
    if (depth == n + 1) {
      res.n_full++;
      if (res.n_full == 1 || node_bound < res.min_full_bound) {
        res.min_full_bound = node_bound;
      }
      if (res.n_full == 1 || node_bound > res.max_full_bound) {
        res.max_full_bound = node_bound;
      }
      if (res.n_full <= max_print) {
        os_out << "FULL STAR: bound=" << node_bound;
        for (int a = 0; a <= n; a++) {
          os_out << " [k=" << a << " w=";
          for (int i = 0; i < n; i++) {
            os_out << (i > 0 ? "," : "") << chosen[a]->w(i);
          }
          os_out << "]";
        }
        os_out << std::endl;
      }
      return;
    }
    // The unassigned facet with the fewest compatible candidates.
    int best_k = -1;
    std::vector<StarCandidate<T> const *> best_list;
    for (int k = 0; k <= n; k++) {
      if (chosen[k] != nullptr) {
        continue;
      }
      std::vector<StarCandidate<T> const *> list;
      for (auto &c : ListCandFacet[k]) {
        if (compatible(c, slot)) {
          list.push_back(&c);
        }
      }
      if (list.empty()) {
        return;
      }
      if (best_k == -1 || list.size() < best_list.size()) {
        best_k = k;
        best_list = list;
      }
    }
    for (auto *c : best_list) {
      chosen[best_k] = c;
      double bnd = bound_config();
      if (bnd >= threshold) {
        node_bound = bnd;
        std::vector<int> slot2 = slot;
        for (auto &pr : c->edges) {
          slot2[pr.first] = pr.second;
        }
        f_rec(slot2, depth + 1);
      }
      chosen[best_k] = nullptr;
      if (!res.complete) {
        return;
      }
    }
  };
  // The root: one candidate per orbit of the stabilizer of the root facet.
  for (auto &c : RootList) {
    if (!compatible(c, slot0)) {
      continue;
    }
    chosen[root_facet] = &c;
    double bnd = bound_config();
    if (bnd >= threshold) {
      std::vector<int> slot2 = slot0;
      for (auto &pr : c.edges) {
        slot2[pr.first] = pr.second;
      }
      f_rec(slot2, 1);
    }
    chosen[root_facet] = nullptr;
    if (!res.complete) {
      break;
    }
  }
  return res;
}

// ---------------------------------------------------------------------------
// Extension of a configuration of unimodular simplices
// ---------------------------------------------------------------------------

/*
  A configuration of unimodular simplices, extended facet by facet. The
  neighbours across a facet are obtained from the list P0 of the apexes of
  the neighbours of the standard simplex sigma_0 across its facet opposite to
  the origin that survive the pair bound: for a simplex N with the facet
  opposite to u_k, the unimodular affine map y -> u_k + y M, M the matrix of
  the rows u_j - u_k (j != k), sends sigma_0 to N and its facet opposite to
  the origin to that facet, and the pair bound is invariant.

  The search fills at each node the open facet with the fewest candidates
  compatible with the parity frame of all the edges of the configuration,
  and prunes with the bound of the whole configuration and the emptiness
  constraints between every simplex and every vertex of the configuration.
  A configuration is refuted when every branch dies before depth_limit
  simplices have been added.
 */
struct ExtendResult {
  size_t n_node;
  size_t n_bound;
  bool refuted;
  bool complete;
  int max_depth;
};

template <typename T, typename Tint>
ExtendResult ExtendConfiguration(std::vector<MyMatrix<T>> const &ListStart,
                                 std::vector<MyVector<T>> const &P0, bool is_pc,
                                 double threshold, int depth_limit,
                                 size_t node_limit, int Vmax_general,
                                 double dstar, std::ostream &os_out,
                                 std::ostream &os) {
  int n = ListStart[0].cols();
  using Pt = std::vector<long>;
  auto to_pt = [&](MyVector<T> const &v) -> Pt {
    Pt p(n);
    for (int i = 0; i < n; i++) {
      p[i] = UniversalScalarConversion<long, T>(v(i));
    }
    return p;
  };
  auto to_vec = [&](Pt const &p) -> MyVector<T> {
    MyVector<T> v(n);
    for (int i = 0; i < n; i++) {
      v(i) = T(p[i]);
    }
    return v;
  };
  // The configuration: simplices as sorted lists of points.
  std::vector<std::vector<Pt>> ListSimp;
  for (auto &Vert : ListStart) {
    std::vector<Pt> simp;
    for (int i = 0; i <= n; i++) {
      simp.push_back(to_pt(GetRow(Vert, i)));
    }
    std::sort(simp.begin(), simp.end());
    ListSimp.push_back(simp);
  }
  ExtendResult res{0, 0, true, true, 0};
  // The parity frame: class index -> canonical edge vector.
  auto edge_key = [&](Pt const &a, Pt const &b) -> std::pair<int, Pt> {
    Pt e(n);
    int cls = 0;
    for (int i = 0; i < n; i++) {
      e[i] = b[i] - a[i];
      if (((e[i] % 2) + 2) % 2 == 1) {
        cls |= (1 << i);
      }
    }
    for (int i = 0; i < n; i++) {
      if (e[i] != 0) {
        if (e[i] < 0) {
          for (auto &x : e) {
            x = -x;
          }
        }
        break;
      }
    }
    return {cls, e};
  };
  auto frame_ok = [&](std::vector<std::vector<Pt>> const &config) -> bool {
    std::map<int, Pt> frame;
    for (auto &simp : config) {
      for (int i = 0; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
          std::pair<int, Pt> pr = edge_key(simp[i], simp[j]);
          if (pr.first == 0) {
            return false;
          }
          auto iter = frame.find(pr.first);
          if (iter == frame.end()) {
            frame[pr.first] = pr.second;
          } else if (iter->second != pr.second) {
            return false;
          }
        }
      }
    }
    return true;
  };
  auto bound = [&](std::vector<std::vector<Pt>> const &config) -> double {
    res.n_bound++;
    std::vector<MyMatrix<T>> ListVert;
    std::set<Pt> SetPt;
    for (auto &simp : config) {
      MyMatrix<T> Vert(n + 1, n);
      for (int i = 0; i <= n; i++) {
        for (int j = 0; j < n; j++) {
          Vert(i, j) = T(simp[i][j]);
        }
        SetPt.insert(simp[i]);
      }
      ListVert.push_back(Vert);
    }
    std::vector<MyVector<T>> ListIneq;
    std::vector<MyMatrix<T>> ListDir;
    for (size_t is = 0; is < config.size(); is++) {
      std::set<Pt> own(config[is].begin(), config[is].end());
      for (auto &p : SetPt) {
        if (own.count(p) == 0) {
          MyVector<T> w = to_vec(p);
          std::optional<MyMatrix<T>> opt = PositiveDirection(ListVert[is], w);
          if (opt) {
            ListIneq.push_back(EmptinessInequality(ListVert[is], w));
            ListDir.push_back(*opt);
          }
        }
      }
    }
    if (is_pc) {
      PairBound pb = BoundConfigPackCov<T, Tint>(ListVert, ListIneq, 3, os);
      return pb.success ? pb.bound : 1e100;
    }
    MyVector<double> x_opt;
    std::optional<std::pair<double, double>> opt =
        SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
    if (!opt) {
      ListIneq.clear();
      ListDir.clear();
      opt = SolveConfiguration(ListVert, ListIneq, ListDir, x_opt, os);
    }
    // Lazy cuts: the lattice points strictly inside a circumsphere at the
    // optimum give valid emptiness constraints; a few rounds. A round whose
    // solve fails keeps the previous bound, which is valid.
    for (int i_cut = 0; i_cut < 5 && opt; i_cut++) {
      std::vector<MyMatrix<double>> ListA = GetSymmetricBasis<double>(n);
      MyMatrix<double> Qd = ZeroMatrix<double>(n, n);
      for (size_t u = 0; u < ListA.size(); u++) {
        Qd += x_opt(u) * ListA[u];
      }
      MyMatrix<T> Q(n, n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Q(i, j) =
              T(std::lround((Qd(i, j) + Qd(j, i)) / 2 * 1e9)) / T(1000000000);
        }
      }
      if (!IsPositiveDefinite(Q, os)) {
        break;
      }
      bool found = false;
      for (auto &Vert : ListVert) {
        std::optional<MyVector<T>> opt_w =
            FindInteriorPoint<T, Tint>(Q, Vert, os);
        if (opt_w) {
          std::optional<MyMatrix<T>> opt_d = PositiveDirection(Vert, *opt_w);
          if (opt_d) {
            ListIneq.push_back(EmptinessInequality(Vert, *opt_w));
            ListDir.push_back(*opt_d);
            found = true;
          }
        }
      }
      if (!found) {
        break;
      }
      MyVector<double> x_new;
      std::optional<std::pair<double, double>> opt_new =
          SolveConfiguration(ListVert, ListIneq, ListDir, x_new, os);
      if (!opt_new) {
        break;
      }
      opt = opt_new;
      x_opt = x_new;
    }
    return opt ? opt->second : 1e100;
  };
  // The open facets: (simplex index, removed vertex index).
  auto open_facets = [&](std::vector<std::vector<Pt>> const &config)
      -> std::vector<std::pair<size_t, int>> {
    std::map<std::vector<Pt>, int> count;
    for (auto &simp : config) {
      for (int k = 0; k <= n; k++) {
        std::vector<Pt> f;
        for (int i = 0; i <= n; i++) {
          if (i != k) {
            f.push_back(simp[i]);
          }
        }
        count[f]++;
      }
    }
    std::vector<std::pair<size_t, int>> ret;
    for (size_t is = 0; is < config.size(); is++) {
      for (int k = 0; k <= n; k++) {
        std::vector<Pt> f;
        for (int i = 0; i <= n; i++) {
          if (i != k) {
            f.push_back(config[is][i]);
          }
        }
        if (count[f] == 1) {
          ret.push_back({is, k});
        }
      }
    }
    return ret;
  };
  // The candidate apexes across the facet of simplex is opposite to vertex k.
  // The general mode (Vmax_general > 0): the apexes of the neighbours of a
  // simplex of any volume across a facet, of volume at most Vmax_general,
  // surviving the pair bound. They are computed once per type of simplex on
  // a representative and transported by a unimodular affine map. The type
  // is the volume V and the set I of the vertices on which the barycentric
  // coordinates of the classes of Z^n modulo the lattice of the vertices are
  // not all integral; the vertices are ordered with I first, which gives an
  // integral map when the stabilizer is Sym(I) x Sym(J), as for the
  // simplices of volume at most 2. When the map is not integral the
  // neighbours are enumerated directly.
  auto ordered_vertices = [&](std::vector<Pt> const &simp)
      -> std::pair<std::pair<int, int>, std::vector<int>> {
    MyMatrix<T> Vert(n + 1, n);
    for (int i = 0; i <= n; i++) {
      for (int j = 0; j < n; j++) {
        Vert(i, j) = T(simp[i][j]);
      }
    }
    SimplexData<T> sd = GetSimplexData(Vert);
    std::vector<bool> in_I(n + 1, false);
    if (sd.V > 1) {
      for (auto &a : CosetRepresentatives(sd)) {
        int64_t sum = 0;
        for (int i = 0; i < n; i++) {
          if (a[i] != 0) {
            in_I[i + 1] = true;
          }
          sum += a[i];
        }
        if (sum % sd.V != 0) {
          in_I[0] = true;
        }
      }
    }
    std::vector<int> order;
    for (int i = 0; i <= n; i++) {
      if (in_I[i]) {
        order.push_back(i);
      }
    }
    int sizI = order.size();
    for (int i = 0; i <= n; i++) {
      if (!in_I[i]) {
        order.push_back(i);
      }
    }
    return {{sd.V, sizI}, order};
  };
  // The apexes across the facet opposite to vertex k of simp, computed
  // directly, filtered by the pair bound, one bound per orbit of the
  // stabilizer of simp.
  auto direct_apexes =
      [&](std::vector<Pt> const &simp) -> std::vector<std::vector<Pt>> {
    MyMatrix<T> Vert(n + 1, n);
    for (int i = 0; i <= n; i++) {
      for (int j = 0; j < n; j++) {
        Vert(i, j) = T(simp[i][j]);
      }
    }
    SimplexData<T> sd = GetSimplexData(Vert);
    std::vector<std::vector<Pt>> ret(n + 1);
    std::optional<double> opt_lambda = InertiaLowerBound(n, sd.V, dstar);
    if (!opt_lambda) {
      return ret;
    }
    double EdgeBound = 4.0 / *opt_lambda;
    std::vector<int> all_facets;
    for (int k = 0; k <= n; k++) {
      all_facets.push_back(k);
    }
    std::vector<NeighborCandidate<T>> ListCand =
        EnumerateNeighborCandidates<T, Tint>(sd, all_facets, EdgeBound,
                                             Vmax_general, os);
    std::vector<AffineMap<T>> ListMap = SimplexStabilizer(sd);
    std::vector<size_t> ListRep = OrbitRepresentatives(ListCand, ListMap);
    std::vector<std::optional<bool>> ListKeep(ListCand.size());
    for (size_t i = 0; i < ListCand.size(); i++) {
      size_t r = ListRep[i];
      if (!ListKeep[r]) {
        PairBound pb =
            is_pc ? BoundPairPackCov<T, Tint>(Vert, ListCand[r], false, 3, os)
                  : BoundPair<T, Tint>(Vert, ListCand[r], false, 0, os);
        ListKeep[r] = !pb.success || pb.bound >= threshold;
      }
      if (*ListKeep[r]) {
        ret[ListCand[i].k].push_back(to_pt(ListCand[i].w));
      }
    }
    return ret;
  };
  struct TypeData {
    std::vector<Pt> rep;
    std::vector<int> order;
    std::vector<std::vector<Pt>> apexes;
  };
  std::map<std::pair<int, int>, TypeData> map_type;
  std::map<std::vector<Pt>, std::vector<std::vector<Pt>>> cache_direct;
  auto general_apexes = [&](std::vector<Pt> const &simp,
                            int k) -> std::vector<Pt> {
    std::pair<std::pair<int, int>, std::vector<int>> pr =
        ordered_vertices(simp);
    auto iter = map_type.find(pr.first);
    if (iter == map_type.end()) {
      TypeData td{simp, pr.second, direct_apexes(simp)};
      iter = map_type.emplace(pr.first, td).first;
    }
    TypeData const &td = iter->second;
    std::vector<int> const &oN = pr.second;
    std::vector<int> const &oR = td.order;
    // The map x -> N_oN[0] + (x - R_oR[0]) L with
    // (N_oN[i] - N_oN[0]) = (R_oR[i] - R_oR[0]) L.
    MyMatrix<T> MR(n, n), MN(n, n);
    for (int i = 1; i <= n; i++) {
      for (int j = 0; j < n; j++) {
        MR(i - 1, j) = T(td.rep[oR[i]][j] - td.rep[oR[0]][j]);
        MN(i - 1, j) = T(simp[oN[i]][j] - simp[oN[0]][j]);
      }
    }
    MyMatrix<T> L = Inverse(MR) * MN;
    if (!IsIntegralMatrix(L)) {
      auto it2 = cache_direct.find(simp);
      if (it2 == cache_direct.end()) {
        it2 = cache_direct.emplace(simp, direct_apexes(simp)).first;
      }
      return it2->second[k];
    }
    int p = 0;
    while (oN[p] != k) {
      p++;
    }
    std::vector<Pt> ret;
    for (auto &x : td.apexes[oR[p]]) {
      Pt w(n);
      for (int j = 0; j < n; j++) {
        T val = T(simp[oN[0]][j]);
        for (int i = 0; i < n; i++) {
          val += T(x[i] - td.rep[oR[0]][i]) * L(i, j);
        }
        w[j] = UniversalScalarConversion<long, T>(val);
      }
      ret.push_back(w);
    }
    return ret;
  };
  auto candidates = [&](std::vector<std::vector<Pt>> const &config, size_t is,
                        int k) -> std::vector<std::vector<Pt>> {
    std::vector<Pt> const &simp = config[is];
    if (Vmax_general > 0) {
      std::set<std::vector<Pt>> existing(config.begin(), config.end());
      std::vector<std::vector<Pt>> ret;
      std::vector<Pt> list_w = general_apexes(simp, k);
      for (auto &w : list_w) {
        std::vector<Pt> nsimp;
        for (int i = 0; i <= n; i++) {
          if (i != k) {
            nsimp.push_back(simp[i]);
          }
        }
        nsimp.push_back(w);
        std::sort(nsimp.begin(), nsimp.end());
        if (existing.count(nsimp) > 0) {
          continue;
        }
        std::vector<std::vector<Pt>> config2 = config;
        config2.push_back(nsimp);
        if (frame_ok(config2)) {
          ret.push_back(nsimp);
        }
      }
      return ret;
    }
    std::vector<Pt> rows;
    for (int i = 0; i <= n; i++) {
      if (i != k) {
        Pt d(n);
        for (int j = 0; j < n; j++) {
          d[j] = simp[i][j] - simp[k][j];
        }
        rows.push_back(d);
      }
    }
    std::set<std::vector<Pt>> existing(config.begin(), config.end());
    std::vector<std::vector<Pt>> ret;
    for (auto &y : P0) {
      Pt w(n);
      for (int j = 0; j < n; j++) {
        long v = simp[k][j];
        for (int i = 0; i < n; i++) {
          v += UniversalScalarConversion<long, T>(y(i)) * rows[i][j];
        }
        w[j] = v;
      }
      std::vector<Pt> nsimp;
      for (int i = 0; i <= n; i++) {
        if (i != k) {
          nsimp.push_back(simp[i]);
        }
      }
      nsimp.push_back(w);
      std::sort(nsimp.begin(), nsimp.end());
      if (existing.count(nsimp) > 0) {
        continue;
      }
      std::vector<std::vector<Pt>> config2 = config;
      config2.push_back(nsimp);
      if (frame_ok(config2)) {
        ret.push_back(nsimp);
      }
    }
    return ret;
  };
  std::function<void(std::vector<std::vector<Pt>> const &, int)> f_rec =
      [&](std::vector<std::vector<Pt>> const &config, int depth) -> void {
    res.n_node++;
    if (res.n_node % 100 == 0) {
      os_out << "progress: n_node=" << res.n_node << " n_bound=" << res.n_bound
             << " depth=" << depth << std::endl;
    }
    if (res.n_node > node_limit) {
      res.complete = false;
      return;
    }
    res.max_depth = std::max(res.max_depth, depth);
    if (depth >= depth_limit) {
      res.refuted = false;
      // The optimum of the configuration reached, for the covering: a
      // candidate form whose actual covering density can be computed.
      if (!is_pc) {
        std::vector<MyMatrix<T>> ListVert;
        for (auto &simp : config) {
          MyMatrix<T> Vert(n + 1, n);
          for (int i = 0; i <= n; i++) {
            for (int j = 0; j < n; j++) {
              Vert(i, j) = T(simp[i][j]);
            }
          }
          ListVert.push_back(Vert);
        }
        MyVector<double> x_opt;
        std::optional<std::pair<double, double>> opt =
            SolveConfiguration(ListVert, {}, {}, x_opt, os);
        if (opt) {
          std::vector<MyMatrix<double>> ListA = GetSymmetricBasis<double>(n);
          MyMatrix<double> Qd = ZeroMatrix<double>(n, n);
          for (size_t u = 0; u < ListA.size(); u++) {
            Qd += x_opt(u) * ListA[u];
          }
          os_out << "OPTIMUM det=" << opt->first << " of the configuration of "
                 << config.size() << " simplices, Gram matrix:\n";
          os_out << n << " " << n << "\n";
          for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
              os_out << " " << std::lround(Qd(i, j) * 1e9) << "/1000000000";
            }
            os_out << "\n";
          }
        }
      }
      return;
    }
    std::vector<std::pair<size_t, int>> list_open = open_facets(config);
    std::vector<std::vector<Pt>> best;
    bool has_best = false;
    for (auto &of : list_open) {
      std::vector<std::vector<Pt>> cand =
          candidates(config, of.first, of.second);
      if (cand.empty()) {
        return;
      }
      if (!has_best || cand.size() < best.size()) {
        best = cand;
        has_best = true;
      }
    }
    if (!has_best) {
      // No open facet: a closed configuration, not refuted.
      res.refuted = false;
      return;
    }
    for (auto &nsimp : best) {
      std::vector<std::vector<Pt>> config2 = config;
      config2.push_back(nsimp);
      if (bound(config2) >= threshold) {
        f_rec(config2, depth + 1);
        if (!res.refuted || !res.complete) {
          return;
        }
      }
    }
  };
  if (!frame_ok(ListSimp) || bound(ListSimp) < threshold) {
    return res;
  }
  f_rec(ListSimp, 0);
  return res;
}

} // namespace local_delaunay

// clang-format off
#endif  // SRC_DELAUNAY_LOCALDELAUNAYCONFIG_H_
// clang-format on
