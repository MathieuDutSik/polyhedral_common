// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_COVERINGPATH_H_
#define SRC_POLYNORM_POLYNORM_COVERINGPATH_H_

// clang-format off
#include "PolyNorm_MinimalCovering.h"
#include <algorithm>
#include <functional>
#include <optional>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_COVERINGPATH
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_COVERINGPATH
#endif

#ifdef TIMINGS
#define TIMINGS_POLYNORM_COVERINGPATH
#endif

/*
  The covering radius along a Minkowski path.

  For two polytopes P_1, P_2 in R^n and K_lambda = (1 - lambda) P_1 +
  lambda P_2, the function g(lambda) = 1 / mu(K_lambda) on [0,1] is
  continuous and piecewise affine, and is computed exactly.

  The formulation. With c_1, ..., c_m the facet normals of P_1 + P_2 and
  h_i(lambda) = (1 - lambda) h_{P_1}(c_i) + lambda h_{P_2}(c_i), the support
  functions being additive, K_lambda = {x : c_i.x <= h_i(lambda)} for all
  lambda in [0,1]. A point p is not covered by the interiors of the
  translates mu K_lambda + z exactly when for every lattice point z some
  i has c_i.(p - z) >= mu h_i(lambda). With q = p / mu and s = 1 / mu this
  reads
      for all z there is i with c_i.q - s c_i.z - h_i(lambda) >= 0,
  which is linear in (q, s, lambda) jointly. Let G be this set. Since
  s -> (q, s, lambda) stays in G when s increases (a point uncovered at
  level mu is uncovered at smaller levels), g(lambda) = min{s : (q, s,
  lambda) in G}, and translating a deepest hole so that 0 is a nearest
  lattice point allows the normalization q in K_lambda, again linear.

  The pieces. A hole system of n+1 pairs (z_k, i_k) with independent rows
  (c_{i_k}, -c_{i_k}.z_k) determines (q, s) affine in lambda, by
      c_{i_k}.q - s c_{i_k}.z_k = h_{i_k}(lambda).
  Where the point p = q/s is a hole, that is where for all z the maximum
  over i of the affine functions c_i.(q - s z) - h_i(lambda) is >= 0, it
  gives the upper bound g <= s. Each z excludes an open interval of
  lambda, so the validity set is closed and computed exactly; the lattice
  points that matter are those with c_i.z > (c_i.q - h_i) / s for all i,
  and these linear fractional bounds are monotone on an interval, so they
  are bounded by their values at the ends.

  The algorithm. A branch and bound over (q, s, lambda), lambda in an
  interval I, minimizes a linear function s - s_H(lambda) over G: at the
  optimum of a node, a lattice point z with q - s z interior to K_lambda is
  branched on, one child (z, i) per facet. The optimum of a node is a
  vertex of its region, whose active rows are the bounds on lambda and
  tight pairs, so a point of G found this way is a hole pinned by its tight
  pairs: they have rank n+1 in (q, s). (The deepest holes found by
  PolyNorm_Covering.h need not be pinned: for a tile every boundary point
  is a deepest hole, and the facet cones of its linear programs are not
  tight pairs.)
  --- g(lambda_0) is the minimum of s with I = {lambda_0}. Every system of
      n+1 tight pairs of the hole found with independent rows gives a
      piece, kept on its validity interval around lambda_0. The minimum U
      of the pieces is an upper bound on g.
  --- On each interval I where U is a single piece s_H, the search starts
      with the bound 0, so that only the nodes of negative value are
      explored. If none leads to a point of G, g = s_H on I, which is
      proved. Otherwise the point found has s < U(lambda*), its pieces are
      added, and the loop goes on.
  --- An interval of [0,1] without any piece is treated by an evaluation at
      its middle.
  There are finitely many hole systems (the pairs come from a bounded
  region since mu is bounded on [0,1]) and each failed certificate adds a
  piece below U at lambda*, hence a new one, so the algorithm terminates, with
  all the pieces and break points rational and every step exact.

  The facet cones used by PolyNorm_Covering.h to make the children disjoint
  are bilinear in (q, lambda), so the children here overlap: the search is
  correct but explores more nodes.

  The lower bound s >= s_lo = 1 / (2 mu_0) keeps the linear programs
  bounded and is never active at a point of G:
  if a translate of a polytope Q fits in t_j P_j, j = 1, 2, both P_j
  containing 0, then a translate of Q fits in max(t_1, t_2) K_lambda, so
  mu(K_lambda) <= max(t_1, t_2) mu(Q) for all lambda. P_1 and P_2 are
  recentered at the isobarycenters of their vertices, which translates
  K_lambda and does not change mu.
 */

template <typename T> struct PolyNormPathData {
  int n;
  int m;
  // The recentered vertices of P_1 and P_2, homogeneous coordinates.
  MyMatrix<T> EXT1;
  MyMatrix<T> EXT2;
  // The facet normals of P_1 + P_2 (rows) and h_i(lambda) = h1_i + lambda
  // d_i.
  MyMatrix<T> C;
  MyVector<T> h1;
  MyVector<T> d;
  T s_lo;
};

// The support function of the polytope with homogeneous vertices EXT.
template <typename T>
T PolyNorm_PathSupport(MyMatrix<T> const &EXT, MyVector<T> const &c) {
  int n = c.size();
  T val(0);
  for (int i = 0; i < EXT.rows(); i++) {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += c(k) * EXT(i, k + 1);
    }
    if (i == 0 || scal > val) {
      val = scal;
    }
  }
  return val;
}

// The vertices of a full dimensional polytope given by points,
// recentered at their isobarycenter.
template <typename T>
MyMatrix<T> PolyNorm_PathRecenter(MyMatrix<T> const &EXT_input,
                                  std::ostream &os) {
  MyMatrix<T> EXT =
      SelectRow(EXT_input, get_distinct_extreme_rays_indices(EXT_input, os));
  MyVector<T> center = PolyNorm_MinCovCenter(EXT);
  for (int i = 0; i < EXT.rows(); i++) {
    for (int k = 0; k < center.size(); k++) {
      EXT(i, k + 1) -= center(k);
    }
  }
  return EXT;
}

// The smallest t such that a translate of Q fits in t P, for the polytope
// P with homogeneous vertices EXT containing 0 in its interior, by the
// fitting linear program of PolyNorm_Covering.h on the normals C.
template <typename T>
T PolyNorm_PathFit(MyMatrix<T> const &EXT, MyMatrix<T> const &C,
                   std::vector<MyVector<T>> const &ListQ, std::ostream &os) {
  int m = C.rows();
  int n = C.cols();
  MyMatrix<T> Lmat(m, n);
  for (int i = 0; i < m; i++) {
    MyVector<T> c = GetMatrixRow(C, i);
    T h = PolyNorm_PathSupport(EXT, c);
    for (int k = 0; k < n; k++) {
      Lmat(i, k) = c(k) / h;
    }
  }
  return PolyNorm_TranslatedFit(Lmat, ListQ, os);
}

template <typename T>
PolyNormPathData<T> PolyNorm_PathBuildData(MyMatrix<T> const &EXT1_input,
                                           MyMatrix<T> const &EXT2_input,
                                           std::ostream &os) {
  PolyNormPathData<T> data;
  data.n = EXT1_input.cols() - 1;
  int n = data.n;
  data.EXT1 = PolyNorm_PathRecenter(EXT1_input, os);
  data.EXT2 = PolyNorm_PathRecenter(EXT2_input, os);
  // The facets of P_1 + P_2, from the sums of vertices.
  int n1 = data.EXT1.rows();
  int n2 = data.EXT2.rows();
  MyMatrix<T> EXTsum(n1 * n2, n + 1);
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      int pos = i1 * n2 + i2;
      EXTsum(pos, 0) = 1;
      for (int k = 0; k < n; k++) {
        EXTsum(pos, k + 1) = data.EXT1(i1, k + 1) + data.EXT2(i2, k + 1);
      }
    }
  }
  MyMatrix<T> FACsum = DirectFacetComputationInequalities_General(
      EXTsum, std::string("cdd"), os);
  data.m = FACsum.rows();
  int m = data.m;
  data.C = MyMatrix<T>(m, n);
  data.h1 = MyVector<T>(m);
  data.d = MyVector<T>(m);
  for (int i = 0; i < m; i++) {
    // b + a.x >= 0, that is c.x <= b with c = -a.
    MyVector<T> c(n);
    for (int k = 0; k < n; k++) {
      c(k) = -FACsum(i, k + 1);
      data.C(i, k) = c(k);
    }
    T h1 = PolyNorm_PathSupport(data.EXT1, c);
    T h2 = PolyNorm_PathSupport(data.EXT2, c);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    if (h1 + h2 != FACsum(i, 0)) {
      std::cerr << "POLYNORM_COVERINGPATH: The support functions are not "
                << "additive on the facet i=" << i << "\n";
      throw TerminalException{1};
    }
#endif
    data.h1(i) = h1;
    data.d(i) = h2 - h1;
  }
  // The upper bound mu_0 on mu(K_lambda) for all lambda.
  std::vector<std::pair<std::vector<MyVector<T>>, T>> ListQmu;
  {
    std::vector<MyVector<T>> ListQ;
    for (int j = 0; j < n; j++) {
      MyVector<T> e = ZeroVector<T>(n);
      e(j) = 1;
      ListQ.push_back(e);
      e(j) = -1;
      ListQ.push_back(e);
    }
    ListQmu.push_back({ListQ, T(n) / T(2)});
  }
  if (n <= 12) {
    std::vector<MyVector<T>> ListQ;
    for (size_t sgn = 0; sgn < (size_t(1) << n); sgn++) {
      MyVector<T> v(n);
      for (int j = 0; j < n; j++) {
        v(j) = ((sgn >> j) & 1) ? T(1) : T(-1);
      }
      ListQ.push_back(v);
    }
    ListQmu.push_back({ListQ, T(1) / T(2)});
  }
  for (int sgn : {1, -1}) {
    std::vector<MyVector<T>> ListQ{ZeroVector<T>(n)};
    for (int j = 0; j < n; j++) {
      MyVector<T> e = ZeroVector<T>(n);
      e(j) = sgn;
      ListQ.push_back(e);
    }
    ListQmu.push_back({ListQ, T(n)});
  }
  std::optional<T> mu0;
  for (auto &ePair : ListQmu) {
    T t1 = PolyNorm_PathFit(data.EXT1, data.C, ePair.first, os);
    T t2 = PolyNorm_PathFit(data.EXT2, data.C, ePair.first, os);
    T bound = std::max(t1, t2) * ePair.second;
    if (!mu0 || bound < *mu0) {
      mu0 = bound;
    }
  }
  // Half the bound, so that s >= s_lo is never active at a point of G.
  data.s_lo = T(1) / (2 * *mu0);
#ifdef DEBUG_POLYNORM_COVERINGPATH
  os << "POLYNORM_COVERINGPATH: n=" << n << " m=" << m << " mu0=" << *mu0
     << "\n";
#endif
  return data;
}

template <typename T>
T PolyNorm_PathH(PolyNormPathData<T> const &data, int i, T const &lambda) {
  return data.h1(i) + lambda * data.d(i);
}

// The vertices of K_lambda.
template <typename T>
MyMatrix<T> PolyNorm_PathVertices(PolyNormPathData<T> const &data,
                                  T const &lambda, std::ostream &os) {
  MyMatrix<T> FAC(data.m, data.n + 1);
  for (int i = 0; i < data.m; i++) {
    FAC(i, 0) = PolyNorm_PathH(data, i, lambda);
    for (int k = 0; k < data.n; k++) {
      FAC(i, k + 1) = -data.C(i, k);
    }
  }
  MyMatrix<T> EXT = PolyNorm_MinCovVertices(FAC, os);
  return SelectRow(EXT, get_distinct_extreme_rays_indices(EXT, os));
}

// A pair of a hole system: the lattice point z and the normal i.
template <typename T, typename Tint> struct PolyNormPathPair {
  MyVector<Tint> z;
  MyVector<T> z_T;
  int i;
};

// A hole at lambda: the point q = s p of G with its tight pairs, those with
// c_i.(q - s z) = h_i(lambda). When it is a vertex of the region of a node
// of the branch and bound below, the tight pairs have rank n+1 in (q, s)
// and pin it.
template <typename T, typename Tint> struct PolyNormPathHole {
  T lambda;
  T s;
  MyVector<T> q;
  std::vector<PolyNormPathPair<T, Tint>> ListTight;
};

// The hole of the point x = (q, s, lambda) of G. Its tight pairs are the
// pairs (z, i) with q - s z on the facet i of K_lambda, and the pairs of
// ExtraPairs (those of the node where x was found) whose row is active:
// q - s z is then on the hyperplane of the facet i but possibly outside
// K_lambda. A system built with them still passes through x and is kept
// only where it is valid.
template <typename T, typename Tint>
PolyNormPathHole<T, Tint>
PolyNorm_PathHole(PolyNormPathData<T> const &data, MyVector<T> const &x,
                  std::vector<PolyNormPathPair<T, Tint>> const &ExtraPairs,
                  std::ostream &os) {
  int n = data.n;
  int m = data.m;
  MyVector<T> q(n);
  for (int k = 0; k < n; k++) {
    q(k) = x(k);
  }
  T s = x(n);
  T lambda = x(n + 1);
  // The lattice points z with q - s z in K_lambda:
  // h_i(lambda) - c_i.q + s c_i.z >= 0.
  MyMatrix<T> FAC(m, n + 1);
  std::vector<T> ListH(m), ListCq(m);
  for (int i = 0; i < m; i++) {
    ListH[i] = PolyNorm_PathH(data, i, lambda);
    T cq(0);
    for (int k = 0; k < n; k++) {
      cq += data.C(i, k) * q(k);
      FAC(i, k + 1) = s * data.C(i, k);
    }
    ListCq[i] = cq;
    FAC(i, 0) = ListH[i] - cq;
  }
  std::vector<PolyNormPathPair<T, Tint>> ListTight;
  for (auto &z : GetListIntegralPoint<T, Tint>(FAC, os)) {
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    for (int i = 0; i < m; i++) {
      T cz(0);
      for (int k = 0; k < n; k++) {
        cz += data.C(i, k) * z_T(k);
      }
      if (ListCq[i] - s * cz == ListH[i]) {
        ListTight.push_back({z, z_T, i});
      }
    }
  }
  for (auto &P : ExtraPairs) {
    bool is_present = false;
    for (auto &Q : ListTight) {
      if (Q.i == P.i && Q.z == P.z) {
        is_present = true;
      }
    }
    T cz(0);
    for (int k = 0; k < n; k++) {
      cz += data.C(P.i, k) * P.z_T(k);
    }
    if (!is_present && ListCq[P.i] - s * cz == ListH[P.i]) {
      ListTight.push_back(P);
    }
  }
#ifdef DEBUG_POLYNORM_COVERINGPATH
  os << "POLYNORM_COVERINGPATH: hole at lambda=" << lambda << " s=" << s
     << " |tight|=" << ListTight.size() << "\n";
#endif
  return {lambda, s, q, ListTight};
}

// The covering radius of K_lambda by PolyNorm_Covering.h, for the checks.
template <typename T, typename Tint, typename Tgroup>
T PolyNorm_PathCoveringRadius(PolyNormPathData<T> const &data, T const &lambda,
                              std::ostream &os) {
  MyMatrix<T> EXT = PolyNorm_PathVertices(data, lambda, os);
  PolyNormData<T, Tint, Tgroup> dataK =
      PolyNorm_BuildDataTrivialGroup<T, Tint, Tgroup>(EXT, os);
  return ComputePolyNormCovering<T, Tint, Tgroup>(dataK, os).mu;
}

// An affine piece s(lambda) = s0 + s1 lambda of g, valid on [alpha, beta].
template <typename T> struct PolyNormPathPiece {
  T s0;
  T s1;
  T alpha;
  T beta;
};

template <typename T>
T PolyNorm_PathPieceValue(PolyNormPathPiece<T> const &P, T const &lambda) {
  return P.s0 + P.s1 * lambda;
}

// The validity interval around lambda0 of the hole (q, s) = (q0, s0) +
// lambda (q1, s1), or nothing when it is not valid at lambda0.
template <typename T, typename Tint>
std::optional<std::pair<T, T>>
PolyNorm_PathValidity(PolyNormPathData<T> const &data, MyVector<T> const &q0,
                      MyVector<T> const &q1, T const &s0, T const &s1,
                      T const &lambda0, std::ostream &os) {
  int n = data.n;
  int m = data.m;
  // J = {lambda in [0,1] : s(lambda) >= s_lo}: a valid hole has
  // s >= g >= s_lo.
  T J_lo(0), J_hi(1);
  if (s1 > 0) {
    T x = (data.s_lo - s0) / s1;
    J_lo = std::max(J_lo, x);
  } else if (s1 < 0) {
    T x = (data.s_lo - s0) / s1;
    J_hi = std::min(J_hi, x);
  } else if (s0 < data.s_lo) {
    return {};
  }
  if (lambda0 < J_lo || lambda0 > J_hi) {
    return {};
  }
  // The lattice points z with c_i.z > r_i(lambda) for all i for some
  // lambda in J, r_i = (c_i.q - h_i) / s being monotone on J.
  auto r_val = [&](int i, T const &lambda) -> T {
    T val = -PolyNorm_PathH(data, i, lambda);
    for (int k = 0; k < n; k++) {
      val += data.C(i, k) * (q0(k) + lambda * q1(k));
    }
    return val / (s0 + lambda * s1);
  };
  MyMatrix<T> FAC(m, n + 1);
  for (int i = 0; i < m; i++) {
    FAC(i, 0) = -std::min(r_val(i, J_lo), r_val(i, J_hi));
    for (int k = 0; k < n; k++) {
      FAC(i, k + 1) = data.C(i, k);
    }
  }
  T alpha = J_lo;
  T beta = J_hi;
  for (auto &z : GetListIntegralPoint<T, Tint>(FAC, os)) {
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    // phi_i(lambda) = A_i + B_i lambda = c_i.(q - s z) - h_i(lambda), and
    // the bad set {phi_i < 0 for all i} is an open interval (lo, hi).
    std::optional<T> lo, hi;
    bool empty = false;
    for (int i = 0; i < m; i++) {
      T cz(0), cq0(0), cq1(0);
      for (int k = 0; k < n; k++) {
        cz += data.C(i, k) * z_T(k);
        cq0 += data.C(i, k) * q0(k);
        cq1 += data.C(i, k) * q1(k);
      }
      T A = cq0 - s0 * cz - data.h1(i);
      T B = cq1 - s1 * cz - data.d(i);
      if (B > 0) {
        T x = -A / B;
        if (!hi || x < *hi) {
          hi = x;
        }
      } else if (B < 0) {
        T x = -A / B;
        if (!lo || x > *lo) {
          lo = x;
        }
      } else if (A >= 0) {
        empty = true;
        break;
      }
    }
    if (empty || (lo && hi && *lo >= *hi)) {
      continue;
    }
    bool below = hi && *hi <= lambda0;
    bool above = lo && *lo >= lambda0;
    if (!below && !above) {
      // lambda0 is in the bad interval.
      return {};
    }
    if (below && *hi > alpha) {
      alpha = *hi;
    }
    if (above && *lo < beta) {
      beta = *lo;
    }
  }
  return std::pair<T, T>{alpha, beta};
}

// The pieces given by the systems of n+1 tight pairs of a hole with
// independent rows, on their validity intervals around the hole.
template <typename T, typename Tint>
std::vector<PolyNormPathPiece<T>>
PolyNorm_PathPieces(PolyNormPathData<T> const &data,
                    PolyNormPathHole<T, Tint> const &hole, std::ostream &os) {
  int n = data.n;
  int n_pair = hole.ListTight.size();
  auto get_row = [&](int u) -> MyVector<T> {
    PolyNormPathPair<T, Tint> const &P = hole.ListTight[u];
    MyVector<T> row(n + 1);
    T cz(0);
    for (int k = 0; k < n; k++) {
      row(k) = data.C(P.i, k);
      cz += data.C(P.i, k) * P.z_T(k);
    }
    row(n) = -cz;
    return row;
  };
  std::vector<PolyNormPathPiece<T>> ListPiece;
  std::set<std::vector<T>> set_system;
  std::vector<int> sel;
  std::function<void(int)> f_rec = [&](int start) -> void {
    if (int(sel.size()) == n + 1) {
      MyMatrix<T> M(n + 1, n + 1);
      MyVector<T> b0(n + 1), b1(n + 1);
      for (int r = 0; r <= n; r++) {
        MyVector<T> row = get_row(sel[r]);
        for (int k = 0; k <= n; k++) {
          M(r, k) = row(k);
        }
        b0(r) = data.h1(hole.ListTight[sel[r]].i);
        b1(r) = data.d(hole.ListTight[sel[r]].i);
      }
      MyMatrix<T> Minv = Inverse(M);
      MyVector<T> x0 = Minv * b0;
      MyVector<T> x1 = Minv * b1;
      std::vector<T> key;
      for (int k = 0; k <= n; k++) {
        key.push_back(x0(k));
        key.push_back(x1(k));
      }
      if (!set_system.insert(key).second) {
        return;
      }
      MyVector<T> q0(n), q1(n);
      for (int k = 0; k < n; k++) {
        q0(k) = x0(k);
        q1(k) = x1(k);
      }
      std::optional<std::pair<T, T>> opt = PolyNorm_PathValidity<T, Tint>(
          data, q0, q1, x0(n), x1(n), hole.lambda, os);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
      if (!opt) {
        std::cerr << "POLYNORM_COVERINGPATH: A system of tight pairs is not "
                  << "valid at its own hole\n";
        throw TerminalException{1};
      }
#endif
      if (opt) {
        ListPiece.push_back({x0(n), x1(n), opt->first, opt->second});
      }
      return;
    }
    for (int u = start; u < n_pair; u++) {
      sel.push_back(u);
      MyMatrix<T> Mpart(sel.size(), n + 1);
      for (size_t r = 0; r < sel.size(); r++) {
        MyVector<T> row = get_row(sel[r]);
        for (int k = 0; k <= n; k++) {
          Mpart(r, k) = row(k);
        }
      }
      if (RankMat(Mpart) == int(sel.size())) {
        f_rec(u + 1);
      }
      sel.pop_back();
    }
  };
  f_rec(0);
#ifdef DEBUG_POLYNORM_COVERINGPATH
  os << "POLYNORM_COVERINGPATH: lambda=" << hole.lambda
     << " |pieces|=" << ListPiece.size() << "\n";
#endif
  return ListPiece;
}

// The branch and bound over (q, s, lambda) in G with lambda in [la, lb]
// that minimizes s - s0 - s1 lambda. Its two uses:
// --- the evaluation of g at lambda_0: la = lb = lambda_0, s0 = s1 = 0, and
//     the exact minimum is computed;
// --- the certificate of g >= s0 + s1 lambda on [la, lb]: the search starts
//     with the bound 0, so that only the nodes of negative value are
//     explored, and stops at the first point of G of negative value.
// The optimum of a node is a vertex of its region, so that the point
// recorded is pinned by its tight pairs.
template <typename T, typename Tint> struct PolyNormPathSolver {
  PolyNormPathData<T> const &data;
  T la;
  T lb;
  T s0;
  T s1;
  bool stop_at_first;
  std::ostream &os;
  int n;
  int m;
  size_t n_node;
  size_t n_lp;
  // The best value found, or the bound to beat, and its point.
  std::optional<T> best;
  std::optional<MyVector<T>> x_best;
  std::vector<PolyNormPathPair<T, Tint>> Node_best;

  PolyNormPathSolver(PolyNormPathData<T> const &_data, T const &_la,
                     T const &_lb, T const &_s0, T const &_s1,
                     bool _stop_at_first, std::ostream &_os)
      : data(_data), la(_la), lb(_lb), s0(_s0), s1(_s1),
        stop_at_first(_stop_at_first), os(_os), n(_data.n), m(_data.m),
        n_node(0), n_lp(0) {
    if (stop_at_first) {
      best = T(0);
    }
  }

  // The rows of a node over the variables (q, s, lambda): b + a.x >= 0.
  MyMatrix<T>
  NodeInequalities(std::vector<PolyNormPathPair<T, Tint>> const &Node) const {
    int n_row = 3 + m + Node.size();
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_row, n + 3);
    int pos = 0;
    // lambda - la >= 0, lb - lambda >= 0, s - s_lo >= 0
    ListIneq(pos, 0) = -la;
    ListIneq(pos, n + 2) = 1;
    pos++;
    ListIneq(pos, 0) = lb;
    ListIneq(pos, n + 2) = -1;
    pos++;
    ListIneq(pos, 0) = -data.s_lo;
    ListIneq(pos, n + 1) = 1;
    pos++;
    // q in K_lambda: h1_i + d_i lambda - c_i.q >= 0
    for (int i = 0; i < m; i++) {
      ListIneq(pos, 0) = data.h1(i);
      for (int k = 0; k < n; k++) {
        ListIneq(pos, k + 1) = -data.C(i, k);
      }
      ListIneq(pos, n + 2) = data.d(i);
      pos++;
    }
    // c_i.q - s c_i.z - h1_i - d_i lambda >= 0
    for (auto &P : Node) {
      T cz(0);
      for (int k = 0; k < n; k++) {
        cz += data.C(P.i, k) * P.z_T(k);
        ListIneq(pos, k + 1) = data.C(P.i, k);
      }
      ListIneq(pos, 0) = -data.h1(P.i);
      ListIneq(pos, n + 1) = -cz;
      ListIneq(pos, n + 2) = -data.d(P.i);
      pos++;
    }
    return ListIneq;
  }

  std::optional<std::pair<MyVector<T>, T>>
  SolveLP(std::vector<PolyNormPathPair<T, Tint>> const &Node) {
    n_lp++;
    MyMatrix<T> ListIneq = NodeInequalities(Node);
    MyVector<T> ToBeMinimized = ZeroVector<T>(n + 3);
    ToBeMinimized(n + 1) = 1;
    ToBeMinimized(n + 2) = -s1;
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.DirectSolution) {
      return {};
    }
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    if (!eSol.DualSolution) {
      std::cerr << "POLYNORM_COVERINGPATH: The node LP is unbounded\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> const &x = *eSol.DirectSolution;
    T value = x(n) - s0 - s1 * x(n + 1);
    return std::pair<MyVector<T>, T>{x, value};
  }

  // Whether the point x of G is pinned: its tight pairs, with the bounds
  // on lambda that are active, have rank n+2 in (q, s, lambda).
  bool IsPinned(MyVector<T> const &x,
                std::vector<PolyNormPathPair<T, Tint>> const &Node) const {
    PolyNormPathHole<T, Tint> hole =
        PolyNorm_PathHole<T, Tint>(data, x, Node, os);
    std::vector<MyVector<T>> ListRow;
    for (auto &P : hole.ListTight) {
      MyVector<T> row(n + 2);
      T cz(0);
      for (int k = 0; k < n; k++) {
        row(k) = data.C(P.i, k);
        cz += data.C(P.i, k) * P.z_T(k);
      }
      row(n) = -cz;
      row(n + 1) = -data.d(P.i);
      ListRow.push_back(row);
    }
    if (x(n + 1) == la || x(n + 1) == lb) {
      MyVector<T> row = ZeroVector<T>(n + 2);
      row(n + 1) = 1;
      ListRow.push_back(row);
    }
    if (ListRow.empty()) {
      return false;
    }
    MyMatrix<T> M(ListRow.size(), n + 2);
    for (size_t r = 0; r < ListRow.size(); r++) {
      for (int k = 0; k < n + 2; k++) {
        M(r, k) = ListRow[r](k);
      }
    }
    return RankMat(M) == n + 2;
  }

  // A vertex of the optimal face of the node: the minimum of a generic
  // linear function over the points of the node of value at most value.
  // It is unique, hence a vertex, and the region of the node being bounded,
  // it exists.
  MyVector<T>
  VertexOfOptimalFace(std::vector<PolyNormPathPair<T, Tint>> const &Node,
                      T const &value) {
    n_lp++;
    MyMatrix<T> ListIneq0 = NodeInequalities(Node);
    int n_row = ListIneq0.rows();
    MyMatrix<T> ListIneq(n_row + 1, n + 3);
    for (int r = 0; r < n_row; r++) {
      for (int k = 0; k < n + 3; k++) {
        ListIneq(r, k) = ListIneq0(r, k);
      }
    }
    // value + s0 + s1 lambda - s >= 0
    for (int k = 0; k < n + 3; k++) {
      ListIneq(n_row, k) = 0;
    }
    ListIneq(n_row, 0) = value + s0;
    ListIneq(n_row, n + 1) = -1;
    ListIneq(n_row, n + 2) = s1;
    // The weights 1/p for distinct primes p on q and lambda.
    std::vector<int> ListPrime{2,  3,  5,  7,  11, 13, 17, 19, 23, 29,
                               31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    MyVector<T> ToBeMinimized = ZeroVector<T>(n + 3);
    for (int k = 0; k <= n; k++) {
      int idx = (k < n) ? k + 1 : n + 2;
      ToBeMinimized(idx) = T(1) / T(ListPrime[k % ListPrime.size()]);
    }
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    if (!eSol.DirectSolution || !eSol.DualSolution) {
      std::cerr << "POLYNORM_COVERINGPATH: The optimal face LP has no "
                << "optimum\n";
      throw TerminalException{1};
    }
#endif
    return *eSol.DirectSolution;
  }

  bool IsPruned(T const &value) const {
    if (stop_at_first && x_best) {
      return true;
    }
    return best && value >= *best;
  }

  void Process(std::vector<PolyNormPathPair<T, Tint>> const &Node,
               MyVector<T> const &x, T const &value) {
    if (IsPruned(value)) {
      return;
    }
    n_node++;
    T s = x(n);
    T lambda = x(n + 1);
    // The lattice points z with q - s z in K_lambda:
    // h_i(lambda) - c_i.q + s c_i.z >= 0.
    MyMatrix<T> FAC(m, n + 1);
    std::vector<T> ListH(m), ListCq(m);
    for (int i = 0; i < m; i++) {
      ListH[i] = PolyNorm_PathH(data, i, lambda);
      T cq(0);
      for (int k = 0; k < n; k++) {
        cq += data.C(i, k) * x(k);
        FAC(i, k + 1) = s * data.C(i, k);
      }
      ListCq[i] = cq;
      FAC(i, 0) = ListH[i] - cq;
    }
    // The most violated lattice point: q - s z interior to K_lambda.
    std::optional<std::pair<T, MyVector<Tint>>> branch;
    for (auto &z : GetListIntegralPoint<T, Tint>(FAC, os)) {
      MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
      std::optional<T> max_phi;
      for (int i = 0; i < m; i++) {
        T cz(0);
        for (int k = 0; k < n; k++) {
          cz += data.C(i, k) * z_T(k);
        }
        T phi = ListCq[i] - s * cz - ListH[i];
        if (!max_phi || phi > *max_phi) {
          max_phi = phi;
        }
      }
      if (*max_phi < 0) {
        if (!branch || *max_phi < branch->first) {
          branch = std::pair<T, MyVector<Tint>>{*max_phi, z};
        }
      }
    }
    if (!branch) {
      // (q, s, lambda) is in G. The point recorded must be a pinned hole:
      // otherwise a vertex of the optimal face of the node is taken, which
      // is either in G and pinned or has a violating lattice point.
      if (!IsPinned(x, Node)) {
        MyVector<T> x_vert = VertexOfOptimalFace(Node, value);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
        T value_vert = x_vert(n) - s0 - s1 * x_vert(n + 1);
        if (value_vert != value || x_vert == x) {
          std::cerr << "POLYNORM_COVERINGPATH: The vertex of the optimal "
                    << "face is not pinned or not of the same value, x="
                    << StringVector(x) << " x_vert=" << StringVector(x_vert)
                    << " value=" << value << " value_vert=" << value_vert
                    << "\n";
          throw TerminalException{1};
        }
#endif
        n_node--;
        Process(Node, x_vert, value);
        return;
      }
      best = value;
      x_best = x;
      Node_best = Node;
#ifdef DEBUG_POLYNORM_COVERINGPATH
      os << "POLYNORM_COVERINGPATH: point of G of value " << value
         << " at lambda=" << lambda << "\n";
#endif
      return;
    }
    MyVector<Tint> const &z = branch->second;
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    struct Child {
      std::vector<PolyNormPathPair<T, Tint>> Node;
      MyVector<T> x;
      T value;
    };
    std::vector<Child> ListChild;
    for (int i = 0; i < m; i++) {
      std::vector<PolyNormPathPair<T, Tint>> NodeNew = Node;
      NodeNew.push_back({z, z_T, i});
      std::optional<std::pair<MyVector<T>, T>> opt = SolveLP(NodeNew);
      if (opt && !IsPruned(opt->second)) {
        ListChild.push_back(
            {std::move(NodeNew), std::move(opt->first), opt->second});
      }
    }
    std::stable_sort(ListChild.begin(), ListChild.end(),
                     [](Child const &a, Child const &b) -> bool {
                       return a.value < b.value;
                     });
    for (auto &eChild : ListChild) {
      Process(eChild.Node, eChild.x, eChild.value);
    }
  }

  void Run() {
    std::vector<PolyNormPathPair<T, Tint>> Node;
    std::optional<std::pair<MyVector<T>, T>> opt = SolveLP(Node);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    if (!opt) {
      std::cerr << "POLYNORM_COVERINGPATH: The root LP is infeasible\n";
      throw TerminalException{1};
    }
#endif
    Process(Node, opt->first, opt->second);
  }
};

template <typename T> struct PolyNormCoveringPath {
  // The certified pieces, consecutive and covering [0,1]:
  // 1 / mu(K_lambda) = s0 + s1 lambda on [alpha, beta].
  std::vector<PolyNormPathPiece<T>> ListPiece;
  size_t n_evaluation;
  size_t n_certificate;
  size_t n_node;
  size_t n_lp;
};

// The minimum of the pieces on [0,1] as consecutive (interval, piece), or
// the middle of an interval that no piece covers.
template <typename T>
std::pair<std::vector<PolyNormPathPiece<T>>, std::optional<T>>
PolyNorm_PathEnvelope(std::vector<PolyNormPathPiece<T>> const &ListPiece) {
  std::set<T> set_break{T(0), T(1)};
  size_t len = ListPiece.size();
  for (size_t u = 0; u < len; u++) {
    PolyNormPathPiece<T> const &P = ListPiece[u];
    set_break.insert(P.alpha);
    set_break.insert(P.beta);
    for (size_t v = u + 1; v < len; v++) {
      PolyNormPathPiece<T> const &Q = ListPiece[v];
      if (P.s1 != Q.s1) {
        T x = (Q.s0 - P.s0) / (P.s1 - Q.s1);
        T lo = std::max(P.alpha, Q.alpha);
        T hi = std::min(P.beta, Q.beta);
        if (lo <= x && x <= hi) {
          set_break.insert(x);
        }
      }
    }
  }
  std::vector<T> ListBreak(set_break.begin(), set_break.end());
  std::vector<PolyNormPathPiece<T>> ListEnv;
  for (size_t k = 0; k + 1 < ListBreak.size(); k++) {
    T const &a = ListBreak[k];
    T const &b = ListBreak[k + 1];
    if (a < 0 || b > 1) {
      continue;
    }
    T mid = (a + b) / 2;
    std::optional<size_t> best;
    for (size_t u = 0; u < len; u++) {
      PolyNormPathPiece<T> const &P = ListPiece[u];
      if (P.alpha <= a && b <= P.beta) {
        if (!best || PolyNorm_PathPieceValue(P, mid) <
                         PolyNorm_PathPieceValue(ListPiece[*best], mid)) {
          best = u;
        }
      }
    }
    if (!best) {
      return {{}, mid};
    }
    PolyNormPathPiece<T> const &P = ListPiece[*best];
    if (!ListEnv.empty() && ListEnv.back().s0 == P.s0 &&
        ListEnv.back().s1 == P.s1) {
      ListEnv.back().beta = b;
    } else {
      ListEnv.push_back({P.s0, P.s1, a, b});
    }
  }
  return {ListEnv, {}};
}

// The exact function 1 / mu((1 - lambda) P_1 + lambda P_2) on [0,1].
template <typename T, typename Tint, typename Tgroup>
PolyNormCoveringPath<T> PolyNorm_ComputeCoveringPath(MyMatrix<T> const &EXT1,
                                                     MyMatrix<T> const &EXT2,
                                                     std::ostream &os) {
#ifdef TIMINGS_POLYNORM_COVERINGPATH
  MicrosecondTime time;
#endif
  PolyNormPathData<T> data = PolyNorm_PathBuildData(EXT1, EXT2, os);
  PolyNormCoveringPath<T> res{{}, 0, 0, 0, 0};
  std::vector<PolyNormPathPiece<T>> ListPiece;
  // The pieces of the hole at the point x of G.
  auto add_hole = [&](MyVector<T> const &x,
                      std::vector<PolyNormPathPair<T, Tint>> const &Node) -> T {
    PolyNormPathHole<T, Tint> hole =
        PolyNorm_PathHole<T, Tint>(data, x, Node, os);
    for (auto &P : PolyNorm_PathPieces<T, Tint>(data, hole, os)) {
      ListPiece.push_back(P);
    }
    return hole.s;
  };
  // The evaluation of g at lambda, with the pieces of its deepest hole.
  auto evaluate = [&](T const &lambda) -> T {
    PolyNormPathSolver<T, Tint> solver(data, lambda, lambda, T(0), T(0), false,
                                       os);
    solver.Run();
    res.n_evaluation++;
    res.n_node += solver.n_node;
    res.n_lp += solver.n_lp;
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    if (!solver.x_best) {
      std::cerr << "POLYNORM_COVERINGPATH: No point of G at lambda=" << lambda
                << "\n";
      throw TerminalException{1};
    }
#endif
    return add_hole(*solver.x_best, solver.Node_best);
  };
  evaluate(T(0));
  evaluate(T(1));
  // The (interval, piece) proved, a sub-interval with the same piece being
  // proved as well.
  std::vector<PolyNormPathPiece<T>> ListCertified;
  auto is_certified = [&](PolyNormPathPiece<T> const &P) -> bool {
    for (auto &Q : ListCertified) {
      if (Q.s0 == P.s0 && Q.s1 == P.s1 && Q.alpha <= P.alpha &&
          P.beta <= Q.beta) {
        return true;
      }
    }
    return false;
  };
  while (true) {
    std::pair<std::vector<PolyNormPathPiece<T>>, std::optional<T>> env =
        PolyNorm_PathEnvelope(ListPiece);
    if (env.second) {
#ifdef DEBUG_POLYNORM_COVERINGPATH
      os << "POLYNORM_COVERINGPATH: no piece at lambda=" << *env.second << "\n";
#endif
      evaluate(*env.second);
      continue;
    }
    std::optional<MyVector<T>> x_bad;
    std::vector<PolyNormPathPair<T, Tint>> Node_bad;
    [[maybe_unused]] std::optional<T> U_bad;
    for (auto &P : env.first) {
      if (is_certified(P)) {
        continue;
      }
      PolyNormPathSolver<T, Tint> solver(data, P.alpha, P.beta, P.s0, P.s1,
                                         true, os);
      solver.Run();
      res.n_certificate++;
      res.n_node += solver.n_node;
      res.n_lp += solver.n_lp;
#ifdef DEBUG_POLYNORM_COVERINGPATH
      os << "POLYNORM_COVERINGPATH: certificate on [" << P.alpha << ","
         << P.beta << "] of s=" << P.s0 << "+" << P.s1
         << " lambda: " << (solver.x_best ? "fails" : "proved")
         << " n_node=" << solver.n_node << "\n";
#endif
      if (solver.x_best) {
        x_bad = solver.x_best;
        Node_bad = solver.Node_best;
        U_bad = PolyNorm_PathPieceValue(P, (*x_bad)(data.n + 1));
        break;
      }
      ListCertified.push_back(P);
    }
    if (!x_bad) {
      res.ListPiece = env.first;
      break;
    }
    // The point of G found has s < U(lambda): its pieces bring the envelope
    // down at lambda, so that the loop makes progress.
    [[maybe_unused]] T s_hole = add_hole(*x_bad, Node_bad);
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
    {
      T lambda_bad = (*x_bad)(data.n + 1);
      std::optional<T> U;
      for (auto &P : ListPiece) {
        if (P.alpha <= lambda_bad && lambda_bad <= P.beta) {
          T val = PolyNorm_PathPieceValue(P, lambda_bad);
          if (!U || val < *U) {
            U = val;
          }
        }
      }
      if (!U || *U > s_hole || *U >= *U_bad) {
        std::cerr << "POLYNORM_COVERINGPATH: The envelope is not brought "
                  << "down at lambda=" << lambda_bad << "\n";
        throw TerminalException{1};
      }
    }
#endif
  }
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH
  // The pieces are consecutive and g is continuous.
  for (size_t u = 0; u + 1 < res.ListPiece.size(); u++) {
    PolyNormPathPiece<T> const &P = res.ListPiece[u];
    PolyNormPathPiece<T> const &Q = res.ListPiece[u + 1];
    if (P.beta != Q.alpha || PolyNorm_PathPieceValue(P, P.beta) !=
                                 PolyNorm_PathPieceValue(Q, Q.alpha)) {
      std::cerr << "POLYNORM_COVERINGPATH: The pieces do not join at " << P.beta
                << "\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef SANITY_CHECK_POLYNORM_COVERINGPATH_EXTENSIVE
  std::cerr << "POLYNORM_COVERINGPATH: SANITY_CHECK_POLYNORM_COVERINGPATH_"
               "EXTENSIVE, evaluating g at the break points and the middles\n";
  for (auto &P : res.ListPiece) {
    T mid = (P.alpha + P.beta) / 2;
    for (T lambda : std::vector<T>{P.alpha, mid, P.beta}) {
      // An independent computation, by PolyNorm_Covering.h.
      T g =
          T(1) / PolyNorm_PathCoveringRadius<T, Tint, Tgroup>(data, lambda, os);
      if (g != PolyNorm_PathPieceValue(P, lambda)) {
        std::cerr << "POLYNORM_COVERINGPATH: At lambda=" << lambda << " g=" << g
                  << " but the piece gives "
                  << PolyNorm_PathPieceValue(P, lambda) << "\n";
        throw TerminalException{1};
      }
    }
  }
  std::cerr << "POLYNORM_COVERINGPATH: SANITY_CHECK_POLYNORM_COVERINGPATH_"
               "EXTENSIVE, done\n";
#endif
#ifdef TIMINGS_POLYNORM_COVERINGPATH
  os << "|POLYNORM_COVERINGPATH: covering path|=" << time << "\n";
#endif
  return res;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_COVERINGPATH_H_
// clang-format on
