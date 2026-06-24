// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_QUANTIZATIONDEFORMATION_H_
#define SRC_DELAUNAY_QUANTIZATIONDEFORMATION_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "LatticeStabEquiCan.h"
#include "QuantizationIntegral.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

// Computation of the quantization constant along the ray Q + t H, with Q a
// positive definite form, H a symmetric form and t >= 0 small.
//
// This file implements the prerequisites (phase 1):
//  1. The common symmetry group of (Q, H): the integral U with U^T Q U = Q and
//     U^T H U = H (compute_qh_symmetry_gens).
//  2. The T-space of the symmetric forms invariant under that group. By
//     construction it contains both Q and H (build_qh_tspace).
//  3. The iso-Delaunay segment: the range t in [0, tmax) for which Q + t H
//     stays positive definite and inside a single iso-Delaunay domain, so that
//     the Delaunay combinatorics is constant along the segment
//     (find_iso_delaunay_segment).
//
// The actual quantization-as-a-function-of-t will (phase 2) sample the existing
// concrete quantization integral at several t inside the segment and rebuild
// the exact rational functions by interpolation.

#ifdef DEBUG
#define DEBUG_QUANTIZATION_DEFORMATION
#endif

// The common symmetry group Stab(Q) cap Stab(H), as a list of matrix
// generators over T. Q must be positive definite (it is used to extract the
// finite vector family).
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>> compute_qh_symmetry_gens(MyMatrix<T> const &Q,
                                                  MyMatrix<T> const &H,
                                                  std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{Q, H};
  std::vector<MyMatrix<Tint>> gens_i =
      ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
  std::vector<MyMatrix<T>> gens_T;
  for (auto &g : gens_i) {
    gens_T.push_back(UniversalMatrixConversion<T, Tint>(g));
  }
  return gens_T;
}

// The T-space of the forms invariant under the symmetry group generators. We
// use Q itself as the positive-definite reference matrix of the space (it is
// invariant under the generators, hence lies in the T-space), so there is no
// need to search for one.
template <typename T, typename Tint, typename Tgroup>
LinSpaceMatrix<T> build_qh_tspace(MyMatrix<T> const &Q,
                                  std::vector<MyMatrix<T>> const &gens_T,
                                  std::ostream &os) {
  int n = Q.rows();
  std::vector<MyMatrix<T>> ListMat = BasisInvariantForm<T>(n, gens_T, os);
  std::vector<MyMatrix<T>> ListComm;
  return BuildLinSpace<T>(Q, ListMat, ListComm);
}

// The Delaunay tesselation of a concrete Gram matrix (no caching).
template <typename T, typename Tint, typename Tgroup>
DelaunayTesselation<T, Tgroup> delaunay_for_gram(MyMatrix<T> const &GramMat,
                                                 std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  return get_delaunay_tessellation_serial<T, Tint, Tgroup>(data, "none", 0, os);
}

template <typename T, typename Tgroup> struct IsoDelaunaySegment {
  // Valid range of t: the segment [0, tmax) (open at tmax, which is the wall).
  T tmax;
  // Whether tmax is a genuine exit of the iso-Delaunay domain (bounded=true) or
  // the ray stays in the domain for all t >= 0 (bounded=false, tmax then a safe
  // interior value).
  bool bounded;
  // The Delaunay combinatorics, constant on the open segment.
  DelaunayTesselation<T, Tgroup> DT;
  // The defining inequalities of the iso-Delaunay domain (T-space coordinates).
  std::vector<FullAdjInfo<T>> ListIneq;
  // The interior form Q + t_probe H whose combinatorics is DT.
  MyMatrix<T> ProbeGram;
};

// Find the iso-Delaunay segment of Q in the direction H. Start at t = t_init
// and halve t while Q + t H is not positive definite, is on a wall of an
// iso-Delaunay domain, or the segment [Q, Q + t H] leaves the domain. Once a
// valid interior probe Q + t H is found, the exact exit point tmax is computed
// from the domain's defining inequalities.
template <typename T, typename Tint, typename Tgroup>
IsoDelaunaySegment<T, Tgroup>
find_iso_delaunay_segment(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &Q,
                          MyMatrix<T> const &H, T const &t_init,
                          std::ostream &os) {
  MyVector<T> cQ = LINSPA_GetVectorOfMatrixExpression(LinSpa, Q);
  MyVector<T> cH = LINSPA_GetVectorOfMatrixExpression(LinSpa, H);
  auto dotprod = [](MyVector<T> const &a, MyVector<T> const &b) -> T {
    T s(0);
    for (int i = 0; i < a.size(); i++) {
      s += a(i) * b(i);
    }
    return s;
  };
  // A defining inequality whose T-space projection is the zero functional does
  // not constrain the segment (this happens for degenerate directions, e.g. H
  // proportional to Q where Q + t H is a pure rescaling): such an inequality
  // must be ignored rather than read as a wall.
  auto is_zero = [](MyVector<T> const &v) -> bool {
    for (int i = 0; i < v.size(); i++) {
      if (v(i) != 0) {
        return false;
      }
    }
    return true;
  };
  T t = t_init;
  T two(2);
  int n_iter = 0;
  int max_iter = 100;
  while (true) {
    n_iter++;
    if (n_iter > max_iter) {
      std::cerr << "QDEF: find_iso_delaunay_segment failed to find a valid "
                   "segment after "
                << max_iter
                << " halvings. The direction H is likely degenerate (e.g. "
                   "proportional to Q) or Q is not generic.\n";
      throw TerminalException{1};
    }
    MyMatrix<T> G = Q + t * H;
    if (!IsPositiveDefinite(G, os)) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: t=" << t << " not positive definite, halving\n";
#endif
      t /= two;
      continue;
    }
    DelaunayTesselation<T, Tgroup> DT =
        delaunay_for_gram<T, Tint, Tgroup>(G, os);
    std::vector<FullAdjInfo<T>> ListIneq =
        ComputeDefiningIneqIsoDelaunayDomain<T, Tgroup>(DT, LinSpa.ListLineMat,
                                                        os);
    MyVector<T> c_probe = cQ + t * cH;
    // Orient every (non-trivial) defining inequality so that the probe, which is
    // strictly interior to its own iso-Delaunay domain by construction,
    // satisfies oe . c_probe > 0. The domain is then the cone {c : oe . c >= 0}.
    // This avoids relying on the sign convention of the returned inequalities.
    std::vector<MyVector<T>> oriented;
    bool probe_generic = true;
    for (auto &fai : ListIneq) {
      if (is_zero(fai.eIneq)) {
        continue;
      }
      T s = dotprod(fai.eIneq, c_probe);
      if (s == 0) {
        probe_generic = false;
        break;
      }
      if (s > 0) {
        oriented.push_back(fai.eIneq);
      } else {
        oriented.push_back(-fai.eIneq);
      }
    }
    if (!probe_generic) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: t=" << t << " probe on a wall (rigid), halving\n";
#endif
      t /= two;
      continue;
    }
    bool q_inside = true;
    for (auto &oe : oriented) {
      if (dotprod(oe, cQ) < 0) {
        q_inside = false;
        break;
      }
    }
    if (!q_inside) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: t=" << t << " segment leaves the domain, halving\n";
#endif
      t /= two;
      continue;
    }
    // The segment [Q, Q + t H] is inside a single iso-Delaunay domain. Compute
    // the exact exit point of the ray from Q in the direction H.
    bool bounded = false;
    T tmax = t;
    for (auto &oe : oriented) {
      T ai = dotprod(oe, cQ);
      T bi = dotprod(oe, cH);
      if (bi < 0) {
        T cand = ai / (-bi);
        if (!bounded || cand < tmax) {
          tmax = cand;
        }
        bounded = true;
      }
    }
    if (!bounded) {
      tmax = t;
    }
    return IsoDelaunaySegment<T, Tgroup>{tmax, bounded, DT, ListIneq, G};
  }
}

// ---------------------------------------------------------------------------
// Equivalence of deformation directions under the symmetry group.
//
// The form Q + t v v^T transforms under a change of basis M as
// M^T (Q + t v v^T) M = Q + t (M^T v)(M^T v)^T, so two vectors v, w give
// isometric deformations iff w = +/- M^T v for some M with M^T Q M = Q.
// ArithmeticAutomorphismGroup returns generators g with g Q g^T = Q, i.e.
// g = M^T; hence the deformation action is v -> g v (g the returned generator),
// NOT v -> g^T v. The invariant v^T Q^{-1} v is preserved by v -> g v (it
// follows from g Q g^T = Q that g^T Q^{-1} g = Q^{-1}) and is a quick necessary
// test. (Contrast with the relevant/lattice vectors in FreeVectors.h, which are
// genuine lattice vectors and transform as c -> g^T c.)
// ---------------------------------------------------------------------------

template <typename Tint>
MyVector<Tint> sign_canonicalize_vector(MyVector<Tint> const &v) {
  for (int i = 0; i < v.size(); i++) {
    if (v(i) != 0) {
      if (v(i) < 0) {
        return MyVector<Tint>(-v);
      }
      return v;
    }
  }
  return v;
}

// The deformation orbit of v0 under v -> g v (and sign), g ranging over the
// generators returned by ArithmeticAutomorphismGroup (g Q g^T = Q).
template <typename Tint>
std::unordered_set<MyVector<Tint>>
orbit_vector_deformation(std::vector<MyMatrix<Tint>> const &gens,
                         MyVector<Tint> const &v0) {
  std::unordered_set<MyVector<Tint>> orb;
  std::vector<MyVector<Tint>> todo;
  MyVector<Tint> c0 = sign_canonicalize_vector(v0);
  orb.insert(c0);
  todo.push_back(c0);
  while (!todo.empty()) {
    MyVector<Tint> v = todo.back();
    todo.pop_back();
    for (auto &U : gens) {
      MyVector<Tint> w = U * v;
      MyVector<Tint> cw = sign_canonicalize_vector(w);
      if (orb.insert(cw).second) {
        todo.push_back(cw);
      }
    }
  }
  return orb;
}

template <typename Tint>
bool vectors_equivalent(std::vector<MyMatrix<Tint>> const &gens,
                        MyVector<Tint> const &v, MyVector<Tint> const &w) {
  std::unordered_set<MyVector<Tint>> orb = orbit_vector_deformation(gens, v);
  return orb.count(sign_canonicalize_vector(w)) > 0;
}

// ---------------------------------------------------------------------------
// Phase 2: the quantization as a function of t, by sampling + exact rational
// interpolation. The combinatorics is constant on the open segment, so the
// concrete quantization integral evaluated at several t inside the segment
// gives samples of the same rational function SecMoment(t), which we rebuild
// exactly.
// ---------------------------------------------------------------------------

// The SecMoment of a concrete Gram matrix (a single sample).
template <typename T, typename Tint, typename Tgroup>
T secmoment_at_gram(MyMatrix<T> const &GramMat, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  DelaunayTesselation<T, Tgroup> DT =
      get_delaunay_tessellation_serial<T, Tint, Tgroup>(data, "none", 0, os);
  QuantizationResult<T> q =
      ComputeQuantizationIntegral<T, Tint, Tgroup>(data, DT, os);
  return q.SecMoment;
}

// Evaluate a polynomial (coefficients c0 + c1 t + ...) at t.
template <typename T> T eval_poly(MyVector<T> const &C, T const &t) {
  T s(0);
  T tp(1);
  for (int i = 0; i < C.size(); i++) {
    s += C(i) * tp;
    tp *= t;
  }
  return s;
}

// Interpolate the polynomial of degree <= deg through the points (t_k, v_k).
template <typename T>
MyVector<T> poly_from_values(std::vector<T> const &ts, std::vector<T> const &vs,
                             int deg) {
  int N = deg + 1;
  MyMatrix<T> Vand(N, N);
  MyVector<T> rhs(N);
  for (int i = 0; i < N; i++) {
    T tp(1);
    for (int j = 0; j < N; j++) {
      Vand(i, j) = tp;
      tp *= ts[i];
    }
    rhs(i) = vs[i];
  }
  MyMatrix<T> VandInv = Inverse(Vand);
  return VandInv * rhs;
}

// det(Q + t H) as an exact polynomial in t (degree <= n).
template <typename T>
MyVector<T> det_polynomial(MyMatrix<T> const &Q, MyMatrix<T> const &H) {
  int n = Q.rows();
  std::vector<T> ts, vs;
  for (int k = 0; k <= n; k++) {
    T tk(k + 1);
    ts.push_back(tk);
    MyMatrix<T> G = Q + tk * H;
    vs.push_back(DeterminantMat(G));
  }
  return poly_from_values(ts, vs, n);
}

// Fit a rational function S(t) = P(t)/D(t) with deg P <= p, deg D <= q and
// D(0) = 1, through the given (t, s) points. Returns (P, D) or nothing if the
// linear system is singular.
template <typename T>
std::optional<std::pair<MyVector<T>, MyVector<T>>>
fit_rational(std::vector<std::pair<T, T>> const &pts, int p, int q) {
  int N = p + 1 + q;
  if (static_cast<int>(pts.size()) < N) {
    return {};
  }
  MyMatrix<T> Eq(N, N);
  MyVector<T> rhs(N);
  for (int i = 0; i < N; i++) {
    T ti = pts[i].first;
    T si = pts[i].second;
    T tp(1);
    for (int j = 0; j <= p; j++) {
      Eq(i, j) = tp;
      tp *= ti;
    }
    T tk = ti;
    for (int k = 1; k <= q; k++) {
      Eq(i, p + k) = -si * tk;
      tk *= ti;
    }
    rhs(i) = si;
  }
  if (RankMat(Eq) < N) {
    return {};
  }
  MyMatrix<T> EqInv = Inverse(Eq);
  MyVector<T> x = EqInv * rhs;
  MyVector<T> P(p + 1);
  MyVector<T> D(q + 1);
  for (int j = 0; j <= p; j++) {
    P(j) = x(j);
  }
  D(0) = T(1);
  for (int k = 1; k <= q; k++) {
    D(k) = x(p + k);
  }
  return std::make_pair(P, D);
}

template <typename T> struct RationalFunc {
  MyVector<T> P;
  MyVector<T> D;
  int degree;
};

// Reconstruct SecMoment(t) = N(t)/den(t) where the denominator den (with
// den(0) = 1) is known exactly (it equals det(Q + t H)/det(Q)) and only the
// numerator polynomial N has to be interpolated. The numerator degree is
// increased until extra held-out samples validate the fit. Returns nothing if
// no polynomial numerator of degree <= max_degree fits (i.e. the denominator is
// not the expected det), so the caller can fall back to a full rational fit.
template <typename T, typename Fsampler>
std::optional<RationalFunc<T>>
reconstruct_secmoment_known_denominator(std::vector<T> const &tpool,
                                        Fsampler sampler, MyVector<T> const &den,
                                        int max_degree, std::ostream &os) {
  std::vector<std::pair<T, T>> cache; // (t, N = SecMoment * den)
  auto get = [&](int i) -> std::pair<T, T> {
    while (static_cast<int>(cache.size()) <= i) {
      T tt = tpool[cache.size()];
      T ss = sampler(tt);
      T nn = ss * eval_poly(den, tt);
      cache.push_back({tt, nn});
    }
    return cache[i];
  };
  int nval = 4;
  for (int d = 0; d <= max_degree; d++) {
    int N = d + 1;
    if (N + nval > static_cast<int>(tpool.size())) {
      break;
    }
    std::vector<T> ts, ns;
    for (int i = 0; i < N; i++) {
      std::pair<T, T> pr = get(i);
      ts.push_back(pr.first);
      ns.push_back(pr.second);
    }
    MyVector<T> P = poly_from_values(ts, ns, d);
    bool ok = true;
    for (int i = N; i < N + nval; i++) {
      std::pair<T, T> pr = get(i);
      if (eval_poly(P, pr.first) != pr.second) {
        ok = false;
        break;
      }
    }
    if (ok) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: numerator (known denominator) converged at degree=" << d
         << "\n";
#endif
      return RationalFunc<T>{P, den, d};
    }
  }
  return {};
}

// Reconstruct SecMoment(t) = P(t)/D(t) from samples produced lazily by sampler,
// at the t-values in tpool. The degree is increased until extra held-out
// samples validate the fit.
template <typename T, typename Fsampler>
RationalFunc<T> reconstruct_secmoment(std::vector<T> const &tpool,
                                      Fsampler sampler, int max_degree,
                                      std::ostream &os) {
  std::vector<std::pair<T, T>> cache;
  auto get = [&](int i) -> std::pair<T, T> {
    while (static_cast<int>(cache.size()) <= i) {
      T tt = tpool[cache.size()];
      T ss = sampler(tt);
      cache.push_back({tt, ss});
    }
    return cache[i];
  };
  int nval = 4;
  for (int d = 1; d <= max_degree; d++) {
    int N = 2 * d + 1;
    if (N + nval > static_cast<int>(tpool.size())) {
      break;
    }
    std::vector<std::pair<T, T>> fitpts;
    for (int i = 0; i < N; i++) {
      fitpts.push_back(get(i));
    }
    std::optional<std::pair<MyVector<T>, MyVector<T>>> opt =
        fit_rational(fitpts, d, d);
    if (!opt) {
      continue;
    }
    MyVector<T> const &P = opt->first;
    MyVector<T> const &D = opt->second;
    bool ok = true;
    for (int i = N; i < N + nval; i++) {
      std::pair<T, T> pr = get(i);
      if (eval_poly(P, pr.first) != pr.second * eval_poly(D, pr.first)) {
        ok = false;
        break;
      }
    }
    if (ok) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: reconstruct_secmoment converged at degree=" << d << "\n";
#endif
      return RationalFunc<T>{P, D, d};
    }
  }
  std::cerr << "QDEF: reconstruct_secmoment failed up to degree=" << max_degree
            << "\n";
  throw TerminalException{1};
}

template <typename T> struct DeformationDerivatives {
  // Exact Taylor data of SecMoment(t) at t = 0: S(0), S'(0), S''(0).
  T S0;
  T S1;
  T S2;
  int secmoment_degree;
  // The reconstructed SecMoment(t) = num(t) / den(t) (den(0) = 1). When den is
  // the constant 1, SecMoment(t) is a polynomial.
  MyVector<T> secmoment_num;
  MyVector<T> secmoment_den;
  // Exact det(Q + t H) = det0 + det1 t + det2 t^2 + ...
  T det0;
  T det1;
  T det2;
  // The dimensionless normalized second moment G(t) = (1/n) det^(-1/n) S(t)
  // and its first two derivatives at 0 (numerical, since det^(-1/n) is
  // irrational).
  double G0;
  double G1;
  double G2;
};

// Compute the Taylor data of SecMoment(t) at 0 from its rational form and the
// derivatives of the normalized quantizer constant G(t) at 0.
template <typename T>
DeformationDerivatives<T> deformation_derivatives(RationalFunc<T> const &S,
                                                  MyVector<T> const &detpoly,
                                                  int n) {
  // Order-2 jets in t (arrays [c0, c1, c2]) over T.
  auto jet_of_poly = [](MyVector<T> const &C) -> std::array<T, 3> {
    std::array<T, 3> j{T(0), T(0), T(0)};
    for (int i = 0; i < C.size() && i < 3; i++) {
      j[i] = C(i);
    }
    return j;
  };
  auto jet_mul = [](std::array<T, 3> const &a,
                    std::array<T, 3> const &b) -> std::array<T, 3> {
    return {a[0] * b[0], a[0] * b[1] + a[1] * b[0],
            a[0] * b[2] + a[1] * b[1] + a[2] * b[0]};
  };
  // Inverse of a jet with non-zero constant term: 1/(d0 + d1 t + d2 t^2).
  auto jet_inv = [](std::array<T, 3> const &d) -> std::array<T, 3> {
    T i0 = T(1) / d[0];
    T i1 = -d[1] * i0 * i0;
    T i2 = (d[1] * d[1] - d[0] * d[2]) * i0 * i0 * i0;
    return {i0, i1, i2};
  };
  std::array<T, 3> Pj = jet_of_poly(S.P);
  std::array<T, 3> Dj = jet_of_poly(S.D);
  std::array<T, 3> Sj = jet_mul(Pj, jet_inv(Dj));
  DeformationDerivatives<T> res;
  res.S0 = Sj[0];
  res.S1 = Sj[1];
  res.S2 = T(2) * Sj[2];
  res.secmoment_degree = S.degree;
  res.secmoment_num = S.P;
  res.secmoment_den = S.D;
  res.det0 = detpoly.size() > 0 ? detpoly(0) : T(0);
  res.det1 = detpoly.size() > 1 ? detpoly(1) : T(0);
  res.det2 = detpoly.size() > 2 ? detpoly(2) : T(0);
  // Numerical part: G(t) = (1/n) det(t)^(-1/n) S(t).
  double a = 1.0 / static_cast<double>(n);
  double c0 = UniversalScalarConversion<double, T>(res.det0);
  double c1 = UniversalScalarConversion<double, T>(res.det1);
  double c2 = UniversalScalarConversion<double, T>(res.det2);
  double s0 = UniversalScalarConversion<double, T>(Sj[0]);
  double s1 = UniversalScalarConversion<double, T>(Sj[1]);
  double s2 = UniversalScalarConversion<double, T>(Sj[2]);
  // jet of det^(-a): c0^(-a) * (1 + r1 t + r2 t^2)^(-a), r1=c1/c0, r2=c2/c0.
  double r1 = c1 / c0;
  double r2 = c2 / c0;
  double base = std::pow(c0, -a);
  double d0 = base;
  double d1 = base * (-a) * r1;
  double d2 = base * ((-a) * r2 + (-a) * (-a - 1.0) / 2.0 * r1 * r1);
  // jet of G = a * det^(-a) * S.
  double g0 = a * (d0 * s0);
  double g1 = a * (d0 * s1 + d1 * s0);
  double g2 = a * (d0 * s2 + d1 * s1 + d2 * s0);
  res.G0 = g0;
  res.G1 = g1;
  res.G2 = 2.0 * g2;
  return res;
}

// Full computation: given Q (positive definite) and H (symmetric), compute the
// Taylor data at t = 0 of SecMoment(t) and of the normalized quantizer G(t)
// along the ray Q + t H.
template <typename T, typename Tint, typename Tgroup>
DeformationDerivatives<T> compute_deformation_derivatives(MyMatrix<T> const &Q,
                                                          MyMatrix<T> const &H,
                                                          std::ostream &os) {
  int n = Q.rows();
  std::vector<MyMatrix<T>> gens_T =
      compute_qh_symmetry_gens<T, Tint, Tgroup>(Q, H, os);
  LinSpaceMatrix<T> LinSpa = build_qh_tspace<T, Tint, Tgroup>(Q, gens_T, os);
  T t_init(1);
  IsoDelaunaySegment<T, Tgroup> seg =
      find_iso_delaunay_segment<T, Tint, Tgroup>(LinSpa, Q, H, t_init, os);
  // Sample strictly inside (0, tmax/2) to stay away from both walls.
  int pool_size = 4 * n + 30;
  std::vector<T> tpool;
  T half_tmax = seg.tmax / T(2);
  for (int k = 1; k <= pool_size; k++) {
    tpool.push_back(half_tmax * T(k) / T(pool_size + 1));
  }
  auto sampler = [&](T const &tt) -> T {
    MyMatrix<T> G = Q + tt * H;
    return secmoment_at_gram<T, Tint, Tgroup>(G, os);
  };
  MyVector<T> detpoly = det_polynomial<T>(Q, H);
  // The denominator of SecMoment(t) is det(Q + t H)/det(Q) (degree <= n, with
  // constant term 1). Interpolate only the numerator polynomial; this halves
  // the number of (expensive) samples. Fall back to a full rational fit if the
  // denominator turns out not to be det.
  T detQ = detpoly(0);
  MyVector<T> den = detpoly / detQ;
  int max_degree = 4 * n;
  std::optional<RationalFunc<T>> opt =
      reconstruct_secmoment_known_denominator<T, decltype(sampler)>(
          tpool, sampler, den, max_degree, os);
  if (opt) {
    return deformation_derivatives<T>(*opt, detpoly, n);
  }
  os << "QDEF: denominator is not det(Q+tH); falling back to full rational "
        "reconstruction\n";
  RationalFunc<T> S =
      reconstruct_secmoment<T, decltype(sampler)>(tpool, sampler, max_degree,
                                                  os);
  return deformation_derivatives<T>(S, detpoly, n);
}

// GAP-readable record of the deformation along Q + t H (exact rationals only;
// the irrational normalized G(t) derivatives are left to the caller).
template <typename T>
void WriteDeformationGAP(std::ostream &os_out,
                        DeformationDerivatives<T> const &der) {
  os_out << "return rec(SecMoment0:=" << der.S0 << ", SecMoment1:=" << der.S1
         << ", SecMoment2:=" << der.S2 << ",\n";
  os_out << "numerator:=" << StringVectorGAP(der.secmoment_num) << ",\n";
  os_out << "denominator:=" << StringVectorGAP(der.secmoment_den) << ",\n";
  os_out << "degree:=" << der.secmoment_degree << ",\n";
  os_out << "det0:=" << der.det0 << ", det1:=" << der.det1
         << ", det2:=" << der.det2 << ");\n";
}

// ---------------------------------------------------------------------------
// Orbit scan: representatives of integer vectors |v_i| <= bound under the
// deformation action of Aut(Q), with the second derivative G''(0) of
// G(Q + t v v^T) for the first K orbits (sorted by the invariant v^T Q^{-1} v).
// ---------------------------------------------------------------------------

template <typename T, typename Tint> struct DeformationOrbitEntry {
  MyVector<Tint> rep;
  long orbit_size;
  T invariant; // v^T Q^{-1} v
  T S0;        // SecMoment(0)
  T S2;        // SecMoment''(0)
  T det0;
  T det1;
  T det2;
  double G2; // G''(0)
};

template <typename T, typename Tint> struct DeformationOrbitResult {
  int total_orbits;
  std::vector<DeformationOrbitEntry<T, Tint>> entries;
};

template <typename T, typename Tint, typename Tgroup>
DeformationOrbitResult<T, Tint>
compute_deformation_orbits(MyMatrix<T> const &Q, int bound, int K,
                           std::ostream &os) {
  int n = Q.rows();
  MyMatrix<T> Qinv = Inverse(Q);
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(Q, os);
  os << "QORB: n=" << n << " |Aut(Q) generators|=" << autom.size() << "\n";
  // Enumerate integer vectors with |v_i| <= bound (antipodally reduced).
  std::vector<MyVector<Tint>> candidates;
  std::unordered_set<MyVector<Tint>> seen_cand;
  MyVector<Tint> v = ZeroVector<Tint>(n);
  int range = 2 * bound + 1;
  long total = 1;
  for (int i = 0; i < n; i++) {
    total *= range;
  }
  for (long code = 0; code < total; code++) {
    long c = code;
    bool is_zero = true;
    for (int i = 0; i < n; i++) {
      int digit = c % range;
      c /= range;
      v(i) = digit - bound;
      if (v(i) != 0) {
        is_zero = false;
      }
    }
    if (is_zero) {
      continue;
    }
    MyVector<Tint> cv = sign_canonicalize_vector(v);
    if (seen_cand.insert(cv).second) {
      candidates.push_back(cv);
    }
  }
  // Classify the candidates into orbits.
  std::unordered_set<MyVector<Tint>> assigned;
  struct OrbitRec {
    MyVector<Tint> rep;
    long orbit_size;
    T invariant;
  };
  std::vector<OrbitRec> orbits;
  for (auto &cand : candidates) {
    if (assigned.count(cand) > 0) {
      continue;
    }
    std::unordered_set<MyVector<Tint>> orb =
        orbit_vector_deformation(autom, cand);
    for (auto &w : orb) {
      assigned.insert(w);
    }
    MyVector<T> cand_T = UniversalVectorConversion<T, Tint>(cand);
    T invariant = cand_T.dot(Qinv * cand_T);
    orbits.push_back(OrbitRec{cand, static_cast<long>(orb.size()), invariant});
  }
  std::sort(orbits.begin(), orbits.end(),
            [](OrbitRec const &a, OrbitRec const &b) -> bool {
              return a.invariant < b.invariant;
            });
  os << "QORB: number of vector orbits (|v_i|<=" << bound
     << ")=" << orbits.size() << "\n";
  DeformationOrbitResult<T, Tint> res;
  res.total_orbits = orbits.size();
  int nb = std::min<int>(K, orbits.size());
  os << "QORB: computing G''(0) for the first " << nb << " orbits\n";
  for (int i = 0; i < nb; i++) {
    OrbitRec const &orb = orbits[i];
    MyVector<T> v_T = UniversalVectorConversion<T, Tint>(orb.rep);
    MyMatrix<T> H = v_T * v_T.transpose();
    DeformationDerivatives<T> der =
        compute_deformation_derivatives<T, Tint, Tgroup>(Q, H, os);
    os << "QORB: orbit " << i << " v=" << StringVectorGAP(orb.rep)
       << " |orbit|=" << orb.orbit_size << " vTQinvV=" << orb.invariant << "\n";
    os << "QORB:   SecMoment''(0)=" << der.S2 << " G''(0)=" << der.G2 << "\n";
    res.entries.push_back(DeformationOrbitEntry<T, Tint>{
        orb.rep, orb.orbit_size, orb.invariant, der.S0, der.S2, der.det0,
        der.det1, der.det2, der.G2});
  }
  return res;
}

template <typename T, typename Tint>
void WriteDeformationOrbitsGAP(std::ostream &os_out,
                              DeformationOrbitResult<T, Tint> const &res) {
  os_out << "return rec(nbVectorOrbit:=" << res.total_orbits << ",\n";
  os_out << "nbEvaluated:=" << res.entries.size() << ",\n";
  os_out << "ListOrbit:=[";
  bool IsFirst = true;
  for (auto &e : res.entries) {
    if (!IsFirst) {
      os_out << ",\n";
    }
    IsFirst = false;
    os_out << "rec(v:=" << StringVectorGAP(e.rep)
           << ", OrbitSize:=" << e.orbit_size << ", vTQinvV:=" << e.invariant
           << ", SecMoment0:=" << e.S0 << ", SecMoment2:=" << e.S2
           << ", det0:=" << e.det0 << ", det1:=" << e.det1
           << ", det2:=" << e.det2 << ", Gpp:=" << e.G2 << ")";
  }
  os_out << "]);\n";
}

// ---------------------------------------------------------------------------
// Hessian of the normalized quantizer constant G at the point Q.
//
// G is scale invariant, so it really lives on the (n(n+1)/2 - 1)-dimensional
// space of forms modulo scaling. We assemble its Hessian as a quadratic form on
// Sym^n (the symmetric n x n matrices): for a direction H, the second
// directional derivative G''(0) along Q + t H equals Hess(H, H) (the path is a
// straight line, so there is no first-order correction). Factoring out the
// universal positive constant (1/n) det(Q)^{-1/n} leaves the *rational*
// quadratic form
//   R(H) = S'' - (2/n)(D'/D) S' - (1/n)(D''/D) S + ((n+1)/n^2)(D'/D)^2 S
// (S, S', S'' the SecMoment Taylor data; D, D', D'' those of det(Q + t H)),
// which has the same signature as Hess. The group Aut(Q) acts on directions by
// H -> g H g^T (g Q g^T = Q) and R is invariant, so a single evaluation R(D)
// yields one linear equation in the unknown coefficients beta_{bb'} for every
// element of the orbit of D. We accumulate equations -- from the orbits of the
// rank-one directions v v^T (v short integer vectors), and, when those do not
// suffice (rank-one directions alone cannot resolve forms that vanish on the
// Veronese cone), from a few rank-two seed directions -- together with the
// radial relations Hess(Q, .) = 0, until the system has full rank. We then
// solve it exactly and read off the signature.
// ---------------------------------------------------------------------------

// The rational value R(H) = G''(0) / ((1/n) det(Q)^{-1/n}).
//
// Beware the differing conventions inside DeformationDerivatives: S0/S1/S2 are
// the *derivatives* S(0), S'(0), S''(0), whereas det0/det1/det2 are the *Taylor
// coefficients* of det(Q + t H) (so det1 = D'(0) but det2 = D''(0)/2). The
// second derivative is therefore D''(0) = 2*det2. (For rank-one directions
// det(Q + t v v^T) is linear in t by the matrix-determinant lemma, so det2 = 0
// and this factor is invisible; it only matters for higher-rank directions.)
template <typename T>
T rational_hessian_value(DeformationDerivatives<T> const &der, int n) {
  T S = der.S0;    // S(0)
  T S1 = der.S1;   // S'(0)
  T S2 = der.S2;   // S''(0)
  T D = der.det0;  // D(0)
  T D1 = der.det1; // D'(0)
  T D2 = T(2) * der.det2; // D''(0) = 2 * [t^2] det(Q + t H)
  T Tn(n);
  T ratio1 = D1 / D; // D'/D
  T ratio2 = D2 / D; // D''/D
  T R = S2 - (T(2) / Tn) * ratio1 * S1 - (T(1) / Tn) * ratio2 * S +
        ((Tn + T(1)) / (Tn * Tn)) * ratio1 * ratio1 * S;
  return R;
}

// The standard basis of Sym^n: pairs (i, j) with i <= j, the basis element
// being E_ii when i == j and E_ij + E_ji when i < j.
inline std::vector<std::pair<int, int>> sym_basis_pairs(int n) {
  std::vector<std::pair<int, int>> pairs;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      pairs.push_back({i, j});
    }
  }
  return pairs;
}

// Coordinates of a symmetric matrix D in the basis sym_basis_pairs(n): the
// coefficient on (i, j) is simply D(i, j).
template <typename T>
MyVector<T> sym_coords(MyMatrix<T> const &D,
                       std::vector<std::pair<int, int>> const &pairs) {
  int N = pairs.size();
  MyVector<T> c(N);
  for (int b = 0; b < N; b++) {
    c(b) = D(pairs[b].first, pairs[b].second);
  }
  return c;
}

// A bounded set of elements of the group generated by gens (BFS from the
// identity, capped). Any subset of genuine group elements gives valid equations
// g D g^T, R(g D g^T) = R(D); we only need enough images to reach full rank.
template <typename Tint>
std::vector<MyMatrix<Tint>>
generate_group_elements(std::vector<MyMatrix<Tint>> const &gens, int n,
                        size_t cap) {
  std::vector<MyMatrix<Tint>> elts;
  std::unordered_set<MyVector<Tint>> seen;
  auto key = [&](MyMatrix<Tint> const &M) -> MyVector<Tint> {
    MyVector<Tint> v(n * n);
    int idx = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v(idx++) = M(i, j);
      }
    }
    return v;
  };
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  elts.push_back(Id);
  seen.insert(key(Id));
  size_t head = 0;
  while (head < elts.size() && elts.size() < cap) {
    MyMatrix<Tint> cur = elts[head++];
    for (auto &g : gens) {
      MyMatrix<Tint> nx = g * cur;
      if (seen.insert(key(nx)).second) {
        elts.push_back(nx);
        if (elts.size() >= cap) {
          break;
        }
      }
    }
  }
  return elts;
}

template <typename T> struct HessianResult {
  int n;
  int N;             // dim Sym^n = n(n+1)/2
  bool solved;       // the linear system had full rank and a solution
  bool with_radial;  // the Hess(Q,.)=0 rows were consistent and used
  int nbEval;        // number of (expensive) R-evaluations
  int rank_reached;  // rank of the assembled system
  MyMatrix<T> beta;  // the N x N Hessian (rational, up to a positive scale)
  T radial_residual; // max |(beta q)_b|, should be 0 at a critical point
  T data_residual;   // max |c^T beta c - R| over all evaluated directions
  int nbPlus;        // signature of the full N x N Hessian
  int nbMinus;
  int nbZero;        // includes the one radial (scaling) zero
};

// Compute the Hessian of G at Q and its signature, by the orbit-equation method
// described above. bound controls the box |v_i| <= bound of integer vectors
// whose rank-one directions seed the equations.
template <typename T, typename Tint, typename Tgroup>
HessianResult<T> compute_hessian_signature(MyMatrix<T> const &Q, int bound,
                                           std::ostream &os) {
  int n = Q.rows();
  MyMatrix<T> Qinv = Inverse(Q);
  std::vector<std::pair<int, int>> pairs = sym_basis_pairs(n);
  int N = pairs.size();
  // Unknown index table: symmetric pairs (b, b'), b <= b'.
  std::vector<std::vector<int>> uidx(N, std::vector<int>(N, -1));
  int Munk = 0;
  for (int b = 0; b < N; b++) {
    for (int bp = b; bp < N; bp++) {
      uidx[b][bp] = Munk;
      uidx[bp][b] = Munk;
      Munk++;
    }
  }
  os << "QHESS: n=" << n << " dim(Sym^n)=" << N << " unknowns=" << Munk << "\n";
  MyVector<T> qc = sym_coords<T>(Q, pairs);
  // The data equation contributed by a direction with coordinates c:
  // sum_{b<=b'} coef * beta_{bb'} = R, coef = c_b^2 (b==b') or 2 c_b c_b'.
  auto row_from_coords = [&](MyVector<T> const &c) -> MyVector<T> {
    MyVector<T> row = ZeroVector<T>(Munk);
    for (int b = 0; b < N; b++) {
      if (c(b) == 0) {
        continue;
      }
      row(uidx[b][b]) += c(b) * c(b);
      for (int bp = b + 1; bp < N; bp++) {
        row(uidx[b][bp]) += T(2) * c(b) * c(bp);
      }
    }
    return row;
  };
  // Radial rows: Hess(Q, H_b) = sum_{b'} q_{b'} beta_{b b'} = 0.
  std::vector<MyVector<T>> radial_rows;
  for (int b0 = 0; b0 < N; b0++) {
    MyVector<T> row = ZeroVector<T>(Munk);
    for (int bp = 0; bp < N; bp++) {
      row(uidx[b0][bp]) += qc(bp);
    }
    radial_rows.push_back(row);
  }
  // We keep only a linearly independent set of equations (at most Munk rows),
  // so the working matrix never grows without bound. A batch of candidate rows
  // is appended and the set is re-reduced to independent rows via
  // TMat_SelectRowCol; the radial rows are seeded first (rhs 0) and kept.
  std::vector<MyVector<T>> indep_rows;
  std::vector<T> indep_rhs;
  auto reduce_indep = [&]() -> void {
    int nb = static_cast<int>(indep_rows.size());
    MyMatrix<T> A(nb, Munk);
    for (int r = 0; r < nb; r++) {
      A.row(r) = indep_rows[r].transpose();
    }
    SelectionRowCol<T> sel = TMat_SelectRowCol(A);
    std::vector<MyVector<T>> nr;
    std::vector<T> nrhs;
    for (int idx : sel.ListRowSelect) {
      nr.push_back(indep_rows[idx]);
      nrhs.push_back(indep_rhs[idx]);
    }
    indep_rows = std::move(nr);
    indep_rhs = std::move(nrhs);
  };
  for (auto &rr : radial_rows) {
    indep_rows.push_back(rr);
    indep_rhs.push_back(T(0));
  }
  reduce_indep();
  // Append a batch of (row, value) candidates, keep an independent set, and
  // return the resulting rank.
  auto add_batch =
      [&](std::vector<std::pair<MyVector<T>, T>> const &cand) -> int {
    for (auto &pr : cand) {
      indep_rows.push_back(pr.first);
      indep_rhs.push_back(pr.second);
    }
    reduce_indep();
    return static_cast<int>(indep_rows.size());
  };
  // The evaluated directions (coordinates, R), replayed at the end to check
  // that the solved Hessian reproduces every measurement.
  std::vector<std::pair<MyVector<T>, T>> evals;
  // Aut(Q) generators and the box of integer vectors, classified into orbits.
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(Q, os);
  os << "QHESS: |Aut(Q) generators|=" << autom.size() << "\n";
  std::vector<MyVector<Tint>> candidates;
  std::unordered_set<MyVector<Tint>> seen_cand;
  int range = 2 * bound + 1;
  long total = 1;
  for (int i = 0; i < n; i++) {
    total *= range;
  }
  for (long code = 0; code < total; code++) {
    long c = code;
    bool is_zero = true;
    MyVector<Tint> v(n);
    for (int i = 0; i < n; i++) {
      int digit = c % range;
      c /= range;
      v(i) = digit - bound;
      if (v(i) != 0) {
        is_zero = false;
      }
    }
    if (is_zero) {
      continue;
    }
    MyVector<Tint> cv = sign_canonicalize_vector(v);
    if (seen_cand.insert(cv).second) {
      candidates.push_back(cv);
    }
  }
  std::unordered_set<MyVector<Tint>> assigned;
  struct OrbitRec {
    MyVector<Tint> rep;
    std::vector<MyVector<Tint>> members;
    T invariant;
  };
  std::vector<OrbitRec> orbits;
  for (auto &cand : candidates) {
    if (assigned.count(cand) > 0) {
      continue;
    }
    std::unordered_set<MyVector<Tint>> orb =
        orbit_vector_deformation(autom, cand);
    std::vector<MyVector<Tint>> members;
    for (auto &w : orb) {
      assigned.insert(w);
      members.push_back(w);
    }
    MyVector<T> cand_T = UniversalVectorConversion<T, Tint>(cand);
    T invariant = cand_T.dot(Qinv * cand_T);
    orbits.push_back(OrbitRec{cand, members, invariant});
  }
  std::sort(orbits.begin(), orbits.end(),
            [](OrbitRec const &a, OrbitRec const &b) -> bool {
              return a.invariant < b.invariant;
            });
  os << "QHESS: " << orbits.size() << " rank-one orbits (|v_i|<=" << bound
     << ")\n";
  int nbEval = 0;
  bool full = false;
  // Phase 1: rank-one orbit seeds, evaluated by increasing invariant. Rank-one
  // directions can only ever resolve the "quartic" part of the Hessian, so the
  // rank saturates below Munk; we stop after a few consecutive orbits add
  // nothing (plateau) rather than evaluating every orbit.
  int prev_rank = static_cast<int>(indep_rows.size());
  int no_gain = 0;
  int const plateau_limit = 3;
  for (auto &orb : orbits) {
    MyVector<T> v_T = UniversalVectorConversion<T, Tint>(orb.rep);
    MyMatrix<T> H = v_T * v_T.transpose();
    DeformationDerivatives<T> der =
        compute_deformation_derivatives<T, Tint, Tgroup>(Q, H, os);
    T R = rational_hessian_value<T>(der, n);
    nbEval++;
    evals.push_back({sym_coords<T>(H, pairs), R});
    std::vector<std::pair<MyVector<T>, T>> cand;
    for (auto &w : orb.members) {
      MyVector<T> w_T = UniversalVectorConversion<T, Tint>(w);
      MyMatrix<T> Dw = w_T * w_T.transpose();
      cand.push_back({row_from_coords(sym_coords<T>(Dw, pairs)), R});
    }
    int rk = add_batch(cand);
    os << "QHESS: rank-one orbit vTQinvV=" << orb.invariant
       << " |orbit|=" << orb.members.size() << " R=" << R << " rank=" << rk
       << "/" << Munk << "\n";
    if (rk == Munk) {
      full = true;
      break;
    }
    if (rk == prev_rank) {
      no_gain++;
    } else {
      no_gain = 0;
      prev_rank = rk;
    }
    if (no_gain >= plateau_limit) {
      os << "QHESS: rank-one plateau at rank=" << rk << "\n";
      break;
    }
  }
  // Phase 2: rank-two seeds E_ij + E_ji, expanded by the group, when rank-one
  // directions cannot resolve the whole quadratic form.
  if (!full) {
    os << "QHESS: rank-one seeds insufficient; adding rank-two seeds\n";
    std::vector<MyMatrix<Tint>> Gelts =
        generate_group_elements<Tint>(autom, n, 800);
    os << "QHESS: |group elements (capped)|=" << Gelts.size() << "\n";
    for (int i = 0; i < n && !full; i++) {
      for (int j = i + 1; j < n && !full; j++) {
        MyMatrix<T> D = ZeroMatrix<T>(n, n);
        D(i, j) = 1;
        D(j, i) = 1;
        DeformationDerivatives<T> der =
            compute_deformation_derivatives<T, Tint, Tgroup>(Q, D, os);
        T R = rational_hessian_value<T>(der, n);
        nbEval++;
        evals.push_back({sym_coords<T>(D, pairs), R});
        std::vector<std::pair<MyVector<T>, T>> cand;
        for (auto &g : Gelts) {
          MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
          MyMatrix<T> Dp = g_T * D * g_T.transpose();
          cand.push_back({row_from_coords(sym_coords<T>(Dp, pairs)), R});
        }
        int rk = add_batch(cand);
        os << "QHESS: rank-two seed (" << i << "," << j << ") R=" << R
           << " rank=" << rk << "/" << Munk << "\n";
        if (rk == Munk) {
          full = true;
        }
      }
    }
  }
  HessianResult<T> res;
  res.n = n;
  res.N = N;
  res.nbEval = nbEval;
  res.rank_reached = static_cast<int>(indep_rows.size());
  res.solved = false;
  res.with_radial = true;
  res.beta = ZeroMatrix<T>(N, N);
  res.radial_residual = T(0);
  res.data_residual = T(0);
  res.nbPlus = 0;
  res.nbMinus = 0;
  res.nbZero = 0;
  if (!full) {
    os << "QHESS: WARNING system underdetermined, rank=" << res.rank_reached
       << "/" << Munk << " after " << nbEval << " evaluations\n";
    return res;
  }
  // The independent set is now a square full-rank system (radial rows included),
  // so it has a unique exact solution.
  MyMatrix<T> A(Munk, Munk);
  MyVector<T> b(Munk);
  for (int r = 0; r < Munk; r++) {
    A.row(r) = indep_rows[r].transpose();
    b(r) = indep_rhs[r];
  }
  MyMatrix<T> AT = A.transpose();
  std::optional<MyVector<T>> opt = SolutionMat(AT, b);
  if (!opt) {
    os << "QHESS: WARNING no exact solution found\n";
    return res;
  }
  MyVector<T> u = *opt;
  for (int b = 0; b < N; b++) {
    for (int bp = b; bp < N; bp++) {
      res.beta(b, bp) = u(uidx[b][bp]);
      res.beta(bp, b) = u(uidx[b][bp]);
    }
  }
  // Radial residual beta q (exactly 0 at a critical point).
  MyVector<T> bq = res.beta * qc;
  for (int b = 0; b < N; b++) {
    T av = bq(b) < 0 ? T(-bq(b)) : bq(b);
    if (av > res.radial_residual) {
      res.radial_residual = av;
    }
  }
  // Data residual: the solved Hessian must reproduce every measured R(D).
  for (auto &pr : evals) {
    MyVector<T> c = pr.first;
    T val = c.dot(res.beta * c);
    T av = val - pr.second;
    if (av < 0) {
      av = -av;
    }
    if (av > res.data_residual) {
      res.data_residual = av;
    }
  }
  if (res.data_residual != T(0) || res.radial_residual != T(0)) {
    os << "QHESS: WARNING residuals nonzero (Q not critical?): radial="
       << res.radial_residual << " data=" << res.data_residual << "\n";
  }
  DiagSymMat<T> diag = DiagonalizeSymmetricMatrix(res.beta, os);
  res.nbPlus = diag.nbPlus;
  res.nbMinus = diag.nbMinus;
  res.nbZero = diag.nbZero;
  res.solved = true;
  int mZero = res.nbZero > 0 ? res.nbZero - 1 : 0;
  os << "QHESS: full Sym^n signature (nbPlus,nbMinus,nbZero)=(" << res.nbPlus
     << "," << res.nbMinus << "," << res.nbZero << ")\n";
  os << "QHESS: Hessian signature on the m=" << (N - 1)
     << " dimensional space (nbPlus,nbMinus,nbZero)=(" << res.nbPlus << ","
     << res.nbMinus << "," << mZero << ")\n";
  os << "QHESS: nbEval=" << nbEval << " radial_residual=" << res.radial_residual
     << " data_residual=" << res.data_residual << "\n";
  return res;
}

template <typename T>
void WriteHessianGAP(std::ostream &os_out, HessianResult<T> const &res) {
  int mZero = res.nbZero > 0 ? res.nbZero - 1 : 0;
  os_out << "return rec(n:=" << res.n << ", dimSpace:=" << res.N
         << ", dimHessian:=" << (res.N - 1) << ",\n";
  os_out << "solved:=" << (res.solved ? "true" : "false")
         << ", usedRadial:=" << (res.with_radial ? "true" : "false")
         << ", nbEval:=" << res.nbEval << ", rank:=" << res.rank_reached
         << ", radialResidual:=" << res.radial_residual
         << ", dataResidual:=" << res.data_residual << ",\n";
  os_out << "signatureFull:=rec(nbPlus:=" << res.nbPlus
         << ", nbMinus:=" << res.nbMinus << ", nbZero:=" << res.nbZero
         << "),\n";
  os_out << "signature:=rec(nbPlus:=" << res.nbPlus
         << ", nbMinus:=" << res.nbMinus << ", nbZero:=" << mZero << "),\n";
  os_out << "Hessian:=" << StringMatrixGAP(res.beta) << ");\n";
}

// clang-format off
#endif  // SRC_DELAUNAY_QUANTIZATIONDEFORMATION_H_
// clang-format on
