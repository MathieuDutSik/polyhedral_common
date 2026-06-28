// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_QUANTIZATIONDEFORMATION_H_
#define SRC_DELAUNAY_QUANTIZATIONDEFORMATION_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "InvariantVectorFamily.h"
#include "LatticeStabEquiCan.h"
#include "QuantizationIntegral.h"
#include "polynomial.h"
#include "jet.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
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

#ifdef SANITY_CHECK
#define SANITY_CHECK_QUANTIZATION_DEFORMATION
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
  return UniversalStdVectorMatrixConversion<T,Tint>(gens_i);
}

// The T-space for the deformation is the 2-dimensional space spanned by Q and H
// -- the ray Q + t H lives in it, and that is the only space the deformation
// explores. The common symmetry generators of (Q, H) are recorded as the point
// stabilizer of the T-space (PtStabGens), NOT used to build a larger invariant
// space. Using the full invariant space is wrong here: a deformation running
// along an iso-Delaunay wall (a cell co-spherical only on the ray) would induce
// an equality in that larger space and be (incorrectly) rejected as non-generic,
// whereas in span{Q, H} such a cell's Voronoi regulator vanishes on both Q and H
// (V . Q = V . H = 0), so no equality is induced and the deformation is handled.
// Q is the positive-definite reference matrix (SuperMat) of the space.
template <typename T, typename Tint, typename Tgroup>
LinSpaceMatrix<T> build_qh_tspace(MyMatrix<T> const &Q, MyMatrix<T> const &H,
                                  std::vector<MyMatrix<T>> const &gens_T,
                                  [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{Q, H};
  std::vector<MyMatrix<T>> ListComm;
  LinSpaceMatrix<T> LinSpa = BuildLinSpace<T>(Q, ListMat, ListComm);
  LinSpa.PtStabGens = gens_T;
  return LinSpa;
}

// The Delaunay tesselation of the concrete Gram matrix GramMat (no caching),
// usable as a generic probe of TestGram's iso-Delaunay domain in the T-space
// LinSpa. Two conditions are enforced during the enumeration and a violation of
// either makes it abort (returning nullopt, so find_iso_delaunay_segment halves
// t and retries):
//  - genericity (f_incorrect): no cell may induce an equality in the T-space
//    (IsDelaunayPolytopeInducingEqualities), i.e. GramMat must be a generic
//    interior point of its iso-Delaunay domain, not sitting on a wall. On a wall
//    the defining inequalities are ill-defined and order-dependent.
//  - acceptability (DataLattice::CommonGramMat): f_adj, which has the
//    adjacencies at hand, checks every cell/neighbour pair with
//    IsDelaunayPairAcceptableForGramMat and aborts as soon as an adjacent apex
//    falls inside a circumsphere under TestGram, i.e. TestGram lies in a
//    different domain than GramMat.
// A non-null result therefore certifies that GramMat is generic and that
// TestGram shares its Delaunay tesselation.
template <typename T, typename Tint, typename Tgroup>
std::optional<DelaunayTesselation<T, Tgroup>>
delaunay_for_gram(MyMatrix<T> const &GramMat, LinSpaceMatrix<T> const &LinSpa,
                  MyMatrix<T> const &TestGram, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  data.CommonGramMat = TestGram;
  auto f_incorrect = [&](Delaunay_Obj<T, Tgroup> const &x) -> bool {
    return IsDelaunayPolytopeInducingEqualities(x.EXT, LinSpa, os);
  };
  int max_runtime_second = 0;
  return EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(
      data, f_incorrect, max_runtime_second);
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
    // Enumerate the Delaunay tesselation of the probe G = Q + t H, aborting if
    // Q is not acceptable for some cell, i.e. the segment [Q, Q + t H] leaves
    // the iso-Delaunay domain. A nullopt therefore means "halve and retry".
    std::optional<DelaunayTesselation<T, Tgroup>> opt_DT =
        delaunay_for_gram<T, Tint, Tgroup>(G, LinSpa, Q, os);
    if (!opt_DT) {
#ifdef DEBUG_QUANTIZATION_DEFORMATION
      os << "QDEF: t=" << t << " segment leaves the domain, halving\n";
#endif
      t /= two;
      continue;
    }
    DelaunayTesselation<T, Tgroup> DT = *opt_DT;
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
      if (IsZeroVector(fai.eIneq)) {
        std::cerr << "QDEF: a defining inequality of the iso-Delaunay domain is "
                     "zero in the T-space; with a generic probe in span{Q,H} "
                     "every entry of ListIneq must be a non-zero wall\n";
        throw TerminalException{1};
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
#ifdef SANITY_CHECK_QUANTIZATION_DEFORMATION
    // delaunay_for_gram already aborted if Q was not acceptable for some cell of
    // the probe's tesselation, so Q must lie in the oriented cone. Re-derive that
    // q_inside condition from the (deduplicated) defining inequalities and fail
    // loudly should the per-cell test and the global cone ever disagree.
    for (auto &oe : oriented) {
      if (dotprod(oe, cQ) < 0) {
        std::cerr << "QDEF: SANITY_CHECK failed at t=" << t
                  << ": Q outside the iso-Delaunay domain although "
                     "delaunay_for_gram accepted every cell\n";
        throw TerminalException{1};
      }
    }
#endif
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


// ---------------------------------------------------------------------------
// Phase 2: the quantization as a function of t, by sampling + exact rational
// interpolation. The combinatorics is constant on the open segment, so the
// concrete quantization integral evaluated at several t inside the segment
// gives samples of the same rational function SecMoment(t), which we rebuild
// exactly.
// ---------------------------------------------------------------------------

// The full quantization result (second moment S, second-moment matrix M, ...)
// of the lattice with the given Gram matrix.
template <typename T, typename Tint, typename Tgroup>
QuantizationResult<T> quant_at_gram(MyMatrix<T> const &GramMat,
                                    std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  DataLattice<T, Tint, Tgroup> data =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, os);
  DelaunayTesselation<T, Tgroup> DT =
      get_delaunay_tessellation_serial<T, Tint, Tgroup>(data, "none", 0, os);
  return ComputeQuantizationIntegral<T, Tint, Tgroup>(data, DT, os);
}

// The SecMoment of a concrete Gram matrix (a single sample). The second moment
// integral is the real object, so this just reads it off quant_at_gram.
template <typename T, typename Tint, typename Tgroup>
T secmoment_at_gram(MyMatrix<T> const &GramMat, std::ostream &os) {
  return quant_at_gram<T, Tint, Tgroup>(GramMat, os).SecMoment;
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

template <typename T> struct DeformationDerivatives {
  // Order-2 Taylor jets at t = 0. Use jet_derivative(., k) for the k-th
  // derivative and .coeffs[k] for the k-th Taylor coefficient.
  // SecMoment(t), exact.
  Jet<T> secmoment;
  // det(Q + t H), exact.
  Jet<T> det;
  // The dimensionless normalized second moment G(t) = (1/n) det^(-1/n) S(t),
  // numerical (det^(-1/n) is irrational).
  Jet<double> G;
  int secmoment_degree;
  // The reconstructed SecMoment(t) = num(t) / den(t) (den(0) = 1). When den is
  // the constant 1, SecMoment(t) is a polynomial.
  MyVector<T> secmoment_num;
  MyVector<T> secmoment_den;
};

// Compute the Taylor data of SecMoment(t) at 0 from its rational form and the
// derivatives of the normalized quantizer constant G(t) at 0.
template <typename T>
DeformationDerivatives<T> deformation_derivatives(RationalFunc<T> const &S,
                                                  MyVector<T> const &detpoly,
                                                  int n) {
  // We need up to the second derivative, so order-2 jets in t suffice. Recall a
  // jet's coefficient of t^k is the k-th Taylor coefficient (the k-th derivative
  // divided by k!), hence the factor 2 below for the second derivatives.
  int order = 2;
  DeformationDerivatives<T> res;
  // Exact jet of SecMoment(t) = P(t) / D(t) at t = 0, over T.
  res.secmoment = jet_from_poly(S.P, order) * jet_inverse(jet_from_poly(S.D, order));
  res.det = jet_from_poly(detpoly, order);
  res.secmoment_degree = S.degree;
  res.secmoment_num = S.P;
  res.secmoment_den = S.D;
  // Numerical part: the normalized quantizer G(t) = (1/n) det(t)^(-1/n) S(t),
  // as a double jet (det^(-1/n) is irrational).
  double a = 1.0 / static_cast<double>(n);
  Jet<double> det_jet = UniversalJetConversion<double, T>(res.det);
  Jet<double> S_jet = UniversalJetConversion<double, T>(res.secmoment);
  res.G = jet_scalar_mult(a, jet_pow(det_jet, -a) * S_jet);
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
  LinSpaceMatrix<T> LinSpa = build_qh_tspace<T, Tint, Tgroup>(Q, H, gens_T, os);
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
  // The denominator of SecMoment(t) is exactly det(Q + t H)/det(Q) (degree <= n,
  // constant term 1), so we only interpolate the numerator polynomial. A nullopt
  // therefore cannot mean "wrong denominator" -- it can only mean the numerator
  // degree exceeds max_degree, i.e. max_degree is too small.
  T detQ = detpoly(0);
  MyVector<T> den = detpoly / detQ;
  int max_degree = 4 * n;
  RationalFunc<T> S = reconstruct_rational_known_denominator<T, decltype(sampler)>(
      tpool, sampler, den, max_degree, os);
  return deformation_derivatives<T>(S, detpoly, n);
}

// GAP-readable record of the deformation along Q + t H (exact rationals only;
// the irrational normalized G(t) derivatives are left to the caller).
template <typename T>
void WriteDeformationGAP(std::ostream &os_out,
                        DeformationDerivatives<T> const &der) {
  os_out << "return rec(SecMoment0:=" << jet_derivative(der.secmoment, 0)
         << ", SecMoment1:=" << jet_derivative(der.secmoment, 1)
         << ", SecMoment2:=" << jet_derivative(der.secmoment, 2) << ",\n";
  os_out << "numerator:=" << StringVectorGAP(der.secmoment_num) << ",\n";
  os_out << "denominator:=" << StringVectorGAP(der.secmoment_den) << ",\n";
  os_out << "degree:=" << der.secmoment_degree << ",\n";
  os_out << "det0:=" << der.det.coeffs[0] << ", det1:=" << der.det.coeffs[1]
         << ", det2:=" << der.det.coeffs[2] << ");\n";
}

// ---------------------------------------------------------------------------
// Hessian of the normalized quantizer constant G at the point Q.
//
// G is scale invariant, so it lives on the (n(n+1)/2 - 1)-dimensional space of
// forms modulo scaling. Its Hessian is a quadratic form on Sym^n; we want the
// full symmetric bilinear form Hess(H1,H2) (to read off its signature), not just
// the diagonal directional second derivatives.
//
// WHY THE SECOND-MOMENT *MATRIX* is used (the initially surprising part).
// The naive route takes the scalar second derivative R(H) = G''(0) along Q + tH
// -- that is the diagonal Hess(H,H) -- for many rank-one directions H = v v^T,
// and polarizes. This cannot work: R(v v^T) is a quartic polynomial in v, so the
// rank-one scalar values, for ALL integer v at once, span only the space of
// quartics (dimension C(n+3,4) = 126 for n=6), far short of the
// C(n(n+1)/2 + 1, 2) = 231 coefficients of the Hessian. The reconstruction
// plateaus (the "147 = 126 + 21" wall) and the off-Veronese part of the Hessian
// is simply invisible to rank-one scalar second derivatives.
//
// The way out is to differentiate the second-moment MATRIX itself,
//   M(Q) = \int_{V_0(Q)} u u^T du   (the SecMomentMat of the quantization),
//   DM[H] = d/dt M(Q + t H) |_{t=0}.
// Unlike the scalar R, the map H |-> DM[H] is *linear*, hence determined by its
// values on a basis of Sym^n -- and a basis of rank-one forms v v^T suffices,
// because each DM[v v^T] is a whole matrix (a full ROW of the Hessian), not one
// number. So rank-one directions DO determine everything; there is no quartic
// plateau.
//
// With M = M(Q), S = SecMoment(Q) = tr(Q M) and A_i = Q^{-1} H_i, the Hessian
// (rational; equal to G''(0) up to the positive factor (1/n) det(Q)^{-1/n}, so
// with the same signature) is the polarization
//
//   R(H1,H2) = <H1, DM[H2]>
//              - (1/n)[ (tr A1) tr(M H2) + (tr A2) tr(M H1) ]
//              + (S/n) tr(A1 A2) + (S/n^2)(tr A1)(tr A2).
//
// For rank-one directions everything simplifies: <v v^T, DM> = v^T DM v,
// tr(M v v^T) = v^T M v, tr A = v^T Q^{-1} v, and
// tr(Q^{-1} v v^T Q^{-1} w w^T) = (v^T Q^{-1} w)^2. Aut(Q) acts by v -> g v and
// DM is equivariant, DM[(g v)(g v)^T] = g^{-T} DM[v v^T] g^{-1}, so one
// matrix-derivative evaluation per Aut(Q)-orbit determines DM on the whole orbit
// (for E6/E7 a single orbit of short vectors already spans Sym^n -> one
// evaluation). Each DM[v v^T] is obtained (in compute_moment_derivative) by
// interpolating M(Q + t v v^T) entrywise, each entry having the known
// denominator det(Q + t v v^T)^2.
// ---------------------------------------------------------------------------

// The rational value R(H) = G''(0) / ((1/n) det(Q)^{-1/n}) from the scalar
// deformation data; used only as an independent cross-check of the Hessian.
template <typename T>
T rational_hessian_value(DeformationDerivatives<T> const &der, int n) {
  T S = jet_derivative(der.secmoment, 0);
  T S1 = jet_derivative(der.secmoment, 1);
  T S2 = jet_derivative(der.secmoment, 2);
  T D = jet_derivative(der.det, 0);
  T D1 = jet_derivative(der.det, 1);
  T D2 = jet_derivative(der.det, 2); // D''(0)
  T Tn(n);
  T r1 = D1 / D;
  T r2 = D2 / D;
  return S2 - (T(2) / Tn) * r1 * S1 - (T(1) / Tn) * r2 * S +
         ((Tn + T(1)) / (Tn * Tn)) * r1 * r1 * S;
}

// DM[B] = d/dt M(Q + t B)|_0, the derivative of the second-moment matrix along
// the ray Q + t B. Same machinery as compute_deformation_derivatives (iso-
// Delaunay segment, sampling, exact interpolation) but interpolating each entry
// of the matrix M(t) instead of the scalar SecMoment(t).
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> compute_moment_derivative(MyMatrix<T> const &Q,
                                      MyMatrix<T> const &B, std::ostream &os) {
  int n = Q.rows();
  std::vector<MyMatrix<T>> gens_T =
      compute_qh_symmetry_gens<T, Tint, Tgroup>(Q, B, os);
  LinSpaceMatrix<T> LinSpa = build_qh_tspace<T, Tint, Tgroup>(Q, B, gens_T, os);
  IsoDelaunaySegment<T, Tgroup> seg =
      find_iso_delaunay_segment<T, Tint, Tgroup>(LinSpa, Q, B, T(1), os);
  int pool_size = 4 * n + 30;
  std::vector<T> tpool;
  T half_tmax = seg.tmax / T(2);
  for (int k = 1; k <= pool_size; k++) {
    tpool.push_back(half_tmax * T(k) / T(pool_size + 1));
  }
  MyVector<T> detpoly = det_polynomial<T>(Q, B);
  T detQ = detpoly(0);
  MyVector<T> den = detpoly / detQ; // det(Q+tB)/det(Q), constant term 1
  // Each entry M_ij(t) of the second-moment matrix has denominator
  // (det(Q+tB)/det(Q))^2: the Voronoi-cell vertices (Delaunay circumcenters,
  // solutions of a linear system in Q+tB) carry denominator det(Q+tB), and the
  // moment u_i u_j integrated over the (unit-volume) cell is quadratic in those
  // vertices, hence det^2. Only the contraction tr((Q+tB) M) = SecMoment(t)
  // collapses back to a single det. So we interpolate each entry with the known
  // denominator den^2 = (det(Q+tB)/det(Q))^2. (Verified exactly on A3/E6/A6/D6.)
  int dl = den.size();
  MyVector<T> den2 = ZeroVector<T>(2 * dl - 1);
  for (int a = 0; a < dl; a++) {
    for (int b = 0; b < dl; b++) {
      den2(a + b) += den(a) * den(b);
    }
  }
  int max_degree = 4 * n;
  // Lazy cache of the (expensive) matrix samples M(Q + t B).
  std::map<T, MyMatrix<T>> mcache;
  auto getM = [&](T const &tt) -> MyMatrix<T> const & {
    auto it = mcache.find(tt);
    if (it == mcache.end()) {
      MyMatrix<T> G = Q + tt * B;
      QuantizationResult<T> q = quant_at_gram<T, Tint, Tgroup>(G, os);
      it = mcache.emplace(tt, q.SecMomentMat).first;
    }
    return it->second;
  };
  MyMatrix<T> DM(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      auto sampler = [&](T const &tt) -> T { return getM(tt)(i, j); };
      RationalFunc<T> Sf =
          reconstruct_rational_known_denominator<T, decltype(sampler)>(
              tpool, sampler, den2, max_degree, os);
      T deriv =
          jet_derivative(deformation_derivatives<T>(Sf, detpoly, n).secmoment, 1);
      DM(i, j) = deriv;
      DM(j, i) = deriv;
    }
  }
  return DM;
}

// The orbit of v0 under v -> g v (sign-canonicalized), together with a group
// element g realizing each member: member = sign_canonicalize(g v0).
template <typename Tint>
std::vector<std::pair<MyVector<Tint>, MyMatrix<Tint>>>
orbit_elements(std::vector<MyMatrix<Tint>> const &gens,
               MyVector<Tint> const &v0) {
  int n = v0.size();
  std::vector<std::pair<MyVector<Tint>, MyMatrix<Tint>>> orb;
  std::unordered_set<MyVector<Tint>> seen;
  MyVector<Tint> c0 = SignCanonicalizeVector(v0);
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  orb.push_back({c0, Id});
  seen.insert(c0);
  size_t head = 0;
  while (head < orb.size()) {
    MyVector<Tint> v = orb[head].first;
    MyMatrix<Tint> g = orb[head].second;
    head++;
    for (auto &U : gens) {
      MyVector<Tint> w = U * v;
      MyVector<Tint> cw = SignCanonicalizeVector(w);
      if (seen.insert(cw).second) {
        orb.push_back({cw, U * g});
      }
    }
  }
  return orb;
}

template <typename T> struct HessianResult {
  int n;
  int N;             // dim Sym^n = n(n+1)/2
  bool solved;
  int nbEval;        // number of (expensive) matrix-derivative evaluations
  int nbBasis;       // number of rank-one basis directions (= N when solved)
  MyMatrix<T> beta;  // the N x N Hessian in the v v^T basis (up to + scale)
  T radial_residual; // max |(beta c)| with sum c_k v_k v_k^T = Q (0 if critical)
  T check_residual;  // |R_predicted - R_measured| on an independent direction
  int nbPlus;
  int nbMinus;
  int nbZero;        // includes the one radial (scaling) zero
};

// Compute the Hessian of G at Q and its signature by the moment-derivative
// method. The rank-one basis of Sym^n is built shell by shell of increasing
// dual norm v^T Q^{-1} v (no bound parameter).
template <typename T, typename Tint, typename Tgroup>
HessianResult<T> compute_hessian_signature(MyMatrix<T> const &Q,
                                           std::ostream &os) {
  int n = Q.rows();
  int N = n * (n + 1) / 2;
  MyMatrix<T> Qinv = Inverse(Q);
  os << "QHESS: n=" << n << " dim(Sym^n)=" << N << "\n";
  HessianResult<T> res;
  res.n = n;
  res.N = N;
  res.solved = false;
  res.nbEval = 0;
  res.nbBasis = 0;
  res.beta = ZeroMatrix<T>(N, N);
  res.radial_residual = T(0);
  res.check_residual = T(0);
  res.nbPlus = 0;
  res.nbMinus = 0;
  res.nbZero = 0;
  // M(Q) and S = SecMoment(Q).
  QuantizationResult<T> q0 = quant_at_gram<T, Tint, Tgroup>(Q, os);
  MyMatrix<T> M0 = q0.SecMomentMat;
  T S = q0.SecMoment;
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(Q, os);
  os << "QHESS: |Aut(Q) generators|=" << autom.size() << "\n";
  // Phase A: a rank-one basis of Sym^n, built shell by shell of increasing dual
  // norm v^T Q^{-1} v -- the deformation invariant that controls the cost and
  // degree of the moment derivative, so the cheapest directions come first. We
  // enumerate integer vectors by the integer dual form Qadj = det(Q) Q^{-1} =
  // adj(Q) (whose norm is det(Q) v^T Q^{-1} v) with a CVPSolver, exactly as in
  // ExtractInvariantVectorFamily, keeping v whenever v v^T raises the rank of
  // the family, until the v v^T span Sym^n. No derivatives, no bound parameter.
  T detQ = DeterminantMat(Q);
  MyMatrix<T> Qadj = detQ * Qinv; // adjugate of Q: integral, positive definite
  CVPSolver<T, Tint> solver(Qadj, os);
  T incr = GetSmallestIncrement(Qadj);
  T norm = incr;
  std::vector<MyVector<Tint>> basis;   // the chosen v_k
  std::vector<MyVector<T>> basis_rows; // SymmetricMatrixToVector(v_k v_k^T)
  std::unordered_set<MyVector<Tint>> seen_cand;
  while (static_cast<int>(basis.size()) < N) {
    std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(norm);
    for (auto &v0 : ListVect) {
      if (static_cast<int>(basis.size()) == N) {
        break;
      }
      MyVector<Tint> v = SignCanonicalizeVector(v0);
      if (!seen_cand.insert(v).second) {
        continue;
      }
      MyVector<T> vT = UniversalVectorConversion<T, Tint>(v);
      MyMatrix<T> B = vT * vT.transpose();
      MyVector<T> row = SymmetricMatrixToVector(B); // length N
      int cur = static_cast<int>(basis_rows.size());
      MyMatrix<T> Test(cur + 1, N);
      for (int r = 0; r < cur; r++) {
        Test.row(r) = basis_rows[r].transpose();
      }
      Test.row(cur) = row.transpose();
      if (RankMat(Test) == cur + 1) {
        basis.push_back(v);
        basis_rows.push_back(row);
      }
    }
    norm += incr;
  }
  res.nbBasis = static_cast<int>(basis.size());
  os << "QHESS: rank-one basis size=" << res.nbBasis << "/" << N
     << " (max dual norm " << (norm - incr) << ")\n";
  // Phase B: DM[v_k v_k^T] for each basis vector, one evaluation per Aut(Q)
  // orbit (others by the equivariance DM[(g v)(g v)^T] = g^{-T} DM[v v^T] g^{-1}).
  std::unordered_map<MyVector<Tint>, std::pair<int, MyMatrix<Tint>>> covered;
  std::vector<MyMatrix<T>> DM_reps;
  std::vector<MyMatrix<T>> DM_basis(N);
  for (int k = 0; k < N; k++) {
    MyVector<Tint> v = basis[k];
    auto it = covered.find(v);
    if (it == covered.end()) {
      MyVector<T> vT = UniversalVectorConversion<T, Tint>(v);
      MyMatrix<T> B = vT * vT.transpose();
      MyMatrix<T> DM = compute_moment_derivative<T, Tint, Tgroup>(Q, B, os);
      res.nbEval++;
      int rep_id = static_cast<int>(DM_reps.size());
      DM_reps.push_back(DM);
      for (auto &pr : orbit_elements<Tint>(autom, v)) {
        covered[pr.first] = {rep_id, pr.second};
      }
      os << "QHESS: orbit rep v=" << StringVectorGAP(v)
         << " vTQinvV=" << vT.dot(Qinv * vT) << " (eval " << res.nbEval << ")\n";
      it = covered.find(v);
    }
    int rep_id = it->second.first;
    MyMatrix<Tint> g = it->second.second; // v = sign_canon(g * rep)
    MyMatrix<T> gT = UniversalMatrixConversion<T, Tint>(g);
    MyMatrix<T> ginv = Inverse(gT);
    DM_basis[k] = ginv.transpose() * DM_reps[rep_id] * ginv;
  }
  // Phase C: assemble the Hessian Gram matrix in the v_k v_k^T basis.
  std::vector<MyVector<T>> vbasis(N);
  for (int k = 0; k < N; k++) {
    vbasis[k] = UniversalVectorConversion<T, Tint>(basis[k]);
  }
  T Tn(n);
  MyMatrix<T> beta(N, N);
  for (int k = 0; k < N; k++) {
    T qk = vbasis[k].dot(Qinv * vbasis[k]);
    T mk = vbasis[k].dot(M0 * vbasis[k]);
    for (int l = 0; l < N; l++) {
      T ql = vbasis[l].dot(Qinv * vbasis[l]);
      T ml = vbasis[l].dot(M0 * vbasis[l]);
      T qkl = vbasis[k].dot(Qinv * vbasis[l]);
      T dm_kl = vbasis[k].dot(DM_basis[l] * vbasis[k]); // <B_k, DM[B_l]>
      T val = dm_kl - (T(1) / Tn) * (qk * ml + ql * mk) +
              (S / Tn) * qkl * qkl + (S / (Tn * Tn)) * qk * ql;
      beta(k, l) = val;
    }
  }
  // Symmetrize (the <B_k,DM[B_l]> term is symmetric in exact arithmetic).
  MyMatrix<T> betaS(N, N);
  for (int k = 0; k < N; k++) {
    for (int l = 0; l < N; l++) {
      betaS(k, l) = (beta(k, l) + beta(l, k)) / T(2);
    }
  }
  res.beta = betaS;
  res.solved = true;
  // Radial check: the coordinates c of Q in the v_k v_k^T basis satisfy
  // beta c = Hess(Q, .) = 0 at a critical point. We solve sum_k c_k vec(B_k) =
  // vec(Q) with vec = SymmetricMatrixToVector (BmatT has row k = vec(B_k)).
  MyMatrix<T> BmatT(N, N);
  for (int k = 0; k < N; k++) {
    BmatT.row(k) = basis_rows[k].transpose();
  }
  MyVector<T> qvec = SymmetricMatrixToVector(Q);
  std::optional<MyVector<T>> optc = SolutionMat(BmatT, qvec);
  if (optc) {
    MyVector<T> bc = betaS * (*optc);
    for (int k = 0; k < N; k++) {
      T av = bc(k) < 0 ? T(-bc(k)) : bc(k);
      if (av > res.radial_residual) {
        res.radial_residual = av;
      }
    }
  }
  // Signature.
  DiagSymMat<T> diag = DiagonalizeSymmetricMatrix(betaS, os);
  res.nbPlus = diag.nbPlus;
  res.nbMinus = diag.nbMinus;
  res.nbZero = diag.nbZero;
  int mZero = res.nbZero > 0 ? res.nbZero - 1 : 0;
  os << "QHESS: full Sym^n signature (nbPlus,nbMinus,nbZero)=(" << res.nbPlus
     << "," << res.nbMinus << "," << res.nbZero << ")\n";
  os << "QHESS: Hessian signature on the m=" << (N - 1)
     << " dimensional space (nbPlus,nbMinus,nbZero)=(" << res.nbPlus << ","
     << res.nbMinus << "," << mZero << ")\n";
  os << "QHESS: nbEval=" << res.nbEval
     << " radial_residual=" << res.radial_residual << "\n";
#ifdef SANITY_CHECK_QUANTIZATION_DEFORMATION
  // Independent cross-check: predict R for a held-out direction H = B_0 + B_1
  // (coordinates e_0 + e_1) via beta, and compare with a direct scalar
  // deformation computation along Q + t H. This costs an extra (expensive)
  // deformation evaluation, so it is only done under SANITY_CHECK; an
  // inconsistency is a hard error.
  if (N >= 2) {
    MyMatrix<T> Ba = vbasis[0] * vbasis[0].transpose();
    MyMatrix<T> Bb = vbasis[1] * vbasis[1].transpose();
    MyMatrix<T> H = Ba + Bb;
    T predicted = betaS(0, 0) + T(2) * betaS(0, 1) + betaS(1, 1);
    DeformationDerivatives<T> der =
        compute_deformation_derivatives<T, Tint, Tgroup>(Q, H, os);
    T measured = rational_hessian_value<T>(der, n);
    T diff = predicted - measured;
    res.check_residual = diff < 0 ? T(-diff) : diff;
    if (res.check_residual != T(0)) {
      std::cerr << "QHESS: cross-check FAILED, the assembled Hessian disagrees "
                   "with the direct deformation: predicted="
                << predicted << " measured=" << measured << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return res;
}

template <typename T>
void WriteHessianGAP(std::ostream &os_out, HessianResult<T> const &res) {
  int mZero = res.nbZero > 0 ? res.nbZero - 1 : 0;
  os_out << "return rec(n:=" << res.n << ", dimSpace:=" << res.N
         << ", dimHessian:=" << (res.N - 1) << ",\n";
  os_out << "solved:=" << (res.solved ? "true" : "false")
         << ", nbEval:=" << res.nbEval << ", nbBasis:=" << res.nbBasis
         << ", radialResidual:=" << res.radial_residual
         << ", checkResidual:=" << res.check_residual << ",\n";
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
