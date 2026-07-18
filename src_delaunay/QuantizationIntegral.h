// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_QUANTIZATIONINTEGRAL_H_
#define SRC_DELAUNAY_QUANTIZATIONINTEGRAL_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "POLY_Fundamental.h"
#include "POLY_lrslib.h"
#include "InvariantVectorFamily.h"
#include "MatrixGroupAverage.h"
#include "jet_number.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <utility>
#include <vector>
// clang-format on

// Port of the GAP function QuantizationIntegral (and its dependencies
// SoftComputation, ListOrbitContainingEXT, SingleListOrbit, OrbitBarycenter,
// OrbitBarycenterSymmetricMatrix, DirectIntegral, ...) from
// GAPpackages/MyPolyhedral/lib/LatticeDelaunays.g
//
// It computes the second moment of the Voronoi cell of a lattice (the
// "quantization integral"). The cell volume is normalized to 1 (the
// computation works in lattice coordinates, so the fundamental cell has
// volume 1).
//
// All the heavy machinery already exists in the C++ codebase:
//  * lattice stabilizer/equivalence with a marked center, via the weight
//    matrix + integral group machinery (mirrors DataLattice.FuncAutom and
//    DataLattice.FuncEquiv).
//  * group theory on the Delaunay vertex permutation groups (Stabilizer_OnSets,
//    RepresentativeAction_OnSets, OrbitFace, OrbitSplittingSet).
//  * the permutation -> affine matrix homomorphism (RepresentVertexPermutation,
//    the GAP PhiPermMat).
//  * triangulation for the direct integral (lrs::GetTriangulation).

#ifdef DEBUG
#define DEBUG_QUANTIZATION_INTEGRAL
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_QUANTIZATION_INTEGRAL
#endif

// n! as a value of the (exact) scalar type T.
template <typename T> T factorial(int const &n) {
  T ret(1);
  for (int i = 2; i <= n; i++)
    ret *= T(i);
  return ret;
}

template <typename T> struct QuantizationResult {
  T SecMoment;
  MyMatrix<T> SecMomentMat;
  T TheVolume;
  MyVector<T> TheBarycenter;
  double NormalizedSecondMoment;
};

// Evaluate a scalar at a concrete value t0 of the deformation parameter. A plain
// scalar is t-independent and returns itself (t0 ignored); a jet collapses to the
// base field via its Horner evaluation. This is what lets direct_integral pick a
// non-degenerate rational point t0 at which to compute the (combinatorial, hence
// t-independent) leaf triangulation, and orient each simplex volume, while the
// integration itself keeps the full jet expansion.
template <typename T> inline T eval_at(T const &x, T const &) { return x; }
template <typename Tb, int N>
inline Tb eval_at(jet<Tb, N> const &x, Tb const &t0) {
  return x.eval(t0);
}

// Valuation (order of the first non-zero Taylor coefficient): 0 for an ordinary
// non-zero scalar, leading_index for a jet. Used only to profile how many leaf
// simplices are degenerate at t = 0 (a simplex volume of valuation d contributes
// to the moment only at order >= d).
template <typename T> inline int scalar_valuation(T const &x) {
  return (x == T(0)) ? 1000 : 0;
}
template <typename Tb, int N> inline int scalar_valuation(jet<Tb, N> const &x) {
  return x.leading_index();
}
#ifdef PROFILE_SIMPLEX_VALUATION
inline long g_simplex_valuation[32] = {0};
inline long g_simplex_total = 0;
#endif

// The computer needs only three things beyond the Delaunay tessellation: the
// dimension n, the invariant vector family SHV (integer vectors, used by the
// marked-center stabilizer/equivalence), and the metric GramMat. It does not
// otherwise depend on the DataLattice (in particular not on its CVP solver), so
// it is templated on a single scalar type T: instantiating it with T = jet<...>
// re-evaluates the whole integral along a ray Q + t H (the tessellation being
// constant along an iso-Delaunay segment), with SHV / DT / GramMat supplied as
// the converted t = 0 data.
template <typename T, typename Tint, typename Tgroup>
struct QuantizationComputer {
  using Telt = typename Tgroup::Telt;
  using TintGrp = typename Tgroup::Tint;

  DelaunayTesselation<T, Tgroup> const &DT;
  std::ostream &os;
  int n;
  MyMatrix<T> SHV;
  MyMatrix<T> GramMat;
  // Inverse of the (fixed) metric, cached: get_basis needs G^{-1} once per cell,
  // but the metric never changes during the recursion.
  MyMatrix<T> GramMatInv;
  // The base field underneath T (T itself for an ordinary scalar, the coefficient
  // type for a jet); the type in which the leaf triangulation is computed.
  using Tscal = std::decay_t<decltype(constant_term(std::declval<T const &>()))>;
  // The deformation value t0 at which the leaf triangulation is evaluated and the
  // simplices oriented. Only meaningful when T is a jet; for a plain scalar it is
  // ignored by eval_at. For jets it must be an interior point of the iso-Delaunay
  // segment (t0 = tmax/2, where the tessellation was pinned down), so that the
  // evaluated polytope is non-degenerate and its triangulation is the t -> 0+ one.
  Tscal t0_triang{};
  // The metric evaluated at t0 (see eval_mat_ref), wrapped in the single-element
  // ListMat the weight-matrix equivalence/automorphism engine expects. Constant
  // during the recursion, so built once in the constructor.
  std::vector<MyMatrix<Tscal>> ListMat_ref;

  // Ordinary constructor: take n and the metric from a DataLattice, and compute
  // the Z-spanning invariant family ExtractInvariantVectorFamilyZbasis(GramMat)
  // -- required for the automorphism group to be the exact lattice stabilizer
  // without the integral refinement. This (LLL-based, integer) computation lives
  // in the constructor rather than in compute() so that compute() never
  // instantiates it: the deformation instantiation (T = jet) uses the direct
  // constructor below and supplies the converted t = 0 family instead.
  QuantizationComputer(DataLattice<T, Tint, Tgroup> &data,
                       DelaunayTesselation<T, Tgroup> const &DT,
                       std::ostream &os)
      : DT(DT), os(os), n(data.n), GramMat(data.solver.GramMat) {
    MyMatrix<Tint> SHV_i =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
    SHV = UniversalMatrixConversion<T, Tint>(SHV_i);
    GramMatInv = Inverse(GramMat);
    ListMat_ref = {eval_mat_ref(GramMat)};
  }

  // Direct constructor: explicit n, SHV, metric. Used for the deformation
  // instantiation (T = jet), where SHV and GramMat are the converted t = 0
  // family and the ray metric Q + t H.
  QuantizationComputer(int n_in, MyMatrix<T> const &SHV_in,
                       MyMatrix<T> const &GramMat_in,
                       DelaunayTesselation<T, Tgroup> const &DT,
                       std::ostream &os, Tscal const &t0_in = Tscal())
      : DT(DT), os(os), n(n_in), SHV(SHV_in), GramMat(GramMat_in),
        t0_triang(t0_in) {
    GramMatInv = Inverse(GramMat);
    ListMat_ref = {eval_mat_ref(GramMat)};
  }

  // -------- small geometric helpers (homogeneous coordinates) --------

  // The center (and radius) of the circumscribing sphere, homogeneous, length
  // n+1 with eCent(0) = 1. The recursion applies these to Delaunay FACES of every
  // dimension (down to single vertices), so they must use the general routine.
  MyVector<T> center_homog(MyMatrix<T> const &EXT) const {
    return CenterDelaunayPolytopeGeneral<T>(GramMat, EXT);
  }

  T square_radius(MyMatrix<T> const &EXT) const {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT);
    return eCP.SquareRadius;
  }

  // Same, but for a FULL Delaunay cell (DT.l_dels[...].EXT), which is always
  // full-dimensional -- so the faster CenterDelaunayPolytope applies.
  MyVector<T> center_homog_fulldim(MyMatrix<T> const &EXT) const {
    return CenterDelaunayPolytope<T>(GramMat, EXT);
  }

  // Isobarycenter, homogeneous, length n+1 (first coordinate equals 1).
  MyVector<T> isobarycenter_homog(MyMatrix<T> const &EXT) const {
    int nbRow = EXT.rows();
    int nbCol = EXT.cols();
    MyVector<T> V = ZeroVector<T>(nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbCol; iCol++)
        V(iCol) += EXT(iRow, iCol);
    T eRow(nbRow);
    return V / eRow;
  }

  // Application of an affine (n+1)x(n+1) matrix on a homogeneous row vector:
  // GAP convention v * M.
  MyVector<T> apply_matrix(MyVector<T> const &v, MyMatrix<T> const &M) const {
    return M.transpose() * v;
  }

  // Convert a group order / index (TintGrp, an mpz_class) to the field T.
  // We go through the decimal string because there is no direct conversion
  // from mpz_class to the boost rational types.
  T tint_grp_to_T(TintGrp const &x) const {
    std::stringstream ss;
    ss << x;
    return ParseScalar<T>(ss.str());
  }

  // -------- the lattice stabilizer / equivalence with a marked center --------
  // These mirror DataLattice.FuncAutom and DataLattice.FuncEquiv. Contrary to
  // the plain Delaunay versions in LatticeDelaunay.h they take an explicit
  // center (needed because the algorithm marks the center of the dual cell),
  // and FuncAutom returns the affine matrix generators (the GAP matrix group),
  // not only a permutation group.

  struct MatrixGroupInfo {
    TintGrp order;
    std::vector<MyMatrix<T>> gens;
  };

  // Build the affine (n+1)x(n+1) matrix [1,x] -> [1, c + (x-c) L] fixing the
  // marked center c.
  MyMatrix<T> affine_from_linear(MyMatrix<T> const &L,
                                 MyVector<T> const &c) const {
    MyMatrix<T> M = ZeroMatrix<T>(n + 1, n + 1);
    M(0, 0) = 1;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        M(i + 1, j + 1) = L(i, j);
    MyVector<T> cL = L.transpose() * c;
    for (int j = 0; j < n; j++)
      M(0, j + 1) = c(j) - cL(j);
    return M;
  }

  MyVector<T> affine_part(MyVector<T> const &v_homog) const {
    MyVector<T> v(n);
    for (int i = 0; i < n; i++)
      v(i) = v_homog(i + 1);
    return v;
  }

  // Evaluate a jet matrix at the non-degenerate reference point t0 = tmax/2,
  // collapsing it to the base field. The lattice stabilizer / equivalence and its
  // realizing integral map are constant along the OPEN iso-Delaunay segment, so
  // computing them at the interior point t0 yields the correct generic-segment
  // group -- which is the t -> 0+ group, NOT the inflated t = 0 wall group -- over
  // a non-degenerate rational configuration. The whole weight-matrix machinery
  // (canonical ordering, GetBasisFromOrdering, TestEquivalence) then runs in the
  // base field and never divides by a zero divisor. For an ordinary scalar T this
  // is the identity (eval_at ignores t0), reproducing the original computation.
  MyMatrix<Tscal> eval_mat_ref(MyMatrix<T> const &M) const {
    MyMatrix<Tscal> R(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++)
      for (int j = 0; j < M.cols(); j++)
        R(i, j) = eval_at(M(i, j), t0_triang);
    return R;
  }

  // FuncAutom(EXT, center): affine automorphisms of EXT fixing the marked
  // center.
  MatrixGroupInfo func_autom_center(MyMatrix<T> const &EXT,
                                    MyVector<T> const &center_h) const {
    using Tidx = typename Telt::Tidx;
    MyVector<T> c = affine_part(center_h);
    MyMatrix<T> EXText = get_reduced_delaunay_shv(EXT, GramMat, SHV, c);
    // Color the augmented point set {vertices - c} U {SHV} so the isometry
    // search cannot MIX the two blocks. This is essential: a 4-dimensional cube
    // face has vertices of norm sqrt(4)/2 = 1, exactly the norm of the SHV
    // vectors +-e_i, so the two families sit on the same sphere. Without the
    // coloring the combinatorial stabilizer contains "accidental" isometries
    // swapping cube corners with +-e_i, whose linear realization has half-integer
    // entries (not in GL_n(Z)); those first appear at n = 8 (index 3) and inflate
    // the orbit multiplicity, over-counting the recursive cell moment by 3. With
    // the coloring the stabilizer is block-preserving, and because SHV is a
    // Z-basis its realization is automatically unimodular -- so we get the exact
    // lattice automorphism group directly, jet-ably, with no integral refinement.
    // The stabilizer (and its realization) is computed at the reference point t0
    // over the base field (see eval_mat_ref): correct generic-segment group, no
    // zero-divisor in the weight-matrix machinery.
    MyMatrix<Tscal> EXText_ref = eval_mat_ref(EXText);
    int nbVert = EXT.rows();
    int nbP = EXText.rows();
    std::vector<Tscal> Vdiag(nbP, Tscal(0));
    for (int i = nbVert; i < nbP; i++)
      Vdiag[i] = Tscal(1);
    std::vector<std::vector<Tidx>> ListGenPerm =
        GetListGenAutomorphism_ListMat_Vdiag<Tscal, Tscal, Tgroup>(
            EXText_ref, ListMat_ref, Vdiag, os);
    std::vector<Telt> LGen;
    for (auto &eList : ListGenPerm)
      LGen.push_back(Telt(eList));
    Tgroup GRP(LGen, nbP);
    MatrixGroupInfo info;
    info.order = GRP.size();
    for (auto &eGen : GRP.GeneratorsOfGroup()) {
#ifdef DEBUG_QUANTIZATION_INTEGRAL
      os << "QI: func_autom_center, Before FindTransformation\n";
#endif
      MyMatrix<Tscal> Lc =
          FindTransformation<Tscal, Telt>(EXText_ref, EXText_ref, eGen);
      MyMatrix<T> L = UniversalMatrixConversion<T, Tscal>(Lc);
#ifdef DEBUG_QUANTIZATION_INTEGRAL
      os << "QI: func_autom_center, After FindTransformation\n";
#endif
#ifdef SANITY_CHECK_QUANTIZATION_INTEGRAL
      // Every stabilizer generator must be realized by a lattice automorphism,
      // i.e. an integral (GL_n(Z)) linear map. A non-integral realization means
      // the isometry group contains a symmetry that does not preserve the
      // lattice (the vertex/SHV coloring is missing or wrong) -- exactly the bug
      // that silently inflated the group and over-counted the cell moment.
      if (!IsIntegralMatrix(L)) {
        std::cerr << "QUANT: func_autom_center produced a NON-INTEGRAL "
                     "automorphism; the isometry group is larger than the "
                     "lattice automorphism group (bad vertex/SHV coloring)\n";
        throw TerminalException{1};
      }
#endif
      info.gens.push_back(affine_from_linear(L, c));
    }
    return info;
  }

  // FuncEquiv(EXT1, EXT2, c1, c2): affine map EXT1 -> EXT2 sending c1 -> c2, or
  // nothing.
  std::optional<MyMatrix<T>>
  func_equiv_center(MyMatrix<T> const &EXT1, MyVector<T> const &c1_h,
                    MyMatrix<T> const &EXT2, MyVector<T> const &c2_h) const {
    using Tidx = typename Telt::Tidx;
    MyVector<T> c1 = affine_part(c1_h);
    MyVector<T> c2 = affine_part(c2_h);
    auto extend = [&](MyMatrix<T> const &L) -> MyMatrix<T> {
      MyMatrix<T> M = ZeroMatrix<T>(n + 1, n + 1);
      M(0, 0) = 1;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          M(i + 1, j + 1) = L(i, j);
      MyVector<T> delta = c2 - L.transpose() * c1;
      for (int i = 0; i < n; i++)
        M(0, i + 1) = delta(i);
      return M;
    };
    MyMatrix<T> EXText1 = get_reduced_delaunay_shv(EXT1, GramMat, SHV, c1);
    MyMatrix<T> EXText2 = get_reduced_delaunay_shv(EXT2, GramMat, SHV, c2);
    // Same coloring as func_autom_center: distinguish cell vertices from SHV so
    // the equivalence cannot mix the two blocks. The single equivalence returned
    // is then realized by a unimodular map (SHV is a Z-basis), so no coset search
    // via the (non-jet-able) LinPolytopeIntegral_Isomorphism is needed.
    // Equivalence (and its realization) at the reference point t0 over the base
    // field (see eval_mat_ref): correct generic-segment equivalence, no
    // zero-divisor in the weight-matrix machinery.
    MyMatrix<Tscal> EXText1_ref = eval_mat_ref(EXText1);
    MyMatrix<Tscal> EXText2_ref = eval_mat_ref(EXText2);
    std::vector<Tscal> Vdiag1(EXText1.rows(), Tscal(0)),
        Vdiag2(EXText2.rows(), Tscal(0));
    for (int i = EXT1.rows(); i < EXText1.rows(); i++)
      Vdiag1[i] = Tscal(1);
    for (int i = EXT2.rows(); i < EXText2.rows(); i++)
      Vdiag2[i] = Tscal(1);
    std::optional<std::vector<Tidx>> eRes =
        TestEquivalence_ListMat_Vdiag<Tscal, Tscal, Tidx>(
            EXText1_ref, ListMat_ref, Vdiag1, EXText2_ref, ListMat_ref, Vdiag2,
            os);
    if (!eRes)
      return {};
    Telt eElt(*eRes);
#ifdef DEBUG_QUANTIZATION_INTEGRAL
    os << "QI: func_equiv_center, Before FindTransformation\n";
#endif
    MyMatrix<Tscal> MatEquivC =
        FindTransformation<Tscal, Telt>(EXText1_ref, EXText2_ref, eElt);
    MyMatrix<T> MatEquiv = UniversalMatrixConversion<T, Tscal>(MatEquivC);
#ifdef DEBUG_QUANTIZATION_INTEGRAL
    os << "QI: func_equiv_center, After FindTransformation\n";
#endif
#ifdef SANITY_CHECK_QUANTIZATION_INTEGRAL
    // The equivalence must be realized by a lattice isomorphism (integral map).
    // A non-integral realization means the coloring failed to keep the vertex
    // and SHV blocks apart and the two cells were wrongly matched.
    if (!IsIntegralMatrix(MatEquiv)) {
      std::cerr << "QUANT: func_equiv_center produced a NON-INTEGRAL "
                   "equivalence; it is not a lattice isomorphism (bad "
                   "vertex/SHV coloring)\n";
      throw TerminalException{1};
    }
#endif
    return extend(MatEquiv);
  }

  // -------- the invariant used as a fast prefilter for equivalence --------
  // A coarser invariant than the GAP DelaunayInvariantLattice is fine: it is
  // only used to skip the equivalence test, never to decide equality on its
  // own.
  struct DelInvariant {
    int nbVert;
    int rank;
    bool iso_eq_center;
    T sqr_radius;
    std::vector<std::pair<T, int>> dist_occ;
    bool operator==(DelInvariant const &o) const {
      return nbVert == o.nbVert && rank == o.rank &&
             iso_eq_center == o.iso_eq_center && sqr_radius == o.sqr_radius &&
             dist_occ == o.dist_occ;
    }
  };

  DelInvariant delaunay_invariant(MyMatrix<T> const &EXT) const {
    DelInvariant inv;
    int nbVert = EXT.rows();
    inv.nbVert = nbVert;
    inv.rank = RankMat(EXT);
    // Center and radius are both needed here, so compute the circumsphere once
    // (center_homog + square_radius would each solve for the center separately).
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT);
    MyVector<T> const &cp = eCP.eCent;
    inv.sqr_radius = eCP.SquareRadius;
    MyVector<T> eIso = isobarycenter_homog(EXT);
    inv.iso_eq_center = (eIso == cp);
    std::map<T, int> map_occ;
    for (int i = 0; i < nbVert; i++) {
      for (int j = i + 1; j < nbVert; j++) {
        MyVector<T> eDiff(n);
        for (int u = 0; u < n; u++)
          eDiff(u) = EXT(i, u + 1) - EXT(j, u + 1);
        MyVector<T> prod = GramMat * eDiff;
        T dist(0);
        for (int u = 0; u < n; u++)
          dist += eDiff(u) * prod(u);
        map_occ[dist] += 1;
      }
    }
    for (auto &kv : map_occ)
      inv.dist_occ.push_back({kv.first, kv.second});
    return inv;
  }

  struct PairInvariant {
    DelInvariant inv_over;
    bool operator==(PairInvariant const &o) const {
      return inv_over == o.inv_over;
    }
  };

  PairInvariant invariant_of_pair([[maybe_unused]] MyMatrix<T> const &EXT,
                                  MyMatrix<T> const &EXTover) const {
    PairInvariant pinv;
    pinv.inv_over = delaunay_invariant(EXTover);
    return pinv;
  }

  // -------- the basis of the dual cell (GAP GetBasis) --------
  // [center] together with a basis of the Gram-orthogonal complement of the
  // affine span of EXT.
  MyMatrix<T> get_basis(MyMatrix<T> const &EXT) const {
    MyVector<T> cp = center_homog(EXT);
    int nbVert = EXT.rows();
    MyMatrix<T> SpaceBas(nbVert, n);
    for (int iVert = 0; iVert < nbVert; iVert++)
      for (int i = 0; i < n; i++)
        SpaceBas(iVert, i) = EXT(iVert, i + 1) - EXT(0, i + 1);
    // The Gram-orthogonal complement of the cell's affine span, { x : x^T G s = 0
    // for all s in span(SpaceBas) }. Rather than NullspaceMat(G * SpaceBas^T) --
    // whose row reduction, over jets, divides by a zero-divisor pivot wherever
    // the t = 0 metric makes the reduced product degenerate -- we use the
    // identity x _|__G span(S)  <=>  G^T x _|__I span(S), i.e.
    //   NSP = NullspaceTrMat(SpaceBas) * G^{-1}.
    // NullspaceTrMat(SpaceBas) is the ordinary (identity) orthogonal complement
    // of the span; SpaceBas is t-independent (differences of the cell's fixed
    // vertices), so this nullspace is over constant jets and never hits a
    // zero divisor. Inverse(GramMat) is safe because the constant term of GramMat
    // is Q, which is invertible. The two together give the same jet-valued
    // Gram-orthogonal complement, division-free.
    MyMatrix<T> Ncomp = NullspaceTrMat(SpaceBas);
    MyMatrix<T> NSP = Ncomp * GramMatInv;
    int dimNSP = NSP.rows();
    MyMatrix<T> TheBasis(1 + dimNSP, n + 1);
    for (int i = 0; i < n + 1; i++)
      TheBasis(0, i) = cp(i);
    for (int iNSP = 0; iNSP < dimNSP; iNSP++) {
      TheBasis(1 + iNSP, 0) = 0;
      for (int i = 0; i < n; i++)
        TheBasis(1 + iNSP, i + 1) = NSP(iNSP, i);
    }
    return TheBasis;
  }

  // -------- integral linear algebra helpers --------

  // TransformIntegral: express INT1 (given in Basis1) in Basis2 (GAP).
  MyMatrix<T> transform_integral(MyMatrix<T> const &INT1,
                                 MyMatrix<T> const &Basis1,
                                 MyMatrix<T> const &Basis2) const {
    int k = Basis1.rows();
    MyMatrix<T> P(k, k);
    for (int i = 0; i < k; i++) {
      MyVector<T> row = GetMatrixRow(Basis1, i);
      std::optional<MyVector<T>> sol = SolutionMat(Basis2, row);
      MyVector<T> c = unfold_opt(sol, "transform_integral: SolutionMat fails");
      for (int j = 0; j < k; j++)
        P(i, j) = c(j);
    }
    T det = DeterminantMat(P);
    if (det < 0)
      det = -det;
    return det * (P.transpose() * INT1 * P);
  }

  // FuncLiftIntegralStd (GAP): lift a (k)x(k) integral to (k+1)x(k+1).
  MyMatrix<T> func_lift_integral_std(MyMatrix<T> const &TheInt) const {
    int nRel = TheInt.rows();
    MyMatrix<T> ret = ZeroMatrix<T>(nRel + 1, nRel + 1);
    T fac2(nRel + 2), fac1(nRel + 1), fac0(nRel);
    for (int i = 0; i < nRel; i++)
      for (int j = 0; j < nRel; j++)
        ret(i + 1, j + 1) = TheInt(i, j) / fac2;
    for (int i = 0; i < nRel; i++) {
      ret(i + 1, 0) = TheInt(i, 0) / fac1;
      ret(0, i + 1) = TheInt(0, i) / fac1;
    }
    ret(0, 0) = TheInt(0, 0) / fac0;
    return ret;
  }

  // -------- orbit barycenter of a symmetric matrix (action M -> P^t M P) -----
  // The average over the finite group is computed without enumerating the orbit
  // by src_group's OrbitBarycenterSymmetricMatrix, which averages over
  // M -> g M g^t; we hand it the transposed generators g = P^t. (The plain
  // vector version OrbitBarycenter is called directly at its single use site,
  // where apply_matrix uses the same v -> g^T v convention.)
  MyMatrix<T>
  orbit_barycenter_symmetric_matrix(MyMatrix<T> const &TheSymMat,
                                    std::vector<MyMatrix<T>> const &gens) const {
    std::vector<MyMatrix<T>> gensT;
    for (auto &eGen : gens)
      gensT.push_back(eGen.transpose());
    return OrbitBarycenterSymmetricMatrix(TheSymMat, gensT);
  }

  // Express a matrix group in a basis (GAP GroupExpressionInTheBasis): for each
  // generator g and basis row b, the coordinates of b*g in the basis.
  std::vector<MyMatrix<T>>
  group_expression_in_basis(MyMatrix<T> const &TheBasis,
                            std::vector<MyMatrix<T>> const &gens) const {
    int k = TheBasis.rows();
    std::vector<MyMatrix<T>> ret;
    for (auto &eGen : gens) {
      MyMatrix<T> M(k, k);
      for (int i = 0; i < k; i++) {
        MyVector<T> row = GetMatrixRow(TheBasis, i);
        MyVector<T> img = apply_matrix(row, eGen);
        std::optional<MyVector<T>> sol = SolutionMat(TheBasis, img);
        MyVector<T> c =
            unfold_opt(sol, "group_expression_in_basis: SolutionMat fails");
        for (int j = 0; j < k; j++)
          M(i, j) = c(j);
      }
      ret.push_back(M);
    }
    return ret;
  }

  // -------- the direct integral over a polytope (GAP DirectIntegralLRS) ------
  // Integral, over the polytope spanned by ListVert, of the degree<=2 monomials
  // expressed in TheBasis coordinates.
  MyMatrix<T> direct_integral(std::vector<MyVector<T>> const &ListVert,
                              MyMatrix<T> const &TheBasis) const {
    int nRel = TheBasis.rows() - 1;
    if (nRel == 0) {
      MyMatrix<T> ret(1, 1);
      ret(0, 0) = 1;
      return ret;
    }
    int nbVert = ListVert.size();
    MyMatrix<T> EXTinBasis(nbVert, nRel + 1);
    for (int iVert = 0; iVert < nbVert; iVert++) {
      std::optional<MyVector<T>> sol = SolutionMat(TheBasis, ListVert[iVert]);
      MyVector<T> c = unfold_opt(sol, "direct_integral: SolutionMat fails");
      for (int j = 0; j < nRel + 1; j++)
        EXTinBasis(iVert, j) = c(j);
    }
    T IntDeg0(0);
    MyVector<T> IntDeg1 = ZeroVector<T>(nRel);
    // Degree-2 integral, symmetric, stored as its upper triangle (a <= b) flat.
    std::vector<T> IntDeg2ut(nRel * (nRel + 1) / 2, T(0));
    // The leaf triangulation is purely combinatorial (any triangulation of the
    // cell integrates the polynomial identically), so it must NOT be run over
    // jets: at t = 0 the split vertices coincide and lrs' fraction-free pivoting
    // would divide by a zero divisor. Instead evaluate the (jet) vertices at the
    // interior point t0, where the polytope is non-degenerate, and triangulate
    // that rational polytope with the ordinary lrs. For a plain scalar T this is
    // the identity (eval_at ignores t0) and reproduces the previous behaviour.
    MyMatrix<Tscal> EXTinBasis_t0(nbVert, nRel + 1);
    for (int iVert = 0; iVert < nbVert; iVert++)
      for (int j = 0; j < nRel + 1; j++)
        EXTinBasis_t0(iVert, j) = eval_at(EXTinBasis(iVert, j), t0_triang);
    // Triangulate the t0-evaluated polytope, keeping each simplex together with
    // its SIGNED determinant there (detT0). lrs produces the determinant
    // essentially for free, so this replaces a per-simplex DeterminantMat.
    std::vector<std::pair<std::vector<int>, Tscal>> trig =
        lrs::GetTriangulationDet(EXTinBasis_t0);
    // The covariance of the uniform distribution on a simplex with n+1 vertices
    // is (1/((n+1)(n+2))) sum_u (v_u - c)(v_u - c)^T; Tnp1 = n+1 is also the
    // barycenter divisor (vertex count). These are the linear factors (n+1),
    // (n+2), NOT factorials -- the overall n! normalization is applied at the
    // very end via factorial(nRel).
    T Tnp1(nRel + 1), Tnp2(nRel + 2);
    T invDenom = T(1) / (Tnp1 * Tnp2); // constant, hoisted out of the loop
    for (auto &[LV, detT0] : trig) {
      int sdim = LV.size();
      MyMatrix<T> EXTsimplex(sdim, nRel + 1);
      for (int u = 0; u < sdim; u++)
        for (int j = 0; j < nRel + 1; j++)
          EXTsimplex(u, j) = EXTinBasis(LV[u], j);
      // Oriented volume: the jet determinant of a split simplex is alpha*t^d + ..
      // (it vanishes at t = 0), so its sign is meaningless as a jet. trig is a
      // valid triangulation at the interior point t0, so the correct orientation
      // is the sign of the EXACT volume there: making every signed jet volume
      // positive at t0 reconstructs the cell integral as an analytic function of t
      // (equal to it on a neighbourhood of t0, hence -- both being polynomials --
      // for all t, including t -> 0^+). That sign must be the EXACT (untruncated)
      // volume at t0, i.e. the determinant of the simplex in the t0-evaluated
      // rational polytope EXTinBasis_t0. Using instead the order-N-truncated jet
      // determinant evaluated at t0 (eval_at(VolSimplex, t0)) can give the WRONG
      // sign at a large t0 -- that truncation error, not any leaf-flip, is what
      // forced the tiny t0 = tmax/1000000. For a plain scalar T, EXTinBasis_t0 ==
      // EXTinBasis, so this reduces to the ordinary sign of |det|.
      T VolSimplex = DeterminantMat(EXTsimplex);
#ifdef PROFILE_SIMPLEX_VALUATION
      {
        int val = scalar_valuation(VolSimplex);
        if (val > 31)
          val = 31;
        g_simplex_valuation[val]++;
        g_simplex_total++;
      }
#endif
#ifdef SANITY_CHECK_QUANTIZATION_INTEGRAL
      // trig triangulates the polytope evaluated at the non-degenerate interior
      // point t0, so every leaf simplex must have non-zero volume there. A zero
      // here means the triangulation is invalid (or t0 is degenerate).
      if (detT0 == Tscal(0)) {
        std::cerr << "QUANTIZATION_INTEGRAL: degenerate leaf simplex (zero volume "
                     "at the interior point t0); the leaf triangulation is "
                     "invalid\n";
        throw TerminalException{1};
      }
#endif
      if (detT0 < Tscal(0))
        VolSimplex = -VolSimplex;
      MyVector<T> bary = ZeroVector<T>(nRel + 1);
      for (int u = 0; u < sdim; u++)
        for (int j = 0; j < nRel + 1; j++)
          bary(j) += EXTsimplex(u, j);
      bary /= Tnp1;
      // Spatial coordinates of the barycenter.
      MyVector<T> barySpace = bary.tail(nRel);
      // Centered vertex coordinates (rows = v_u - c, spatial part only).
      MyMatrix<T> Diff(sdim, nRel);
      for (int u = 0; u < sdim; u++)
        for (int a = 0; a < nRel; a++)
          Diff(u, a) = EXTsimplex(u, a + 1) - barySpace(a);
      // Second-moment contribution of this simplex, accumulated on the upper
      // triangle only (the matrix is symmetric). Forming the covariance
      // Diff^T Diff replaces the previous contraction of each of the ~nRel^2 unit
      // matrices H against every centered vertex (a dense H*eDiff per H and per
      // vertex, O(nRel^4*sdim)); for the uniform distribution on the simplex the
      // second moment about the origin is
      //   M2(a,b) = Vol * ( c_a c_b + (Diff^T Diff)(a,b) / ((nRel+1)(nRel+2)) ).
      int pos = 0;
      for (int a = 0; a < nRel; a++)
        for (int b = a; b < nRel; b++) {
          T s(0);
          for (int u = 0; u < sdim; u++)
            s += Diff(u, a) * Diff(u, b);
          IntDeg2ut[pos] +=
              VolSimplex * (barySpace(a) * barySpace(b) + s * invDenom);
          pos++;
        }
      IntDeg0 += VolSimplex;
      for (int a = 0; a < nRel; a++)
        IntDeg1(a) += VolSimplex * barySpace(a);
    }
    MyMatrix<T> ret = ZeroMatrix<T>(nRel + 1, nRel + 1);
    {
      int pos = 0;
      for (int a = 0; a < nRel; a++)
        for (int b = a; b < nRel; b++) {
          ret(a + 1, b + 1) = IntDeg2ut[pos];
          ret(b + 1, a + 1) = IntDeg2ut[pos];
          pos++;
        }
    }
    for (int i = 0; i < nRel; i++) {
      ret(0, i + 1) = IntDeg1(i);
      ret(i + 1, 0) = IntDeg1(i);
    }
    ret(0, 0) = IntDeg0;
    return ret / factorial<T>(nRel);
  }

  // -------- the recursive structure --------

  struct DelRecord {
    int iOrb;
    MyMatrix<T> eMat; // affine (n+1)x(n+1)
    Face eSet;        // subset of the vertices of Delaunay iOrb
  };

  struct AddiInfo {
    Tgroup Stab;                      // stabilizer of eSet in Delaunay iOrb
    std::vector<Face> ListOrbitRelFacet;
  };

  struct Symbol {
    MyMatrix<T> EXT;
    MatrixGroupInfo StabEXT;
    TintGrp DStabEXT_order;
    int TheLevel;
    DelRecord TheDel;
    PairInvariant eInvariant;
    bool has_invariant = false;
  };

  struct SoftComp {
    MyMatrix<T> EXT;
    std::vector<DelRecord> ListRecord;
    MyMatrix<T> TheBasis;
    DelInvariant eInvariant;
    std::vector<AddiInfo> ListAdditionalInfo;
    TintGrp NumberIncident;
    MyVector<T> TheBarycenter;
    MyMatrix<T> TheIntegral;
  };

  std::vector<SoftComp> ListOrbitIntegrals; // the bank

  // GAP SoftComputation.
  SoftComp soft_computation(MyMatrix<T> const &EXT,
                            MatrixGroupInfo const &StabEXT,
                            DelRecord const &OneRecord) {
    SoftComp sc;
    sc.EXT = EXT;
    sc.ListRecord.push_back(OneRecord);
    std::vector<int> ListStatus; // 0 = NO, 1 = YES
    ListStatus.push_back(0);
    sc.ListAdditionalInfo.resize(1);
    TintGrp OrdStabEXT = StabEXT.order;
    TintGrp NumberIncident(0);
    MyVector<T> SumElement = ZeroVector<T>(n + 1);
    auto FuncInsertRecord = [&](DelRecord const &rec) -> void {
      for (auto &existing : sc.ListRecord) {
        if (existing.iOrb == rec.iOrb) {
          Tgroup const &TheStab = DT.l_dels[rec.iOrb].GRP;
          std::optional<Telt> g =
              TheStab.RepresentativeAction_OnSets(existing.eSet, rec.eSet);
          if (g)
            return;
        }
      }
      ListStatus.push_back(0);
      sc.ListRecord.push_back(rec);
      sc.ListAdditionalInfo.push_back(AddiInfo());
    };
    while (true) {
      bool IsFinished = true;
      for (size_t iRecord = 0; iRecord < sc.ListRecord.size(); iRecord++) {
        if (ListStatus[iRecord] == 0) {
          IsFinished = false;
          DelRecord eRecord = sc.ListRecord[iRecord];
          MyMatrix<T> const &TheEXT = DT.l_dels[eRecord.iOrb].EXT;
          std::vector<Delaunay_AdjO<T>> const &TheAdjacencies =
              DT.l_dels[eRecord.iOrb].ListAdj;
          Tgroup const &TheStab = DT.l_dels[eRecord.iOrb].GRP;
          Tgroup DStabEXT = TheStab.Stabilizer_OnSets(eRecord.eSet);
          TintGrp index = OrdStabEXT / DStabEXT.size();
          NumberIncident += index;
          MyVector<T> TheCenter =
              apply_matrix(center_homog_fulldim(TheEXT), eRecord.eMat);
          MyVector<T> SingleInv = OrbitBarycenter(TheCenter, StabEXT.gens);
          T indexT = tint_grp_to_T(index);
          SumElement += indexT * SingleInv;
          std::vector<Face> ListOrbitRelFacet;
          std::vector<Telt> StabGens = TheStab.GeneratorsOfGroup();
          for (auto &eAdjacency : TheAdjacencies) {
            vectface fullorbit = OrbitFace(eAdjacency.eInc, StabGens);
            vectface reps = OrbitSplittingSet(fullorbit, DStabEXT);
            for (auto &fInc : reps) {
              if (is_subset(eRecord.eSet, fInc)) {
                ListOrbitRelFacet.push_back(fInc);
                std::optional<Telt> optElt =
                    TheStab.RepresentativeAction_OnSets(eAdjacency.eInc, fInc);
                Telt eElt = unfold_opt(
                    optElt, "QI: soft_computation: RepresentativeAction fails");
                MyMatrix<T> const &TheEXT2 = DT.l_dels[eAdjacency.iOrb].EXT;
                MyMatrix<T> PhiMat =
                    RepresentVertexPermutation<T, Telt>(TheEXT, TheEXT, eElt);
                MyMatrix<T> eG2 = eAdjacency.eBigMat * PhiMat * eRecord.eMat;
                MyMatrix<T> Mimg = TheEXT2 * eG2;
                ContainerMatrix<T> cont(Mimg);
                Face eSet2(TheEXT2.rows());
                int nbE = EXT.rows();
                for (int iE = 0; iE < nbE; iE++) {
                  MyVector<T> V = GetMatrixRow(EXT, iE);
                  std::optional<size_t> opt = cont.GetIdx_v(V);
                  size_t pos = unfold_opt(
                      opt, "QI: soft_computation: vertex not found in adjacent");
                  eSet2[pos] = 1;
                }
                FuncInsertRecord({eAdjacency.iOrb, eG2, eSet2});
              }
            }
          }
          ListStatus[iRecord] = 1;
          AddiInfo ai;
          ai.Stab = DStabEXT;
          ai.ListOrbitRelFacet = ListOrbitRelFacet;
          sc.ListAdditionalInfo[iRecord] = ai;
        }
      }
      if (IsFinished)
        break;
    }
    sc.TheBasis = get_basis(EXT);
    sc.eInvariant = delaunay_invariant(EXT);
    sc.NumberIncident = NumberIncident;
    T NumberIncidentT = tint_grp_to_T(NumberIncident);
    sc.TheBarycenter = SumElement / NumberIncidentT;
    return sc;
  }

  // GAP VoronoiPolytopeListVertices.
  std::vector<MyVector<T>>
  voronoi_polytope_list_vertices(std::vector<MyMatrix<T>> const &StabGens,
                                 SoftComp const &sc) const {
    std::vector<MyVector<T>> ListVert;
    std::set<MyVector<T>> seen;
    for (auto &eSoft : sc.ListRecord) {
      MyMatrix<T> const &TheEXT = DT.l_dels[eSoft.iOrb].EXT;
      MyVector<T> eCenter = apply_matrix(center_homog_fulldim(TheEXT), eSoft.eMat);
      std::vector<MyVector<T>> orb =
          OrbitComputation_vector(StabGens, eCenter);
      for (auto &v : orb)
        if (seen.insert(v).second)
          ListVert.push_back(v);
    }
    return ListVert;
  }

  std::vector<MyVector<T>>
  OrbitComputation_vector(std::vector<MyMatrix<T>> const &gens,
                          MyVector<T> const &v) const {
    auto f_prod = [](MyVector<T> const &x,
                     MyMatrix<T> const &M) -> MyVector<T> {
      return M.transpose() * x;
    };
    return OrbitComputation<MyMatrix<T>, MyVector<T>, decltype(f_prod)>(
        gens, v, f_prod, os);
  }

  // GAP SingleListOrbit.
  std::vector<Face> single_list_orbit(MyMatrix<T> const &TheEXT,
                                      Face const &eSet,
                                      AddiInfo const &ai) const {
    int nbVert = TheEXT.rows();
    std::vector<Telt> StabGens = ai.Stab.GeneratorsOfGroup();
    // Orbits of the complement of eSet under ai.Stab.
    std::vector<int> complement;
    for (int i = 0; i < nbVert; i++)
      if (eSet[i] == 0)
        complement.push_back(i);
    std::vector<int> orbit_of(nbVert, -1);
    std::vector<std::vector<int>> O;
    for (int pt : complement) {
      if (orbit_of[pt] != -1)
        continue;
      int idxOrb = O.size();
      std::vector<int> orb;
      std::vector<int> todo{pt};
      orbit_of[pt] = idxOrb;
      orb.push_back(pt);
      while (!todo.empty()) {
        int cur = todo.back();
        todo.pop_back();
        for (auto &eGen : StabGens) {
          int img = OnPoints(cur, eGen);
          if (orbit_of[img] == -1) {
            orbit_of[img] = idxOrb;
            orb.push_back(img);
            todo.push_back(img);
          }
        }
      }
      O.push_back(orb);
    }
    // RelFAC: the orbit (under ai.Stab) of every representative facet.
    std::vector<Face> RelFAC;
    for (auto &eRep : ai.ListOrbitRelFacet) {
      vectface orb = OrbitFace(eRep, StabGens);
      for (auto &f : orb)
        RelFAC.push_back(f);
    }
    int nbFAC = RelFAC.size();
    MyMatrix<T> FAC(nbFAC, n + 1);
    for (int x = 0; x < nbFAC; x++) {
      MyVector<T> ineq = FindFacetInequality(TheEXT, RelFAC[x]);
      for (int j = 0; j < n + 1; j++)
        FAC(x, j) = ineq(j);
    }
    int RNK = (nbFAC == 0) ? 0 : RankMat(FAC);
    std::vector<int> ListVertStatus(O.size(), 1);
    std::vector<Face> ListOrbit;
    for (size_t iOrb = 0; iOrb < O.size(); iOrb++) {
      if (ListVertStatus[iOrb] == 1) {
        int firstVert = O[iOrb][0];
        std::vector<int> ListSel;
        for (int x = 0; x < nbFAC; x++)
          if (RelFAC[x][firstVert] == 1)
            ListSel.push_back(x);
        ListVertStatus[iOrb] = 0;
        if (ListSel.size() > 0) {
          MyMatrix<T> FACinc(ListSel.size(), n + 1);
          for (size_t u = 0; u < ListSel.size(); u++)
            for (int j = 0; j < n + 1; j++)
              FACinc(u, j) = FAC(ListSel[u], j);
          if (RankMat(FACinc) == RNK - 1) {
            // VertSet = Intersection over the selected facets.
            Face VertSet(nbVert);
            for (int i = 0; i < nbVert; i++)
              VertSet[i] = 1;
            for (int sel : ListSel)
              for (int i = 0; i < nbVert; i++)
                if (RelFAC[sel][i] == 0)
                  VertSet[i] = 0;
            // Mark off all complement vertices of VertSet \ eSet.
            for (auto idx = VertSet.find_first(); idx != Face::npos;
                 idx = VertSet.find_next(idx)) {
              if (eSet[idx] == 0) {
                int oi = orbit_of[idx];
                if (oi != -1)
                  ListVertStatus[oi] = 0;
              }
            }
            Face newFace = VertSet;
            for (auto idx = eSet.find_first(); idx != Face::npos;
                 idx = eSet.find_next(idx))
              newFace[idx] = 1;
            ListOrbit.push_back(newFace);
          }
        }
      }
    }
    return ListOrbit;
  }

  // GAP ListOrbitContainingEXT.
  std::vector<Symbol> list_orbit_containing_EXT(Symbol const &TheSymbol,
                                                SoftComp const &sc) {
    MyMatrix<T> const &EXT = TheSymbol.EXT;
    int TheLevel = TheSymbol.TheLevel;
    std::vector<Symbol> ListOut;
    T DelVal(1);
    MyVector<T> isoEXT = isobarycenter_homog(EXT);
    auto FuncInsert = [&](DelRecord const &rec) -> void {
      MyMatrix<T> const &TheEXT1 = DT.l_dels[rec.iOrb].EXT;
      std::vector<int> LV = FaceToVector<int>(rec.eSet);
      MyMatrix<T> sub(LV.size(), n + 1);
      for (size_t u = 0; u < LV.size(); u++)
        for (int j = 0; j < n + 1; j++)
          sub(u, j) = TheEXT1(LV[u], j);
      MyMatrix<T> EXT1 = sub * rec.eMat;
      PairInvariant eInvariant = invariant_of_pair(EXT, EXT1);
      MyVector<T> isoEXT1 = isobarycenter_homog(EXT1);
      MyVector<T> eIso1 = (isoEXT1 + DelVal * isoEXT) / (DelVal + T(1));
      for (auto &eSymbol2 : ListOut) {
        MyMatrix<T> const &EXT2 = eSymbol2.EXT;
        MyVector<T> isoEXT2 = isobarycenter_homog(EXT2);
        MyVector<T> eIso2 = (isoEXT2 + DelVal * isoEXT) / (DelVal + T(1));
        if (eInvariant == eSymbol2.eInvariant) {
          if (func_equiv_center(EXT, eIso1, EXT, eIso2))
            return;
        }
      }
      MatrixGroupInfo DStabInfo = func_autom_center(EXT1, eIso1);
      MatrixGroupInfo StabInfo = func_autom_center(EXT1, isoEXT1);
      Symbol sym;
      sym.EXT = EXT1;
      sym.StabEXT = StabInfo;
      sym.DStabEXT_order = DStabInfo.order;
      sym.TheLevel = TheLevel + 1;
      sym.TheDel = rec;
      sym.eInvariant = eInvariant;
      sym.has_invariant = true;
      ListOut.push_back(sym);
    };
    for (size_t iRecord = 0; iRecord < sc.ListRecord.size(); iRecord++) {
      DelRecord const &eRecord = sc.ListRecord[iRecord];
      AddiInfo const &ai = sc.ListAdditionalInfo[iRecord];
      MyMatrix<T> const &TheEXT = DT.l_dels[eRecord.iOrb].EXT;
      std::vector<Face> ListOrbGenerated =
          single_list_orbit(TheEXT, eRecord.eSet, ai);
      for (auto &eOrb : ListOrbGenerated)
        FuncInsert({eRecord.iOrb, eRecord.eMat, eOrb});
    }
    return ListOut;
  }

  // GAP FuncCheckInBank.
  std::optional<MyMatrix<T>> func_check_in_bank(SoftComp const &sc1) {
    MyVector<T> eIso1 = sc1.TheBarycenter;
    MyMatrix<T> const &TheBasis1 = sc1.TheBasis;
    for (auto &sc2 : ListOrbitIntegrals) {
      if (sc2.eInvariant == sc1.eInvariant) {
        MyVector<T> eIso2 = sc2.TheBarycenter;
        std::optional<MyMatrix<T>> eEquiv =
            func_equiv_center(sc2.EXT, eIso2, sc1.EXT, eIso1);
        if (eEquiv) {
          MyMatrix<T> ImageTheBasis2 = sc2.TheBasis * (*eEquiv);
          return transform_integral(sc2.TheIntegral, ImageTheBasis2, TheBasis1);
        }
      }
    }
    return {};
  }

  // GAP FuncRespawn heuristic.
  bool func_respawn(TintGrp const &OrdGrp, TintGrp const &NBV,
                    int TheLevel) const {
    TintGrp v20(20), v130(130), v70(70), v100(100), v110(110), v1000(1000),
        v10000(10000);
    if (NBV < v20)
      return false;
    if (NBV > v130)
      return true;
    if (OrdGrp > v10000)
      return true;
    if (NBV > v100 && OrdGrp > v1000)
      return true;
    if (TheLevel == 2)
      return false;
    if (NBV < v70)
      return false;
    if (NBV > v110 && OrdGrp > v100 && TheLevel < 2)
      return true;
    return false;
  }

  // GAP __RecursiveIntegralEvaluation.
  SoftComp recursive_integral_evaluation(Symbol const &TheSymbol) {
    MyMatrix<T> const &EXT = TheSymbol.EXT;
    MatrixGroupInfo const &StabEXT = TheSymbol.StabEXT;
    SoftComp sc = soft_computation(EXT, StabEXT, TheSymbol.TheDel);
#ifdef DEBUG_QUANTIZATION_INTEGRAL
    os << "QUANT: Level=" << TheSymbol.TheLevel << " |EXT|=" << EXT.rows()
       << " |Stab|=" << StabEXT.order << " |Dels|=" << sc.NumberIncident
       << "\n";
#endif
    std::optional<MyMatrix<T>> testCheck = func_check_in_bank(sc);
    if (testCheck) {
      sc.TheIntegral = *testCheck;
      return sc;
    }
    bool testRespawn =
        func_respawn(StabEXT.order, sc.NumberIncident, TheSymbol.TheLevel);
    MyMatrix<T> TheIntegral;
    if (testRespawn) {
      int k = sc.TheBasis.rows();
      TheIntegral = ZeroMatrix<T>(k, k);
      std::vector<Symbol> ListOrbit = list_orbit_containing_EXT(TheSymbol, sc);
      MyVector<T> Ccenter = sc.TheBarycenter;
      std::vector<MyMatrix<T>> gens_in_basis =
          group_expression_in_basis(sc.TheBasis, StabEXT.gens);
      for (auto &eOrbRecord : ListOrbit) {
        TintGrp mult = StabEXT.order / eOrbRecord.DStabEXT_order;
        T multT = tint_grp_to_T(mult);
        SoftComp eOrbComp = recursive_integral_evaluation(eOrbRecord);
        MyMatrix<T> const &eOrbCompBasis = eOrbComp.TheBasis;
        int m = eOrbCompBasis.rows();
        // NewBasis = [Ccenter, b0 - Ccenter, b1, ..., b_{m-1}].
        MyMatrix<T> NewBasis(m + 1, n + 1);
        for (int j = 0; j < n + 1; j++) {
          NewBasis(0, j) = Ccenter(j);
          NewBasis(1, j) = eOrbCompBasis(0, j) - Ccenter(j);
        }
        for (int u = 1; u < m; u++)
          for (int j = 0; j < n + 1; j++)
            NewBasis(u + 1, j) = eOrbCompBasis(u, j);
        MyMatrix<T> lifted = func_lift_integral_std(eOrbComp.TheIntegral);
        MyMatrix<T> TheIntegralOrb =
            transform_integral(lifted, NewBasis, sc.TheBasis);
        MyMatrix<T> sym =
            orbit_barycenter_symmetric_matrix(TheIntegralOrb, gens_in_basis);
        TheIntegral += multT * sym;
      }
    } else {
      std::vector<MyVector<T>> ListVert =
          voronoi_polytope_list_vertices(StabEXT.gens, sc);
      TheIntegral = direct_integral(ListVert, sc.TheBasis);
    }
    sc.TheIntegral = TheIntegral;
    ListOrbitIntegrals.push_back(sc);
    return sc;
  }

  // GAP InitialPair: seed from a vertex of the first Delaunay. The vertex of
  // the lattice we integrate around is the origin; every Delaunay vertex is a
  // lattice point, so a translation suffices and no random walk is needed.
  Symbol initial_pair() {
    MyVector<T> TheVert = ZeroVector<T>(n + 1);
    TheVert(0) = 1;
    int iOrb = 0;
    MyMatrix<T> const &TheEXT = DT.l_dels[iOrb].EXT;
    // Translation mapping vertex 0 of TheEXT to the origin.
    MyMatrix<T> eMat = IdentityMat<T>(n + 1);
    for (int i = 0; i < n; i++)
      eMat(0, i + 1) = TheVert(i + 1) - TheEXT(0, i + 1);
    MyMatrix<T> EXT(1, n + 1);
    for (int j = 0; j < n + 1; j++)
      EXT(0, j) = TheVert(j);
    Face eSet(TheEXT.rows());
    eSet[0] = 1;
    Symbol sym;
    sym.EXT = EXT;
    sym.StabEXT = func_autom_center(EXT, TheVert);
    sym.DStabEXT_order = sym.StabEXT.order;
    sym.TheLevel = 0;
    sym.TheDel = {iOrb, eMat, eSet};
    return sym;
  }

  // GAP QuantizationIntegral.
  QuantizationResult<T> compute() {
    // SHV (the Z-spanning invariant family used by the marked-center stabilizer/
    // equivalence) is filled by the constructor -- see the note there on why the
    // LLL-based computation is kept out of compute().
    Symbol TheInit = initial_pair();
    SoftComp TheRes = recursive_integral_evaluation(TheInit);
    MyMatrix<T> Id = IdentityMat<T>(n + 1);
    MyMatrix<T> TheInt = transform_integral(TheRes.TheIntegral, TheRes.TheBasis, Id);
    T SecMoment(0);
    MyMatrix<T> SecMomentMat = ZeroMatrix<T>(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        SecMomentMat(i, j) = TheInt(i + 1, j + 1);
        SecMoment += GramMat(i, j) * TheInt(i + 1, j + 1);
      }
    }
    T TheVol = TheInt(0, 0);
    MyVector<T> TheBarycenter(n + 1);
    for (int j = 0; j < n + 1; j++)
      TheBarycenter(j) = TheInt(0, j) / TheVol;
    QuantizationResult<T> res;
    res.SecMoment = SecMoment;
    res.SecMomentMat = SecMomentMat;
    res.TheVolume = TheVol;
    res.TheBarycenter = TheBarycenter;
    // Normalized second moment (Conway-Sloane): G = (1/n) det^(-1/n) SecMoment.
    double det = UniversalScalarConversion<double, T>(DeterminantMat(GramMat));
    double secmom = UniversalScalarConversion<double, T>(SecMoment);
    res.NormalizedSecondMoment =
        (1.0 / n) * std::pow(det, -1.0 / n) * secmom;
    return res;
  }
};

template <typename T, typename Tint, typename Tgroup>
QuantizationResult<T>
ComputeQuantizationIntegral(DataLattice<T, Tint, Tgroup> &data,
                            DelaunayTesselation<T, Tgroup> const &DT,
                            std::ostream &os) {
  QuantizationComputer<T, Tint, Tgroup> comp(data, DT, os);
  return comp.compute();
}

// Reinterpret a Delaunay tessellation (whose EXT / adjacency matrices are
// integer, hence constant along an iso-Delaunay segment) over another scalar
// type -- e.g. T -> jet<T, N> to evaluate the integral along a ray Q + t H.
template <typename Tout, typename Tin, typename Tgroup>
DelaunayTesselation<Tout, Tgroup>
ConvertTesselationScalar(DelaunayTesselation<Tin, Tgroup> const &DT) {
  DelaunayTesselation<Tout, Tgroup> out;
  for (auto &e : DT.l_dels) {
    Delaunay_Entry<Tout, Tgroup> eo;
    eo.EXT = UniversalMatrixConversion<Tout, Tin>(e.EXT);
    eo.GRP = e.GRP;
    for (auto &adj : e.ListAdj) {
      Delaunay_AdjO<Tout> ao;
      ao.eInc = adj.eInc;
      ao.eBigMat = UniversalMatrixConversion<Tout, Tin>(adj.eBigMat);
      ao.iOrb = adj.iOrb;
      eo.ListAdj.push_back(std::move(ao));
    }
    out.l_dels.push_back(std::move(eo));
  }
  return out;
}

// Compile/feasibility probe: instantiate the quantization computer over jets.
// Given the t = 0 tessellation, invariant family and a jet Gram matrix, run the
// whole integral over jet<T, N> and return the second-moment matrix as jets.
template <typename T, int N, typename Tint, typename Tgroup>
MyMatrix<jet<T, N>>
QuantizationSecMomentMatJet(DelaunayTesselation<T, Tgroup> const &DT_base,
                            MyMatrix<T> const &SHV_base,
                            MyMatrix<jet<T, N>> const &GramMat_jet, int n,
                            T const &t0, std::ostream &os) {
  using Tj = jet<T, N>;
  DelaunayTesselation<Tj, Tgroup> DT_jet =
      ConvertTesselationScalar<Tj, T, Tgroup>(DT_base);
  MyMatrix<Tj> SHV_jet = UniversalMatrixConversion<Tj, T>(SHV_base);
  // t0: the interior point of the iso-Delaunay segment at which the leaf
  // triangulation is evaluated (see direct_integral). Tscal of jet<T,N> is T.
  QuantizationComputer<Tj, Tint, Tgroup> comp(n, SHV_jet, GramMat_jet, DT_jet,
                                              os, t0);
  MyMatrix<jet<T, N>> res = comp.compute().SecMomentMat;
#ifdef PROFILE_SIMPLEX_VALUATION
  os << "SIMPLEX_VALUATION total=" << g_simplex_total << "\n";
  for (int v = 0; v < 8; v++)
    if (g_simplex_valuation[v] > 0)
      os << "  valuation " << v << ": " << g_simplex_valuation[v] << " ("
         << (100.0 * g_simplex_valuation[v] / g_simplex_total) << "%)\n";
#endif
  return res;
}

template <typename T> std::string normg_string(double val) {
  std::stringstream ss;
  ss << std::setprecision(15) << val;
  return ss.str();
}

template <typename T>
void WriteQuantizationGAP(std::ostream &os_out,
                          QuantizationResult<T> const &res) {
  os_out << "rec(SecMoment:=" << res.SecMoment << ",\n";
  os_out << "SecMomentMat:=" << StringMatrixGAP(res.SecMomentMat) << ",\n";
  os_out << "TheVolume:=" << res.TheVolume << ",\n";
  os_out << "TheBarycenter:=" << StringVectorGAP(res.TheBarycenter) << ",\n";
  os_out << "NormalizedSecondMoment:="
         << normg_string<T>(res.NormalizedSecondMoment) << ")";
}

template <typename T>
void WriteQuantizationPYTHON(std::ostream &os_out,
                             QuantizationResult<T> const &res) {
  os_out << "{\"SecMoment\":" << res.SecMoment;
  os_out << ", \"SecMomentMat\":" << StringMatrixPYTHON(res.SecMomentMat);
  os_out << ", \"TheVolume\":" << res.TheVolume;
  os_out << ", \"TheBarycenter\":" << StringVectorPYTHON(res.TheBarycenter);
  os_out << ", \"NormalizedSecondMoment\":"
         << normg_string<T>(res.NormalizedSecondMoment) << "}";
}

// Isotropy (extremal / "white" quantizer) test. By the Zamir-Feder theorem the
// optimal lattice quantizer is white: the second-moment matrix M of the Voronoi
// cell is proportional to the inverse Gram matrix, equivalently
//   GramMat * M = Lambda * I
// for a scalar Lambda. We test this exactly: with Lambda = trace(GramMat*M)/n,
// the lattice is isotropic iff the defect DefectMat = GramMat*M - Lambda*I is
// the zero matrix. The Voronoi-cell volume normalization only rescales M, Lambda
// and DefectMat by a common factor, so it does not affect the exact test.
template <typename T> struct IsotropyResult {
  MyMatrix<T> SecMomentMat;          // M = int_{V_0} u_i u_j (volume = TheVolume)
  MyMatrix<T> GramTimesSecMomentMat; // GramMat * M
  T Lambda;                          // trace(GramMat * M) / n
  MyMatrix<T> DefectMat;             // GramMat * M - Lambda * I
  bool IsIsotropic;                  // DefectMat == 0
};

template <typename T, typename Tint, typename Tgroup>
IsotropyResult<T> ComputeIsotropy(DataLattice<T, Tint, Tgroup> &data,
                                  DelaunayTesselation<T, Tgroup> const &DT,
                                  MyMatrix<T> const &GramMat, std::ostream &os) {
  QuantizationResult<T> qres =
      ComputeQuantizationIntegral<T, Tint, Tgroup>(data, DT, os);
  int n = GramMat.rows();
  MyMatrix<T> M = qres.SecMomentMat;
  MyMatrix<T> Prod = GramMat * M;
  T tr(0);
  for (int i = 0; i < n; i++) {
    tr += Prod(i, i);
  }
  T Lambda = tr / T(n);
  MyMatrix<T> DefectMat = Prod - Lambda * IdentityMat<T>(n);
  bool IsIsotropic = true;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (DefectMat(i, j) != T(0)) {
        IsIsotropic = false;
      }
    }
  }
  IsotropyResult<T> res;
  res.SecMomentMat = M;
  res.GramTimesSecMomentMat = Prod;
  res.Lambda = Lambda;
  res.DefectMat = DefectMat;
  res.IsIsotropic = IsIsotropic;
  return res;
}

template <typename T>
void WriteIsotropyGAP(std::ostream &os_out, IsotropyResult<T> const &res) {
  os_out << "return rec(SecMomentMat:=" << StringMatrixGAP(res.SecMomentMat)
         << ",\n";
  os_out << "GramTimesSecMomentMat:="
         << StringMatrixGAP(res.GramTimesSecMomentMat) << ",\n";
  os_out << "Lambda:=" << res.Lambda << ",\n";
  os_out << "DefectMat:=" << StringMatrixGAP(res.DefectMat) << ",\n";
  os_out << "IsIsotropic:=" << (res.IsIsotropic ? "true" : "false") << ");\n";
}

// clang-format off
#endif  // SRC_DELAUNAY_QUANTIZATIONINTEGRAL_H_
// clang-format on
