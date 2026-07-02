// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_QUANTIZATIONINTEGRAL_H_
#define SRC_DELAUNAY_QUANTIZATIONINTEGRAL_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "POLY_Fundamental.h"
#include "POLY_lrslib.h"
#include "InvariantVectorFamily.h"
#include "MatrixGroupAverage.h"
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

template <typename T, typename Tint, typename Tgroup>
struct QuantizationComputer {
  using Telt = typename Tgroup::Telt;
  using TintGrp = typename Tgroup::Tint;

  DataLattice<T, Tint, Tgroup> &data;
  DelaunayTesselation<T, Tgroup> const &DT;
  std::ostream &os;
  int n;
  MyMatrix<T> const &GramMat;

  QuantizationComputer(DataLattice<T, Tint, Tgroup> &data,
                       DelaunayTesselation<T, Tgroup> const &DT,
                       std::ostream &os)
      : data(data), DT(DT), os(os), n(data.n), GramMat(data.solver.GramMat) {}

  // -------- small geometric helpers (homogeneous coordinates) --------

  // The center (and radius) of the circumscribing sphere, homogeneous, length
  // n+1 with eCent(0) = 1.
  MyVector<T> center_homog(MyMatrix<T> const &EXT) const {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT);
    return eCP.eCent;
  }

  T square_radius(MyMatrix<T> const &EXT) const {
    CP<T> eCP = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT);
    return eCP.SquareRadius;
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

  // The lcm of the denominators of the entries (GAP ListFactors).
  T list_factors(MyVector<T> const &v) const {
    MyVector<T> scaled = RemoveFractionVector(v);
    int len = v.size();
    for (int i = 0; i < len; i++) {
      if (v(i) != 0) {
        T quot = scaled(i) / v(i);
        if (quot < 0)
          quot = -quot;
        return quot;
      }
    }
    return T(1);
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

  // FuncAutom(EXT, center): affine automorphisms of EXT fixing the marked
  // center.
  MatrixGroupInfo func_autom_center(MyMatrix<T> const &EXT,
                                    MyVector<T> const &center_h) const {
    using Tgr = GraphListAdj;
    using Tidx_value = int16_t;
    MyVector<T> c = affine_part(center_h);
    MyMatrix<T> EXText = get_reduced_delaunay_shv(EXT, GramMat, data.SHV, c);
    WeightMatrix<true, T, Tidx_value> WMat =
        GetSimpleWeightMatrix<T, Tidx_value>(EXText, GramMat, os);
    Tgroup GRPisom =
        GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
    MyMatrix<T> EXTextInt = RemoveFractionMatrix(EXText);
    Tgroup GRPlatt = LinPolytopeIntegral_Stabilizer(EXTextInt, GRPisom, os);
    MatrixGroupInfo info;
    info.order = GRPlatt.size();
    for (auto &eGen : GRPlatt.GeneratorsOfGroup()) {
      MyMatrix<T> L = FindTransformation<T, Telt>(EXText, EXText, eGen);
      info.gens.push_back(affine_from_linear(L, c));
    }
    return info;
  }

  // FuncEquiv(EXT1, EXT2, c1, c2): affine map EXT1 -> EXT2 sending c1 -> c2, or
  // nothing.
  std::optional<MyMatrix<T>>
  func_equiv_center(MyMatrix<T> const &EXT1, MyVector<T> const &c1_h,
                    MyMatrix<T> const &EXT2, MyVector<T> const &c2_h) const {
    using Tgr = GraphListAdj;
    using Tidx_value = int16_t;
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
    MyMatrix<T> EXText1 = get_reduced_delaunay_shv(EXT1, GramMat, data.SHV, c1);
    MyMatrix<T> EXText2 = get_reduced_delaunay_shv(EXT2, GramMat, data.SHV, c2);
    WeightMatrix<true, T, Tidx_value> WMat1 =
        GetSimpleWeightMatrix<T, Tidx_value>(EXText1, GramMat, os);
    WeightMatrix<true, T, Tidx_value> WMat2 =
        GetSimpleWeightMatrix<T, Tidx_value>(EXText2, GramMat, os);
    std::optional<Telt> eRes =
        TestEquivalenceWeightMatrix<T, Telt, Tidx_value>(WMat1, WMat2, os);
    if (!eRes)
      return {};
    Telt const &eElt = *eRes;
    MyMatrix<T> MatEquiv = FindTransformation<T, Telt>(EXText1, EXText2, eElt);
    if (IsIntegralMatrix(MatEquiv))
      return extend(MatEquiv);
    Tgroup GRPisom1 =
        GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat1, os);
    MyMatrix<T> EXTextInt1 = RemoveFractionMatrix(EXText1);
    MyMatrix<T> EXTextInt2 = RemoveFractionMatrix(EXText2);
    std::optional<MyMatrix<T>> opt =
        LinPolytopeIntegral_Isomorphism<T, Tgroup>(EXTextInt1, EXTextInt2,
                                                   GRPisom1, eElt, os);
    if (!opt)
      return {};
    return extend(*opt);
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
    T lfactors;
    std::vector<std::pair<T, int>> dist_occ;
    bool operator==(DelInvariant const &o) const {
      return nbVert == o.nbVert && rank == o.rank &&
             iso_eq_center == o.iso_eq_center && sqr_radius == o.sqr_radius &&
             lfactors == o.lfactors && dist_occ == o.dist_occ;
    }
  };

  DelInvariant delaunay_invariant(MyMatrix<T> const &EXT) const {
    DelInvariant inv;
    int nbVert = EXT.rows();
    inv.nbVert = nbVert;
    inv.rank = RankMat(EXT);
    MyVector<T> cp = center_homog(EXT);
    inv.sqr_radius = square_radius(EXT);
    inv.lfactors = list_factors(cp);
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
    std::vector<T> list_iso_center;
    bool operator==(PairInvariant const &o) const {
      return inv_over == o.inv_over && list_iso_center == o.list_iso_center;
    }
  };

  PairInvariant invariant_of_pair(MyMatrix<T> const &EXT,
                                  MyMatrix<T> const &EXTover) const {
    PairInvariant pinv;
    pinv.inv_over = delaunay_invariant(EXTover);
    MyVector<T> isoOver = isobarycenter_homog(EXTover);
    MyVector<T> isoEXT = isobarycenter_homog(EXT);
    for (int DelVal : {3, 5, 7, 11}) {
      T dv(DelVal);
      MyVector<T> eIso = (isoOver + dv * isoEXT) / (dv + T(1));
      pinv.list_iso_center.push_back(list_factors(eIso));
    }
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
    MyMatrix<T> SpaceBasInt = RemoveFractionMatrix(SpaceBas);
    MyMatrix<T> Prod = GramMat * SpaceBasInt.transpose();
    MyMatrix<T> NSP = NullspaceIntMat(Prod);
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
    // List of unit Gram matrices and their coefficients.
    std::vector<MyMatrix<T>> ListGramMat;
    std::vector<T> ListCoef;
    for (int i = 0; i < nRel; i++) {
      MyMatrix<T> H = ZeroMatrix<T>(nRel, nRel);
      H(i, i) = 1;
      ListGramMat.push_back(H);
      ListCoef.push_back(T(1));
    }
    for (int i = 0; i < nRel - 1; i++) {
      for (int j = i + 1; j < nRel; j++) {
        MyMatrix<T> H = ZeroMatrix<T>(nRel, nRel);
        H(i, j) = 1;
        H(j, i) = 1;
        ListGramMat.push_back(H);
        ListCoef.push_back(T(1) / T(2));
      }
    }
    T IntDeg0(0);
    MyVector<T> IntDeg1 = ZeroVector<T>(nRel);
    MyMatrix<T> IntDeg2 = ZeroMatrix<T>(nRel, nRel);
    vectface trig = lrs::GetTriangulation(EXTinBasis);
    // The covariance of the uniform distribution on a simplex with n+1 vertices
    // is (1/((n+1)(n+2))) sum_u (v_u - c)(v_u - c)^T; Tnp1 = n+1 is also the
    // barycenter divisor (vertex count). These are the linear factors (n+1),
    // (n+2), NOT factorials -- the overall n! normalization is applied at the
    // very end via factorial(nRel).
    T Tnp1(nRel + 1), Tnp2(nRel + 2);
    for (auto &eSimplex : trig) {
      std::vector<int> LV = FaceToVector<int>(eSimplex);
      int sdim = LV.size();
      MyMatrix<T> EXTsimplex(sdim, nRel + 1);
      for (int u = 0; u < sdim; u++)
        for (int j = 0; j < nRel + 1; j++)
          EXTsimplex(u, j) = EXTinBasis(LV[u], j);
      T VolSimplex = DeterminantMat(EXTsimplex);
      if (VolSimplex < 0)
        VolSimplex = -VolSimplex;
      MyVector<T> bary = ZeroVector<T>(nRel + 1);
      for (int u = 0; u < sdim; u++)
        for (int j = 0; j < nRel + 1; j++)
          bary(j) += EXTsimplex(u, j);
      bary /= Tnp1;
      // Spatial coordinates of the barycenter, and of the origin relative to it.
      MyVector<T> barySpace = bary.tail(nRel);
      MyVector<T> eDiff0 = -barySpace;
      // Centered vertex coordinates (rows = v_u - c, spatial part only).
      MyMatrix<T> Diff(sdim, nRel);
      for (int u = 0; u < sdim; u++)
        for (int a = 0; a < nRel; a++)
          Diff(u, a) = EXTsimplex(u, a + 1) - barySpace(a);
      for (size_t iMat = 0; iMat < ListCoef.size(); iMat++) {
        MyMatrix<T> const &H = ListGramMat[iMat];
        // sum_u (v_u - c)^T H (v_u - c) over the simplex vertices.
        T TheSum(0);
        for (int u = 0; u < sdim; u++) {
          MyVector<T> eDiff = GetMatrixRow(Diff, u);
          TheSum += eDiff.dot(H * eDiff);
        }
        // c^T H c (value of the quadratic form at the barycenter).
        T quad0 = eDiff0.dot(H * eDiff0);
        T TheInt = quad0 + TheSum / (Tnp1 * Tnp2);
        T fact = VolSimplex * ListCoef[iMat] * TheInt;
        IntDeg2 += fact * H;
      }
      IntDeg0 += VolSimplex;
      for (int a = 0; a < nRel; a++)
        IntDeg1(a) += VolSimplex * barySpace(a);
    }
    MyMatrix<T> ret = ZeroMatrix<T>(nRel + 1, nRel + 1);
    for (int i = 0; i < nRel; i++)
      for (int j = 0; j < nRel; j++)
        ret(i + 1, j + 1) = IntDeg2(i, j);
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
              apply_matrix(center_homog(TheEXT), eRecord.eMat);
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
                    optElt, "soft_computation: RepresentativeAction fails");
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
                      opt, "soft_computation: vertex not found in adjacent");
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
      MyVector<T> eCenter = apply_matrix(center_homog(TheEXT), eSoft.eMat);
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
    // The marked-center face stabilizer/equivalence needs an
    // automorphism-invariant spanning vector family (the Delaunay enumeration
    // itself works with full-dimensional polytopes, so it leaves SHV empty).
    if (data.SHV.rows() == 0) {
      MyMatrix<Tint> SHV_i =
          ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
      data.SHV = UniversalMatrixConversion<T, Tint>(SHV_i);
    }
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
