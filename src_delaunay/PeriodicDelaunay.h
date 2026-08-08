// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_PERIODICDELAUNAY_H_
#define SRC_DELAUNAY_PERIODICDELAUNAY_H_

// clang-format off
#include "PeriodicStructures.h"
#include "GRP_GroupFct.h"
#include "LatticeDelaunay.h"
#include "Shvec_exact.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERIODIC_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERIODIC_DELAUNAY
#endif

/*
  Closest vectors of a periodic point set Z^n + {c_1, ..., c_m} to a point.

  The points are expressed by their numerators over the denominator N of
  the point set: the point of numerator u is u / N, and u is congruent
  modulo N to the numerator of the coset it belongs to. Keeping the scaled
  form means the vertex matrices of the periodic Delaunay polytopes are
  integral, exactly like in the lattice case.
 */
template <typename T, typename Tint> struct PeriodicCVPResult {
  T TheNorm;
  MyMatrix<Tint> ListVectScaled;
};

/*
  The closest points of Z^n + c_k to v are c_k plus the closest points of
  Z^n to v - c_k, so the periodic problem is m lattice problems whose
  results are merged on the smallest norm. That merge is the only periodic
  part: the enumeration itself is the existing CVP solver.
 */
template <typename T, typename Tint>
PeriodicCVPResult<T, Tint>
PeriodicClosestVectors(CVPSolver<T, Tint> const &solver,
                       PeriodicPointSet<Tint> const &pps,
                       MyVector<T> const &eV) {
  int n = eV.size();
  int n_coset = pps.cosets_num.rows();
  T N_T = UniversalScalarConversion<T, Tint>(pps.N);
  std::optional<T> min_norm;
  std::vector<MyVector<Tint>> l_vect;
  for (int k = 0; k < n_coset; k++) {
    MyVector<T> eV_shift(n);
    for (int j = 0; j < n; j++) {
      T num = UniversalScalarConversion<T, Tint>(pps.cosets_num(k, j));
      eV_shift(j) = eV(j) - num / N_T;
    }
    resultCVP<T, Tint> res = solver.nearest_vectors(eV_shift);
    if (min_norm && res.TheNorm > *min_norm) {
      continue;
    }
    if (!min_norm || res.TheNorm < *min_norm) {
      min_norm = res.TheNorm;
      l_vect.clear();
    }
    int n_vect = res.ListVect.rows();
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      // The point is c_k + z, of numerator cosets_num(k) + N z.
      MyVector<Tint> u(n);
      for (int j = 0; j < n; j++) {
        u(j) = pps.cosets_num(k, j) + pps.N * res.ListVect(i_vect, j);
      }
      l_vect.push_back(std::move(u));
    }
  }
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (!min_norm) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicClosestVectors: the point set "
                 "has no coset at all\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> ListVectScaled = MatrixFromVectorFamily(l_vect);
  return {*min_norm, std::move(ListVectScaled)};
}

/*
  Whether an affine transformation preserves a periodic point set, and the
  permutation it induces on the cosets when it does.

  With p = c_k + z and the transformation p -> p A + w / d, the image is
  (c_k A + w / d) + z A. It stays in the point set for every z if and only
  if A is unimodular (the lattice Z^n is preserved) and c_k A + w / d is a
  coset modulo Z^n. Multiplied by the denominator N of the point set, that
  last condition is the integral congruence

      cosets_num(k) A + N w / d = cosets_num(sigma(k))   modulo N

  so the whole test is integer arithmetic, provided N w / d is integral --
  which is the invariant of the periodic setting, every transformation of
  the point set having its translation in (1 / N) Z^n.
 */
template <typename Tring>
std::optional<std::vector<size_t>>
PeriodicCosetPermutation(PeriodicPointSet<Tring> const &pps,
                         PeriodicAffineTransform<Tring> const &x) {
  int n = pps.cosets_num.cols();
  int n_coset = pps.cosets_num.rows();
  Tring det = DeterminantMat(x.A);
  if (det != 1 && det != -1) {
    return {};
  }
  // The scaled translation N w / d, integral exactly when the translation
  // lies in (1 / N) Z^n.
  MyVector<Tring> w_scal(n);
  for (int j = 0; j < n; j++) {
    Tring prod = pps.N * x.w(j);
    Tring res = ResInt(prod, x.d);
    if (res != 0) {
      return {};
    }
    w_scal(j) = QuoInt(prod, x.d);
  }
  std::vector<size_t> sigma(n_coset);
  Face f_hit(n_coset);
  for (int k = 0; k < n_coset; k++) {
    MyVector<Tring> img(n);
    for (int j = 0; j < n; j++) {
      Tring eSum = w_scal(j);
      for (int i = 0; i < n; i++) {
        AddMul(eSum, pps.cosets_num(k, i), x.A(i, j));
      }
      img(j) = eSum;
    }
    std::optional<size_t> opt = GetCosetIndex(pps, img);
    if (!opt) {
      return {};
    }
    sigma[k] = *opt;
    f_hit[*opt] = 1;
  }
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  // The induced map on the cosets is injective, being the restriction of a
  // bijection of the point set, so all the cosets must be reached.
  if (static_cast<int>(f_hit.count()) != n_coset) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicCosetPermutation: the induced "
                 "map on the cosets is not a bijection\n";
    throw TerminalException{1};
  }
#endif
  return sigma;
}

template <typename Tring>
bool IsPeriodicPointSetPreserved(PeriodicPointSet<Tring> const &pps,
                                 PeriodicAffineTransform<Tring> const &x) {
  return PeriodicCosetPermutation(pps, x).has_value();
}

/*
  The action of the transformations on the classes of (Z / N)^n. The
  condition of preserving the point set lives there: a transformation
  preserves it exactly when the induced permutation of the classes maps the
  cosets onto the cosets. Shared by the subgroup computation and by the
  equivalence.
 */
template <typename Tring> struct PeriodicClassAction {
  PeriodicPointSet<Tring> const &pps;
  int n;
  size_t N_s;
  size_t n_class;
  PeriodicClassAction(PeriodicPointSet<Tring> const &_pps)
      : pps(_pps), n(_pps.cosets_num.cols()),
        N_s(UniversalScalarConversion<size_t, Tring>(_pps.N)), n_class(1) {
    for (int i = 0; i < n; i++) {
      n_class *= N_s;
    }
  }
  // The classes are indexed in the mixed radix order.
  size_t index(MyVector<Tring> const &u) const {
    size_t idx = 0;
    for (int i = n - 1; i >= 0; i--) {
      Tring res = ResInt(u(i), pps.N);
      idx = idx * N_s + UniversalScalarConversion<size_t, Tring>(res);
    }
    return idx;
  }
  MyVector<Tring> get_vector(size_t idx) const {
    MyVector<Tring> u(n);
    for (int i = 0; i < n; i++) {
      u(i) = UniversalScalarConversion<Tring, size_t>(idx % N_s);
      idx /= N_s;
    }
    return u;
  }
  // The permutation of the classes induced by a transformation, which is
  // defined as soon as the linear part is unimodular and the translation
  // lies in (1 / N) Z^n.
  std::vector<size_t>
  permutation(PeriodicAffineTransform<Tring> const &x) const {
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
    Tring det = DeterminantMat(x.A);
    if (det != 1 && det != -1) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicClassAction: the linear part "
                   "is not unimodular\n";
      throw TerminalException{1};
    }
    for (int j = 0; j < n; j++) {
      Tring prod = pps.N * x.w(j);
      if (ResInt(prod, x.d) != 0) {
        std::cerr << "PERIODIC_DELAUNAY: PeriodicClassAction: the translation "
                     "is outside of (1/N) Z^n\n";
        throw TerminalException{1};
      }
    }
#endif
    MyVector<Tring> w_scal(n);
    for (int j = 0; j < n; j++) {
      Tring prod = pps.N * x.w(j);
      w_scal(j) = QuoInt(prod, x.d);
    }
    std::vector<size_t> ret(n_class);
    for (size_t i_class = 0; i_class < n_class; i_class++) {
      MyVector<Tring> u = get_vector(i_class);
      MyVector<Tring> img(n);
      for (int j = 0; j < n; j++) {
        Tring eSum = w_scal(j);
        for (int i = 0; i < n; i++) {
          AddMul(eSum, u(i), x.A(i, j));
        }
        img(j) = eSum;
      }
      ret[i_class] = index(img);
    }
    return ret;
  }
  // The classes that are cosets, over an index space where the classes
  // start at shift.
  Face coset_face(size_t shift, size_t n_tot) const {
    Face f(n_tot);
    int n_coset = pps.cosets_num.rows();
    for (int k = 0; k < n_coset; k++) {
      MyVector<Tring> u = GetMatrixRow(pps.cosets_num, k);
      f[shift + index(u)] = 1;
    }
    return f;
  }
};

/*
  A group of transformations made to act on the disjoint union of the
  polytope vertices and the classes of (Z / N)^n. Preserving the point set
  is not a condition that can be checked generator by generator, the
  failing elements not forming a complement, but on this union it becomes a
  condition on sets: the subgroup preserving the point set is the setwise
  stabilizer of the cosets, and the search for an equivalence preserving it
  is a RepresentativeAction on the same sets.

  f_trans maps a permutation of the vertices to the transformation it
  realizes.
 */
template <typename Tring, typename Tgroup> struct PeriodicUnionGroup {
  Tgroup GRPbig;
  Face f_coset;
  Face f_vert;
};

template <typename Tring, typename Tgroup, typename Fget_trans>
PeriodicUnionGroup<Tring, Tgroup>
BuildPeriodicUnionGroup(PeriodicClassAction<Tring> const &act,
                        Tgroup const &GRP, Fget_trans f_trans) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  size_t n_vert = GRP.n_act();
  size_t n_tot = n_vert + act.n_class;
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (n_tot > static_cast<size_t>(std::numeric_limits<Tidx>::max())) {
    std::cerr << "PERIODIC_DELAUNAY: BuildPeriodicUnionGroup: the N^n classes "
                 "overflow the permutation index type (n_class="
              << act.n_class << ")\n";
    throw TerminalException{1};
  }
#endif
  std::vector<Telt> ListGenBig;
  for (auto &eGen : GRP.SmallGeneratingSet()) {
    std::vector<size_t> cperm = act.permutation(f_trans(eGen));
    std::vector<Tidx> eList(n_tot);
    for (size_t i = 0; i < n_vert; i++) {
      eList[i] = static_cast<Tidx>(OnPoints(static_cast<Tidx>(i), eGen));
    }
    for (size_t i_class = 0; i_class < act.n_class; i_class++) {
      eList[n_vert + i_class] = static_cast<Tidx>(n_vert + cperm[i_class]);
    }
    ListGenBig.push_back(Telt(eList));
  }
  Face f_vert(n_tot);
  for (size_t i = 0; i < n_vert; i++) {
    f_vert[i] = 1;
  }
  return {Tgroup(ListGenBig, static_cast<Tidx>(n_tot)),
          act.coset_face(n_vert, n_tot), std::move(f_vert)};
}

// The subgroup of the transformations preserving the periodic point set:
// the setwise stabilizer of the cosets on the union, restricted back to
// the vertices. The symmetries permuting the cosets are kept, which is the
// point of doing it this way rather than by colouring the vertices with
// their coset.
template <typename Tring, typename Tgroup, typename Fget_trans>
Tgroup PeriodicCosetPreservingSubgroup(PeriodicPointSet<Tring> const &pps,
                                       Tgroup const &GRP, Fget_trans f_trans,
                                       [[maybe_unused]] std::ostream &os) {
  PeriodicClassAction<Tring> act(pps);
  PeriodicUnionGroup<Tring, Tgroup> pug =
      BuildPeriodicUnionGroup<Tring, Tgroup, Fget_trans>(act, GRP, f_trans);
  Tgroup GRPstab = pug.GRPbig.Stabilizer_OnSets(pug.f_coset);
  return ReducedGroupActionFace(GRPstab, pug.f_vert);
}

/*
  Norm of the difference between the point of numerator u and the point eV,
  for the quadratic form GramMat.
 */
template <typename T, typename Tint>
T PeriodicNormDiff(MyMatrix<T> const &GramMat, Tint const &N,
                   MyVector<Tint> const &u, MyVector<T> const &eV) {
  int n = eV.size();
  T N_T = UniversalScalarConversion<T, Tint>(N);
  MyVector<T> eDiff(n);
  for (int j = 0; j < n; j++) {
    T num = UniversalScalarConversion<T, Tint>(u(j));
    eDiff(j) = num / N_T - eV(j);
  }
  return EvaluationQuadForm<T, T>(GramMat, eDiff);
}

/*
  Everything a Delaunay polytope of a periodic point set needs for its
  isometry computations, built once and used by both the stabilizer and
  the equivalence: the vertices in homogeneous unscaled form, their
  circumcenter, the extended set of the recentered vertices followed by
  the short vectors, the face of the vertex block, and the weight matrix
  of the extended set with the two blocks separated by their diagonal
  value.

  Two facts of the periodic setting shape what follows:
  --- the short vectors are required to generate Z^n, so an isometry
  preserving that block permutes vectors generating the lattice and its
  linear part is automatically in GL_n(Z): the integral stabilizer
  computation of the lattice case is not needed here,
  --- the circumcenter is mapped to the circumcenter, so the affine map of
  a permutation is x -> (x - c1) A + c2, whose translation c2 - c1 A lies
  in (1 / N) Z^n: the precondition of the class action holds by
  construction.
 */
template <typename T, typename Tidx_value> struct PeriodicDelaunayIsoData {
  MyMatrix<T> EXT_T;
  MyVector<T> Cent;
  MyMatrix<T> EXText;
  // The recentered vertices alone: the permutations handed around act on
  // the vertex block, so this is what determines the linear part.
  MyMatrix<T> EXTvert;
  Face eFace;
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat;
};

template <typename T, typename Tring, typename Tidx, typename Tidx_value>
PeriodicDelaunayIsoData<T, Tidx_value>
BuildPeriodicDelaunayIsoData(MyMatrix<T> const &GramMat,
                             PeriodicPointSet<Tring> const &pps,
                             MyMatrix<T> const &SHV,
                             MyMatrix<Tring> const &EXT_scaled,
                             std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n = GramMat.rows();
  int nbVert = EXT_scaled.rows();
  T N_T = UniversalScalarConversion<T, Tring>(pps.N);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (RankMat(SHV) != n) {
    std::cerr << "PERIODIC_DELAUNAY: BuildPeriodicDelaunayIsoData: the short "
                 "vectors do not span the space, so the linear parts are not "
                 "forced to be integral\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> EXT_T(nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    EXT_T(iVert, 0) = T(1);
    for (int i = 0; i < n; i++) {
      EXT_T(iVert, i + 1) =
          UniversalScalarConversion<T, Tring>(EXT_scaled(iVert, i)) / N_T;
    }
  }
  MyVector<T> Cent = get_reduced_center(GramMat, EXT_T);
  MyMatrix<T> EXText = get_reduced_delaunay_shv(EXT_T, GramMat, SHV, Cent);
  std::vector<MyMatrix<T>> ListMat{GramMat};
  std::vector<T> Vdiag = get_vdiag_delaunay_shv<T>(EXT_T, SHV);
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListMat_Vdiag<T, Tfield, Tidx, Tidx_value>(
          EXText, ListMat, Vdiag, os);
  Face eFace = get_face_delaunay_shv(EXT_T, SHV);
  MyMatrix<T> EXTvert(nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      EXTvert(iVert, i) = EXText(iVert, i);
    }
  }
  return {std::move(EXT_T),   std::move(Cent),  std::move(EXText),
          std::move(EXTvert), std::move(eFace), std::move(WMat)};
}

// The transformation realized by a permutation of the extended set of
// data1 onto the extended set of data2: the linear part from the
// recentered vertices, the translation from the circumcenters.
template <typename T, typename Tring, typename Telt, typename Tidx_value>
PeriodicAffineTransform<Tring>
PeriodicTransformFromPerm(PeriodicDelaunayIsoData<T, Tidx_value> const &data1,
                          PeriodicDelaunayIsoData<T, Tidx_value> const &data2,
                          Telt const &eElt, Tring const &N) {
  int n = data1.Cent.size();
  MyMatrix<T> A_T =
      FindTransformation<T, Telt>(data1.EXTvert, data2.EXTvert, eElt);
  MyMatrix<Tring> A = UniversalMatrixConversion<Tring, T>(A_T);
  T N_T = UniversalScalarConversion<T, Tring>(N);
  MyVector<Tring> w(n);
  for (int j = 0; j < n; j++) {
    T eSum = data2.Cent(j);
    for (int i = 0; i < n; i++) {
      eSum -= data1.Cent(i) * A_T(i, j);
    }
    w(j) = UniversalScalarConversion<Tring, T>(eSum * N_T);
  }
  return {std::move(A), std::move(w), N};
}

// The isometries of the extended set that preserve the vertex block,
// acting on the vertices. Their linear parts are integral, the short
// vectors generating the lattice.
template <typename T, typename Tgroup, typename Tidx_value>
Tgroup PeriodicDelaunayIsometryGroup(
    PeriodicDelaunayIsoData<T, Tidx_value> const &data, std::ostream &os) {
  using Tgr = GraphListAdj;
  Tgroup PreGRPisom =
      GetStabilizerWeightMatrix<std::vector<T>, Tgr, Tgroup, Tidx_value>(
          data.WMat, os);
  Tgroup GRPisomExt = PreGRPisom.Stabilizer_OnSets(data.eFace);
  return ReducedGroupActionFace(GRPisomExt, data.eFace);
}

/*
  The stabilizer of a Delaunay polytope of a periodic point set: the
  isometries of the extended set, restricted to those preserving the point
  set.
 */
template <typename T, typename Tring, typename Tgroup>
Tgroup PeriodicDelaunay_Stabilizer(MyMatrix<T> const &GramMat,
                                   PeriodicPointSet<Tring> const &pps,
                                   MyMatrix<T> const &SHV,
                                   MyMatrix<Tring> const &EXT_scaled,
                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tidx_value = int16_t;
  PeriodicDelaunayIsoData<T, Tidx_value> data =
      BuildPeriodicDelaunayIsoData<T, Tring, Tidx, Tidx_value>(
          GramMat, pps, SHV, EXT_scaled, os);
  Tgroup GRPisom = PeriodicDelaunayIsometryGroup<T, Tgroup, Tidx_value>(data, os);
  auto f_trans = [&](Telt const &eElt) -> PeriodicAffineTransform<Tring> {
    return PeriodicTransformFromPerm<T, Tring, Telt, Tidx_value>(data, data,
                                                                 eElt, pps.N);
  };
  return PeriodicCosetPreservingSubgroup<Tring, Tgroup, decltype(f_trans)>(
      pps, GRPisom, f_trans, os);
}

/*
  An equivalence between two Delaunay polytopes of a periodic point set,
  as a transformation preserving the point set and mapping the first onto
  the second.

  The weight matrices give one isometry m0 of the extended sets, and the
  isometries mapping the first onto the second are exactly the h * m0 for
  h in the isometry group of the first. Requiring that h * m0 preserve the
  point set reads, on the classes, h(S) = m0^-1(S) with S the set of the
  cosets, so the search is a RepresentativeAction on the union group of
  the first polytope. Testing m0 alone would miss the equivalences, which
  is what makes this more than a single check.
 */
template <typename T, typename Tring, typename Tgroup>
std::optional<PeriodicAffineTransform<Tring>>
PeriodicDelaunay_TestEquivalence(MyMatrix<T> const &GramMat,
                                 PeriodicPointSet<Tring> const &pps,
                                 MyMatrix<T> const &SHV,
                                 MyMatrix<Tring> const &EXT1_scaled,
                                 MyMatrix<Tring> const &EXT2_scaled,
                                 std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tidx_value = int16_t;
  PeriodicDelaunayIsoData<T, Tidx_value> data1 =
      BuildPeriodicDelaunayIsoData<T, Tring, Tidx, Tidx_value>(
          GramMat, pps, SHV, EXT1_scaled, os);
  PeriodicDelaunayIsoData<T, Tidx_value> data2 =
      BuildPeriodicDelaunayIsoData<T, Tring, Tidx, Tidx_value>(
          GramMat, pps, SHV, EXT2_scaled, os);
  std::optional<Telt> eRes =
      TestEquivalenceWeightMatrix<std::vector<T>, Telt, Tidx_value>(
          data1.WMat, data2.WMat, os);
  if (!eRes) {
    return {};
  }
  PeriodicAffineTransform<Tring> m0 =
      PeriodicTransformFromPerm<T, Tring, Telt, Tidx_value>(data1, data2,
                                                            *eRes, pps.N);
  if (IsPeriodicPointSetPreserved(pps, m0)) {
    return m0;
  }
  // The other isometries of the pair are the h * m0, so the point set is
  // preserved by one of them exactly when some h sends the cosets to their
  // preimage by m0.
  PeriodicClassAction<Tring> act(pps);
  Tgroup GRPisom1 = PeriodicDelaunayIsometryGroup<T, Tgroup, Tidx_value>(data1, os);
  auto f_trans1 = [&](Telt const &eElt) -> PeriodicAffineTransform<Tring> {
    return PeriodicTransformFromPerm<T, Tring, Telt, Tidx_value>(data1, data1,
                                                                 eElt, pps.N);
  };
  PeriodicUnionGroup<Tring, Tgroup> pug =
      BuildPeriodicUnionGroup<Tring, Tgroup, decltype(f_trans1)>(act, GRPisom1,
                                                                 f_trans1);
  size_t n_vert = GRPisom1.n_act();
  size_t n_tot = n_vert + act.n_class;
  std::vector<size_t> cperm_m0 = act.permutation(m0);
  Face f_target(n_tot);
  for (size_t i_class = 0; i_class < act.n_class; i_class++) {
    if (pug.f_coset[n_vert + cperm_m0[i_class]] == 1) {
      f_target[n_vert + i_class] = 1;
    }
  }
  std::optional<Telt> opt =
      pug.GRPbig.RepresentativeAction_OnSets(pug.f_coset, f_target);
  if (!opt) {
    return {};
  }
  std::vector<Tidx> eList_vert(n_vert);
  for (size_t i = 0; i < n_vert; i++) {
    eList_vert[i] = static_cast<Tidx>(OnPoints(static_cast<Tidx>(i), *opt));
  }
  PeriodicAffineTransform<Tring> h = f_trans1(Telt(eList_vert));
  PeriodicAffineTransform<Tring> ret = h * m0;
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (!IsPeriodicPointSetPreserved(pps, ret)) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicDelaunay_TestEquivalence: the "
                 "equivalence found does not preserve the point set\n";
    throw TerminalException{1};
  }
#endif
  return ret;
}

// clang-format off
#endif  // SRC_DELAUNAY_PERIODICDELAUNAY_H_
// clang-format on
