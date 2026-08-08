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
  The subgroup of the transformations preserving the periodic point set.

  Preserving the point set is not a condition that can be tested generator
  by generator: the elements failing it are not a subgroup complement, so
  filtering the generators would not give the subgroup. It is however the
  setwise stabilizer of the cosets under the action on the finite group
  (Z / N)^n induced by the transformations, so it is obtained by letting
  the group act on the disjoint union of the polytope vertices and the N^n
  classes, taking the setwise stabilizer of the classes that are cosets,
  and restricting the result back to the vertices.

  f_trans maps a permutation of the vertices to the transformation it
  realizes. Every generator must have a unimodular linear part and a
  translation in (1 / N) Z^n, so that its action on the classes is defined;
  that holds once the group has been restricted to the transformations
  preserving the lattice.
 */
template <typename Tring, typename Tgroup, typename Fget_trans>
Tgroup PeriodicCosetPreservingSubgroup(PeriodicPointSet<Tring> const &pps,
                                       Tgroup const &GRP, Fget_trans f_trans,
                                       [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = pps.cosets_num.cols();
  int n_coset = pps.cosets_num.rows();
  size_t n_vert = GRP.n_act();
  size_t N_s = UniversalScalarConversion<size_t, Tring>(pps.N);
  // The classes of (Z / N)^n in the mixed radix order, which is the index
  // used for the second block of the permutations.
  size_t n_class = 1;
  for (int i = 0; i < n; i++) {
    n_class *= N_s;
  }
  size_t n_tot = n_vert + n_class;
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (n_tot > static_cast<size_t>(std::numeric_limits<Tidx>::max())) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicCosetPreservingSubgroup: the "
                 "N^n classes overflow the permutation index type (n_class="
              << n_class << ")\n";
    throw TerminalException{1};
  }
#endif
  auto class_index = [&](MyVector<Tring> const &u) -> size_t {
    size_t idx = 0;
    for (int i = n - 1; i >= 0; i--) {
      Tring res = ResInt(u(i), pps.N);
      idx = idx * N_s + UniversalScalarConversion<size_t, Tring>(res);
    }
    return idx;
  };
  auto class_vector = [&](size_t idx) -> MyVector<Tring> {
    MyVector<Tring> u(n);
    for (int i = 0; i < n; i++) {
      u(i) = UniversalScalarConversion<Tring, size_t>(idx % N_s);
      idx /= N_s;
    }
    return u;
  };
  std::vector<Telt> ListGenBig;
  for (auto &eGen : GRP.SmallGeneratingSet()) {
    PeriodicAffineTransform<Tring> x = f_trans(eGen);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
    // The action on the classes is defined only for the transformations
    // preserving the lattice: the group has to be restricted to those
    // before being handed over.
    Tring det = DeterminantMat(x.A);
    if (det != 1 && det != -1) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicCosetPreservingSubgroup: a "
                   "generator has a non unimodular linear part\n";
      throw TerminalException{1};
    }
    for (int j = 0; j < n; j++) {
      Tring prod = pps.N * x.w(j);
      if (ResInt(prod, x.d) != 0) {
        std::cerr << "PERIODIC_DELAUNAY: PeriodicCosetPreservingSubgroup: a "
                     "generator has a translation outside of (1/N) Z^n\n";
        throw TerminalException{1};
      }
    }
#endif
    MyVector<Tring> w_scal(n);
    for (int j = 0; j < n; j++) {
      Tring prod = pps.N * x.w(j);
      w_scal(j) = QuoInt(prod, x.d);
    }
    std::vector<Tidx> eList(n_tot);
    for (size_t i = 0; i < n_vert; i++) {
      eList[i] = static_cast<Tidx>(OnPoints(static_cast<Tidx>(i), eGen));
    }
    for (size_t i_class = 0; i_class < n_class; i_class++) {
      MyVector<Tring> u = class_vector(i_class);
      MyVector<Tring> img(n);
      for (int j = 0; j < n; j++) {
        Tring eSum = w_scal(j);
        for (int i = 0; i < n; i++) {
          AddMul(eSum, u(i), x.A(i, j));
        }
        img(j) = eSum;
      }
      eList[n_vert + i_class] =
          static_cast<Tidx>(n_vert + class_index(img));
    }
    ListGenBig.push_back(Telt(eList));
  }
  Tgroup GRPbig(ListGenBig, static_cast<Tidx>(n_tot));
  Face f_coset(n_tot);
  for (int k = 0; k < n_coset; k++) {
    MyVector<Tring> u = GetMatrixRow(pps.cosets_num, k);
    f_coset[n_vert + class_index(u)] = 1;
  }
  Tgroup GRPstab = GRPbig.Stabilizer_OnSets(f_coset);
  Face f_vert(n_tot);
  for (size_t i = 0; i < n_vert; i++) {
    f_vert[i] = 1;
  }
  return ReducedGroupActionFace(GRPstab, f_vert);
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
  The stabilizer of a Delaunay polytope of a periodic point set.

  The scheme is the one of PolytopeGen_StabilizerKernel: the vertices,
  recentered at the circumcenter, are extended by the short vectors, the
  two blocks are separated by their diagonal value, and the isometries of
  the result are computed. Two facts make the rest simpler than in the
  lattice case:
  --- the short vectors are required to generate Z^n, so an isometry
  preserving that block permutes vectors generating the lattice and its
  linear part is automatically in GL_n(Z): the integral stabilizer
  computation of the lattice case is not needed,
  --- the circumcenter is fixed by every symmetry of the polytope, so the
  affine map is x -> (x - c) A + c and the translation is c - c A, which
  lies in (1/N) Z^n with c: the precondition of the coset restriction
  holds by construction.
  What remains is the restriction to the transformations preserving the
  point set, which is PeriodicCosetPreservingSubgroup.

  EXT_scaled holds the numerators of the vertices over the denominator N.
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
  using Tgr = GraphListAdj;
  using Tfield = typename overlying_field<T>::field_type;
  int n = GramMat.rows();
  int nbVert = EXT_scaled.rows();
  T N_T = UniversalScalarConversion<T, Tring>(pps.N);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  // The short vectors have to generate Z^n for the linear parts to be
  // integral without further work.
  if (RankMat(SHV) != n || !IsIntegralMatrix(SHV) ||
      Tring(1) != UniversalScalarConversion<Tring, T>(
                      T(DeterminantMat(GetZbasis(SHV))))) {
    if (RankMat(SHV) != n) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicDelaunay_Stabilizer: the short "
                   "vectors do not span the space\n";
      throw TerminalException{1};
    }
  }
#endif
  // The vertices in the homogeneous unscaled form used by the geometry.
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
  Tgroup PreGRPisom =
      GetStabilizerWeightMatrix<std::vector<T>, Tgr, Tgroup, Tidx_value>(WMat,
                                                                        os);
  Tgroup GRPisomExt = PreGRPisom.Stabilizer_OnSets(eFace);
  Tgroup GRPisom = ReducedGroupActionFace(GRPisomExt, eFace);
  // The transformation realized by a permutation of the vertices: the
  // linear part from the recentered vertices, the translation from the
  // circumcenter being fixed.
  MyMatrix<T> EXTvert(nbVert, n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int i = 0; i < n; i++) {
      EXTvert(iVert, i) = EXT_T(iVert, i + 1) - Cent(i);
    }
  }
  auto f_trans = [&](Telt const &eElt) -> PeriodicAffineTransform<Tring> {
    MyMatrix<T> A_T = FindTransformation<T, Telt>(EXTvert, EXTvert, eElt);
    MyMatrix<Tring> A = UniversalMatrixConversion<Tring, T>(A_T);
    // Translation c - c A, scaled by N so as to be integral.
    MyVector<Tring> w(n);
    for (int j = 0; j < n; j++) {
      T eSum = Cent(j);
      for (int i = 0; i < n; i++) {
        eSum -= Cent(i) * A_T(i, j);
      }
      w(j) = UniversalScalarConversion<Tring, T>(eSum * N_T);
    }
    return {std::move(A), std::move(w), pps.N};
  };
  return PeriodicCosetPreservingSubgroup<Tring, Tgroup, decltype(f_trans)>(
      pps, GRPisom, f_trans, os);
}

// clang-format off
#endif  // SRC_DELAUNAY_PERIODICDELAUNAY_H_
// clang-format on
