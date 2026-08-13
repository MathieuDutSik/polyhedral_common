// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_PERIODICDELAUNAY_H_
#define SRC_DELAUNAY_PERIODICDELAUNAY_H_

// clang-format off
#include "PeriodicStructures.h"
#include "GRP_GroupFct.h"
#include "FundamentalDelaunay.h"
#include "IsoDelaunayDomains.h"
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
  Closest vectors of a periodic point set to a point.

  The point set is D Z^n + {c_0 = 0, c_1, ..., c_m} with the c_i integral,
  so its points ARE integral vectors -- those congruent to one of the c_i
  modulo D -- and there is a single frame throughout, no scaled and
  unscaled version of anything. The vertex matrices of the periodic
  Delaunay polytopes are integral for the same reason as in the lattice
  case, and a transformation preserving the set is an ordinary integral
  affine matrix.
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
  T D_T = UniversalScalarConversion<T, Tint>(pps.N);
  std::optional<T> min_norm;
  std::vector<MyVector<Tint>> l_vect;
  for (int k = 0; k < n_coset; k++) {
    // The points of the set on the coset c_k are the c_k + D y with y in
    // Z^n, so the solver of Z^n is asked about (eV - c_k) / D and the
    // distance it returns is the one sought divided by D^2.
    MyVector<T> eV_shift(n);
    for (int j = 0; j < n; j++) {
      T num = UniversalScalarConversion<T, Tint>(pps.cosets_num(k, j));
      eV_shift(j) = (eV(j) - num) / D_T;
    }
    resultCVP<T, Tint> res = solver.nearest_vectors(eV_shift);
    T norm = res.TheNorm * D_T * D_T;
    if (min_norm && norm > *min_norm) {
      continue;
    }
    if (!min_norm || norm < *min_norm) {
      min_norm = norm;
      l_vect.clear();
    }
    int n_vect = res.ListVect.rows();
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
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
  A solver over a periodic point set, with the interface the Delaunay
  geometry asks for: the member GramMat and the method nearest_vectors.
  With it, FindDelaunayPolytope and FindAdjacentDelaunayPolytope apply to
  a periodic point set without a line of their code being specific to it.

  Everything is expressed in the coordinates scaled by the denominator N
  of the point set, in which the points are the integral vectors
  congruent to a coset modulo N: the vertex matrices stay integral, as
  the geometry requires. Scaling the points by N and keeping the form
  amounts to scaling the form by N^2, which leaves the Delaunay
  decomposition unchanged, so the only thing to respect is that the
  distances returned here are the scaled ones -- the geometry compares
  them to circumradii computed from the same scaled coordinates.
 */
template <typename T, typename Tint> struct PeriodicCVPSolver {
  MyMatrix<T> GramMat;
  CVPSolver<T, Tint> solver;
  PeriodicPointSet<Tint> pps;
  int n;
  T N_T;
  PeriodicCVPSolver(MyMatrix<T> const &_GramMat,
                    PeriodicPointSet<Tint> const &_pps, std::ostream &_os)
      : GramMat(_GramMat), solver(_GramMat, _os), pps(_pps),
        n(_GramMat.rows()),
        N_T(UniversalScalarConversion<T, Tint>(_pps.N)) {
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
    // The geometry needs a point of the set at the origin, which
    // PeriodicPointSetFromRational arranges by translating the set.
    MyVector<Tint> zero = ZeroVector<Tint>(n);
    if (!GetCosetIndex(pps, zero)) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicCVPSolver: the point set does "
                   "not contain the origin, it has to be translated so that "
                   "one of its cosets is the origin\n";
      throw TerminalException{1};
    }
#endif
  }
  // In scaled coordinates the period lattice is N Z^n, so the moves
  // between points of the set are the unit vectors scaled by N.
  std::vector<MyVector<T>> get_seed_differences() const {
    std::vector<MyVector<T>> ret;
    for (int i = 0; i < n; i++) {
      MyVector<T> V1 = ZeroVector<T>(n);
      V1(i) = N_T;
      ret.push_back(V1);
      MyVector<T> V2 = ZeroVector<T>(n);
      V2(i) = -N_T;
      ret.push_back(V2);
    }
    return ret;
  }
  // The closest points of the point set to eV, one single frame.
  resultCVP<T, Tint> nearest_vectors(MyVector<T> const &eV) const {
    PeriodicCVPResult<T, Tint> res = PeriodicClosestVectors(solver, pps, eV);
    return {res.TheNorm, res.ListVectScaled};
  }
};

// The moves between points of a periodic point set, in scaled
// coordinates: the period lattice is N Z^n there, so the moves of the
// lattice have to be scaled by N as well.
template <typename Tint>
MyMatrix<Tint> PeriodicScaledGraverBasis(MyMatrix<Tint> const &ShvGraverBasis,
                                         Tint const &N) {
  int n_row = ShvGraverBasis.rows();
  int n_col = ShvGraverBasis.cols();
  MyMatrix<Tint> ret(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      ret(i, j) = N * ShvGraverBasis(i, j);
    }
  }
  return ret;
}

/*
  Whether an affine transformation preserves a periodic point set, and the
  permutation it induces on the cosets when it does.

  The point set being D Z^n + {c_0 = 0, c_1, ..., c_m} with the c_i
  integral, a transformation of it is an ordinary integral affine matrix M
  of order n+1 acting on the right of the rows (1, x): its linear part is
  the lower right block A and its translation the first row w. The image
  of p = c_k + D z is (c_k A + w) + D (z A), so the set is preserved
  exactly when A is unimodular, which is what keeps D Z^n, and c_k A + w
  is a coset modulo D. Everything is integer arithmetic.
 */
template <typename Tring>
struct PeriodicAffineParts {
  MyMatrix<Tring> A;
  MyVector<Tring> w;
};

// The linear part and the translation of an integral affine matrix, or
// nothing when it is not integral -- which simply means the candidate is
// not a transformation of the point set.
template <typename Tring, typename T>
std::optional<PeriodicAffineParts<Tring>>
GetPeriodicAffineParts(MyMatrix<T> const &M) {
  int n = M.rows() - 1;
  if (!IsIntegralMatrix(M)) {
    return {};
  }
  MyMatrix<Tring> A(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A(i, j) = UniversalScalarConversion<Tring, T>(M(i + 1, j + 1));
    }
  }
  MyVector<Tring> w(n);
  for (int j = 0; j < n; j++) {
    w(j) = UniversalScalarConversion<Tring, T>(M(0, j + 1));
  }
  return PeriodicAffineParts<Tring>{std::move(A), std::move(w)};
}

template <typename Tring, typename T>
std::optional<std::vector<size_t>>
PeriodicCosetPermutation(PeriodicPointSet<Tring> const &pps,
                         MyMatrix<T> const &M) {
  int n = pps.cosets_num.cols();
  int n_coset = pps.cosets_num.rows();
  std::optional<PeriodicAffineParts<Tring>> opt_parts =
      GetPeriodicAffineParts<Tring, T>(M);
  if (!opt_parts) {
    return {};
  }
  MyMatrix<Tring> const &A = opt_parts->A;
  MyVector<Tring> const &w = opt_parts->w;
  Tring det = DeterminantMat(A);
  if (det != 1 && det != -1) {
    return {};
  }
  std::vector<size_t> sigma(n_coset);
  Face f_hit(n_coset);
  for (int k = 0; k < n_coset; k++) {
    MyVector<Tring> img(n);
    for (int j = 0; j < n; j++) {
      Tring eSum = w(j);
      for (int i = 0; i < n; i++) {
        AddMul(eSum, pps.cosets_num(k, i), A(i, j));
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

template <typename Tring, typename T>
bool IsPeriodicPointSetPreserved(PeriodicPointSet<Tring> const &pps,
                                 MyMatrix<T> const &M) {
  return PeriodicCosetPermutation<Tring, T>(pps, M).has_value();
}

/*
  The affine extension of a linear part, if there is one: the smallest
  translation making x -> x A + w preserve the point set. The image of
  c_0 = 0 is the translation itself, so w has to be congruent to one of
  the cosets and only n_coset candidates exist, not D^n. Note that a
  linear part can fail on its own (it moves a coset off the coset set)
  and be adequate all the same, -Id on {0, c} being the basic example:
  the translation by c brings the image cosets {0, -c} back onto {0, c}.

  The linear part acts the same in the scaled and the unscaled frame, so
  the caller does not need to know about the scaling; the matrix returned
  is homogeneous in the scaled frame like every periodic transformation.
 */
template <typename Tring, typename T>
std::optional<MyMatrix<T>>
PeriodicAffineExtension(PeriodicPointSet<Tring> const &pps,
                        MyMatrix<T> const &LinPart) {
  int n = LinPart.rows();
  int n_coset = pps.cosets_num.rows();
  for (int k = 0; k < n_coset; k++) {
    MyMatrix<T> H = IdentityMat<T>(n + 1);
    for (int j = 0; j < n; j++) {
      H(0, j + 1) = UniversalScalarConversion<T, Tring>(pps.cosets_num(k, j));
      for (int i = 0; i < n; i++) {
        H(i + 1, j + 1) = LinPart(i, j);
      }
    }
    if (IsPeriodicPointSetPreserved<Tring, T>(pps, H)) {
      return H;
    }
  }
  return {};
}

/*
  The action of the transformations on the classes of (Z / D)^n. The
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
  // The permutation of the classes induced by a transformation, defined
  // as soon as it is one of the point set.
  template <typename T>
  std::vector<size_t> permutation(MyMatrix<T> const &M) const {
    std::optional<PeriodicAffineParts<Tring>> opt_parts =
        GetPeriodicAffineParts<Tring, T>(M);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
    // The action on the classes is defined for the transformations of the
    // point set only, so the group has to be restricted to those before
    // being handed over.
    if (!opt_parts) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicClassAction: the matrix is not "
                   "integral, so it is not a transformation of the set\n";
      throw TerminalException{1};
    }
    Tring det = DeterminantMat(opt_parts->A);
    if (det != 1 && det != -1) {
      std::cerr << "PERIODIC_DELAUNAY: PeriodicClassAction: the linear part "
                   "is not unimodular\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tring> const &A = opt_parts->A;
    MyVector<Tring> const &w = opt_parts->w;
    std::vector<size_t> ret(n_class);
    for (size_t i_class = 0; i_class < n_class; i_class++) {
      MyVector<Tring> u = get_vector(i_class);
      MyVector<Tring> img(n);
      for (int j = 0; j < n; j++) {
        Tring eSum = w(j);
        for (int i = 0; i < n; i++) {
          AddMul(eSum, u(i), A(i, j));
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
T PeriodicNormDiff(MyMatrix<T> const &GramMat, MyVector<Tint> const &u,
                   MyVector<T> const &eV) {
  int n = eV.size();
  MyVector<T> eDiff(n);
  for (int j = 0; j < n; j++) {
    eDiff(j) = UniversalScalarConversion<T, Tint>(u(j)) - eV(j);
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
  // The same vertices in the scaled frame, in which the point set is
  // D Z^n + {c_i} with the c_i integral. The transformations are read off
  // from those, so that they come out integral.
  MyMatrix<T> EXT_scaled_T;
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
  if (EXT_scaled.cols() != n + 1) {
    std::cerr << "PERIODIC_DELAUNAY: BuildPeriodicDelaunayIsoData: the vertex "
                 "matrix should be homogeneous, a leading column of 1 "
                 "followed by the scaled coordinates\n";
    throw TerminalException{1};
  }
#endif
  // EXT_scaled is homogeneous, as the Delaunay geometry produces it: the
  // leading column of 1 is kept and only the coordinates are unscaled.
  MyMatrix<T> EXT_T(nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    EXT_T(iVert, 0) = T(1);
    for (int i = 0; i < n; i++) {
      EXT_T(iVert, i + 1) =
          UniversalScalarConversion<T, Tring>(EXT_scaled(iVert, i + 1)) / N_T;
    }
  }
  MyMatrix<T> EXT_scaled_T = UniversalMatrixConversion<T, Tring>(EXT_scaled);
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
  return {std::move(EXT_T),  std::move(EXT_scaled_T), std::move(Cent),
          std::move(EXText), std::move(EXTvert),      std::move(eFace),
          std::move(WMat)};
}

// A permutation of the extended set restricted to the vertex block. The
// two blocks are separated by their diagonal value in the weight matrix,
// so a vertex goes to a vertex and the restriction is a permutation.
template <typename Telt>
Telt RestrictPermToVertexBlock(Telt const &eElt, int const &nbVert) {
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> eList(nbVert);
  for (int i = 0; i < nbVert; i++) {
    Tidx img = OnPoints(static_cast<Tidx>(i), eElt);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
    if (static_cast<int>(img) >= nbVert) {
      std::cerr << "PERIODIC_DELAUNAY: RestrictPermToVertexBlock: a vertex is "
                   "mapped outside of the vertex block\n";
      throw TerminalException{1};
    }
#endif
    eList[i] = img;
  }
  return Telt(eList);
}

/*
  The transformation realized by a permutation of the vertices of data1
  onto the vertices of data2, as an integral affine matrix.

  It is read off the HOMOGENEOUS vertex matrices, whose rows (1, v) are
  integral in this representation, rather than through the circumcentres:
  the affine matrix then comes out of one FindTransformation and its
  integrality is automatic, the extended set containing short vectors
  which generate the period lattice D Z^n -- preserving them forces the
  linear part to be integral, and the translation is the image of an
  integral vertex minus its linear image. The circumcentres are still
  needed to recentre the vertices for the weight matrix, which is what
  makes the isometries linear, but nothing here uses them.

  The permutation acts on the VERTEX block, so it is the homogeneous
  vertex matrices that are passed, not the extended ones.
 */
template <typename T, typename Telt, typename Tidx_value>
MyMatrix<T>
PeriodicTransformFromPerm(PeriodicDelaunayIsoData<T, Tidx_value> const &data1,
                          PeriodicDelaunayIsoData<T, Tidx_value> const &data2,
                          Telt const &eElt) {
  return FindTransformation<T, Telt>(data1.EXT_scaled_T, data2.EXT_scaled_T,
                                    eElt);
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
  auto f_trans = [&](Telt const &eElt) -> MyMatrix<T> {
    return PeriodicTransformFromPerm<T, Telt, Tidx_value>(data, data, eElt);
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
std::optional<MyMatrix<T>>
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
  // The equivalence is found on the extended sets, so it has the degree
  // of those; the transformation is read off the vertices, so it is
  // restricted to the vertex block first. The colouring of the two blocks
  // makes the restriction well defined, a vertex going to a vertex.
  Telt eRes_vert = RestrictPermToVertexBlock<Telt>(*eRes, data1.EXTvert.rows());
  MyMatrix<T> m0 =
      PeriodicTransformFromPerm<T, Telt, Tidx_value>(data1, data2, eRes_vert);
  if (IsPeriodicPointSetPreserved(pps, m0)) {
    return m0;
  }
  // The other isometries of the pair are the h * m0, so the point set is
  // preserved by one of them exactly when some h sends the cosets to their
  // preimage by m0.
  PeriodicClassAction<Tring> act(pps);
  Tgroup GRPisom1 = PeriodicDelaunayIsometryGroup<T, Tgroup, Tidx_value>(data1, os);
  auto f_trans1 = [&](Telt const &eElt) -> MyMatrix<T> {
    return PeriodicTransformFromPerm<T, Telt, Tidx_value>(data1, data1, eElt);
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
  MyMatrix<T> h = f_trans1(Telt(eList_vert));
  MyMatrix<T> ret = h * m0;
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (!IsPeriodicPointSetPreserved<Tring, T>(pps, ret)) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicDelaunay_TestEquivalence: the "
                 "equivalence found does not preserve the point set\n";
    throw TerminalException{1};
  }
#endif
  return ret;
}

/*
  The stabilizer and the equivalence of Gram matrices under the
  transformations that also preserve the periodic point set, which is what
  the enumeration of the periodic iso-Delaunay domains recognizes its
  domains with. Both are the lattice computations with one extra adequacy
  condition on the matrices: the linear part has to admit an affine
  extension preserving the point set. The condition is subgroup-defining
  (the point-set-preserving affine transformations form a group and their
  linear parts are its image), as the kernels require.

  The generators of LinSpa.PtStabGens are handed to the kernels as
  known-good elements, so they have to be adequate themselves: that is an
  invariant of the input T-space, checked here, and the programs reading a
  T-space together with a point set have to validate it unconditionally.
 */
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
template <typename T, typename Tring>
void CheckPtStabGensPreserveCosets(LinSpaceMatrix<T> const &LinSpa,
                                   PeriodicPointSet<Tring> const &pps) {
  for (auto &eGen : LinSpa.PtStabGens) {
    if (!PeriodicAffineExtension<Tring, T>(pps, eGen)) {
      std::cerr << "PERIODIC_DELAUNAY: a generator of PtStabGens admits no "
                   "affine extension preserving the point set. The point "
                   "stabilizer of the T-space has to stabilize the cosets "
                   "for the periodic computations to be correct. eGen=\n";
      WriteMatrix(std::cerr, eGen);
      throw TerminalException{1};
    }
  }
}
#endif

template <typename T, typename Tring, typename Tgroup>
Result_ComputeStabilizer_SHV<T, Tgroup> LINSPA_ComputeStabilizer_SHV_Periodic(
    LinSpaceMatrix<T> const &LinSpa, PeriodicPointSet<Tring> const &pps,
    MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
    std::optional<MyMatrix<T>> const &CommonGramMat, std::ostream &os) {
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  CheckPtStabGensPreserveCosets(LinSpa, pps);
#endif
  auto f_extra = [&](MyMatrix<T> const &x) -> bool {
    return PeriodicAffineExtension<Tring, T>(pps, x).has_value();
  };
  return LINSPA_ComputeStabilizer_SHV_Kernel<T, Tgroup, decltype(f_extra)>(
      LinSpa, eMat, SHV_T, CommonGramMat, f_extra, os);
}

template <typename T, typename Tring, typename Tgroup>
std::optional<MyMatrix<T>> LINSPA_TestEquivalenceGramMatrix_SHV_Periodic(
    LinSpaceMatrix<T> const &LinSpa, PeriodicPointSet<Tring> const &pps,
    MyMatrix<T> const &eMat1, MyMatrix<T> const &eMat2,
    MyMatrix<T> const &SHV1_T, MyMatrix<T> const &SHV2_T,
    std::optional<MyMatrix<T>> const &CommonGramMat, std::ostream &os) {
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  CheckPtStabGensPreserveCosets(LinSpa, pps);
#endif
  auto f_extra = [&](MyMatrix<T> const &x) -> bool {
    return PeriodicAffineExtension<Tring, T>(pps, x).has_value();
  };
  return LINSPA_TestEquivalenceGramMatrix_SHV_Kernel<T, Tgroup,
                                                     decltype(f_extra)>(
      LinSpa, eMat1, eMat2, SHV1_T, SHV2_T, CommonGramMat, f_extra, os);
}

// The vertices of a cell back in the unscaled rational coordinates, the
// leading column of 1 being kept.
template <typename T, typename Tring>
MyMatrix<T> PeriodicUnscaledVertices(MyMatrix<Tring> const &EXT_scaled,
                                     Tring const &N) {
  int nbVert = EXT_scaled.rows();
  int n = EXT_scaled.cols() - 1;
  T N_T = UniversalScalarConversion<T, Tring>(N);
  MyMatrix<T> EXT_T(nbVert, n + 1);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    EXT_T(iVert, 0) = T(1);
    for (int i = 0; i < n; i++) {
      EXT_T(iVert, i + 1) =
          UniversalScalarConversion<T, Tring>(EXT_scaled(iVert, i + 1)) / N_T;
    }
  }
  return EXT_T;
}

/*
  An invariant of a cell of a periodic point set, used to sort the orbits
  before the equivalence is tried. Beyond what the lattice invariant uses
  -- the number of vertices and the multiset of the distances between them
  -- the cosets the vertices belong to contribute through the multiset of
  their counts, which is invariant since a symmetry permutes the cosets.
 */
template <typename T, typename Tring>
size_t ComputeInvariantPeriodicDelaunay(MyMatrix<T> const &GramMat,
                                        PeriodicPointSet<Tring> const &pps,
                                        size_t const &seed,
                                        MyMatrix<Tring> const &EXT_scaled) {
  int nbVert = EXT_scaled.rows();
  int n = GramMat.rows();
  MyMatrix<T> EXT_T = PeriodicUnscaledVertices<T, Tring>(EXT_scaled, pps.N);
  std::map<T, size_t> map_dist;
  MyVector<T> eDiff(n);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int jVert = iVert + 1; jVert < nbVert; jVert++) {
      for (int i = 0; i < n; i++) {
        eDiff(i) = EXT_T(iVert, i + 1) - EXT_T(jVert, i + 1);
      }
      map_dist[EvaluationQuadForm<T, T>(GramMat, eDiff)] += 1;
    }
  }
  // How many vertices sit on each coset, as a multiset: a symmetry may
  // permute the cosets, so only the counts are invariant, not which coset
  // carries which count.
  std::map<size_t, size_t> map_coset;
  {
    std::vector<size_t> count(pps.cosets_num.rows(), 0);
    for (int iVert = 0; iVert < nbVert; iVert++) {
      MyVector<Tring> u(n);
      for (int i = 0; i < n; i++) {
        u(i) = EXT_scaled(iVert, i + 1);
      }
      std::optional<size_t> opt = GetCosetIndex(pps, u);
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
      if (!opt) {
        std::cerr << "PERIODIC_DELAUNAY: ComputeInvariantPeriodicDelaunay: a "
                     "vertex is not a point of the point set\n";
        throw TerminalException{1};
      }
#endif
      count[*opt] += 1;
    }
    for (auto &val : count) {
      map_coset[val] += 1;
    }
  }
  size_t hash = seed;
  auto combine_hash = [](size_t &seed_i, size_t new_hash) -> void {
    seed_i ^= new_hash + 0x9e3779b8 + (seed_i << 6) + (seed_i >> 2);
  };
  combine_hash(hash, static_cast<size_t>(nbVert));
  for (auto &kv : map_dist) {
    combine_hash(hash, std::hash<T>()(kv.first));
    combine_hash(hash, kv.second);
  }
  for (auto &kv : map_coset) {
    combine_hash(hash, kv.first);
    combine_hash(hash, kv.second);
  }
  return hash;
}

/*
  The data of a periodic Delaunay enumeration: the form, the point set,
  the short vectors used to rigidify the isometry computations, the
  solver, and the moves of the local improvement of the adjacency. The
  short vectors have to generate Z^n, and the Graver moves are those of
  the point set, hence scaled by N.
 */
template <typename T, typename Tring, typename Tgroup>
struct PeriodicDataDelaunay {
  MyMatrix<T> GramMat;
  PeriodicPointSet<Tring> pps;
  MyMatrix<T> SHV;
  PeriodicCVPSolver<T, Tring> solver;
  MyMatrix<Tring> ShvGraverBasisScaled;
  RecordDualDescOperation<T, Tgroup> rddo;
};

template <typename T, typename Tring, typename Tgroup>
PeriodicDataDelaunay<T, Tring, Tgroup> GetPeriodicDataDelaunay(
    MyMatrix<T> const &GramMat, PeriodicPointSet<Tring> const &pps,
    MyMatrix<T> const &SHV, MyMatrix<Tring> const &ShvGraverBasis,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, std::ostream &os) {
  return {GramMat,
          pps,
          SHV,
          PeriodicCVPSolver<T, Tring>(GramMat, pps, os),
          PeriodicScaledGraverBasis(ShvGraverBasis, pps.N),
          RecordDualDescOperation<T, Tgroup>(AllArr, os)};
}

// A cell of the tessellation, its vertices homogeneous and scaled.
template <typename Tring, typename Tgroup> struct PeriodicDelaunay_Obj {
  MyMatrix<Tring> EXT;
  Tgroup GRP;
};

template <typename Tring> struct PeriodicDelaunay_AdjI {
  Face eInc;
  MyMatrix<Tring> EXT;
};

template <typename Tring> struct PeriodicDelaunay_AdjO {
  Face eInc;
  MyMatrix<Tring> eBigMat;
};

/*
  The enumeration of the orbits of Delaunay cells of a periodic point set,
  in the form the shared adjacency scheme consumes. Every piece it is made
  of is the periodic counterpart of the lattice one and nothing of the
  scheme itself is periodic: the initial cell and the adjacent cells come
  from the geometry templated on the solver, the orbits of facets from the
  usual dual description, and the recognition of an already found cell
  from the periodic equivalence.
 */
template <typename T, typename Tring, typename Tgroup>
struct PeriodicDataDelaunayFunc {
  PeriodicDataDelaunay<T, Tring, Tgroup> &data;
  using Tobj = PeriodicDelaunay_Obj<Tring, Tgroup>;
  using TadjI = PeriodicDelaunay_AdjI<Tring>;
  using TadjO = PeriodicDelaunay_AdjO<Tring>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
    MyMatrix<Tring> EXT =
        FindDelaunayPolytope<T, Tring>(data.solver, data.rddo.os);
    return {std::move(EXT), {}};
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantPeriodicDelaunay<T, Tring>(data.GramMat, data.pps,
                                                      seed, x.EXT);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<T>> opt =
        PeriodicDelaunay_TestEquivalence<T, Tring, Tgroup>(
            data.GramMat, data.pps, data.SHV, y.EXT, x.EXT, data.rddo.os);
    if (!opt) {
      return {};
    }
    TadjO ret{y.eInc, UniversalMatrixConversion<Tring, T>(*opt)};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    Tobj x_ret{x.EXT, {}};
    int n = data.GramMat.rows();
    MyMatrix<Tring> eBigMat = IdentityMat<Tring>(n + 1);
    TadjO ret{x.eInc, std::move(eBigMat)};
    return {std::move(x_ret), std::move(ret)};
  }
  std::optional<std::vector<TadjI>> f_adj(Tobj &x) {
    std::ostream &os = data.rddo.os;
    x.GRP = PeriodicDelaunay_Stabilizer<T, Tring, Tgroup>(
        data.GramMat, data.pps, data.SHV, x.EXT, os);
    MyMatrix<T> EXT_T = PeriodicUnscaledVertices<T, Tring>(x.EXT, data.pps.N);
    vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, x.GRP, data.rddo);
    // The facets are described on the scaled vertices, which the geometry
    // of the adjacency works with.
    MyMatrix<T> EXT_scaled_T = UniversalMatrixConversion<T, Tring>(x.EXT);
    SubsetRankOneSolver<T> ext_solver(EXT_scaled_T);
    std::vector<TadjI> ListAdj;
    for (auto &eOrbB : TheOutput) {
      MyMatrix<Tring> EXTadj = FindAdjacentDelaunayPolytope<T, Tring>(
          data.solver, data.ShvGraverBasisScaled, EXT_scaled_T, ext_solver,
          eOrbB, os);
      ListAdj.push_back({eOrbB, std::move(EXTadj)});
    }
    return ListAdj;
  }
  Tobj f_adji_obj(TadjI const &x) { return {x.EXT, {}}; }
  size_t f_complexity(Tobj const &x) { return x.EXT.rows(); }
};

/*
  The tessellation, in the form the iso-Delaunay machinery consumes, of
  the orbits the enumeration produced.

  Everything is in the scaled frame: the vertices are carried over as they
  are and the transformations are read through their acting matrices, so
  the products below are between matching frames. The inequalities come
  out multiplied by N^2 with respect to the unscaled ones, and
  ComputeDelaunayAdjIneq canonicalizes them anyway, so the defining system
  of the iso-Delaunay domain is the same.

  The form has to be GENERIC for the point set, that is its Delaunay cells
  have to be simplices: a cell with more vertices is co-spherical only on
  a wall of the iso-Delaunay domain, where the inequalities of the
  tessellation are not defined. BuildVoronoiIneqPreComputeChecked says so
  when the sanity checks are on.
 */
template <typename T, typename Tring, typename Tgroup>
DelaunayTesselationIneq<T, Tgroup>
BuildPeriodicDelaunayTesselationIneq(
    std::vector<DatabaseEntry_Serial<PeriodicDelaunay_Obj<Tring, Tgroup>,
                                     PeriodicDelaunay_AdjO<Tring>>> const
        &l_ent,
    [[maybe_unused]] PeriodicPointSet<Tring> const &pps,
    std::vector<std::vector<T>> const &ListGram, std::ostream &os) {
  int n_del = l_ent.size();
  std::vector<Delaunay_EntryIneq<T, Tgroup>> l_dels(n_del);
  for (int i_del = 0; i_del < n_del; i_del++) {
    MyMatrix<Tring> const &EXT = l_ent[i_del].x.EXT;
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tring>(EXT);
    VoronoiInequalityPreComput<T> vipc =
        BuildVoronoiIneqPreComputeChecked<T>(EXT_T, ListGram, os);
    ContainerMatrix<T> cont(EXT_T);
    std::vector<Delaunay_AdjIneqO<T>> ListAdj;
    for (auto &eAdj : l_ent[i_del].ListAdj) {
      MyMatrix<T> EXT_target_T =
          UniversalMatrixConversion<T, Tring>(l_ent[eAdj.iOrb].x.EXT);
      MyMatrix<T> M = UniversalMatrixConversion<T, Tring>(eAdj.x.eBigMat);
      MyVector<T> eIneq =
          ComputeDelaunayAdjIneq(vipc, cont, EXT_target_T, M, ListGram, os);
      ListAdj.push_back({eAdj.x.eInc, eAdj.x.eBigMat, eAdj.iOrb,
                         std::move(eIneq)});
    }
    l_dels[i_del] = {EXT, l_ent[i_del].x.GRP, std::move(ListAdj)};
  }
  return {std::move(l_dels)};
}

// clang-format off
#endif  // SRC_DELAUNAY_PERIODICDELAUNAY_H_
// clang-format on
