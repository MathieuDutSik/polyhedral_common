// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_FUNCTIONS_H_
#define SRC_LATT_TSPACE_FUNCTIONS_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "Shvec_exact.h"
#include "Positivity.h"
#include "SHORT_Realizability.h"
#include "boost_serialization.h"
#include <set>
#include <vector>
#include <unordered_map>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_TSPACE_FUNCTIONS
#endif

#ifdef DISABLE_DEBUG_TSPACE_FUNCTIONS
#undef DEBUG_TSPACE_FUNCTIONS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_TSPACE_FUNCTIONS
#endif

/*
  By a T-space, we mean a vector space which intersect the cone of positive
  definite matrices. We are interested here only on those positive definite
  elements.
  ---
  The motivation for this kind of construction is for finding good packing,
  covering, etc. The idea is that the smaller dimensionality allows for a
  smaller number of combinatorial possibilities.
  ---
  However, if the total number of combinatorial possibilities is expected to be
  smaller, the combinatorial complexity on programming and thinking is
  definitely higher. The thing that are considered in this section:
  * Finding the integral saturation (useful for integral separation)
  * The generation of the T-spaces (see other include)
  * The computation of stabilizer and equivalence of positive definite matrices
  in the T-space.
 */

/*
  Aspects of equivalence issues.
  * The point stabilizer of the T-space X, that is:
  PtStab(X) = {g in GL_n(Z) s.t. g M g^T = M for all M in X}
  * The global stabilizer of the T-space X, that is
  GlStab(X) = {g in GL_n(Z) s.t. g M g^T \in X for all M in X}
  -----
  Facts:
  * The group PtStab(X) is finite.
  * The group GlStab(X) is finitely generated.
  * If X = {M s.t. g M g^T = M for all g in G0} for G0 a finite
  subgroup of GL_n(Z) then GlStab(X) is the normalizer of G0.
  * If a saturating integral basis is chosen for X then the
  expression of GlStab(X) is integral in that basis.
  * The group PtStab(X) is easy to compute: Take a basis of the
  T-space, a supermat, an invariant set for the supermat and then
  the combined scalar product for the basis.
  * The group GlStab(X) is much harder to compute. A generating
  set is usually a byproduct of the polyhedral tessellation like
  coming from the perfect form or the Delaunay polytopes. This is
  the Opgenorth algorithm.
  ----
  Usually, we want to compute with the full group GlStab(X).
  However, for some application for example to number theory we
  want to compute with the group GL_n(Z[t]) with t a generating
  element. If the degree of t is d then expressed in a good basis
  the matrix is in GL_{nd}(Z) and the matrix of GL_{n}(Z[t]) are
  characterized as the ones commuting with the multiplication by t.
  This is especially convenient since we can nicely characterize
  commuting operation [1, Section 3.3].
  Another example is if a group G stabilize a subspace. Then the
  normalizer of G will also stabilize the subspace.
  The stabilization of a subspace and the fact that there is a
  canonically associated scalar product allows to define an orthognal
  and so a projector to commute with.
  ----
  The big problem that we have is given A in X, to compute the
  stabilizer in GlStab(X) and given two matrices A, B in X, check
  if there is an element g in GlStab(X) such that g A g^T = B.
  Both problems are clearly finite: We can compute the stabilizer
  in GL_n(Z) of a matrix A and then keep the elements that
  stabilize X.
  Can we do better?
  * We can use cosets with the point stabilizer, but it is never that
  large and if the index is too high, that can sink us.
  * The preserved subspace trick should help us in some cases.
  * For case like number ring with us not looking at GL_n(Z[t]) what
  can we actually do? Can we describe some normalizer or something?
  For gaussian or eisenstein, we have a conjugation operation. Could we
  generalize that?
  * We could map the problem into the matrix space, but the transformation
  would be matricial in S^n. How could we characterize the ones coming as
  X -> P X P^T ? See the linear preserver problem literature about that.

  References:
  [1]: David Bremner, Mathieu Dutour Sikiric, Dmitrii V. Pasechnik,
  Thomas Rehn and Achill Schuermann, Computing symmetry groups of polyhedra
  LMS J. Comput. Math. 17 (1) (2014) 565–581, doi:10.1112/S1461157014000400
  -----
  The question of Bravais space had been forgotten. That could get us some
  improvement and so we should put it.
  The idea is that if a T-space is Bravais, then we can take the equivalence
  of the inner form and their stabilizer. So, in that case there is no need
  to check for the stabilization of the T-space. But we need to find the
  reference.
  -----
  The question of characteristic vector set ought to be considered as we
  have this in the GAP code.
 */


/*
  The matrix of pairwise scalar products.
  It is useful for the self-duality, but not only.
 */
template<typename T>
struct PairwiseScalarInfo {
  MyMatrix<T> PairwiseScalarInv;
};

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, PairwiseScalarInfo<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("PairwiseScalarInv", val.PairwiseScalarInv);
}
}  // namespace boost::serialization

/*
  The structure of a T-space (see the papers by Mathieu Dutour Sikiric,
  Frank Vallentin and Achill Schuermann).

  We have the following:
  * A SuperMat is a positive definite matrix in the T-space.
  * The ListMat is a basis of matrices of the T-space.
  * The ListLineMat is the line matrix for effective computation.
  * The ListComm is a set of matrices which are not necessarily of
    determinant 1. The group of symmetries has to commutes with the
    elements in ListComm.
 */
template <typename T> struct LinSpaceMatrix {
  int n;
  // Whether it is a Bravais space or not
  bool isBravais;
  // The positive definite matrix.
  MyMatrix<T> SuperMat;
  // The basis of the T-space
  std::vector<MyMatrix<T>> ListMat;
  // The basis but expressed as line matrices
  std::vector<std::vector<T>> ListLineMat;
  // The basis expressed as a big matrix, useful for membership questions
  MyMatrix<T> ListMatAsBigMat;
  // The list of matrices with which the global stabilizing elements must
  // commute
  std::vector<MyMatrix<T>> ListComm;
  // The list of preserved subspaces by the element of the global stabilizer
  std::vector<MyMatrix<T>> ListSubspaces;
  // The point stabilizer of the T-space
  std::vector<MyMatrix<T>> PtStabGens;
  // The list of spanning elements. If empty, does not apply. This is used
  // for the T-spaces coming from
  std::vector<MyMatrix<T>> l_spanning_elements;
  // The list of outer elements used for stabilizer / isomorphism. For
  // example for quadratic spaces, the conjugation might be put.
  std::vector<MyMatrix<T>> l_outer_elements;
  // Pairwise scalar info.
  PairwiseScalarInfo<T> pairwise_scalar_info;
  // Whether the T-space is known to be self-dual
  bool is_self_dual;
};

// Equality check for LinSpaceMatrix
template <typename T>
bool LinSpaceMatrixEqual(LinSpaceMatrix<T> const &a, LinSpaceMatrix<T> const &b) {
  if (a.n != b.n)
    return false;
  if (a.isBravais != b.isBravais)
    return false;
  if (a.is_self_dual != b.is_self_dual)
    return false;
  if (a.SuperMat != b.SuperMat)
    return false;
  if (a.ListMat != b.ListMat)
    return false;
  if (a.ListLineMat != b.ListLineMat)
    return false;
  if (a.ListMatAsBigMat != b.ListMatAsBigMat)
    return false;
  if (a.ListComm != b.ListComm)
    return false;
  if (a.ListSubspaces != b.ListSubspaces)
    return false;
  if (a.PtStabGens != b.PtStabGens)
    return false;
  if (a.l_spanning_elements != b.l_spanning_elements)
    return false;
  if (a.l_outer_elements != b.l_outer_elements)
    return false;
  if (a.pairwise_scalar_info.PairwiseScalarInv !=
      b.pairwise_scalar_info.PairwiseScalarInv)
    return false;
  return true;
}

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, LinSpaceMatrix<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("n", val.n);
  ar &make_nvp("isBravais", val.isBravais);
  ar &make_nvp("SuperMat", val.SuperMat);
  ar &make_nvp("ListMat", val.ListMat);
  ar &make_nvp("ListLineMat", val.ListLineMat);
  ar &make_nvp("ListMatAsBigMat", val.ListMatAsBigMat);
  ar &make_nvp("ListComm", val.ListComm);
  ar &make_nvp("ListSubspaces", val.ListSubspaces);
  ar &make_nvp("PtStabGens", val.PtStabGens);
  ar &make_nvp("l_spanning_elements", val.l_spanning_elements);
  ar &make_nvp("l_outer_elements", val.l_outer_elements);
  ar &make_nvp("pairwise_scalar_info", val.pairwise_scalar_info);
  ar &make_nvp("is_self_dual", val.is_self_dual);
}
}  // namespace boost::serialization

template <typename T>
MyMatrix<T> GetListMatAsBigMat(std::vector<MyMatrix<T>> const &ListMat) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "TSPACE: We have n_mat = 0\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  int sym_dim = (n * (n + 1)) / 2;
  MyMatrix<T> BigMat(n_mat, sym_dim);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyVector<T> V = SymmetricMatrixToVector(ListMat[i_mat]);
    AssignMatrixRow(BigMat, i_mat, V);
  }
  return BigMat;
}


/*
  This information is used for self-duality but also for
  representing evaluation functions in a matrix space.
 */
template <typename T>
PairwiseScalarInfo<T> get_pairwise_scalar_info(std::vector<MyMatrix<T>> const& ListMat, MyMatrix<T> const& SuperMat) {
  int n_mat = ListMat.size();
  MyMatrix<T> PairwiseScalar(n_mat, n_mat);
  MyMatrix<T> InvSuperMat = Inverse(SuperMat);
  std::vector<MyMatrix<T>> l_mat;
  for (auto & M1: ListMat) {
    MyMatrix<T> M2 = M1 * InvSuperMat;
    l_mat.push_back(M2);
  }
  for (int i=0; i<n_mat; i++) {
    for (int j=i; j<n_mat; j++) {
      MyMatrix<T> prod = l_mat[i] * l_mat[j];
      T scal(0);
      for (int i=0; i<prod.rows(); i++) {
        scal += prod(i,i);
      }
      PairwiseScalar(i, j) = scal;
      PairwiseScalar(j, i) = scal;
    }
  }
  MyMatrix<T> PairwiseScalarInv = Inverse(PairwiseScalar);
  PairwiseScalarInfo<T> psi{PairwiseScalarInv};
  return psi;
}



template <typename T>
LinSpaceMatrix<T> BuildLinSpace(MyMatrix<T> const &SuperMat,
                                std::vector<MyMatrix<T>> const &ListMat,
                                std::vector<MyMatrix<T>> const &ListComm) {
  int n = SuperMat.rows();
  std::vector<std::vector<T>> ListLineMat;
  for (auto &eMat : ListMat) {
    std::vector<T> eV = GetLineVector(eMat);
    ListLineMat.push_back(eV);
  }
  MyMatrix<T> BigMat = GetListMatAsBigMat(ListMat);
  std::vector<MyMatrix<T>> ListSubspaces;
  MyMatrix<T> eGen = -IdentityMat<T>(n);
  std::vector<MyMatrix<T>> PtStab{eGen};
  // It may be actually a Bravais space, but setting up to false avoids
  // potential problems.
  bool isBravais = false;
  std::vector<MyMatrix<T>> l_spanning_elements;
  std::vector<MyMatrix<T>> l_outer_elements;
  PairwiseScalarInfo<T> pairwise_scalar_info = get_pairwise_scalar_info(ListMat, SuperMat);
  bool is_self_dual = false;
  return {n,      isBravais, SuperMat,      ListMat, ListLineMat,
          BigMat, ListComm,  ListSubspaces, PtStab,
          l_spanning_elements, l_outer_elements,
          pairwise_scalar_info, is_self_dual};
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ComputePointStabilizerTspace(MyMatrix<T> const &SuperMat,
                             std::vector<MyMatrix<T>> const &ListMat,
                             std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tfield = T;
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(SuperMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  std::vector<T> Vdiag(SHV_T.rows(), 0);
  std::vector<std::vector<Tidx>> ListGenPerm =
      GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tgroup>(SHV_T, ListMat,
                                                              Vdiag, os);
  std::vector<MyMatrix<Tint>> ListGenMatr;
  for (auto &eList : ListGenPerm) {
    MyMatrix<T> eMatr_T = FindTransformation_vect(SHV_T, SHV_T, eList);
    MyMatrix<Tint> eMatr = UniversalMatrixConversion<Tint, T>(eMatr_T);
    ListGenMatr.push_back(eMatr);
  }
  return ListGenMatr;
}


/*
  The Frobenius scalar product is defined as Tr(AB^T)
 */
template <typename T>
T frobenius_inner(MyMatrix<T> const& M1, MyMatrix<T> const& M2) {
  int n_row = M1.rows();
  int n_col = M1.cols();
  T sum(0);
  for (int i_row=0; i_row<n_row; i_row++) {
    for (int i_col=0; i_col<n_col; i_col++) {
      sum += M1(i_row,i_col) * M2(i_row,i_col);
    }
  }
  return sum;
}

template <typename T>
MyMatrix<T> GetRandomPositiveDefinite(LinSpaceMatrix<T> const &LinSpa,
                                      int const &N, std::ostream &os) {
  int n = LinSpa.n;
  MyMatrix<T> TheMat = ZeroMatrix<T>(n, n);
  for (auto &eMat : LinSpa.ListMat) {
    int coef = random() % (2 * N + 1) - N;
    TheMat += coef * eMat;
  }
  while (true) {
    if (IsPositiveDefinite(TheMat, os)) {
      return TheMat;
    }
    TheMat += LinSpa.SuperMat;
  }
}

template <typename T, typename Tint, typename Tgroup>
bool IsSymmetryGroupCorrect(MyMatrix<T> const &GramMat,
                            LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row, 0);
  std::vector<MyMatrix<T>> ListMat = {GramMat};
  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tgroup>(SHV_T, ListMat,
                                                              Vdiag, os);
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
        FindTransformationGeneral_vect(SHV_T, SHV_T, eList);
    if (!opt) {
      std::cerr << "TSPACE: Failed to find the matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &M_T = *opt;
    if (!IsIntegralMatrix(M_T)) {
      std::cerr << "TSPACE: Bug: The matrix should be integral\n";
      throw TerminalException{1};
    }
    for (auto &eMat : LinSpa.ListMat) {
      MyMatrix<T> eMatImg = M_T * eMat * M_T.transpose();
      if (eMatImg != eMat) {
        return false;
      }
    }
  }
  return true;
}

template <typename T, typename Tint>
bool IsBravaisSpace(int n, std::vector<MyMatrix<T>> const &ListMat,
                    std::vector<MyMatrix<Tint>> const &ListGen,
                    [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrix<T>> ListGen_T;
  for (auto &eGen : ListGen) {
    MyMatrix<T> eGen_T = UniversalMatrixConversion<T, Tint>(eGen);
    ListGen_T.push_back(eGen_T);
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: IsBravaisSpace: |ListGen|=" << ListGen.size() << "\n";
  os << "TSPACE: IsBravaisSpace: |ListMat|=" << ListMat.size() << "\n";
  for (auto &eGen : ListGen) {
    for (auto &eMat : ListMat) {
      MyMatrix<T> eMatImg = eGen * eMat * eGen.transpose();
      if (eMatImg != eMat) {
        std::cerr << "TSPACE: IsBravaisSpace: |eMat|=" << eMat.rows() << " / "
                  << eMat.cols() << "\n";
        std::cerr << "TSPACE: IsBravaisSpace: |eGen|=" << eGen.rows() << " / "
                  << eGen.cols() << "\n";
        std::cerr << "TSPACE: IsBravaisSpace: eMat should equal to eMatImg\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  std::vector<MyMatrix<T>> BasisInv = BasisInvariantForm(n, ListGen_T, os);
  MyMatrix<T> Big_ListMat = GetListMatAsBigMat(ListMat);
  MyMatrix<T> Big_BasisInv = GetListMatAsBigMat(BasisInv);
  //
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: IsBravaisSpace: |Big_ListMat|=" << Big_ListMat.rows() << " / "
     << Big_ListMat.cols() << "\n";
  os << "TSPACE: IsBravaisSpace: |Big_BasisInv|=" << Big_BasisInv.rows()
     << " / " << Big_BasisInv.cols() << "\n";
  if (!IsSubspaceContained(Big_ListMat, Big_BasisInv)) {
    std::cerr << "TSPACE: The elements of ListMat are not in the invariant "
                 "space which "
                 "is not acceptable\n";
    throw TerminalException{1};
  }
#endif
  //
  return IsSubspaceContained(Big_BasisInv, Big_ListMat);
}

/*
  We want to find symmetries of the polytope that are pointwise stabilizer of
  the T-space. Two behaviors that we want to avoid:
  * Not globally stabilizing the T-space.
  * Having an action on the T-space that is non-trivial.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T>
GetRandomPositiveDefiniteNoNontrivialSymm(LinSpaceMatrix<T> const &LinSpa,
                                          int const &N, std::ostream &os) {
  int N_work = N;
  while (true) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: GetRandomPositiveDefiniteNoNontrivialSymm: Before "
          "GetRandomPositiveDefinite\n";
#endif
    MyMatrix<T> TheMat = GetRandomPositiveDefinite(LinSpa, N_work, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: GetRandomPositiveDefiniteNoNontrivialSymm: After "
          "GetRandomPositiveDefinite\n";
    os << "TSPACE: GetRandomPositiveDefiniteNoNontrivialSymm: TestMat=\n";
    WriteMatrix(os, TheMat);
#endif
    bool test = IsSymmetryGroupCorrect<T, Tint, Tgroup>(TheMat, LinSpa, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: GetRandomPositiveDefiniteNoNontrivialSymm: test=" << test
       << "\n";
#endif
    if (test) {
#ifdef DEBUG_TSPACE_FUNCTIONS
      os << "TSPACE: GetRandomPositiveDefiniteNoNontrivialSymm: return with "
            "N_work="
         << N_work << "\n";
#endif
      return TheMat;
    }
    N_work += 1;
  }
}

template <typename T>
MyMatrix<T> LINSPA_GetMatrixInTspace(LinSpaceMatrix<T> const &LinSpa,
                                     MyVector<T> const &eVect) {
  return GetMatrixFromBasis(LinSpa.ListMat, eVect);
}

template <typename T>
MyVector<T> LINSPA_GetVectorOfMatrixExpression(LinSpaceMatrix<T> const &LinSpa,
                                               MyMatrix<T> const &eMat) {
  MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
  int dimSymm = eMatVect.size();
  int dimLinSpa = LinSpa.ListMat.size();
  MyMatrix<T> TotalMatrix(dimLinSpa, dimSymm);
  for (int iLinSpa = 0; iLinSpa < dimLinSpa; iLinSpa++) {
    MyVector<T> V = SymmetricMatrixToVector(LinSpa.ListMat[iLinSpa]);
    AssignMatrixRow(TotalMatrix, iLinSpa, V);
  }
  std::optional<MyVector<T>> RecSol = SolutionMat(TotalMatrix, eMatVect);
  MyVector<T> V = unfold_opt(RecSol, "Failure in SolutionMat");
  return V;
}

template <typename T>
MyMatrix<T> matrix_in_t_space(MyMatrix<T> const &Pmat,
                              LinSpaceMatrix<T> const &LinSpa) {
  int dimSpace = LinSpa.ListMat.size();
  MyMatrix<T> MatSpace(dimSpace, dimSpace);
  int pos = 0;
  for (auto &eMatSp : LinSpa.ListMat) {
    MyMatrix<T> eMatSpImg = Pmat * eMatSp * Pmat.transpose();
    MyVector<T> eMatSpImg_V = SymmetricMatrixToVector(eMatSpImg);
    std::optional<MyVector<T>> opt =
        SolutionMat(LinSpa.ListMatAsBigMat, eMatSpImg_V);
    MyVector<T> eV = unfold_opt(opt, "matrix row");
    AssignMatrixRow(MatSpace, pos, eV);
    pos += 1;
  }
  return MatSpace;
}

template <typename T> void reset_paperwork(LinSpaceMatrix<T> &LinSpa) {
  LinSpa.ListLineMat.clear();
  for (auto &eMat : LinSpa.ListMat) {
    std::vector<T> eV = GetLineVector(eMat);
    LinSpa.ListLineMat.push_back(eV);
  }
  int n_mat = LinSpa.ListMat.size();
  if (n_mat > 0) {
    LinSpa.ListMatAsBigMat = GetListMatAsBigMat(LinSpa.ListMat);
    LinSpa.n = LinSpa.ListMat[0].rows();
  } else {
    std::cerr << "TSPACE: We have 0 matrices for ListMat\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tint, typename Tgroup>
void reset_pt_stab_gens(LinSpaceMatrix<T> &LinSpa, std::ostream &os) {
  std::vector<MyMatrix<Tint>> ListGens =
      ComputePointStabilizerTspace<T, Tint, Tgroup>(LinSpa.SuperMat,
                                                    LinSpa.ListMat, os);
  LinSpa.PtStabGens =
    UniversalStdVectorMatrixConversion<T,Tint>(ListGens);
}

template <typename T, typename Tint, typename Tgroup>
LinSpaceMatrix<T> BuildLinSpaceMatrix(std::vector<MyMatrix<T>> const &ListMat,
                                      std::ostream &os) {
  LinSpaceMatrix<T> LinSpa;
  LinSpa.ListMat = ListMat;
  reset_paperwork(LinSpa);
  std::optional<MyMatrix<T>> opt =
    GetOnePositiveDefiniteMatrix<T, Tint>(LinSpa.ListMat, os);
  if (opt) {
    LinSpa.SuperMat = *opt;
  } else {
    std::cerr << "TSPACE: Failed to find a positive definite matrix\n";
    throw TerminalException{1};
  }
  // ListComm not set as it cannot be guessed.
  // Same for ListSubspaces.
  reset_pt_stab_gens<T, Tint, Tgroup>(LinSpa, os);
  LinSpa.isBravais =
      IsBravaisSpace(LinSpa.n, LinSpa.ListMat, LinSpa.PtStabGens, os);
  return LinSpa;
}

template<typename T, typename Tint>
MyMatrix<Tint> vector_family_saturation(MyMatrix<Tint> const& SHV, std::vector<MyMatrix<T>> const& l_gen) {
  std::vector<MyVector<Tint>> l_vect;
  std::unordered_set<MyVector<Tint>> set_vect;
  auto f_insert=[&](MyVector<Tint> const& v) -> void {
    if (!set_vect.contains(v)) {
      set_vect.insert(v);
      l_vect.push_back(v);
    }
  };
  for (int u=0; u<SHV.rows(); u++) {
    MyVector<Tint> V = GetMatrixRow(SHV, u);
    f_insert(V);
  }
  size_t start = 0;
  while(true) {
    size_t len = l_vect.size();
    for (size_t u=start; u<len; u++) {
      MyVector<Tint> v1 = l_vect[u];
      for (auto & eGen_T: l_gen) {
        MyMatrix<Tint> eGen = UniversalMatrixConversion<Tint,T>(eGen_T);
        MyVector<Tint> v2 = eGen.transpose() * v1;
        f_insert(v2);
      }
    }
    start = len;
    if (start == l_vect.size()) {
      break;
    }
  }
  return MatrixFromVectorFamily(l_vect);
}

// clang-format off
#endif  // SRC_LATT_TSPACE_FUNCTIONS_H_
// clang-format on
