// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_FUNCTIONS_H_
#define SRC_LATT_TSPACE_FUNCTIONS_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include "Positivity.h"
#include "Tspace_ListMatTransform.h"
#include "Tspace_ListMatSaturation.h"
#include "SHORT_Realizability.h"
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

template<typename T>
struct SelfDualInfo {
  MyMatrix<T> PairwiseScalar;
  MyMatrix<T> PairwiseScalarInv;
};

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
  // Whether the T-space is known to be self-dual
  std::optional<SelfDualInfo<T>> self_dual_info;
};

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
  std::optional<SelfDualInfo<T>> self_dual_info;
  return {n,      isBravais, SuperMat,      ListMat, ListLineMat,
          BigMat, ListComm,  ListSubspaces, PtStab,
          l_spanning_elements, self_dual_info};
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
  Find one positive definite matrix in the space assuming that one exists.
  If one of the matrix of the basis is positive definite then the first one
  is provided.
  ---
  Maybe we should use semidefinite programming? Is that absolutely
  needed? Or just slower?
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveDefiniteMatrix_ListV(std::vector<MyMatrix<T>> const &ListMat,
                                   std::vector<MyVector<Tint>> const& ListV_init,
                                   std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPFCT: GetOnePositiveDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> const &eMat = ListMat[i_mat];
    if (IsPositiveDefinite(eMat, os)) {
      return eMat;
    }
  }
  std::vector<MyVector<Tint>> ListV = ListV_init;
#ifdef DEBUG_TSPACE_FUNCTIONS
  int iter=0;
#endif
  while (true) {
    int n_vect = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
    iter += 1;
    os << "TSPFCT: GetOnePositiveDefiniteMatrix iter=" << iter << " n_vect=" << n_vect << "\n";
#endif
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_vect, 1 + n_mat);
    MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_mat);
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      ListIneq(i_vect, 0) = -1;
      MyVector<Tint> V = ListV[i_vect];
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      for (int i_mat = 0; i_mat < n_mat; i_mat++) {
        T val = EvaluationQuadForm(ListMat[i_mat], V_T);
        ListIneq(i_vect, 1 + i_mat) = val;
        ToBeMinimized(1 + i_mat) += val;
      }
    }
    //
    // Solving the linear program
    //
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.PrimalDefined || !eSol.DualDefined) {
      return {};
    }
    MyMatrix<T> TrySuperMat = ZeroMatrix<T>(n, n);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      TrySuperMat += eSol.DirectSolution(i_mat) * ListMat[i_mat];
    }
    if (IsPositiveDefinite(TrySuperMat, os)) {
      return TrySuperMat;
    }
    //
    // Failed, trying to find a vector
    //
    int rnk = RankMat(TrySuperMat);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPFCT: GetOnePositiveDefiniteMatrix get_one_vect, rnk=" << rnk << " n=" << n << "\n";
    os << "TSPFCT: TrySuperMat=\n";
    WriteMatrix(os, TrySuperMat);
#endif
    if (rnk < n) {
      std::vector<MyVector<Tint>> list_v = GetShortVectorDegenerate<T, Tint>(TrySuperMat, os);
      for (auto & eV: list_v) {
        ListV.push_back(eV);
      }
    } else {
      T CritNorm(0);
      bool StrictIneq = false;
      MyVector<Tint> V = GetShortIntegralVector<T, Tint>(TrySuperMat, CritNorm, StrictIneq, os);
      ListV.push_back(V);
    }
  }
}

template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveDefiniteMatrix(std::vector<MyMatrix<T>> const &ListMat,
                             std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  std::vector<MyVector<Tint>> ListV_init = get_initial_vector_test_v<Tint>(n, {}, os);
  return GetOnePositiveDefiniteMatrix_ListV(ListMat, ListV_init, os);
}

/*
  Test group finiteness.
  An integral matrix group G is finite if and only if it preserves
  a positive definite quadratic form.
  Proof:
  * If G is finite then take the identity and the average under G. The average
    is positive definite and is preserved.
  * If G preserves a positive definite quadratic form Q, then G is a subgroup
    of the stabilizer of Q. The stabilizer of Q is finite so that follows.

  If a matrix group is finite then it preserves a positive definite quadratic form
  ---
  We only need to determine finiteness, so no need to care too much about issues
  of transpose and similar.
 */
template<typename T, typename Tint>
bool test_finiteness_group(std::vector<MyMatrix<Tint>> const& ListGens, std::ostream &os) {
  if (ListGens.size() == 0) {
    // No generators, so the group is trivial.
    return true;
  }
  int n = ListGens[0].rows();
  std::vector<MyMatrix<T>> ListGens_T = UniversalStdVectorMatrixConversion<T,Tint>(ListGens);
  std::vector<MyMatrix<T>> ListMat = BasisInvariantForm<T>(n, ListGens_T, os);
  if (ListMat.size() == 0) {
    // No invariant matrix
    return false;
  }
  std::optional<MyMatrix<T>> opt = GetOnePositiveDefiniteMatrix<T,Tint>(ListMat, os);
  if (opt) {
    return true;
  } else {
    return false;
  }
}

/*
  Tries to find a semi-definite matrix.
  Unfortunately, we tend to go into infinite loops.
  And that seems to be an intrinsic problem since
  S_{n,>0} is not a polyhedral cone.
  ---
  But it works sometimes.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveSemiDefiniteMatrix_ListV(std::vector<MyMatrix<T>> const &ListMat,
                                       std::vector<MyVector<Tint>> const& ListV_init,
                                       std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> const &eMat = ListMat[i_mat];
    if (IsPositiveSemiDefinite(eMat, os)) {
      // A useful ansatz
      return eMat;
    }
  }
  std::vector<MyVector<Tint>> ListV = ListV_init;
  int n_vect_init = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
  int iter=0;
#endif
  while (true) {
    int n_vect = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
    iter += 1;
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix iter=" << iter << " n_vect=" << n_vect << "\n";
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_vect + 1, 1 + n_mat);
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      MyVector<Tint> const& V = ListV[i_vect];
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      for (int i_mat = 0; i_mat < n_mat; i_mat++) {
        T val = EvaluationQuadForm(ListMat[i_mat], V_T);
        ListIneq(i_vect, 1 + i_mat) = val;
      }
    }
    ListIneq(n_vect, 0) = -1;
    MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_mat);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      T sum(0);
      for (int i_vect=0; i_vect<n_vect_init; i_vect++) {
        sum += ListIneq(i_vect, 1 + i_mat);
      }
      ListIneq(n_vect, 1 + i_mat) = sum;
      ToBeMinimized(1 + i_mat) = sum;
    }
    //
    // Solving the linear program
    //
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.PrimalDefined || !eSol.DualDefined) {
      return {};
    }
    MyMatrix<T> TrySuperMat = ZeroMatrix<T>(n, n);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      TrySuperMat += eSol.DirectSolution(i_mat) * ListMat[i_mat];
    }
    if (IsPositiveSemiDefinite(TrySuperMat, os)) {
      return TrySuperMat;
    }
    //
    // Failed, trying to find a vector
    //
    MyMatrix<T> NSP_T = NullspaceIntMat(TrySuperMat);
    MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
    MyMatrix<Tint> Compl = SubspaceCompletionInt<Tint>(NSP, n);
    MyMatrix<T> Compl_T = UniversalMatrixConversion<T,Tint>(Compl);
    MyMatrix<T> M = Compl_T * TrySuperMat * Compl_T.transpose();
    T CritNorm(0);
    bool StrictIneq = true;
    MyVector<Tint> V1 = GetShortIntegralVector<T, Tint>(M, CritNorm, StrictIneq, os);
    MyVector<Tint> V2 = Compl.transpose() * V1;
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix V2=" << StringVectorGAP(V2) << "\n";
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix TrySuperMat=\n";
    WriteMatrixGAP(os, TrySuperMat);
    os << "\n";
#endif
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
    T norm = EvaluationQuadForm(TrySuperMat, V2);
    if (norm >= 0) {
      std::cerr << "TSPFCT: The norm should be negative since it should be a counter example\n";
      throw TerminalException{1};
    }
#endif
    ListV.push_back(V2);
  }
}

template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveSemiDefiniteMatrix(std::vector<MyMatrix<T>> const &ListMat,
                                 std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  std::vector<MyVector<Tint>> ListV_init = get_initial_vector_test_v<Tint>(n, {}, os);
  return GetOnePositiveSemiDefiniteMatrix_ListV(ListMat, ListV_init, os);
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
void set_self_dual_info(LinSpaceMatrix<T> & LinSpa) {
  int n = LinSpa.n;
  int n_mat = LinSpa.ListMat.size();
  MyMatrix<T> PairwiseScalar(n_mat, n_mat);
  MyMatrix<T> InvSuperMat = Inverse(LinSpa.SuperMat);
  std::vector<MyMatrix<T>> l_mat;
  for (auto & M1: LinSpa.ListMat) {
    MyMatrix<T> M2 = M1 * InvSuperMat;
    l_mat.push_back(M2);
  }
  for (int i=0; i<n_mat; i++) {
    for (int j=i; j<n_mat; j++) {
      MyMatrix<T> prod = l_mat[i] * l_mat[j];
      T scal(0);
      for (int i=0; i<n; i++) {
        scal += prod(i,i);
      }
      PairwiseScalar(i, j) = scal;
      PairwiseScalar(j, i) = scal;
    }
  }
  MyMatrix<T> PairwiseScalarInv = Inverse(PairwiseScalar);
  SelfDualInfo<T> self_dual_info{PairwiseScalar, PairwiseScalarInv};
  LinSpa.self_dual_info = self_dual_info;
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

/*
  Compute the orthogonal projector from the subspace and the gram matrix
 */
template <typename T>
MyMatrix<T> GetOrthogonalProjectorMatrix(MyMatrix<T> const &eG,
                                         MyMatrix<T> const &subspace) {
  int n = eG.rows();
  int dim_space = subspace.rows();
  MyMatrix<T> prod = subspace * eG;
  MyMatrix<T> OrthSpace = NullspaceTrMat(prod);
  MyMatrix<T> FullBasis = Concatenate(subspace, OrthSpace);
  MyMatrix<T> InvFullBasis = Inverse(FullBasis);
  MyMatrix<T> DiagMat = ZeroMatrix<T>(n, n);
  for (int i = 0; i < dim_space; i++) {
    DiagMat(i, i) = 1;
  }
  MyMatrix<T> RetMat = InvFullBasis * DiagMat * FullBasis;
  return RetMat;
}

/*
  Computes the set of matrices to be used for the product preservation
 */
template <typename T>
std::vector<MyMatrix<T>>
GetFamilyDiscMatrices(MyMatrix<T> const &eG,
                      std::vector<MyMatrix<T>> const &ListComm,
                      std::vector<MyMatrix<T>> const &ListSubspace) {
  std::vector<MyMatrix<T>> ListDisc{eG};
  for (auto &eComm : ListComm) {
    MyMatrix<T> prod = eComm * eG;
    ListDisc.emplace_back(std::move(prod));
  }
  for (auto &subspace : ListSubspace) {
    MyMatrix<T> ProjMat = GetOrthogonalProjectorMatrix(eG, subspace);
    MyMatrix<T> prod = ProjMat * eG;
    ListDisc.emplace_back(std::move(prod));
  }
  return ListDisc;
}

template <typename T>
bool is_stab_space(MyMatrix<T> const &Pmat, LinSpaceMatrix<T> const &LinSpa) {
  for (auto &eMatSp : LinSpa.ListMat) {
    MyMatrix<T> eMatSpImg = Pmat * eMatSp * Pmat.transpose();
    MyVector<T> eMatSpImg_V = SymmetricMatrixToVector(eMatSpImg);
    std::optional<MyVector<T>> opt =
        SolutionMat(LinSpa.ListMatAsBigMat, eMatSpImg_V);
    if (!opt) {
      return false;
    }
  }
  return true;
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

template <typename T, typename Telt>
MyMatrix<T> get_mat_from_shv_perm(Telt const &elt, MyMatrix<T> const &SHV_T,
                                  [[maybe_unused]] MyMatrix<T> const &eMat) {
  std::optional<MyMatrix<T>> opt = FindTransformationGeneral(SHV_T, SHV_T, elt);
  MyMatrix<T> Pmat = unfold_opt(opt, "Failed to get transformation");
#ifdef DEBUG_TSPACE_FUNCTIONS
  if (!IsIntegralMatrix(Pmat)) {
    std::cerr << "TSPACE: The matrix TransMat should be integral\n";
    throw TerminalException{1};
  }
  MyMatrix<T> eMatImg = Pmat * eMat * Pmat.transpose();
  if (eMatImg != eMat) {
    std::cerr << "TSPACE: The matrix TransMat does not preserve eMat\n";
    throw TerminalException{1};
  }
#endif
  return Pmat;
}

template <typename T, typename Telt>
std::optional<MyMatrix<T>>
is_corr_and_solve(Telt const &elt, MyMatrix<T> const &SHV_T,
                  MyMatrix<T> const &eMat, LinSpaceMatrix<T> const &LinSpa) {
  MyMatrix<T> Pmat = get_mat_from_shv_perm(elt, SHV_T, eMat);
  if (!is_stab_space(Pmat, LinSpa)) {
    return {};
  }
  return Pmat;
}

template <typename T, typename Telt> struct PermutationBuilder {
public:
  using Tidx = typename Telt::Tidx;
  int n_row;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, Tidx> MapV;
  PermutationBuilder(MyMatrix<T> const &SHV_T) {
    n_row = SHV_T.rows();
    for (int i_row = 0; i_row < n_row; i_row++) {
      MyVector<T> eV = GetMatrixRow(SHV_T, i_row);
      Tidx i_row_idx = static_cast<Tidx>(i_row);
      ListV.push_back(eV);
      MapV[eV] = i_row_idx;
    }
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
    if (ListV.size() != MapV.size()) {
      std::cerr << "TSPACE: PermutationBuilder, we have duplication |ListV|="
                << ListV.size() << " |MapV|=" << MapV.size() << "\n";
      throw TerminalException{1};
    }
#endif
  }
  Telt get_permutation(MyMatrix<T> const &M,
                       [[maybe_unused]] std::ostream &os) {
    std::vector<Tidx> eList(n_row);
    for (int i_row = 0; i_row < n_row; i_row++) {
      MyVector<T> Vimg = M.transpose() * ListV[i_row];
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
      if (MapV.count(Vimg) == 0) {
        std::cerr << "TSPACE: MapV should contain i_row=" << i_row
                  << " Vimg=" << StringVectorGAP(Vimg) << "\n";
        std::cerr << "TSPACE: M=\n";
        WriteMatrix(std::cerr, M);
        for (int x_row = 0; x_row < n_row; x_row++) {
          std::cerr << "TSPACE: x_row=" << i_row
                    << " V=" << StringVectorGAP(ListV[x_row]) << "\n";
        }
        throw TerminalException{1};
      }
#endif
      Tidx pos = MapV.at(Vimg);
      eList[i_row] = pos;
    }
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: get_permutation, n_row=" << n_row << " eList=[";
    for (int i_row = 0; i_row < n_row; i_row++) {
      int val = eList[i_row];
      os << " " << val;
    }
    os << " ]\n";
#endif
    return Telt(std::move(eList));
  }
};

template <typename T, typename Telt>
Telt get_elt_from_matrix(MyMatrix<T> const &mat, MyMatrix<T> const &SHV_T, std::ostream &os) {
  PermutationBuilder<T, Telt> builder(SHV_T);
  return builder.get_permutation(mat, os);
}



template <typename T, typename Telt>
std::vector<Telt> get_list_elt_from_list_matrices(std::vector<MyMatrix<T>> const &l_matr,
                                                  MyMatrix<T> const &SHV_T,
                                                  std::ostream &os) {
  PermutationBuilder<T, Telt> builder(SHV_T);
  std::vector<Telt> LGenGlobStab_perm;
#ifdef DEBUG_TSPACE_FUNCTIONS
  size_t pos = 0;
#endif
  for (auto &eMatr : l_matr) {
    Telt ePerm = builder.get_permutation(eMatr, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: get_list_elt_from_list_matrices pos=" << pos
       << " ePerm=" << ePerm << "\n";
    pos += 1;
#endif
    LGenGlobStab_perm.emplace_back(std::move(ePerm));
  }
  return LGenGlobStab_perm;
}

template <typename T, typename Tgroup>
Tgroup get_perm_group_from_list_matrices(std::vector<MyMatrix<T>> const &l_matr,
                                         MyMatrix<T> const &SHV_T,
                                         std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n_row = SHV_T.rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: get_perm_group_from_list_matrices n_row=" << n_row << "\n";
#endif
  std::vector<Telt> LGenGlobStab_perm = get_list_elt_from_list_matrices<T,Telt>(l_matr, SHV_T, os);
  return Tgroup(LGenGlobStab_perm, n_row);
}

template <typename T, typename Tgroup> struct Result_ComputeStabilizer_SHV {
  using Telt = typename Tgroup::Telt;
  std::optional<std::vector<MyMatrix<T>>> l_gens;
  std::optional<std::pair<std::vector<Telt>, Tgroup>> perms_and_group;
  std::vector<MyMatrix<T>>
  get_list_matrix(MyMatrix<T> const &SHV_T, MyMatrix<T> const &eMat,
                  LinSpaceMatrix<T> const &LinSpa,
                  [[maybe_unused]] std::ostream &os) const {
    if (l_gens) {
      return *l_gens;
    }
    if (perms_and_group) {
      std::vector<MyMatrix<T>> LGenGlobStab_matr;
      std::vector<Telt> const &LGenGlobStab_perm = perms_and_group->first;
      for (auto &eGen : LGenGlobStab_perm) {
        std::optional<MyMatrix<T>> opt =
            is_corr_and_solve(eGen, SHV_T, eMat, LinSpa);
        MyMatrix<T> eGenMatr = unfold_opt(opt, "Failed to unfold");
        LGenGlobStab_matr.push_back(eGenMatr);
      }
      return LGenGlobStab_matr;
    }
    std::cerr << "We did not generate l_gens / perms_and_group\n";
    throw TerminalException{1};
  }
  std::vector<Telt> get_list_perms(MyMatrix<T> const &SHV_T,
                                   std::ostream &os) const {
    if (l_gens) {
      PermutationBuilder<T, Telt> builder(SHV_T);
      std::vector<MyMatrix<T>> const &l_matr = *l_gens;
      std::vector<Telt> LGenGlobStab_perm;
      for (auto &eGen : l_matr) {
        Telt ePerm = builder.get_permutation(eGen, os);
        LGenGlobStab_perm.emplace_back(std::move(ePerm));
      }
      return LGenGlobStab_perm;
    }
    if (perms_and_group) {
      return perms_and_group->first;
    }
    std::cerr << "We did not generate l_gens / perms_and_group\n";
    throw TerminalException{1};
  }
  Tgroup get_perm_group(MyMatrix<T> const &SHV_T, std::ostream &os) const {
    if (l_gens) {
      std::vector<MyMatrix<T>> const &l_matr = *l_gens;
      return get_perm_group_from_list_matrices<T, Tgroup>(l_matr, SHV_T, os);
    }
    if (perms_and_group) {
      return perms_and_group->second;
    }
    std::cerr << "We did not generate l_gens / perms_and_group\n";
    throw TerminalException{1};
  }
};

template <typename T, typename Tgroup>
Result_ComputeStabilizer_SHV<T, Tgroup>
get_from_gens(std::vector<MyMatrix<T>> const &l_gens) {
  return {l_gens, {}};
}

template <typename T, typename Tgroup>
Result_ComputeStabilizer_SHV<T, Tgroup> get_from_perms_and_group(
    std::pair<std::vector<typename Tgroup::Telt>, Tgroup> const &pair) {
  return {{}, pair};
}

/*
  For a positive definite matrix in the T-space, we compute the group
  of transformations that preserves:
  * The positive definite quadratic form
  * The T-space itself

  We incorporate all the ideas that we have put forward:
  * Use of the ListComm
  * Use of the subspaces to build commutting projectors which are stabilized.
  * If all else fails, use of the cosets.
 */
template <typename T, typename Tgroup>
Result_ComputeStabilizer_SHV<T, Tgroup>
LINSPA_ComputeStabilizer_SHV(LinSpaceMatrix<T> const &LinSpa,
                             MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using LeftCosets = typename Tgroup::LeftCosets;
  using Tfield = T;
  int n_row = SHV_T.rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer_SHV n_row=" << n_row << "\n";
#endif
  std::vector<T> Vdiag(n_row, 0);
  std::vector<MyMatrix<T>> ListMat =
      GetFamilyDiscMatrices(eMat, LinSpa.ListComm, LinSpa.ListSubspaces);

  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tgroup>(SHV_T, ListMat,
                                                              Vdiag, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer_SHV |ListGen|=" << ListGen.size()
     << "\n";
#endif
  //
  // Try the direct strategy and hopes to be lucky
  //
  auto get_generators = [&]() -> std::optional<std::vector<MyMatrix<T>>> {
    std::vector<MyMatrix<T>> ListTransMat;
    for (auto &eGen : ListGen) {
      Telt elt(eGen);
      std::optional<MyMatrix<T>> opt =
          is_corr_and_solve(elt, SHV_T, eMat, LinSpa);
      if (opt) {
        ListTransMat.push_back(*opt);
      } else {
        return {};
      }
    }
    return ListTransMat;
  };
  std::optional<std::vector<MyMatrix<T>>> opt = get_generators();
  if (opt) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: LINSPA_ComputeStabilizer_SHV success of the direct "
          "approach\n";
#endif
    return get_from_gens<T, Tgroup>(*opt);
  }
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets for
  // resolving that.
  //
  std::vector<Telt> LGenPerm;
  for (auto &eList : ListGen) {
    Telt ePerm(eList);
    LGenPerm.push_back(ePerm);
  }
  Tgroup FullGRP(LGenPerm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer_SHV |FullGRP|=" << FullGRP.size()
     << "\n";
#endif
  PermutationBuilder<T, Telt> builder(SHV_T);
  std::vector<Telt> LGenGlobStab_perm;
  for (auto &eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder.get_permutation(eGen, os);
    LGenGlobStab_perm.push_back(ePerm);
  }
  Tgroup GRPsub(LGenGlobStab_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer_SHV |GRPsub|=" << GRPsub.size()
     << "\n";
#endif
  auto try_upgrade = [&]() -> std::optional<Telt> {
    // We can use either the left or right cosets.
    // This is because the right thing to use is the double cosets.
    // However, we do not have the formalism for having iterator
    // over the double cosets. We build all of them.
    //
    // For the left/right cosets we have efficient iterators
    // and that is why we use them here. We cannot afford at all
    // to enumerate all the double cosets because we will have
    // some scenario where FullGRP is big, GRPsub very small and
    // that would mean enumerating all the elements of the group.
    //
    // Left  transversals are g H
    // Right transversals are H g
    LeftCosets rc = FullGRP.left_cosets(GRPsub);
    for (auto &eCosReprPerm : rc) {
      std::optional<MyMatrix<T>> opt =
          is_corr_and_solve(eCosReprPerm, SHV_T, eMat, LinSpa);
      if (opt) {
        // We have this problem that the first cosets is not necessarily the one
        // of GRPsub and that the coset of the GRPsub is also not necessarily
        // the identity.
        if (!GRPsub.isin(eCosReprPerm)) {
#ifdef DEBUG_TSPACE_FUNCTIONS
          os << "TSPACE: LINSPA_ComputeStabilizer_SHV Finding a new "
                "eCosReprPerm\n";
#endif
          return eCosReprPerm;
        }
      }
    }
    return {};
  };
  while (true) {
    std::optional<Telt> opt = try_upgrade();
    if (opt) {
      // Found another stabilizing element, upgrading the group and retry.
      LGenGlobStab_perm.push_back(*opt);
      GRPsub = Tgroup(LGenGlobStab_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
      os << "TSPACE: LINSPA_ComputeStabilizer_SHV Now |GRPsub|="
         << GRPsub.size() << "\n";
#endif
    } else {
      break;
    }
  }
  std::pair<std::vector<typename Tgroup::Telt>, Tgroup> pair{
      std::move(LGenGlobStab_perm), std::move(GRPsub)};
  return get_from_perms_and_group<T, Tgroup>(pair);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
LINSPA_ComputeStabilizer(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &eMat, std::ostream &os) {
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  Result_ComputeStabilizer_SHV<T, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<T, Tgroup>(LinSpa, eMat, SHV_T, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer, we have result\n";
#endif
  return result.get_list_matrix(SHV_T, eMat, LinSpa, os);
}

template <typename T>
size_t LINSPA_Invariant_SHV(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                            MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
                            std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, 0);
  std::vector<MyMatrix<T>> ListMat =
      GetFamilyDiscMatrices(eMat, LinSpa.ListComm, LinSpa.ListSubspaces);
  return GetInvariant_ListMat_Vdiag<T, Tfield>(seed, SHV_T, ListMat, Vdiag, os);
}

template <typename T, typename Tint>
size_t LINSPA_Invariant(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &eMat, std::ostream &os) {
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  return LINSPA_Invariant_SHV<T>(seed, LinSpa, eMat, SHV_T, os);
}

/*
  For two positive definite matrices M1 find if it exists a transformation P
  such that
  * P M1 P^T = M2
  * P LinSpa.ListMat P^T  image is LinSpa.ListMat
*/
template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LINSPA_TestEquivalenceGramMatrix_SHV(
    LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &eMat1,
    MyMatrix<T> const &eMat2, MyMatrix<T> const &SHV1_T,
    MyMatrix<T> const &SHV2_T, std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using LeftCosets = typename Tgroup::LeftCosets;
  //  using Tfield = T;
  if (SHV1_T.rows() != SHV2_T.rows()) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, Exiting here at |SHV1| <> |SHV2|\n";
#endif
    return {};
  }
  int n_row = SHV1_T.rows();
  std::vector<T> Vdiag1(n_row, 0), Vdiag2(n_row, 0);
  std::vector<MyMatrix<T>> ListMat1 =
      GetFamilyDiscMatrices(eMat1, LinSpa.ListComm, LinSpa.ListSubspaces);
  std::vector<MyMatrix<T>> ListMat2 =
      GetFamilyDiscMatrices(eMat2, LinSpa.ListComm, LinSpa.ListSubspaces);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Equiv, |ListComm|=" << LinSpa.ListComm.size()
     << " |ListSubspaces|=" << LinSpa.ListSubspaces.size() << "\n";
  os << "TSPACE: Equiv, |ListMat1|=" << ListMat1.size()
     << " |ListMat1|=" << ListMat1.size() << "\n";
#endif
  std::optional<std::vector<Tidx>> opt1 =
      TestEquivalence_ListMat_Vdiag<T, Tfield, Tidx>(
          SHV1_T, ListMat1, Vdiag1, SHV2_T, ListMat2, Vdiag2, os);
  if (!opt1) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, Exiting here at opt1\n";
#endif
    return {};
  }
  //
  // Building one equivalence that should preserve the basics and which we can
  // work on.
  //
  Telt eltEquiv(*opt1);
  Telt eltInv = Inverse(eltEquiv);
  std::optional<MyMatrix<T>> opt2 =
      FindTransformationGeneral(SHV2_T, SHV1_T, eltInv);
  MyMatrix<T> OneEquiv = unfold_opt(opt2, "Failed to get transformation");
#ifdef DEBUG_TSPACE_FUNCTIONS
  if (!IsIntegralMatrix(OneEquiv)) {
    std::cerr << "TSPACE: Equiv, The matrix TransMat should be integral\n";
    throw TerminalException{1};
  }
  MyMatrix<T> eMat1_img = OneEquiv * eMat1 * OneEquiv.transpose();
  if (eMat1_img != eMat2) {
    std::cerr << "TSPACE: Equiv, The matrix TransMat does not preserve eMat\n";
    throw TerminalException{1};
  }
#endif
  if (is_stab_space(OneEquiv, LinSpa)) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, Direct approach success, no need to go further\n";
#endif
    return OneEquiv;
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Equiv, Direct approach failure, computing stabilizer and "
        "iterating\n";
#endif
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets for
  // resolving that.
  //
  std::vector<std::vector<Tidx>> ListGen1 =
      GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tgroup>(SHV1_T, ListMat1,
                                                              Vdiag1, os);
  std::vector<Telt> LGenPerm1;
  for (auto &eList1 : ListGen1) {
    Telt ePerm1(eList1);
    LGenPerm1.push_back(ePerm1);
  }
  Tgroup FullGRP1(LGenPerm1, n_row);
  PermutationBuilder<T, Telt> builder1(SHV1_T);
  std::vector<Telt> LGenGlobStab1_perm;
  for (auto &eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder1.get_permutation(eGen, os);
    LGenGlobStab1_perm.push_back(ePerm);
  }
  Tgroup GRPsub1(LGenGlobStab1_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Equiv, n_row=" << n_row << " |FullGRP1|=" << FullGRP1.size()
     << " |GRPsub1|=" << GRPsub1.size() << "\n";
#endif
#ifdef DEBUG_TSPACE_FUNCTIONS_DISABLE
  auto f_get_group_size = [&]() -> size_t {
    size_t n_elt = 0;
    for (auto &elt : FullGRP1) {
      MyMatrix<T> eMatr = get_mat_from_shv_perm(elt, SHV1_T, eMat1);
      if (is_stab_space(eMatr, LinSpa)) {
        n_elt += 1;
      }
    }
    return n_elt;
  };
  size_t n_elt = f_get_group_size();
  os << "TSPACE: Equiv, |FullGRP1|=" << FullGRP1.size() << " n_elt=" << n_elt
     << " |GRPsub1|=" << GRPsub1.size() << "\n";
  size_t pos_equiv_grp = 0;
  os << "TSPACE: Equiv(" << pos_equiv_grp << "), |GRPsub1|=" << GRPsub1.size()
     << "\n";
#endif
#ifdef SANITY_CHECK_EXTENSIVE_TSPACE_FUNCTIONS
  auto f_exhaustive = [&]() -> std::optional<MyMatrix<T>> {
    os << "TSPACE: Starting extensive check over |FullGRP1|=" << FullGRP1.order() << "\n";
    for (auto &elt : FullGRP1) {
      MyMatrix<T> eMatr = get_mat_from_shv_perm(elt, SHV1_T, eMat1);
      MyMatrix<T> eProd_T = OneEquiv * eMatr;
      if (is_stab_space(eProd_T, LinSpa)) {
        return eProd_T;
      }
    }
    return {};
  };
#endif
  struct PartSol {
    std::optional<Telt> new_gen;
    std::optional<MyMatrix<T>> sol;
  };
  auto try_solution = [&]() -> PartSol {
  // The left cosets are the cosets such as G = \cup_c c H
  // We need to use left cosets for the computation
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, |FullGRP1|=" << FullGRP1.size()
       << " |GRPsub1|=" << GRPsub1.size() << "\n";
#endif
    LeftCosets rc = FullGRP1.left_cosets(GRPsub1);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, we have rc\n";
#endif
    for (auto &eCosReprPerm : rc) {
      MyMatrix<T> eCosReprMatr =
          get_mat_from_shv_perm(eCosReprPerm, SHV1_T, eMat1);
      MyMatrix<T> eProd = OneEquiv * eCosReprMatr;
      if (is_stab_space(eProd, LinSpa)) {
#ifdef DEBUG_TSPACE_FUNCTIONS
        os << "TSPACE: Equiv, We found an equivalence\n";
#endif
#ifdef DEBUG_TSPACE_FUNCTIONS
        MyMatrix<T> eMat1_img = eProd * eMat1 * eProd.transpose();
        if (eMat1_img != eMat2) {
          std::cerr << "TSPACE: Equiv, we do not have an equivalence\n";
          throw TerminalException{1};
        }
#endif
        return {{}, eProd};
      }
      if (is_stab_space(eCosReprMatr, LinSpa)) {
        // We have this problem that the first cosets is not necessarily the one
        // of GRPsub and that the coset of the GRPsub is also not necessarily
        // the identity.
        if (!GRPsub1.isin(eCosReprPerm)) {
#ifdef DEBUG_TSPACE_FUNCTIONS
          os << "TSPACE: Equiv, We found some new generator\n";
#endif
          return {eCosReprPerm, {}};
        }
      }
    }
    return {{}, {}};
  };
  while (true) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, before try_solution\n";
#endif
    PartSol p_sol = try_solution();
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, after try_solution\n";
#endif
    if (p_sol.sol) {
      MyMatrix<T> const &Pmat_T = *p_sol.sol;
#ifdef SANITY_CHECK_EXTENSIVE_TSPACE_FUNCTIONS
      if (!f_exhaustive()) {
        std::cerr << "TSPACE: We found equiv with one method but the "
                     "exhaustive does not\n";
        throw TerminalException{1};
      }
#endif
#ifdef DEBUG_TSPACE_FUNCTIONS
      os << "TSPACE: Equiv, before returning Pmat\n";
#endif
      return Pmat_T;
    }
    if (!p_sol.new_gen) {
#ifdef SANITY_CHECK_EXTENSIVE_TSPACE_FUNCTIONS
      if (f_exhaustive()) {
        std::cerr << "TSPACE: We found non-equiv with one method but the "
                     "exhaustive does\n";
        throw TerminalException{1};
      }
#endif
#ifdef DEBUG_TSPACE_FUNCTIONS
      os << "TSPACE: Equiv, before returning None\n";
#endif
      return {};
    }
    LGenGlobStab1_perm.push_back(*p_sol.new_gen);
    GRPsub1 = Tgroup(LGenGlobStab1_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS_DISABLE
    pos_equiv_grp += 1;
    os << "TSPACE: Equiv(" << pos_equiv_grp << "), |GRPsub1|=" << GRPsub1.size()
       << "\n";
#endif
  }
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LINSPA_TestEquivalenceGramMatrix(LinSpaceMatrix<T> const &LinSpa,
                                 MyMatrix<T> const &eMat1,
                                 MyMatrix<T> const &eMat2, std::ostream &os) {
  T det1 = DeterminantMat(eMat1);
  T det2 = DeterminantMat(eMat2);
  if (det1 != det2) {
    return {};
  }
  MyMatrix<Tint> SHV1 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat1, os);
  MyMatrix<Tint> SHV2 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat2, os);
  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(SHV2);
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(LinSpa, eMat1, eMat2,
                                                      SHV1_T, SHV2_T, os);
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &M_T = *opt;
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
  return M;
}

/*
  Compute an invariant of the gram matrix
  so that it does not change after arithmetic
  equivalence.
  However, it is not invariant under scaling.
*/
template <typename T, typename Tint>
size_t GetInvariantGramShortest(MyMatrix<T> const &eGram,
                                MyMatrix<Tint> const &SHV, size_t const &seed,
                                [[maybe_unused]] std::ostream &os) {
  T eDet = DeterminantMat(eGram);
  int nbVect = SHV.rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: eDet=" << eDet << " nbVect=" << nbVect << "\n";
#endif
  std::vector<MyVector<T>> ListV;
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<Tint> V = GetMatrixRow(SHV, iVect);
    MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
    ListV.push_back(V_T);
  }
  std::map<T, size_t> map_diag, map_off_diag;
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<T> Vprod = eGram * ListV[iVect];
    T eNorm = Vprod.dot(ListV[iVect]);
    map_diag[eNorm] += 1;
    for (int jVect = iVect + 1; jVect < nbVect; jVect++) {
      T eScal = Vprod.dot(ListV[jVect]);
      map_off_diag[eScal] += 1;
    }
  }
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t hash_ret = seed;
  size_t hash_det = std::hash<T>()(eDet);
  combine_hash(hash_ret, hash_det);
  for (auto &kv : map_diag) {
    size_t hash_T = std::hash<T>()(kv.first);
    combine_hash(hash_ret, hash_T);
    combine_hash(hash_ret, kv.second);
  }
  for (auto &kv : map_off_diag) {
    size_t hash_T = std::hash<T>()(kv.first);
    combine_hash(hash_ret, hash_T);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
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

// clang-format off
#endif  // SRC_LATT_TSPACE_FUNCTIONS_H_
// clang-format on
