// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_FUNCTIONS_H_
#define SRC_LATT_TSPACE_FUNCTIONS_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "ShortestUniversal.h"
#include "PolytopeEquiStab.h"
#include "Positivity.h"
#include "SHORT_ShortestConfig.h"
#include <set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_TSPACE_FUNCTIONS
#endif

/*
  By a T-space, we mean a vector space which intersect the cone of positive definite
  matrices. We are interested here only on those positive definite elements.
  ---
  The motivation for this kind of construction is for finding good packing, covering, etc.
  The idea is that the smaller dimensionality allows for a smaller number of
  combinatorial possibilities.
  ---
  However, if the total number of combinatorial possibilities is expected to be smaller,
  the combinatorial complexity on programming and thinking is definitely higher.
  The thing that are considered in this section:
  * Finding the integral saturation (useful for integral separation)
  * The generation of the T-spaces (see other include)
  * The computation of stabilizer and equivalence of positive definite matrices in
  the T-space.
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
  // The list of matrices with which the global stabilizing elements must commute
  std::vector<MyMatrix<T>> ListComm;
  // The list of preserved subspaces by the element of the global stabilizer
  std::vector<MyMatrix<T>> ListSubspaces;
  // The point stabilizer of the T-space
  std::vector<MyMatrix<T>> PtStabGens;
};

template<typename T>
MyMatrix<T> GetListMatAsBigMat(std::vector<MyMatrix<T>> const& ListMat) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "We have n_mat = 0\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  int sym_dim = (n*(n+1)) / 2;
  MyMatrix<T> BigMat(n_mat, sym_dim);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    MyVector<T> V = SymmetricMatrixToVector(ListMat[i_mat]);
    AssignMatrixRow(BigMat, i_mat, V);
  }
  return BigMat;
}

//
// We search for the set of matrices satisfying g M g^T = M for all g in ListGen
//
template <typename T>
std::vector<MyMatrix<T>>
BasisInvariantForm(int const &n, std::vector<MyMatrix<T>> const &ListGen) {
  std::vector<std::vector<int>> ListCoeff;
  MyMatrix<int> ListCoeffRev(n,n);
  int pos = 0;
  for (int iLin = 0; iLin < n; iLin++) {
    for (int iCol = 0; iCol <= iLin; iCol++) {
      ListCoeff.push_back({iLin, iCol});
      ListCoeffRev(iLin, iCol) = pos;
      ListCoeffRev(iCol, iLin) = pos;
      pos++;
    }
  }
  int nbCoeff = pos;
  auto FuncPos = [&](int const &i, int const &j) -> int {
    return ListCoeffRev(i, j);
  };
  std::vector<MyVector<T>> ListEquations;
  for (auto &eGen : ListGen)
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++) {
        MyVector<T> TheEquation = ZeroVector<T>(nbCoeff);
        TheEquation(FuncPos(i, j)) += 1;
        for (int k = 0; k < n; k++)
          for (int l = 0; l < n; l++) {
            int pos = FuncPos(k, l);
            TheEquation(pos) += -eGen(i, k) * eGen(j, l);
          }
        ListEquations.push_back(TheEquation);
      }
  MyMatrix<T> MatEquations = MatrixFromVectorFamily(ListEquations);
  std::cerr << "Before call to NullspaceTrMat nbGen=" << ListGen.size()
            << " MatEquations(rows/cols)=" << MatEquations.rows() << " / "
            << MatEquations.cols() << "\n";
  MyMatrix<T> NSP = NullspaceTrMat(MatEquations);
  std::cerr << "After call to NullspaceTrMat NSP.rows=" << NSP.rows()
            << " cols=" << NSP.cols() << "\n";
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    MyMatrix<T> eMat(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        eMat(i, j) = eRow(FuncPos(i, j));
    TheBasis[iDim] = eMat;
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  for (auto & eBasis : TheBasis) {
    for (auto & eGen : ListGen) {
      MyMatrix<T> eProd = eGen * eBasis * eGen.transpose();
      if (eProd != eBasis) {
        std::cerr << "We have eProd <> eBasis, so a linear algebra bug\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return TheBasis;
}



template<typename T>
LinSpaceMatrix<T> BuildLinSpace(MyMatrix<T> const& SuperMat, std::vector<MyMatrix<T>> const& ListMat, std::vector<MyMatrix<T>> const& ListComm) {
  int n = SuperMat.rows();
  std::vector<std::vector<T>> ListLineMat;
  for (auto & eMat : ListMat) {
    std::vector<T> eV = GetLineVector(eMat);
    ListLineMat.push_back(eV);
  }
  MyMatrix<T> BigMat = GetListMatAsBigMat(ListMat);
  std::vector<MyMatrix<T>> ListSubspaces;
  MyMatrix<T> eGen = -IdentityMat<T>(n);
  std::vector<MyMatrix<T>> PtStab{eGen};
  // It may be actually a Bravais space, but setting up to false avoids potential problems.
  bool isBravais = false;
  return {n, isBravais, SuperMat, ListMat, ListLineMat, BigMat, ListComm, ListSubspaces, PtStab};
}

template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> ComputePointStabilizerTspace(MyMatrix<T> const& SuperMat, std::vector<MyMatrix<T>> const& ListMat, std::ostream & os) {
  using Tidx = uint32_t;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(SuperMat);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  std::vector<T> Vdiag(SHV_T.rows(), 0);
  std::vector<std::vector<Tidx>> ListGenPerm =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx>(SHV_T, ListMat, Vdiag, os);
  std::vector<MyMatrix<Tint>> ListGenMatr;
  for (auto & eList : ListGenPerm) {
    MyMatrix<T> eMatr_T = FindTransformation_vect(SHV_T, SHV_T, eList);
    MyMatrix<Tint> eMatr = UniversalMatrixConversion<Tint,T>(eMatr_T);
    ListGenMatr.push_back(eMatr);
  }
  return ListGenMatr;
}

template<typename T, typename Tint>
MyMatrix<T> GetOnePositiveDefiniteMatrix(std::vector<MyMatrix<T>> const& ListMat, std::ostream& os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "The number of matrices is 0 so we cannot build a positive definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  std::vector<MyVector<Tint>> ListV;
  for (int i=0; i<n; i++) {
    MyVector<Tint> V = ZeroVector<Tint>(n);
    V(i) = 1;
    ListV.push_back(V);
  }
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      for (int sign=0; sign<2; sign++) {
        int signB = -1 + 2*sign;
        MyVector<Tint> V = ZeroVector<Tint>(n);
        V(i) = 1;
        V(j) = signB;
        ListV.push_back(V);
      }
    }
  }
  while(true) {
    int n_vect = ListV.size();
    MyMatrix<T> ListIneq(n_vect, 1 + n_mat);
    MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_mat);
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      ListIneq(i_vect,0) = -1;
      MyVector<Tint> V = ListV[i_vect];
      MyVector<T> V_T = UniversalVectorConversion<T,Tint>(V);
      for (int i_mat=0; i_mat<n_mat; i_mat++) {
        T val = EvaluationQuadForm(ListMat[i_mat], V_T);
        ListIneq(i_vect, 1 + i_mat) = val;
        ToBeMinimized(1+i_mat) += val;
      }
    }
    //
    // Solving the linear program
    //
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.PrimalDefined || !eSol.DualDefined) {
      std::cerr << "The LpSolution dual and primal solutions should be defined\n";
      throw TerminalException{1};
    }
    MyMatrix<T> TrySuperMat = ZeroMatrix<T>(n, n);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      TrySuperMat += eSol.DirectSolution(i_mat) * ListMat[i_mat];
    }
    if (IsPositiveDefinite(TrySuperMat)) {
      return TrySuperMat;
    }
    //
    // Failed, trying to find a vector
    //
    auto get_one_vect=[&]() -> MyVector<Tint> {
      T CritNorm = 0;
      if (RankMat(TrySuperMat) < n) {
        return GetShortVectorDegenerate<T, Tint>(TrySuperMat, CritNorm);
      } else {
        bool NeedNonZero = true;
        bool StrictIneq = true;
        return GetShortVector_unlimited_float<Tint, T>(TrySuperMat, CritNorm, StrictIneq, NeedNonZero);
      }
    };
    MyVector<Tint> V = get_one_vect();
    ListV.push_back(V);
  }
}

template<typename T>
MyMatrix<T> GetRandomPositiveDefinite(LinSpaceMatrix<T> const& LinSpa) {
  int n = LinSpa.n;
  MyMatrix<T> TheMat = ZeroMatrix<T>(n, n);
  int N = 2;
  for (auto & eMat : LinSpa.ListMat) {
    int coef = random() % (2 * N + 1) - N;
    TheMat += coef * eMat;
  }
  while(true) {
    if (IsPositiveDefinite(TheMat)) {
      return TheMat;
    }
    TheMat += LinSpa.SuperMat;
  }
}

template<typename T, typename Tint>
bool IsSymmetryGroupCorrect(MyMatrix<T> const& GramMat, LinSpaceMatrix<T> const& LinSpa, std::ostream & os) {
  using Tidx = uint32_t;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row,0);
  std::vector<MyMatrix<T>> ListMat = {GramMat};
  std::vector<std::vector<Tidx>> ListGen =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx>(SHV_T, ListMat, Vdiag, os);
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
      FindTransformationGeneral_vect(SHV_T, SHV_T, eList);
    if (!opt) {
      std::cerr << "Failed to find the matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &M_T = *opt;
    if (!IsIntegralMatrix(M_T)) {
      std::cerr << "Bug: The matrix should be integral\n";
      throw TerminalException{1};
    }
    for (auto & eMat : LinSpa.ListMat) {
      MyMatrix<T> eMatImg = M_T * eMat * M_T.transpose();
      if (eMatImg != eMat) {
        return false;
      }
    }
  }
  return true;
}



template<typename T, typename Tint>
bool IsBravaisSpace(int n, std::vector<MyMatrix<T>> const& ListMat, std::vector<MyMatrix<Tint>> const& ListGen) {
  std::vector<MyMatrix<T>> ListGen_T;
  for (auto & eGen : ListGen) {
    MyMatrix<T> eGen_T = UniversalMatrixConversion<T,Tint>(eGen);
    ListGen_T.push_back(eGen_T);
  }
  std::vector<MyMatrix<T>> BasisInv = BasisInvariantForm(n, ListGen_T);
  MyMatrix<T> Big_ListMat = GetListMatAsBigMat(ListMat);
  MyMatrix<T> Big_BasisInv = GetListMatAsBigMat(BasisInv);
  //
#ifdef DEBUG_TSPACE_FUNCTIONS
  if (!IsSubspaceContained(Big_ListMat, Big_BasisInv)) {
    std::cerr << "The elements of ListMat are not in the invariant space which is not acceptable\n";
    throw TerminalException{1};
  }
#endif
  //
  return IsSubspaceContained(Big_BasisInv, Big_ListMat);
}


/*
  We want to find symmetries of the polytope that are pointwise stabilizer of the T-space.
  Two behaviors that we want to avoid:
  * Not globally stabilizing the T-space.
  * Having an action on the T-space that is non-trivial.
 */
template<typename T, typename Tint>
MyMatrix<T> GetRandomPositiveDefiniteNoNontrialSymm(LinSpaceMatrix<T> const& LinSpa, std::ostream & os) {
  while(true) {
    MyMatrix<T> TheMat = GetRandomPositiveDefinite(LinSpa);
    bool test = IsSymmetryGroupCorrect<T,Tint>(TheMat, LinSpa, os);
    if (test) {
      return TheMat;
    }
  }
}

/*
  We need the space of matrices to be spanning the integral saturation so that the elements
  of GlStab are integral.
 */
template<typename T>
std::vector<MyMatrix<T>> IntegralSaturationSpace(std::vector<MyMatrix<T>> const& ListMat) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "We have n_mat=0\n";
    std::cerr << "The code could work with n_mat=0 but we are not sure it makes sense\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  int sym_dim = (n+1) *n / 2;
  if (n_mat == sym_dim) {
    std::cerr << "We have n_mat=" << n_mat << " equal to sym_dim=" << sym_dim << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> BigMat(n_mat, sym_dim);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    int pos = 0;
    for (int i=0; i<n; i++) {
      for (int j=i; j<n; j++) {
        BigMat(i_mat, pos) = ListMat[i_mat](i,j);
        pos++;
      }
    }
  }
  MyMatrix<T> NSP1 = NullspaceIntTrMat(BigMat);
  MyMatrix<T> BigMat_renorm = NullspaceIntTrMat(NSP1);
  if (BigMat_renorm.rows() != n_mat) {
    std::cerr << "Incoherence in the dimensions\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListMatRet;
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    MyMatrix<T> eMat(n, n);
    int pos = 0;
    for (int i=0; i<n; i++) {
      for (int j=i; j<n; j++) {
        T val = BigMat_renorm(i_mat, pos);
        eMat(i, j) = val;
        eMat(j, i) = val;
      }
    }
    ListMatRet.push_back(eMat);
  }
  return ListMatRet;
}

template <typename T>
MyMatrix<T> GetMatrixFromBasis(std::vector<MyMatrix<T>> const &ListMat,
                               MyVector<T> const &eVect) {
  int n = ListMat[0].rows();
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  int nbMat = ListMat.size();
  int siz = eVect.size();
  if (siz != nbMat) {
    std::cerr << "Error in GetMatrixFromBasis\n";
    std::cerr << "  siz=" << siz << "\n";
    std::cerr << "nbMat=" << nbMat << "\n";
    std::cerr << "But they should be both equal\n";
    throw TerminalException{1};
  }
  for (int iMat = 0; iMat < nbMat; iMat++)
    RetMat += eVect(iMat) * ListMat[iMat];
  return RetMat;
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
template<typename T>
MyMatrix<T> GetOrthogonalProjectorMatrix(MyMatrix<T> const& eG, MyMatrix<T> const& subspace) {
  int n = eG.rows();
  int dim_space = subspace.rows();
  MyMatrix<T> prod = subspace * eG;
  MyMatrix<T> OrthSpace = NullspaceTrMat(prod);
  MyMatrix<T> FullBasis = Concatenate(subspace, OrthSpace);
  MyMatrix<T> InvFullBasis = Inverse(FullBasis);
  MyMatrix<T> DiagMat = ZeroMatrix<T>(n, n);
  for (int i=0; i<dim_space; i++) {
    DiagMat(i,i) = 1;
  }
  MyMatrix<T> RetMat = InvFullBasis * DiagMat * FullBasis;
  return RetMat;
}

/*
  Computes the set of matrices to be used for the product preservation
 */
template<typename T>
std::vector<MyMatrix<T>> GetFamilyDiscMatrices(MyMatrix<T> const& eG, std::vector<MyMatrix<T>> const& ListComm, std::vector<MyMatrix<T>> const& ListSubspace) {
  std::vector<MyMatrix<T>> ListDisc {eG};
  for (auto & eComm : ListComm) {
    MyMatrix<T> prod = eComm * eG;
    ListDisc.push_back(prod);
  }
  for (auto & subspace : ListSubspace) {
    MyMatrix<T> ProjMat = GetOrthogonalProjectorMatrix(eG, subspace);
    MyMatrix<T> prod = ProjMat * eG;
    ListDisc.push_back(prod);
  }
  return ListDisc;
}


template<typename T>
bool is_stab_space(MyMatrix<T> const& Pmat, LinSpaceMatrix<T> const &LinSpa) {
  for (auto & eMatSp : LinSpa.ListMat) {
    MyMatrix<T> eMatSpImg = Pmat * eMatSp * Pmat.transpose();
    MyVector<T> eMatSpImg_V = SymmetricMatrixToVector(eMatSpImg);
    std::optional<MyVector<T>> opt = SolutionMat(LinSpa.ListMatAsBigMat, eMatSpImg_V);
    if (!opt) {
      return false;
    }
  }
  return true;
}

template<typename T, typename Telt>
MyMatrix<T> get_mat_from_shv_perm(Telt const& elt, MyMatrix<T> const& SHV_T, [[maybe_unused]] MyMatrix<T> const& eMat) {
  std::optional<MyMatrix<T>> opt =
    FindTransformationGeneral(SHV_T, SHV_T, elt);
  MyMatrix<T> Pmat = unfold_opt(opt, "Failed to get transformation");
#ifdef DEBUG_TSPACE_FUNCTIONS
  if (!IsIntegralMatrix(Pmat)) {
    std::cerr << "The matrix TransMat should be integral\n";
    throw TerminalException{1};
  }
  MyMatrix<T> eMatImg = Pmat * eMat * Pmat.transpose();
  if (eMatImg != eMat) {
    std::cerr << "The matrix TransMat does not preserve eMat\n";
    throw TerminalException{1};
  }
#endif
  return Pmat;
}


template<typename T, typename Telt>
std::optional<MyMatrix<T>> is_corr_and_solve(Telt const& elt, MyMatrix<T> const& SHV_T, MyMatrix<T> const& eMat, LinSpaceMatrix<T> const &LinSpa) {
  MyMatrix<T> Pmat = get_mat_from_shv_perm(elt, SHV_T, eMat);
  if (!is_stab_space(Pmat, LinSpa)) {
    return {};
  }
  return Pmat;
}

template<typename T, typename Telt>
struct PermutationBuilder {
public:
  using Tidx = typename Telt::Tidx;
  int n_row;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, Tidx> MapV;
  PermutationBuilder(MyMatrix<T> const& SHV_T) {
    n_row = SHV_T.rows();
    for (int i_row=0; i_row<n_row; i_row++) {
      MyVector<T> eV = GetMatrixRow(SHV_T, i_row);
      Tidx i_row_idx = static_cast<Tidx>(i_row);
      ListV.push_back(eV);
      MapV[eV] = i_row_idx;
    }
  }
  Telt get_permutation(MyMatrix<T> const& M) {
    std::vector<Tidx> eList(n_row);
    for (int i_row=0; i_row<n_row; i_row++) {
      MyVector<T> Vimg = M.transpose() * ListV[i_row];
      Tidx pos = MapV.at(Vimg);
      eList[i_row] = pos;
    }
    return Telt(eList);
  }
};

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
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>> LINSPA_ComputeStabilizer(LinSpaceMatrix<T> const &LinSpa,
                                                  MyMatrix<T> const& eMat,
                                                  std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  //  using TintGroup = typename Tgroup::Tint;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row,0);
  std::vector<MyMatrix<T>> ListMat = GetFamilyDiscMatrices(eMat, LinSpa.ListComm, LinSpa.ListSubspaces);

  std::vector<std::vector<Tidx>> ListGen =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx>(SHV_T, ListMat, Vdiag, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer |ListGen|=" << ListGen.size() << "\n";
#endif
  //
  // Try the direct strategy and hopes to be lucky
  //
  auto get_generators=[&]() -> std::optional<std::vector<MyMatrix<T>>> {
    std::vector<MyMatrix<T>> ListTransMat;
    for (auto & eGen : ListGen) {
      Telt elt(eGen);
      std::optional<MyMatrix<T>> opt = is_corr_and_solve(elt, SHV_T, eMat, LinSpa);
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
    return *opt;
  }
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets for resolving that.
  //
  std::vector<Telt> LGenPerm;
  for (auto & eList : ListGen) {
    Telt ePerm(eList);
    LGenPerm.push_back(ePerm);
  }
  Tgroup FullGRP(LGenPerm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer |FullGRP|=" << FullGRP.size() << "\n";
#endif
  std::vector<Telt> LGenPermPtWiseStab;
  PermutationBuilder<T, Telt> builder(SHV_T);
  std::vector<Telt> LGenGlobStab_perm;
  for (auto & eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder.get_permutation(eGen);
    LGenGlobStab_perm.push_back(ePerm);
  }
  Tgroup GRPsub(LGenGlobStab_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer |GRPsub|=" << GRPsub.size() << "\n";
#endif
  auto try_upgrade=[&]() -> std::optional<Telt> {
    // Not sure if left or right cosets.
    // Left  transversals are g H
    // Right transversals are H g
    std::vector<Telt> ListCos = FullGRP.LeftTransversal_Direct(GRPsub);
    for (auto & eCosReprPerm : ListCos) {
      if (!eCosReprPerm.isIdentity()) {
        std::optional<MyMatrix<T>> opt = is_corr_and_solve(eCosReprPerm, SHV_T, eMat, LinSpa);
        if (opt) {
#ifdef DEBUG_TSPACE_FUNCTIONS
          os << "TSPACE: LINSPA_ComputeStabilizer Finding a new eCosReprPerm\n";
#endif
          return eCosReprPerm;
        }
      }
    }
    return {};
  };
  while(true) {
    std::optional<Telt> opt = try_upgrade();
    if (opt) {
      // Found another stabilizing element, upgrading the group and retry.
      LGenGlobStab_perm.push_back(*opt);
      GRPsub = Tgroup(LGenGlobStab_perm, n_row);
#ifdef DEBUG_TSPACE_FUNCTIONS
      os << "TSPACE: LINSPA_ComputeStabilizer Now |GRPsub|=" << GRPsub.size() << "\n";
#endif
    } else {
      break;
    }
  }
  std::vector<MyMatrix<T>> LGenGlobStab_matr;
  for (auto & eGen : LGenGlobStab_perm) {
    std::optional<MyMatrix<T>> opt = is_corr_and_solve(eGen, SHV_T, eMat, LinSpa);
    MyMatrix<T> eGenMatr = unfold_opt(opt, "Failed to unfold");
    LGenGlobStab_matr.push_back(eGenMatr);
  }
  return LGenGlobStab_matr;
}

/*
  For two positive definite matrices M1 find if it exists a transformation P such that
  * P M1 P^T = M2
  * P LinSpa.ListMat P^T  image is LinSpa.ListMat
*/
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> LINSPA_TestEquivalenceGramMatrix(LinSpaceMatrix<T> const &LinSpa,
                                                               MyMatrix<T> const& eMat1,
                                                               MyMatrix<T> const& eMat2,
                                                               std::ostream & os) {
  using Tfield = typename overlying_field<T>::field_type;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  //  using Tfield = T;
  MyMatrix<Tint> SHV1 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat1);
  MyMatrix<Tint> SHV2 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat2);
  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T,Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T,Tint>(SHV2);
  if (SHV1_T.rows() != SHV2_T.rows()) {
    return {};
  }
  int n_row = SHV1_T.rows();
  std::vector<T> Vdiag1(n_row, 0), Vdiag2(n_row, 0);
  std::vector<MyMatrix<T>> ListMat1 = GetFamilyDiscMatrices(eMat1, LinSpa.ListComm, LinSpa.ListSubspaces);
  std::vector<MyMatrix<T>> ListMat2 = GetFamilyDiscMatrices(eMat2, LinSpa.ListComm, LinSpa.ListSubspaces);
  std::optional<std::vector<Tidx>> opt1 = TestEquivalence_ListMat_Vdiag<T, Tfield, Tidx>(SHV1_T, ListMat1, Vdiag1, SHV2_T, ListMat2, Vdiag2, os);
  if (!opt1) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Exiting here at opt1\n";
#endif
    return {};
  }
  //
  // Building one equivalence that should preserve the basics and which we can work on.
  //
  Telt eltEquiv(*opt1);
  Telt eltInv = Inverse(eltEquiv);
  std::optional<MyMatrix<T>> opt2 = FindTransformationGeneral(SHV2_T, SHV1_T, eltInv);
  MyMatrix<T> OneEquiv = unfold_opt(opt2, "Failed to get transformation");
#ifdef DEBUG_TSPACE_FUNCTIONS
  if (!IsIntegralMatrix(OneEquiv)) {
    std::cerr << "The matrix TransMat should be integral\n";
    throw TerminalException{1};
  }
  MyMatrix<T> eMat1_img = OneEquiv * eMat1 * OneEquiv.transpose();
  if (eMat1_img != eMat2) {
    std::cerr << "The matrix TransMat does not preserve eMat\n";
    throw TerminalException{1};
  }
#endif
  if (is_stab_space(OneEquiv, LinSpa)) {
    MyMatrix<Tint> OneEquiv_tint = UniversalMatrixConversion<Tint,T>(OneEquiv);
    return OneEquiv_tint;
  }
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets for resolving that.
  //
  std::vector<std::vector<Tidx>> ListGen1 =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx>(SHV1_T, ListMat1, Vdiag1, os);
  std::vector<Telt> LGenPerm1;
  for (auto & eList1 : ListGen1) {
    Telt ePerm1(eList1);
    LGenPerm1.push_back(ePerm1);
  }
  Tgroup FullGRP1(LGenPerm1, n_row);
  std::vector<Telt> LGenPermPtWiseStab1;
  PermutationBuilder<T, Telt> builder1(SHV1_T);
  std::vector<Telt> LGenGlobStab1_perm;
  for (auto & eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder1.get_permutation(eGen);
    LGenGlobStab1_perm.push_back(ePerm);
  }
  Tgroup GRPsub1(LGenGlobStab1_perm, n_row);
  struct PartSol {
    std::optional<Telt> new_gen;
    std::optional<MyMatrix<T>> sol;
  };
  auto try_solution=[&]() -> PartSol {
    // Not sure if left or right cosets.
    std::vector<Telt> ListCos = FullGRP1.LeftTransversal_Direct(GRPsub1);
    for (auto & eCosReprPerm : ListCos) {
      if (!eCosReprPerm.isIdentity()) {
        MyMatrix<T> eCosReprMatr = get_mat_from_shv_perm(eCosReprPerm, SHV1_T, eMat1);
        MyMatrix<T> eProd = eCosReprMatr * OneEquiv;
        if (is_stab_space(eProd, LinSpa)) {
          return {{}, eProd};
        }
        if (is_stab_space(eCosReprMatr, LinSpa)) {
          return {eCosReprPerm, {}};
        }
      }
    }
    return {{},{}};
  };
  while(true) {
    PartSol p_sol = try_solution();
    if (p_sol.sol) {
      MyMatrix<T> const& Pmat_T = *p_sol.sol;
      MyMatrix<Tint> Pmat = UniversalMatrixConversion<Tint,T>(Pmat_T);
      return Pmat;
    }
    if (!p_sol.new_gen) {
      return {};
    }
    LGenGlobStab1_perm.push_back(*p_sol.new_gen);
    GRPsub1 = Tgroup(LGenGlobStab1_perm, n_row);
  }
}



/*
  Compute an invariant of the gram matrix
  so that it does not change after arithmetic
  equivalence.
  However, it is not invariant under scaling.
*/
template <typename T, typename Tint>
size_t GetInvariantGramShortest(MyMatrix<T> const &eGram,
                                MyMatrix<Tint> const &SHV,
                                size_t const& seed,
                                [[maybe_unused]] std::ostream & os) {
  T eDet = DeterminantMat(eGram);
  int nbVect = SHV.rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: eDet=" << eDet << " nbVect=" << nbVect << "\n";
#endif
  std::vector<MyVector<T>> ListV;
  for (int iVect=0; iVect<nbVect; iVect++) {
    MyVector<Tint> V = GetMatrixRow(SHV, iVect);
    MyVector<T> V_T = UniversalVectorConversion<T,Tint>(V);
    ListV.push_back(V_T);
  }
  std::map<T,size_t> map_diag, map_off_diag;
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<T> Vprod = eGram * ListV[iVect];
    T eNorm = Vprod.dot(ListV[iVect]);
    map_diag[eNorm] += 1;
    for (int jVect=iVect+1; jVect<nbVect; jVect++) {
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
  for (auto & kv : map_diag) {
    size_t hash_T = std::hash<T>()(kv.first);
    combine_hash(hash_ret, hash_T);
    combine_hash(hash_ret, kv.second);
  }
  for (auto & kv : map_off_diag) {
    size_t hash_T = std::hash<T>()(kv.first);
    combine_hash(hash_ret, hash_T);
    combine_hash(hash_ret, kv.second);
  }
  return hash_ret;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_FUNCTIONS_H_
// clang-format on
