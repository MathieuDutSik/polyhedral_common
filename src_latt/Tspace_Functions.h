// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_FUNCTIONS_H_
#define SRC_LATT_TSPACE_FUNCTIONS_H_

// clang-format off
#include "POLY_PolytopeFct.h"
#include "ShortestUniversal.h"
#include "Temp_PolytopeEquiStab.h"
#include "Temp_Positivity.h"
#include "SHORT_ShortestConfig.h"
#include <set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_TSPACE_GENERAL
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
LinSpaceMatrix<T> BuildLinSpace(MyMatrix<T> const& SuperMat, std::vector<MyMatrix<T>> const& ListMat, std::vector<MyMatrix<T>> const& ListComm) {
  int n = SuperMat.rows();
  std::vector<std::vector<T>> ListLineMat;
  for (auto & eMat : ListMat) {
    std::vector<T> eV = GetLineVector(eMat);
    ListLineMat.push_back(eV);
  }
  std::vector<MyMatrix<T>> ListSubspaces;
  MyMatrix<T> eGen = -IdentityMat<T>(n);
  std::vector<MyMatrix<T>> PtStab{eGen};
  return {n, SuperMat, ListMat, ListLineMat, ListComm, ListSubspaces, PtStab};
}

template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> ComputePointStabilizerTspace(MyMatrix<T> const& SuperMat, std::vector<MyMatrix<T>> const& ListMat, std::ostream & os) {
  using Tidx = uint32_t;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(SuperMat);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  std::vector<T> Vdiag(SHV_T.rows(), 0);
  bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGenPerm =
    GetListGenAutomorphism_ListMat_Vdiag<T, T, Tidx, use_scheme>(SHV_T, ListMat, Vdiag, os);
  std::vector<MyMatrix<Tint>> ListGenMatr;
  for (auto & eList : ListGenPerm) {
    MyMatrix<T> eMatr_T = FindTransformation_vect(SHV_T, SHV_T, eList);
    MyMatrix<Tint> eMatr = UniversalMatrixConversion<Tint,T>(eMatr_T);
    ListGenMatr.push_back(eMatr);
  }
  return ListGenMatr;
}

template<typename T, typename Tint>
MyMatrix<T> GetOnePositiveDefiniteMatrix(std::vector<MyMatrix<T>> const& ListMat) {
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
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized);
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
    ListV.push_back();
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
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row,0);
  std::vector<MyMatrix<T>> ListMat = {GramMat};
  const bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGen =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(SHV_T, ListMat, Vdiag, os);
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
      FindMatrixTransformationTest(SHV_T, SHV_T, eList);
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
    bool test = IsSymmetryGroupCorrect(TheMat, LinSpa, os);
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
  For a positive definite matrix in the T-space, we compute the group
  of transformations that presrves:
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
  using TintGroup = typename Tgroup::Tint;
  using Tfield = T;
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  int n_row = SHV.rows();
  std::vector<T> Vdiag(n_row,0);
  std::vector<MyMatrix<T>> ListMat = {eMat};
  const bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGen =
    GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(SHV_T, ListMat, Vdiag, os);
  //
  // Try the direct strategy and hopes to be lucky
  //
  auto is_corr_and_solve=[&](Telt const& elt) -> std::optional<MyMatrix<T>> {
    std::optional<MyMatrix<T>> opt = FindMatrixTransformationTest(SHV_T, SHV_T, elt);
    MyMatrix<T> TransMat = unfold_opt(opt, "Failed to get transformation");
    for (auto & eMat : LinSpa.ListMat) {
      MyMatrix<T> eMatImg = eGen * eMat * eGen.transpose();
      MyVector<T> eMatImg_V = SymmetricMatrixToVector(eMatImg);
      std::optional<MyVector<T>> opt = SolutionMat(LinSpa.ListMatAsBigMat, eMatImg_V);
      if (!opt) {
        return {};
      }
    }
    return TransMat;
  };
  auto get_generators=[&]() -> std::optional<std::vector<MyMatrix<T>>> {
    std::vector<MyMatrix<T>> ListTransMat;
    for (auto & eGen : ListGen) {
      Telt elt(eGen);
      std::optional<MyMatrix<T>> opt = is_corr_and_solve(elt);
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
  std::vector<Telt> LGenPermPtWiseStab;
  std::vector<MyVector<T>> ListV;
  std::unordered_map<MyVector<T>, size_t> MapV;
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV_T, i_row);
    ListV.push_back(V);
    MapV[V] = i_row;
  }
  std::vector<Telt> LGenGlobStab_perm;
  for (auto & eGen : LinSpa.PtStabGens) {
    std::vector<Tidx> eList(n_row);
    for (int i_row=0; i_row<n_row; i_row++) {
      MyVector<T> Vimg = eGen.transpose() * ListV[i_row];
      int pos = MapV.at(Vimg);
      eList[i_row] = pos;
    }
    Telt ePerm(eList);
    LGenGlobStab_perm.push_back(ePerm);
  }
  Tgroup GRPsub(LGenGlobStab, n_row);
  auto try_upgrade=[&]() -> std::optional<Telt> {
    // Not sure if left or right cosets.
    std::vector<Telt> ListCos = FullGRP.LeftTransversal_Direct(GRPsub);
    for (auto & eGen : ListCos) {
      // Not sure if the is_identity works
      if (!elt.is_identity()) {
        std::optional<MyMatrix<T>> opt = is_corr_and_solve(elt);
        if (opt) {
          return *opt;
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
      GRPsub = Tgroup(LGenGlobStab_perm, n);
    } else {
      break;
    }
  }
  std::vector<MyMatrix<T>> LGenGlobStab_matr;
  for (auto & eGen : LGenGlobStab_perm) {
    std::optional<MyMatrix<T>> opt = is_corr_and_solve();
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
std::optional<MyMatrix<T>> LINSPA_TestEquivalenceGramMatrix(LinSpaceMatrix<T> const &LinSpa,
                                                            MyMatrix<T> const& eMat1,
                                                            MyMatrix<T> const& eMat1,
                                                            std::ostream & os) {
  using Tidx = uint32_t;
  MyMatrix<Tint> SHV1 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat1);
  MyMatrix<Tint> SHV2 = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat2);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
  std::vector<T> Vdiag(SHV_T.rows(), 0);
  bool use_scheme = true;
  std::vector<std::vector<Tidx>> ListGenPerm =
    GetListGenAutomorphism_ListMat_Vdiag<T, T, Tidx, use_scheme>(SHV_T, ListMat, Vdiag, os);


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
#ifdef DEBUG_TSPACE_GENERAL
  os << "eDet=" << eDet << " nbVect=" << nbVect << "\n";
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
