// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_STABEQUIINV_H_
#define SRC_LATT_TSPACE_STABEQUIINV_H_

/*
  Functions for computing stabilizer, equivalence and invariant
  of matrices in a T-space.
  ---
  Unfortunately, there is no canonical form in full generality.
 */

// clang-format off
#include "Tspace_Fundamental.h"
#include "POLY_Fundamental.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include "Positivity.h"
#include "subgroup_algorithms.h"
#include <set>
#include <vector>
#include <unordered_map>
#include <map>
// clang-format on

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
  The member of the family attached to a preserved subspace. The projector
  inverts a basis, so it is rational even for an integral Gram matrix. Over a
  ring it is therefore computed on the overlying field and the primitive
  integral matrix proportional to the result is returned.

  That rescaling is harmless for the three uses of the family. For the
  invariant and the stabilizer, a permutation preserves the rescaled member
  exactly when it preserves the original one. For the equivalence, the first
  member of the family is the Gram matrix itself, which is carried exactly:
  a transformation matching the families maps eG1 to eG2, hence maps the
  projector of one to the projector of the other, so the two rescalings
  agree and nothing can be matched that was not matched before.
 */
template <typename T>
MyMatrix<T> GetSubspaceDiscMatrix(MyMatrix<T> const &eG,
                                  MyMatrix<T> const &subspace) {
  if constexpr (is_ring_field<T>::value) {
    MyMatrix<T> ProjMat = GetOrthogonalProjectorMatrix(eG, subspace);
    return ProjMat * eG;
  } else {
    using Tfield = typename overlying_field<T>::field_type;
    MyMatrix<Tfield> eG_F = UniversalMatrixConversion<Tfield, T>(eG);
    MyMatrix<Tfield> subspace_F = UniversalMatrixConversion<Tfield, T>(subspace);
    MyMatrix<Tfield> ProjMat = GetOrthogonalProjectorMatrix(eG_F, subspace_F);
    MyMatrix<Tfield> prod = ProjMat * eG_F;
    return RemoveFractionMatrixPlusCoeffRing(prod).TheMat;
  }
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
    ListDisc.emplace_back(GetSubspaceDiscMatrix(eG, subspace));
  }
  return ListDisc;
}

template <typename T>
bool is_stab_space(MyMatrix<T> const &Pmat, LinSpaceMatrix<T> const &LinSpa) {
  // The membership of the image in the space is a rational question even for
  // an integral space, so the solution comes back over the overlying field
  // (which is T itself when T is already a field).
  using Tfield = typename overlying_field<T>::field_type;
  for (auto &eMatSp : LinSpa.ListMat) {
    MyMatrix<T> eMatSpImg = Pmat * eMatSp * Pmat.transpose();
    MyVector<T> eMatSpImg_V = SymmetricMatrixToVector(eMatSpImg);
    std::optional<MyVector<Tfield>> opt =
	SolutionMat(LinSpa.ListMatAsBigMat, eMatSpImg_V);
    if (!opt) {
      return false;
    }
  }
  return true;
}


/*
  The matrix realizing a permutation of the vector family. With a family
  spanning the lattice over Z the matrix is integral and the solve always
  succeeds over the ring; with a full-rank family the matrix can be
  rational, in which case the ring solve reports it as non-representable
  and the empty optional is returned. The caller decides whether that is a
  rejected candidate or a reason to go to the integral machinery.
 */
template <typename T, typename Telt>
std::optional<MyMatrix<T>>
get_mat_from_shv_perm(Telt const &elt, MyMatrix<T> const &SHV_T,
                      [[maybe_unused]] MyMatrix<T> const &eMat) {
  std::optional<MyMatrix<T>> opt = FindTransformationGeneral(SHV_T, SHV_T, elt);
  if (!opt) {
    return {};
  }
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
  MyMatrix<T> const &Pmat = *opt;
  MyMatrix<T> eMatImg = Pmat * eMat * Pmat.transpose();
  if (eMatImg != eMat) {
    std::cerr << "TSPACE: The matrix TransMat does not preserve eMat\n";
    throw TerminalException{1};
  }
#endif
  return opt;
}

template <typename T, typename Telt>
std::optional<MyMatrix<T>>
is_corr_and_solve(Telt const &elt, MyMatrix<T> const &SHV_T,
                  MyMatrix<T> const &eMat, LinSpaceMatrix<T> const &LinSpa) {
  std::optional<MyMatrix<T>> opt = get_mat_from_shv_perm(elt, SHV_T, eMat);
  if (!opt) {
    return {};
  }
  if (!is_stab_space(*opt, LinSpa)) {
    return {};
  }
  return opt;
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
      if (!MapV.contains(Vimg)) {
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
/*
  The f_extra predicate is an extra adequacy condition on the solved
  matrices, on top of the stabilization of the T-space. It has to be
  subgroup-defining: satisfied by the identity and closed under product
  and inverse, since generator adequacy is taken to imply group adequacy.
  The trivial instantiation is the lattice case; the periodic point sets
  add the preservation of their cosets.
 */
template <typename T, typename Tgroup, typename Fextra>
Result_ComputeStabilizer_SHV<T, Tgroup>
LINSPA_ComputeStabilizer_SHV_Kernel(LinSpaceMatrix<T> const &LinSpa,
                             MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
                             std::optional<MyMatrix<T>> const &CommonGramMat,
                             Fextra const &f_extra, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n_row = SHV_T.rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer_SHV n_row=" << n_row << "\n";
#endif
#ifdef TIMINGS_TSPACE_FUNCTIONS
  MicrosecondTime time_st;
#endif
  std::vector<T> Vdiag(n_row, T(0));
  std::vector<MyMatrix<T>> ListMat =
      GetFamilyDiscMatrices(eMat, LinSpa.ListComm, LinSpa.ListSubspaces);
  // When a common Gram matrix is imposed on all the domains, the stabilizer
  // must preserve it. Adding it to the list restricts the automorphisms to
  // the subgroup that also fixes CommonGramMat.
  if (CommonGramMat) {
    ListMat.push_back(*CommonGramMat);
  }

#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Stab, disc_matrices|=" << time_st << "\n";
#endif
  //
  // The generators of the integral symmetry group of the configuration.
  // The direct strategy: when all of them satisfy the T-space and f_extra
  // conditions, the integral group is the answer, since the conditions are
  // subgroup-defining.
  //
  std::vector<MyMatrix<T>> LGenInt = GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(
      SHV_T, ListMat, Vdiag, os);
  bool all_gens_ok = true;
  for (auto &eGen : LGenInt) {
    if (!is_stab_space(eGen, LinSpa) || !f_extra(eGen)) {
      all_gens_ok = false;
      break;
    }
  }
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Stab, int_automorphism ok=" << all_gens_ok
     << "|=" << time_st << "\n";
#endif
  if (all_gens_ok) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: LINSPA_ComputeStabilizer_SHV success of the direct "
          "approach\n";
#endif
    return get_from_gens<T, Tgroup>(LGenInt);
  }
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets for
  // resolving that. Everything in the search is integral: the big group is
  // the integral symmetry group and the pt-wise stabilizer is inside it.
  //
  PermutationBuilder<T, Telt> builder(SHV_T);
  std::vector<Telt> LGenPerm_big;
  for (auto &eGen : LGenInt) {
    LGenPerm_big.push_back(builder.get_permutation(eGen, os));
  }
  std::vector<Telt> LGenPerm_sma;
  for (auto &eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder.get_permutation(eGen, os);
    LGenPerm_sma.emplace_back(std::move(ePerm));
  }
  auto f_correct=[&](Telt const& x) -> bool {
    std::optional<MyMatrix<T>> opt =
      is_corr_and_solve(x, SHV_T, eMat, LinSpa);
    return opt.has_value() && f_extra(*opt);
  };

  std::pair<std::vector<Telt>, Tgroup> pair = get_intermediate_group<Tgroup,decltype(f_correct)>(n_row,
                                                                                                 LGenPerm_sma,
                                                                                                 LGenPerm_big,
                                                                                                 f_correct, os);
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Stab, coset_search|=" << time_st << "\n";
#endif
  return get_from_perms_and_group<T, Tgroup>(pair);
}

template <typename T, typename Tgroup>
Result_ComputeStabilizer_SHV<T, Tgroup>
LINSPA_ComputeStabilizer_SHV(LinSpaceMatrix<T> const &LinSpa,
                             MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
                             std::optional<MyMatrix<T>> const &CommonGramMat,
                             std::ostream &os) {
  auto f_extra=[]([[maybe_unused]] MyMatrix<T> const& x) -> bool {
    return true;
  };
  return LINSPA_ComputeStabilizer_SHV_Kernel<T, Tgroup, decltype(f_extra)>(
      LinSpa, eMat, SHV_T, CommonGramMat, f_extra, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
LINSPA_ComputeStabilizer(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &eMat,
                         std::optional<MyMatrix<T>> const &CommonGramMat,
                         std::ostream &os) {
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<T, Tint>(eMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  Result_ComputeStabilizer_SHV<T, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<T, Tgroup>(LinSpa, eMat, SHV_T,
                                              CommonGramMat, os);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: LINSPA_ComputeStabilizer, we have result\n";
#endif
  return result.get_list_matrix(SHV_T, eMat, LinSpa, os);
}

template <typename T>
size_t LINSPA_Invariant_SHV(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                            MyMatrix<T> const &eMat, MyMatrix<T> const &SHV_T,
                            std::optional<MyMatrix<T>> const &CommonGramMat,
                            std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, T(0));
  std::vector<MyMatrix<T>> ListMat =
      GetFamilyDiscMatrices(eMat, LinSpa.ListComm, LinSpa.ListSubspaces);
  // When a common Gram matrix is imposed on all the domains, the relevant
  // equivalences must preserve it. Adding it to the list of matrices makes
  // the invariant compatible with that restricted equivalence.
  if (CommonGramMat) {
    ListMat.push_back(*CommonGramMat);
  }
  return GetInvariant_ListMat_Vdiag<T, Tfield>(seed, SHV_T, ListMat, Vdiag, os);
}

template <typename T, typename Tint>
size_t LINSPA_Invariant(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &eMat, std::ostream &os) {
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyFullRank<T, Tint>(eMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  std::optional<MyMatrix<T>> CommonGramMat;
  return LINSPA_Invariant_SHV<T>(seed, LinSpa, eMat, SHV_T, CommonGramMat, os);
}


/*
  For two positive definite matrices M1 find if it exists a transformation P
  such that
  * P M1 P^T = M2
  * P LinSpa.ListMat P^T  image is LinSpa.ListMat
*/
/*
  Same f_extra contract as LINSPA_ComputeStabilizer_SHV_Kernel: an extra
  subgroup-defining adequacy condition on the matrices. An equivalence is
  returned only if it satisfies it; the search over the automorphisms of
  eMat1 continues past the candidates that do not.
 */
template <typename T, typename Tgroup, typename Fextra>
std::optional<MyMatrix<T>> LINSPA_TestEquivalenceGramMatrix_SHV_Kernel(
    LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &eMat1,
    MyMatrix<T> const &eMat2, MyMatrix<T> const &SHV1_T,
    MyMatrix<T> const &SHV2_T,
    std::optional<MyMatrix<T>> const &CommonGramMat, Fextra const &f_extra,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
  int nbCol = SHV1_T.cols();
  int rnk1 = RankMat(SHV1_T);
  int rnk2 = RankMat(SHV2_T);
  if (nbCol != rnk1 || nbCol != rnk2) {
    std::cerr << "TSPACE: Equiv, rnk1=" << rnk1 << " rnk2=" << rnk2 << " nbCol=" << nbCol << "\n";
    std::cerr << "TSPACE: SHV1_T and SHV2_T should be of full rank\n";
    throw TerminalException{1};
  }
#endif
  if (SHV1_T.rows() != SHV2_T.rows()) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, Exiting here at |SHV1| <> |SHV2|\n";
#endif
    return {};
  }
  int n_row = SHV1_T.rows();
#ifdef TIMINGS_TSPACE_FUNCTIONS
  MicrosecondTime time_eq;
#endif
  std::vector<T> Vdiag1(n_row, T(0)), Vdiag2(n_row, T(0));
  std::vector<MyMatrix<T>> ListMat1 =
      GetFamilyDiscMatrices(eMat1, LinSpa.ListComm, LinSpa.ListSubspaces);
  std::vector<MyMatrix<T>> ListMat2 =
      GetFamilyDiscMatrices(eMat2, LinSpa.ListComm, LinSpa.ListSubspaces);
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Equiv, disc_matrices|=" << time_eq << "\n";
#endif
  // When a common Gram matrix is imposed on all the domains, the equivalence
  // has to preserve it. Adding it identically to both lists forces the
  // matching permutation (and thus the recovered transformation) to satisfy
  // P * CommonGramMat * P^T = CommonGramMat.
  if (CommonGramMat) {
    ListMat1.push_back(*CommonGramMat);
    ListMat2.push_back(*CommonGramMat);
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Equiv, |ListComm|=" << LinSpa.ListComm.size()
     << " |ListSubspaces|=" << LinSpa.ListSubspaces.size() << "\n";
  os << "TSPACE: Equiv, |ListMat1|=" << ListMat1.size()
     << " |ListMat2|=" << ListMat2.size() << "\n";
#endif
  std::optional<MyMatrix<T>> optEquivInt =
      TestIntEquivalence_ListMat_Vdiag<T, Tgroup>(
          SHV1_T, ListMat1, Vdiag1, SHV2_T, ListMat2, Vdiag2, os);
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Equiv, int_equivalence n_row=" << n_row
     << " n_mat=" << ListMat1.size() << " dim=" << SHV1_T.cols()
     << " found=" << optEquivInt.has_value() << "|=" << time_eq << "\n";
#endif
  if (!optEquivInt) {
    return {};
  }
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
  {
    MyMatrix<T> const &EquivInt = *optEquivInt;
    MyMatrix<T> eMat1_img = EquivInt * eMat1 * EquivInt.transpose();
    if (eMat1_img != eMat2) {
      std::cerr << "TSPACE: Equiv, the matrix does not map eMat1 to eMat2\n";
      throw TerminalException{1};
    }
    if (CommonGramMat) {
      MyMatrix<T> CommonImg = EquivInt * (*CommonGramMat) * EquivInt.transpose();
      if (CommonImg != *CommonGramMat) {
        std::cerr << "TSPACE: Equiv, the matrix does not preserve "
                     "CommonGramMat\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  auto f_is_ok = [&](MyMatrix<T> const &x) -> bool {
    return is_stab_space(x, LinSpa) && f_extra(x);
  };
  bool direct_ok = f_is_ok(*optEquivInt);
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Equiv, is_stab_space ok=" << direct_ok << "|=" << time_eq
     << "\n";
#endif
  if (direct_ok) {
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPACE: Equiv, Direct approach success, no need to go further\n";
#endif
    return *optEquivInt;
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Equiv, Direct approach failure, computing stabilizer and "
        "iterating\n";
#endif
  //
  // The direct approach failed, let us use the pt-wise-stab and the cosets
  // for resolving that. The big group is the integral symmetry group of the
  // first configuration: composing the integral equivalence with an integral
  // symmetry enumerates exactly the integral equivalences.
  //
  std::vector<MyMatrix<T>> LGenInt = GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(
      SHV1_T, ListMat1, Vdiag1, os);
  PermutationBuilder<T, Telt> builder1(SHV1_T);
  std::vector<Telt> LGenPerm_big;
  for (auto &eGen : LGenInt) {
    LGenPerm_big.push_back(builder1.get_permutation(eGen, os));
  }
  // The group elements that preserve the whole T-space, will preserve
  // the matrix considered here. And so it is a subgroup of the full group.
  std::vector<Telt> LGenPerm_sma;
  for (auto &eGen : LinSpa.PtStabGens) {
    Telt ePerm = builder1.get_permutation(eGen, os);
    LGenPerm_sma.push_back(ePerm);
  }
  auto f_get_out=[&](Telt const& x) -> MyMatrix<T> {
    std::optional<MyMatrix<T>> opt_mat =
        get_mat_from_shv_perm(x, SHV1_T, eMat1);
    return unfold_opt(opt_mat, "the coset representative should be integral");
  };
  std::optional<MyMatrix<T>> result = get_intermediate_equivalence<MyMatrix<T>,Tgroup,decltype(f_get_out),decltype(f_is_ok)>(n_row,
                                                                                                                             LGenPerm_sma,
                                                                                                                             LGenPerm_big,
                                                                                                                             *optEquivInt,
                                                                                                                             f_get_out,
                                                                                                                             f_is_ok,
                                                                                                                             os);
#ifdef TIMINGS_TSPACE_FUNCTIONS
  os << "|TSPACE: Equiv, coset_search|=" << time_eq << "\n";
#endif
#ifdef SANITY_CHECK_EXTENSIVE_TSPACE_FUNCTIONS
  if (result) {
    MyMatrix<T> const& eMat_T = *result;
    if (!is_stab_space(eMat_T, LinSpa)) {
      std::cerr << "TSPACE: The matrix is not stabilizing the space\n";
      throw TerminalException{1};
    }
    MyMatrix<T> eMat1_img = eMat_T * eMat1 * eMat_T.transpose();
    if (eMat1_img != eMat2) {
      std::cerr << "TSPACE: Equiv, we do not have an equivalence\n";
      throw TerminalException{1};
    }
    if (CommonGramMat) {
      MyMatrix<T> CommonImg = eMat_T * (*CommonGramMat) * eMat_T.transpose();
      if (CommonImg != *CommonGramMat) {
        std::cerr << "TSPACE: Equiv, equivalence does not preserve "
                     "CommonGramMat\n";
        throw TerminalException{1};
      }
    }
  } else {
    // Exhaustive check over the integral group anchored at the integral
    // equivalence: those products are exactly the integral equivalences.
    Tgroup FullGRP1(LGenPerm_big, n_row);
    for (auto &elt : FullGRP1) {
      std::optional<MyMatrix<T>> opt_mat =
          get_mat_from_shv_perm(elt, SHV1_T, eMat1);
      MyMatrix<T> eMatr =
          unfold_opt(opt_mat, "the integral group element should solve");
      MyMatrix<T> eProd_T = *optEquivInt * eMatr;
      if (f_is_ok(eProd_T)) {
        std::cerr << "TSPACE: We found an equivalence when we do not expect any\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return result;
}

template <typename T, typename Tgroup>
std::optional<MyMatrix<T>> LINSPA_TestEquivalenceGramMatrix_SHV(
    LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &eMat1,
    MyMatrix<T> const &eMat2, MyMatrix<T> const &SHV1_T,
    MyMatrix<T> const &SHV2_T,
    std::optional<MyMatrix<T>> const &CommonGramMat, std::ostream &os) {
  auto f_extra=[]([[maybe_unused]] MyMatrix<T> const& x) -> bool {
    return true;
  };
  return LINSPA_TestEquivalenceGramMatrix_SHV_Kernel<T, Tgroup,
                                                     decltype(f_extra)>(
      LinSpa, eMat1, eMat2, SHV1_T, SHV2_T, CommonGramMat, f_extra, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LINSPA_TestEquivalenceGramMatrix(LinSpaceMatrix<T> const &LinSpa,
                                 MyMatrix<T> const &eMat1,
                                 MyMatrix<T> const &eMat2, std::ostream &os) {
  MyMatrix<Tint> SHV1 = ExtractInvariantVectorFamilyFullRank<T, Tint>(eMat1, os);
  MyMatrix<Tint> SHV2 = ExtractInvariantVectorFamilyFullRank<T, Tint>(eMat2, os);
  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(SHV2);
  std::optional<MyMatrix<T>> CommonGramMat;
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(
          LinSpa, eMat1, eMat2, SHV1_T, SHV2_T, CommonGramMat, os);
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &M_T = *opt;
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
  return M;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_STABEQUIINV_H_
// clang-format on
