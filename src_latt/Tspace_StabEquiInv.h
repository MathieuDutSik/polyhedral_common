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
#include "POLY_Fundamental.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include "Positivity.h"
#include "Tspace_Functions.h"
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

// clang-format off
#endif  // SRC_LATT_TSPACE_STABEQUIINV_H_
// clang-format on
