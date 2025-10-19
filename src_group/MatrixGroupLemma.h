// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPLEMMA_H_
#define SRC_GROUP_MATRIXGROUPLEMMA_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "ClassicLLL.h"
#include "Timings.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_LEMMA
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP_LEMMA
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_GROUP_LEMMA
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP_LEMMA
#endif

/*
  This is about The Schreier Lemma, finding pre-image of subgroups.
  The problem is to make it as efficient as possible.
  That is to reduce the number of generators needed.
 */

template<typename T, typename Tgroup>
std::vector<MyMatrix<T>> PreImageSubgroupOneStep(std::vector<MyMatrix<T>> const& ListMatr, std::vector<typename Tgroup::Telt> const& ListPerm, MyMatrix<T> const& id_matr, Tgroup const& eGRP, std::ostream& os) {
#ifdef TIMINGS_MATRIX_GROUP_LEMMA
  MicrosecondTime time;
#endif
#ifdef DEBUG_MATRIX_GROUP_LEMMA
  os << "MATGRPBAS: PreImageSubgroupOneStep, |eGRP|=" << eGRP.size() << " |ListPerm|=" << ListPerm.size() << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatr1 =
    permutalib::PreImageSubgroup<Tgroup, MyMatrix<T>>(ListMatr, ListPerm, id_matr, eGRP);
#ifdef TIMINGS_MATRIX_GROUP_LEMMA
  os << "|MATGRPBAS: PreImageSubgroupOneStep, permutalib::PreImageSubgroup|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_LEMMA
  os << "MATGRPBAS: PreImageSubgroupOneStep, comp(ListMatr1)=" << compute_complexity_listmat(ListMatr1) << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatr2 = ExhaustiveReductionComplexityGroupMatrix<T>(ListMatr1, os);
#ifdef TIMINGS_MATRIX_GROUP_LEMMA
  os << "|MATGRPBAS: PreImageSubgroupOneStep, ExhaustiveReductionComplexityGroupMatrix|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_LEMMA
  os << "MATGRPBAS: PreImageSubgroupOneStep, comp(ListMatr2)=" << compute_complexity_listmat(ListMatr2) << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP_LEMMA_DISABLE
  CheckGroupEquality<T,Tgroup>(ListMatr1, ListMatr2, os);
#endif
  return ListMatr2;
}

template<typename T, typename Tgroup>
std::vector<MyMatrix<T>> PreImageSubgroup(std::vector<MyMatrix<T>> const& ListMatr, std::vector<typename Tgroup::Telt> const& ListPerm, std::function<typename Tgroup::Telt(MyMatrix<T> const&)> f_get_perm, MyMatrix<T> const& id_matr, Tgroup const& eGRP, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Telt id_perm = eGRP.get_identity();
  Tidx len = id_perm.size();
  Tgroup GRPbig(ListPerm, len);
  if (GRPbig.size() == eGRP.size()) {
    return ListMatr;
  }
#ifdef TRACK_INFO_MATRIX_GROUP_LEMMA
  WriteGroupFile("GRPbig", GRPbig);
  WriteGroupFile("GRPsub", eGRP);
  WriteGroupFileGAP("GRPbig_gap", GRPbig);
  WriteGroupFileGAP("GRPsub_gap", eGRP);
#endif
#ifdef DEBUG_MATRIX_GROUP_LEMMA
  os << "MATGRPBAS: PreImageSubgroup, beginning\n";
#endif
  std::vector<Tgroup> l_grp = GRPbig.GetAscendingChainSubgroup(eGRP);
  size_t len_stab = l_grp.size() - 1;
#ifdef DEBUG_MATRIX_GROUP_LEMMA
  os << "MATGRPBAS: PreImageSubgroup, len_stab=" << len_stab << "\n";
  for (size_t iGRP=0; iGRP<=len_stab; iGRP++) {
    os << "MATGRPBAS: PreImageSubgroup, iGRP=" << iGRP << "/" << len_stab << " |eGRP|=" << l_grp[iGRP].size() << "\n";
  }
#endif
  std::vector<MyMatrix<T>> LGenMatr = ListMatr;
  std::vector<Telt> LGenPerm = ListPerm;
  for (size_t u=0; u<len_stab; u++) {
    size_t idx = len_stab - 1 - u;
#ifdef DEBUG_MATRIX_GROUP_LEMMA
    os << "MATGRPBAS: PreImageSubgroup, len_stab=" << len_stab << " u=" << u << " idx=" << idx << "\n";
#endif
    LGenMatr = PreImageSubgroupOneStep<T,Tgroup>(LGenMatr, LGenPerm, id_matr, l_grp[idx], os);
    if (idx > 0) {
      LGenPerm.clear();
      for (auto & eMatr: LGenMatr) {
        Telt ePerm = f_get_perm(eMatr);
        LGenPerm.push_back(ePerm);
      }
    }
  }
  return LGenMatr;
}



// clang-format off
#endif  // SRC_GROUP_MATRIXGROUPLEMMA_H_
// clang-format on
