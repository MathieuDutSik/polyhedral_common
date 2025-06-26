// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPBASIC_H_
#define SRC_GROUP_MATRIXGROUPBASIC_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixMod.h"
#include "ClassicLLL.h"
#include "Timings.h"
#include "TestGroup.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <map>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_BASIC
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP_BASIC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_GROUP_BASIC
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP_BASIC
#endif

template <typename T>
void write_matrix_group(std::vector<MyMatrix<T>> const& list_mat, std::string const& context) {
  if (list_mat.size() == 0) {
    return;
  }
  int dim = list_mat[0].rows();
  std::string Prefix = "MatrixGroup_" + context + "_dim" + std::to_string(dim) + "_idx";
  std::string FileOut = FindAvailableFileFromPrefix(Prefix);
  WriteListMatrixFile(FileOut, list_mat);
}

template<typename T>
std::string compute_complexity_matrix(MyMatrix<T> const& mat) {
  int n = mat.rows();
  T ell1(0);
  T ellinfinity(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = mat(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
      if (abs_val > ellinfinity) {
        ellinfinity = abs_val;
      }
    }
  }
  return "(ell1=" + std::to_string(ell1) + ", ellinf=" + std::to_string(ellinfinity) + ")";
}

template<typename T>
std::string compute_complexity_listmat(std::vector<MyMatrix<T>> const& list_mat) {
  if (list_mat.size() == 0) {
    return "zero generators";
  }
  int n = list_mat[0].rows();
  size_t n_mat = list_mat.size();
  T ell1_global(0);
  T ellinfinite_global(0);
  for (auto & e_mat: list_mat) {
    T ell1(0);
    T ellinfinity(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        T val = e_mat(i,j);
        T abs_val = T_abs(val);
        ell1 += abs_val;
        if (abs_val > ellinfinity) {
          ellinfinity = abs_val;
        }
      }
    }
    ell1_global += ell1;
    if (ellinfinity > ellinfinite_global) {
      ellinfinite_global = ellinfinity;
    }
  }
  return "(n_gen=" + std::to_string(n_mat) + ", ell1_global=" + std::to_string(ell1_global) + ", ellinfinity=" + std::to_string(ellinfinite_global) + ")";
}


template<bool always_equal>
std::string compute_complexity_listseq(std::vector<permutalib::SequenceType<always_equal>> const& list_seq) {
  size_t n_seq = list_seq.size();
  size_t ell1_global = 0;
  size_t ellinfinite_global = 0;
  for (auto & seq: list_seq) {
    std::vector<int64_t> const& ListIdx = seq.getVect();
    size_t len = ListIdx.size();
    ell1_global += len;
    if (ellinfinite_global < len) {
      ellinfinite_global = len;
    }
  }
  return "(n_seq=" + std::to_string(n_seq) + ", ell1_global=" + std::to_string(ell1_global) + ", ellinfinity=" + std::to_string(ellinfinite_global) + ")";
}




template <typename T>
size_t GetRationalInvariant(std::vector<MyMatrix<T>> const &ListGen) {
  std::set<T> set_den;
  for (auto &eGen : ListGen) {
    if (!IsIntegralMatrix(eGen)) {
      int n = eGen.rows();
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          T eDen = GetDenominator(eGen(i, j));
          set_den.insert(eDen);
        }
      }
    }
  }
  std::set<T> primes;
  for (auto &eDen : set_den) {
    std::map<T, size_t> l_primes = FactorsIntMap(eDen);
    for (auto &kv : l_primes) {
      primes.insert(kv.first);
    }
  }
  T prod(1);
  for (auto &p : primes) {
    prod *= p;
  }
  return std::hash<T>()(prod);
}






template <typename T, typename Thelper>
T L1normMatrixGroup(Thelper const &helper,
                    std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  T sum(0);
  for (auto &eMat : ListMatr) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(eMat(i, j));
  }
  return sum;
}

template <typename T, typename Tint, typename Thelper>
std::pair<std::vector<MyMatrix<T>>, MyMatrix<Tint>>
LLLMatrixGroupReduction(Thelper const &helper,
                        std::vector<MyMatrix<T>> const &ListMatr) {
  int n = helper.n;
  MyMatrix<T> PosDefMat = IdentityMat<T>(n);
  for (auto &eMat : ListMatr) {
    if (!IsIdentity(eMat)) {
      MyMatrix<T> eProd = eMat * eMat.transpose();
      PosDefMat += eProd;
    }
  }
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(PosDefMat);
  MyMatrix<Tint> const &Pmat = pair.Pmat;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> PmatInv_T = Inverse(Pmat_T);
  std::vector<MyMatrix<T>> ListMatrNew;
  for (auto &eMat : ListMatr) {
    MyMatrix<T> eMatNew = Pmat_T * eMat * PmatInv_T;
    ListMatrNew.emplace_back(std::move(eMatNew));
  }
  return {std::move(ListMatrNew), Pmat};
}

template <typename T> T LinearSpace_GetDivisor(MyMatrix<T> const &TheSpace) {
  T TheDet = T_abs(DeterminantMat(TheSpace));
  T eDiv(1);
  int n = TheSpace.rows();
  RecSolutionIntMat<T> eCan(TheSpace);
  while (true) {
    MyMatrix<T> M = eDiv * IdentityMat<T>(n);
    bool test = eCan.is_containing_m(M);
    if (test) {
      return eDiv;
    }
#ifdef SANITY_CHECK_MATRIX_GROUP_BASIC
    if (eDiv > TheDet) {
      std::cerr << "eDiv=" << eDiv << " TheDet=" << TheDet << "\n";
      std::cerr << "TheSpace=\n";
      WriteMatrix(std::cerr, TheSpace);
      std::cerr << "Clear error in LinearSpace_GetDivisor\n";
      throw TerminalException{1};
    }
#endif
    eDiv += 1;
  }
}

template <typename T>
MyMatrix<T>
MatrixIntegral_GetInvariantSpace(int const &n,
                                 std::vector<MyMatrix<T>> const &LGen,
                                 [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrix<T>> LGenTot;
  for (auto &eGen : LGen) {
    LGenTot.push_back(eGen);
    LGenTot.push_back(Inverse(eGen));
  }
  LGenTot.push_back(IdentityMat<T>(n));
  MyMatrix<T> TheSpace = IdentityMat<T>(n);
  T TheDet(1);
#ifdef DEBUG_MATRIX_GROUP_BASIC
  size_t iter = 0;
#endif
  while (true) {
    std::vector<MyVector<T>> ConcatSpace;
    for (auto &eGen : LGenTot) {
      MyMatrix<T> TheSpaceImg = TheSpace * eGen;
      for (int i = 0; i < n; i++) {
        ConcatSpace.push_back(GetMatrixRow(TheSpaceImg, i));
      }
    }
    MyMatrix<T> NewSpace1 = GetZbasis(MatrixFromVectorFamily(ConcatSpace));
    // The LLL reduction appears quite efficient
    MyMatrix<T> NewSpace = SublatticeBasisReduction(NewSpace1);
    T NewDet = T_abs(DeterminantMat(NewSpace));
    if (NewDet == TheDet) {
#ifdef DEBUG_MATRIX_GROUP_BASIC
      os << "MAT_GRP: MatrixIntegral_GetInvariantSpace, NewSpace=\n";
      WriteMatrix(os, NewSpace);
      os << "MAT_GRP: MatrixIntegral_GetInvariantSpace, returning after n_iter=" << iter << " TheDet=" << TheDet << "\n";
#endif
      return NewSpace;
    }
    TheSpace = NewSpace;
    TheDet = NewDet;
#ifdef DEBUG_MATRIX_GROUP_BASIC
    os << "MAT_GRP: MatrixIntegral_GetInvariantSpace, iter=" << iter << " TheDet=" << TheDet << "\n";
    iter += 1;
#endif
  }
}



template<typename T, typename Tgroup>
std::vector<MyMatrix<T>> PreImageSubgroupOneStep(std::vector<MyMatrix<T>> const& ListMatr, std::vector<typename Tgroup::Telt> const& ListPerm, MyMatrix<T> const& id_matr, Tgroup const& eGRP, std::ostream& os) {
  using Tseq = permutalib::SequenceType<false>;
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  MicrosecondTime time;
#endif
  std::vector<Tseq> ListSeq;
  std::vector<MyMatrix<T>> ListMatrInv;
  for (size_t i_elt=0; i_elt<ListPerm.size(); i_elt++) {
    std::vector<int64_t> ListIdx{int64_t(i_elt) + 1};
    Tseq seq(ListIdx);
    ListSeq.push_back(seq);
    MyMatrix<T> eMatrInv = Inverse(ListMatr[i_elt]);
    ListMatrInv.push_back(eMatrInv);
  }
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MAT_GRP: PreImageSubgroupOneStep, ListSeq / ListMatrInv|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroupOneStep, |eGRP|=" << eGRP.size() << " |ListPerm|=" << ListPerm.size() << "\n";
#endif
  Tseq id_seq;
  std::vector<Tseq> ListSeq_sub =
    permutalib::PreImageSubgroup<Tgroup, Tseq>(ListSeq, ListPerm, id_seq, eGRP);
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MAT_GRP: PreImageSubgroupOneStep, permutalib::PreImageSubgroup|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroupOneStep, comp(ListSeq_sub)=" << compute_complexity_listseq(ListSeq_sub) << "\n";
#endif
  std::vector<Tseq> ListSeq_sub_red = ExhaustiveReductionComplexitySequences(ListSeq_sub, os);
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MAT_GRP: PreImageSubgroupOneStep, ExhaustiveReductionComplexitySequences|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroupOneStep, comp(ListSeq_sub_red)=" << compute_complexity_listseq(ListSeq_sub_red) << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatr_sub;
  for (auto & eSeq : ListSeq_sub_red) {
    MyMatrix<T> eMatr = id_matr;
    for (auto & eVal : eSeq.getVect()) {
      if (eVal > 0) {
        size_t pos = eVal - 1;
        eMatr *= ListMatr[pos];
      } else {
        size_t fVal = - eVal;
        size_t pos = fVal - 1;
        eMatr *= ListMatrInv[pos];
      }
    }
    ListMatr_sub.push_back(eMatr);
  }
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MAT_GRP: PreImageSubgroupOneStep, ListMatr_sub|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroupOneStep, comp(ListMatr_sub)=" << compute_complexity_listmat(ListMatr_sub) << "\n";
#endif
#ifdef TRACK_INFO_MATRIX_GROUP_BASIC
  write_matrix_group(ListGen1, "PreImageSubgroupOneStep");
#endif
  std::vector<MyMatrix<T>> ListMatr_ret = ExhaustiveReductionComplexityGroupMatrix<T>(ListMatr_sub, os);
#ifdef TIMINGS_MATRIX_GROUP_BASIC
  os << "|MAT_GRP: PreImageSubgroupOneStep, ExhaustiveReductionComplexityGroupMatrix|=" << time << "\n";
#endif
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroupOneStep, comp(ListMatr_ret)=" << compute_complexity_listmat(ListMatr_ret) << "\n";
#endif
#ifdef SANITY_CHECK_MATRIX_GROUP_BASIC
  CheckGroupEquality<T,Tgroup>(ListMatr_ret, ListMatr_sub, os);
#endif
  return ListMatr_ret;
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
  std::vector<Tgroup> l_grp = GRPbig.GetAscendingChainSubgroup(eGRP);
  size_t len_stab = l_grp.size() - 1;
#ifdef DEBUG_MATRIX_GROUP_BASIC
  os << "MAT_GRP: PreImageSubgroup, len_stab=" << len_stab << "\n";
  for (size_t iGRP=0; iGRP<=len_stab; iGRP++) {
    os << "MAT_GRP: PreImageSubgroup, iGRP=" << iGRP << "/" << len_stab << " |eGRP|=" << l_grp[iGRP].size() << "\n";
  }
#endif
  std::vector<MyMatrix<T>> LGenMatr = ListMatr;
  std::vector<Telt> LGenPerm = ListPerm;
  for (size_t u=0; u<len_stab; u++) {
    size_t idx = len_stab - 1 - u;
#ifdef DEBUG_MATRIX_GROUP_BASIC
    os << "MAT_GRP: PreImageSubgroup, len_stab=" << len_stab << " u=" << u << " idx=" << idx << "\n";
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
#endif  // SRC_GROUP_MATRIXGROUPBASIC_H_
// clang-format on
