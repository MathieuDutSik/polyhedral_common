// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_POLYTOPEEQUISTAB_H_
#define SRC_GROUP_POLYTOPEEQUISTAB_H_

/*
  The code for working with polytope, canonical ordering of
  their vertices stabilizer.

  This is based on WeightMatrix and WeightMatrixSpecified
  codes.

  Other specific interest:
  ---Computing the integral stabilizers and computing the
  integral equivalence.
  ---Computing the integral subgroup rightcoset computations.
  ---For centrally symmetric polytope, implement the
  absolute value trick.

 */

// clang-format off
#include "Timings.h"
#include "WeightMatrix.h"
#include "WeightMatrixLimited.h"
#include "WeightMatrixSpecified.h"
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYTOPE_EQUI_STAB
#define DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
#define DEBUG_POLYTOPE_EQUI_STAB_TRACK_ISOMORPHISM
#endif

#ifdef TIMINGS
#define TIMINGS_POLYTOPE_EQUI_STAB
#define TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_THRESHOLD_SCHEME
#define SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_STAB
#define SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_CANONIC
#endif


static const size_t THRESHOLD_USE_SUBSET_SCHEME_STAB = 1000;
static const size_t THRESHOLD_USE_SUBSET_SCHEME_CANONIC = 1000;

static const size_t THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC = 2;

//
// Equivalence of subsets and stabilizer of a WeightMatrix
//

template <typename T, typename Telt, typename Tidx_value>
std::optional<Telt>
TestEquivalenceSubset(WeightMatrix<true, T, Tidx_value> const &WMat,
                      Face const &f1, Face const &f2, std::ostream &os) {
  using Tidx = typename Telt::Tidx;
  size_t siz = WMat.GetWeightSize();
  size_t n = WMat.rows();
  auto g = [&](Face const &f, size_t iRow, size_t iCol) -> int {
    if (iRow < n && iCol < n)
      return WMat.GetValue(iRow, iCol);
    if (iRow == n && iCol == n)
      return siz + 2;
    if (iRow == n) {
      // Thus iCol < n.
      if (f[iCol] == 0)
        return siz;
      else
        return siz + 1;
    }
    // Last case: Necessarily we have iCol == n && iRow < n
    if (f[iRow] == 0)
      return siz;
    else
      return siz + 1;
  };
  WeightMatrix<true, int, Tidx_value> WMat1(
      n + 1, [&](size_t iRow, size_t iCol) -> int { return g(f1, iRow, iCol); },
      os);
  WeightMatrix<true, int, Tidx_value> WMat2(
      n + 1, [&](size_t iRow, size_t iCol) -> int { return g(f2, iRow, iCol); },
      os);
  std::optional<Telt> test =
      TestEquivalenceWeightMatrix_norenorm_perm<int, Telt>(WMat1, WMat2, os);
  if (!test)
    return {};
  Telt const& eElt = *test;
  std::vector<Tidx> eList(n);
  for (size_t i = 0; i < n; i++) {
    Tidx eVal = OnPoints(i, eElt);
    eList[i] = eVal;
  }
  Telt fElt(std::move(eList));
  return fElt;
}

template <typename T, typename Tgroup, typename Tidx_value>
Tgroup StabilizerSubset(WeightMatrix<true, T, Tidx_value> const &WMat,
                        Face const &f, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphListAdj;
  size_t siz = WMat.GetWeightSize();
  size_t n = WMat.rows();
  auto g = [&](size_t iRow, size_t iCol) -> int {
    if (iRow < n && iCol < n)
      return WMat.GetValue(iRow, iCol);
    if (iRow == n && iCol == n)
      return siz + 2;
    if (iRow == n) {
      if (f[iCol] == 0)
        return siz;
      else
        return siz + 1;
    }
    // Last case: Necessarily we have iCol == n && iRow < n
    if (f[iRow] == 0)
      return siz;
    else
      return siz + 1;
  };
  WeightMatrix<true, int, Tidx_value> WMatW(n + 1, g, os);
  Tgroup GRP = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMatW, os);
  std::vector<Telt> ListPerm;
  for (auto &ePerm : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(n);
    for (size_t i = 0; i < n; i++)
      eList[i] = OnPoints(i, ePerm);
    ListPerm.push_back(Telt(eList));
  }
  return Tgroup(ListPerm, n);
}

//
// The basic of EXT functionality
//

template <typename T> MyMatrix<T> Kernel_GetQmatrix(MyMatrix<T> const &TheEXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  size_t nbRow = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  MyMatrix<T> QMat(nbCol, nbCol);
  for (size_t iCol = 0; iCol < nbCol; iCol++) {
    for (size_t jCol = iCol; jCol < nbCol; jCol++) {
      T eSum(0);
      for (size_t iRow = 0; iRow < nbRow; iRow++)
        eSum += TheEXT(iRow, jCol) * TheEXT(iRow, iCol);
      QMat(iCol, jCol) = eSum;
      QMat(jCol, iCol) = eSum;
    }
  }
  return Inverse_destroy(QMat);
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
GetQmatrix(MyMatrix<T> const &TheEXT, [[maybe_unused]] std::ostream &os) {
  return Kernel_GetQmatrix(TheEXT);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
GetQmatrix(MyMatrix<T> const &TheEXT, [[maybe_unused]] std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif

  MyMatrix<Tfield> TheEXT_F = UniversalMatrixConversion<Tfield, T>(TheEXT);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: UniversalMatrixConversion1|=" << time << "\n";
#endif

  MyMatrix<Tfield> Q_F = Kernel_GetQmatrix(TheEXT_F);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: Kernel_GetQmatrix|=" << time << "\n";
#endif

  MyMatrix<Tfield> Q_F_red = RemoveFractionMatrix(Q_F);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: RemoveFractionMatrix|=" << time << "\n";
#endif

  MyMatrix<T> RetMat = UniversalMatrixConversion<T, Tfield>(Q_F_red);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: UniversalMatrixConversion2|=" << time << "\n";
#endif
  return RetMat;
}

template <typename T>
MyMatrix<T> GetQmatrix_NotFullRank(MyMatrix<T> const &TheEXT,
                                   std::ostream &os) {
  std::vector<int> eList = ColumnReductionSet(TheEXT);
  int len = eList.size();
  int n_cols = TheEXT.cols();
  MyMatrix<T> TheEXT_red = SelectColumn(TheEXT, eList);
  MyMatrix<T> Qmat_red = GetQmatrix(TheEXT_red, os);
  MyMatrix<T> Qmat = ZeroMatrix<T>(n_cols, n_cols);
  for (int i = 0; i < len; i++) {
    int i2 = eList[i];
    for (int j = 0; j < len; j++) {
      int j2 = eList[j];
      Qmat(i2, j2) = Qmat_red(i, j);
    }
  }
  return Qmat;
}

template <typename T, typename Tidx, typename Treturn, typename F>
Treturn FCT_EXT_Qinput(MyMatrix<T> const &TheEXT, MyMatrix<T> const &Qinput,
                       F f) {
  size_t nbRow = TheEXT.rows();
  size_t max_val = std::numeric_limits<Tidx>::max();
  if (nbRow > max_val) {
    std::cerr << "PES: Error in FCT_EXT_Qinput due to too small coefficient range\n";
    std::cerr << "PES: nbRow=" << nbRow
              << " std::numeric_limits<Tidx>::max()=" << max_val << "\n";
    throw TerminalException{1};
  }
  size_t nbCol = TheEXT.cols();
  using Tfield = typename overlying_field<T>::field_type;
  MyVector<T> V(nbCol);
  MyVector<T> Vtr(nbCol);
  // Functions for computing the weighted matrix entries.
  auto f1 = [&](size_t iRow) -> void {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eSum(0);
      for (size_t jCol = 0; jCol < nbCol; jCol++)
        eSum += Qinput(iCol, jCol) * TheEXT(iRow, jCol);
      V(iCol) = eSum;
    }
  };
  auto f2 = [&](size_t jRow) -> T {
    T eSum(0);
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      eSum += V(iCol) * TheEXT(jRow, iCol);
    return eSum;
  };
  auto f1tr = [&](size_t iRow) -> void {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eSum(0);
      for (size_t jCol = 0; jCol < nbCol; jCol++)
        eSum += Qinput(iCol, jCol) * TheEXT(iRow, jCol);
      Vtr(iCol) = eSum;
    }
  };
  auto f2tr = [&](size_t jRow) -> T {
    T eSum(0);
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      eSum += Vtr(iCol) * TheEXT(jRow, iCol);
    return eSum;
  };
  // Preemptive check that the subset is adequate
  auto f3 = [&](std::vector<Tidx> const &Vsubset) -> bool {
    return IsSubsetFullRank<T, Tfield, Tidx>(TheEXT, Vsubset);
  };
  // Extension of the partial automorphism
  auto f4 = [&](const std::vector<Tidx> &Vsubset, const std::vector<Tidx> &Vin,
                const std::vector<std::vector<Tidx>> &ListBlocks)
      -> DataMapping<Tidx> {
    std::optional<MyMatrix<Tfield>> test1 =
        FindMatrixTransformationTest_Subset<T, Tfield, Tidx>(TheEXT, Vsubset,
                                                             Vin);
    Face block_status(ListBlocks.size());
    if (!test1)
      return {false, block_status, {}};
    return RepresentVertexPermutationTest_Blocks<T, Tfield, Tidx>(
        TheEXT, *test1, Vsubset, Vin, ListBlocks);
  };
  // Extension of the partial canonicalization
  auto f5 = [&](std::vector<Tidx> const &Vsubset,
                std::vector<Tidx> const &PartOrd) -> std::vector<Tidx> {
    return ExtendPartialCanonicalization<T, Tfield, Tidx>(TheEXT, Vsubset,
                                                          PartOrd);
  };
  bool is_symm = IsSymmetricMatrix(Qinput);
  return f(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, is_symm);
}

template <typename T, typename Tidx, typename Treturn, typename F>
Treturn FCT_EXT_Qinv(MyMatrix<T> const &TheEXT, F f, std::ostream &os) {
  MyMatrix<T> Qmat = GetQmatrix(TheEXT, os);
  return FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(TheEXT, Qmat, f);
}

template <typename T, typename Tfield, typename Tidx>
std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>>
IsomorphismFromCanonicReord(const MyMatrix<T> &EXT1, const MyMatrix<T> &EXT2,
                            const std::vector<Tidx> &CanonicReord1,
                            const std::vector<Tidx> &CanonicReord2) {
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  size_t nbRow = EXT1.rows();
  // Building the combinatorial equivalence
  std::vector<Tidx> ListIdx(nbRow);
  for (size_t idx = 0; idx < nbRow; idx++)
    ListIdx[CanonicReord1[idx]] = CanonicReord2[idx];
  // Building the matrix equivalence
  MyMatrix<Tfield> Basis1 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT1, CanonicReord1);
  MyMatrix<Tfield> Basis2 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT2, CanonicReord2);
  MyMatrix<Tfield> P = Inverse(Basis1) * Basis2;
  // Now testing the obtained mappings
  bool test = CheckEquivalence(EXT1, EXT2, ListIdx, P);
  if (!test) {
    // We fail the polytope equivalence
    return {};
  }
  std::pair<std::vector<Tidx>, MyMatrix<Tfield>> IsoInfo{ListIdx, P};
  return IsoInfo;
}

template <typename T, typename Tfield, typename Tidx>
std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>>
IsomorphismFromCanonicReord_GramMat(const MyMatrix<T> &EXT1,
                                    const MyMatrix<T> &GramMat1,
                                    const MyMatrix<T> &EXT2,
                                    const MyMatrix<T> &GramMat2,
                                    const std::vector<Tidx> &CanonicReord1,
                                    const std::vector<Tidx> &CanonicReord2,
                                    [[maybe_unused]] std::ostream &os) {
  size_t nbRow = EXT1.rows();
  // Building the combinatorial equivalence
  std::vector<Tidx> ListIdx(nbRow);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: |EXT1|=" << EXT1.rows() << " |EXT2|=" << EXT2.rows() << "\n";
  os << "PES: |CanonicRedord1|=" << CanonicReord1.size()
     << " |CanonicReord2|=" << CanonicReord2.size() << "\n";
  os << "PES: |GramMat1|=" << GramMat1.rows() << " / " << GramMat1.cols()
     << "\n";
  os << "PES: |GramMat2|=" << GramMat2.rows() << " / " << GramMat2.cols()
     << "\n";
  if (nbRow != CanonicReord1.size()) {
    std::cerr << "PES: nbRow=" << nbRow
              << " |CanonicReord1|=" << CanonicReord1.size() << "\n";
    std::cerr << "PES: CanonicReord1 should be of length nbRow\n";
    throw TerminalException{1};
  }
  if (nbRow != CanonicReord2.size()) {
    std::cerr << "PES: nbRow=" << nbRow
              << " |CanonicReord2|=" << CanonicReord2.size() << "\n";
    std::cerr << "PES: CanonicReord2 should be of length nbRow\n";
    throw TerminalException{1};
  }
#endif
  for (size_t idx = 0; idx < nbRow; idx++) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
    size_t pos = static_cast<size_t>(CanonicReord1[idx]);
    if (pos >= nbRow) {
      std::cerr << "PES: pos=" << pos << " nbRow=" << nbRow << "\n";
      std::cerr << "PES: CanonicReord1 entry is above the tange\n";
      throw TerminalException{1};
    }
#endif
    ListIdx[CanonicReord1[idx]] = CanonicReord2[idx];
  }
  // Building the matrix equivalence
  MyMatrix<Tfield> Basis1 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT1, CanonicReord1);
  MyMatrix<Tfield> Basis2 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT2, CanonicReord2);
  MyMatrix<Tfield> P = Inverse(Basis1) * Basis2;
  // Now testing the obtained mappings
  bool test = CheckEquivalence(EXT1, EXT2, ListIdx, P);
  if (!test) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB_TRACK_ISOMORPHISM
    os << "PES: We fail the polytope equivalence\n";
#endif
    return {};
  }
  // We check that EXT1 P = EXT2
  // EXT1 GramMat1 EXT1^T = EXT2 GramMat2 EXT2^T
  // EXT1 GramMat1 EXT1^T = EXT1 P GramMat2 P^T EXT1^T
  // So, GramMat1 = P GramMat2 P^T
  MyMatrix<Tfield> GramMat1_Tfield =
      UniversalMatrixConversion<Tfield, T>(GramMat1);
  MyMatrix<Tfield> GramMat2_Tfield =
      UniversalMatrixConversion<Tfield, T>(GramMat2);
  MyMatrix<Tfield> eProd = P * GramMat2_Tfield * P.transpose();
  if (eProd != GramMat1_Tfield) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB_TRACK_ISOMORPHISM
    os << "PES: We fail the Gram test\n";
#endif
    return {};
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB_TRACK_ISOMORPHISM
  os << "PES: We find an isomorphism\n";
#endif
  std::pair<std::vector<Tidx>, MyMatrix<Tfield>> IsoInfo{ListIdx, P};
  return IsoInfo;
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetSimpleWeightMatrix(MyMatrix<T> const &TheEXT, MyMatrix<T> const &Qinput,
                      std::ostream &os) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2,
               [[maybe_unused]] auto f1tr, [[maybe_unused]] auto f2tr, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5,
               [[maybe_unused]] bool is_symm) -> Treturn {
    return WeightMatrix<true, T, Tidx_value>(nbRow, f1, f2, os);
  };
  //
  size_t n_rows = TheEXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max())) {
    using Tidx = uint8_t;
    return FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(TheEXT, Qinput, f);
  }
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max())) {
    using Tidx = uint16_t;
    return FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(TheEXT, Qinput, f);
  }
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max())) {
    using Tidx = uint32_t;
    return FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(TheEXT, Qinput, f);
  }
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max())) {
    using Tidx = uint64_t;
    return FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(TheEXT, Qinput, f);
  }
#endif
  std::cerr << "PES: Failed to find matching numeric in GetSimpleWeightMatrix\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetWeightMatrix(MyMatrix<T> const &TheEXT,
                                                  std::ostream &os) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f1tr,
               [[maybe_unused]] auto f2tr,
               [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5,
               [[maybe_unused]] bool is_symm) -> Treturn {
    return WeightMatrix<true, T, Tidx_value>(nbRow, f1, f2, os);
  };
  //
  size_t n_rows = TheEXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max())) {
    using Tidx = uint8_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
  }
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max())) {
    using Tidx = uint16_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
  }
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max())) {
    using Tidx = uint32_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
  }
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max())) {
    using Tidx = uint64_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
  }
#endif
  std::cerr << "PES: Failed to find matching numeric in GetWeightMatrix\n";
  throw TerminalException{1};
}

template <typename T>
WeightMatrixLimited<true, T> GetWeightMatrixLimited(MyMatrix<T> const &TheEXT,
                                                    size_t max_offdiag,
                                                    std::ostream &os) {
  using Treturn = WeightMatrixLimited<true, T>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f1tr, [[maybe_unused]] auto f2tr,
               [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5,
               [[maybe_unused]] bool is_symm) -> Treturn {
    return WeightMatrixLimited<true, T>(nbRow, f1, f2, max_offdiag, os);
  };
  using Tidx = size_t;
  return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
}

template <typename Tvalue, typename Tgroup, typename Tidx_value, typename F1,
          typename F2, typename F1tr, typename F2tr, typename F3, typename F4>
std::vector<std::vector<typename Tgroup::Telt::Tidx>> f_for_stab(size_t nbRow, F1 f1, F2 f2, F1tr f1tr, F2tr f2tr, F3 f3,
                                                                 F4 f4, bool is_symm,
                                                                 std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: f_for_stab: nbRow=" << nbRow << " threshold=" << THRESHOLD_USE_SUBSET_SCHEME_STAB
     << "\n";
#endif
  auto f_heuristic=[&]() -> std::vector<std::vector<Tidx>> {
    if (is_symm) {
      return GetStabilizerWeightMatrix_Heuristic<Tvalue, Tidx, true>(nbRow, f1, f2, f1tr, f2tr, f3, f4, os);
    } else {
      return GetStabilizerWeightMatrix_Heuristic<Tvalue, Tidx, false>(nbRow, f1, f2, f1tr, f2tr, f3, f4, os);
    }
  };
  auto f_kernel=[&]() -> std::vector<std::vector<Tidx>> {
    if (is_symm) {
      WeightMatrix<true, Tvalue, Tidx_value> WMat(nbRow, f1, f2, os);
      return GetStabilizerWeightMatrix_Kernel<Tvalue, Tgr, Tidx, Tidx_value,
                                              true>(WMat, os);
    } else {
      WeightMatrix<false, Tvalue, Tidx_value> WMat(nbRow, f1, f2, os);
      return GetStabilizerWeightMatrix_Kernel<Tvalue, Tgr, Tidx, Tidx_value,
                                              false>(WMat, os);
    }
  };
#ifdef SANITY_CHECK_THRESHOLD_SCHEME
  std::vector<Tgroup> LGrp;
  for (int i=0; i<2; i++) {
    auto iife=[&]() -> std::vector<std::vector<Tidx>> {
      if (i == 0) {
        return f_heuristic();
      } else {
        return f_kernel();
      }
    };
    std::vector<std::vector<Tidx>> LGen1 = iife();
    std::vector<Telt> LGen2;
    Tidx n_act = 0;
    for (auto & eGen1 : LGen1) {
      n_act = eGen1.size();
      Telt eGen2(eGen1);
      LGen2.push_back(eGen2);
    }
    Tgroup eGrp(LGen2, n_act);
    LGrp.push_back(eGrp);
  }
  if (LGrp[0] != LGrp[1]) {
    std::cerr << "PES: The groups computed from the different scheme are not coherent\n";
    throw TerminalException{1};
  }
#endif

  if (nbRow > THRESHOLD_USE_SUBSET_SCHEME_STAB) {
    return f_heuristic();
  } else {
    return f_kernel();
  }
}

template <typename T, typename Tgroup, typename Tidx_value>
Tgroup LinPolytope_Automorphism_GramMat_Tidx_value(MyMatrix<T> const &EXT,
                                                   MyMatrix<T> const &GramMat,
                                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  HumanTime time;
#endif
  using Treturn = std::vector<std::vector<Tidx>>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f1tr, auto f2tr, auto f3, auto f4,
               [[maybe_unused]] auto f5, bool is_symm) -> Treturn {
    return f_for_stab<T, Tgroup, Tidx_value, decltype(f1), decltype(f2), decltype(f1tr), decltype(f2tr),
                      decltype(f3), decltype(f4)>(nbRow, f1, f2, f1tr, f2tr, f3, f4,
                                                  is_symm, os);
  };
  Treturn ListGen =
      FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(EXT, GramMat, f);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: LinPolytope_Aut : FCT_EXT_Qinput|=" << time << "\n";
#endif
  std::vector<Telt> LGen;
  for (auto &eList : ListGen) {
    LGen.push_back(Telt(eList));
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: LinPolytope_Aut : LGen|=" << time << "\n";
#endif
  return Tgroup(LGen, nbRow);
}

template <typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism_GramMat(MyMatrix<T> const &EXT,
                                        MyMatrix<T> const &GramMat,
                                        std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, Tgroup, uint8_t>(
        EXT, GramMat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, Tgroup, uint16_t>(
        EXT, GramMat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, Tgroup, uint32_t>(
        EXT, GramMat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, Tgroup, uint64_t>(
        EXT, GramMat, os);
  }
  std::cerr << "PES: Failed to find a matching type\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism(MyMatrix<T> const &EXT, std::ostream &os) {
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  return LinPolytope_Automorphism_GramMat<T, Tgroup>(EXTred, Qmat, os);
}

template <typename Tvalue, typename Tidx, typename Tidx_value, typename F1,
          typename F2, typename F1tr, typename F2tr, typename F3, typename F4, typename F5>
std::vector<Tidx> f_for_canonic(size_t nbRow, F1 f1, F2 f2, F1tr f1tr, F2tr f2tr, F3 f3, F4 f4, F5 f5,
                                bool is_symm, size_t threshold, std::ostream &os) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: f_for_canonic: nbRow=" << nbRow << " threshold=" << threshold << " is_symm=" << is_symm << "\n";
#endif
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  if (nbRow > threshold) {
    if (is_symm) {
      return GetGroupCanonicalizationVector_Heuristic<Tvalue, Tidx, true>(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, os)
        .first;
    } else {
      return GetGroupCanonicalizationVector_Heuristic<Tvalue, Tidx, false>(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, os)
        .first;
    }
  } else {
    if (is_symm) {
      WeightMatrix<true, Tvalue, Tidx_value> WMat(nbRow, f1, f2, os);
      WMat.ReorderingSetWeight();
      return GetGroupCanonicalizationVector_Kernel<Tvalue, Tgr, Tidx,
                                                   Tidx_value, true>(WMat, os)
          .first;
    } else {
      WeightMatrix<false, Tvalue, Tidx_value> WMat(nbRow, f1, f2, os);
      WMat.ReorderingSetWeight();
      return GetGroupCanonicalizationVector_Kernel<Tvalue, Tgr, Tidx,
                                                   Tidx_value, false>(WMat, os)
          .first;
    }
  }
}

template <typename T, typename Tidx, typename Tidx_value>
std::vector<Tidx> LinPolytope_CanonicOrdering_GramMat_Tidx_value(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat, size_t threshold, std::ostream &os) {
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif

  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f1tr, auto f2tr, auto f3, auto f4, auto f5,
               bool is_symm) -> Treturn {
    return f_for_canonic<T, Tidx, Tidx_value, decltype(f1), decltype(f2),
                         decltype(f1tr), decltype(f2tr),
                         decltype(f3), decltype(f4), decltype(f5)>(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, is_symm, threshold, os);
  };
  std::vector<Tidx> CanonicOrd =
      FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(EXT, GramMat, f);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: FCT_EXT_Qinput|=" << time << "\n";
#endif
  return CanonicOrd;
}

template <typename T, typename Tidx>
std::vector<Tidx> LinPolytope_CanonicOrdering_GramMat(MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat, size_t threshold, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, uint8_t>(EXT, GramMat, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, uint16_t>(EXT, GramMat, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, uint32_t>(EXT, GramMat, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, uint64_t>(EXT, GramMat, threshold, os);
  }
  std::cerr << "PES: Failed to find a match for Tidx_value\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx>
std::vector<Tidx> LinPolytope_CanonicOrdering(MyMatrix<T> const &EXT,
                                              size_t threshold, std::ostream &os) {
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  return LinPolytope_CanonicOrdering_GramMat<T, Tidx>(EXTred, Qmat, threshold, os);
}

template <typename T, typename Tidx>
MyMatrix<T> LinPolytope_CanonicForm_Tidx(MyMatrix<T> const &EXT,
                                         size_t threshold, std::ostream &os) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  std::vector<Tidx> CanonicOrd = LinPolytope_CanonicOrdering<T, Tidx>(EXT, threshold, os);
  MyMatrix<T> EXTreord(n_rows, n_cols);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    size_t j_row = CanonicOrd[i_row];
    EXTreord.row(i_row) = EXT.row(j_row);
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: CanonicOrdering + EXTreord|=" << time << "\n";
#endif

  MyMatrix<T> RedMat = CanonicalizeOrderedMatrix(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: CanonicalizeOrderedMatrix|=" << time << "\n";
#endif
  return RedMat;
}

template <typename T>
MyMatrix<T> LinPolytope_CanonicForm(MyMatrix<T> const &EXT, size_t threshold, std::ostream &os) {
  size_t n_rows = EXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, uint8_t>(EXT, threshold, os);
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, uint16_t>(EXT, threshold, os);
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, uint32_t>(EXT, threshold, os);
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, uint64_t>(EXT, threshold, os);
#endif
  std::cerr << "PES: LinPolytope_CanonicForm : Failed to find matching numeric\n";
  throw TerminalException{1};
}

template <typename T>
void check_iso_info_coherence(std::optional<T> const& IsoInfo1, std::optional<T> const& IsoInfo2, std::string const& context) {
  if (!IsoInfo1 && IsoInfo2) {
    std::cerr << "PES: IsoInfo1 fails to find isomorphism, but IsoInfo2 does\n";
    std::cerr << "PES: context=" << context << "\n";
    throw TerminalException{1};
  }
  if (IsoInfo1 && !IsoInfo2) {
    std::cerr << "PES: IsoInfo1 finds isomorphism but IsoInfo2 does not\n";
    std::cerr << "PES: context=" << context << "\n";
    throw TerminalException{1};
  }
}



template <typename T, typename Tidx>
std::optional<std::vector<Tidx>>
LinPolytope_Isomorphism(const MyMatrix<T> &EXT1, const MyMatrix<T> &EXT2,
                        std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  auto f_eval=[&](size_t threshold) -> std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> {
    std::vector<Tidx> CanonicReord1 = LinPolytope_CanonicOrdering<T, Tidx>(EXT1, threshold, os);
    std::vector<Tidx> CanonicReord2 = LinPolytope_CanonicOrdering<T, Tidx>(EXT2, threshold, os);
    return IsomorphismFromCanonicReord<T, Tfield, Tidx>(EXT1, EXT2, CanonicReord1, CanonicReord2);
  };
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo = f_eval(THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
#ifdef SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_CANONIC
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo_B = f_eval(THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
  check_iso_info_coherence(IsoInfo, IsoInfo_B, "LinPolytope_Isomorphism");
#endif
  if (!IsoInfo)
    return {};
  return IsoInfo->first;
}

template <typename T, typename Tidx>
std::optional<std::vector<Tidx>> LinPolytope_Isomorphism_GramMat(
    const MyMatrix<T> &EXT1, const MyMatrix<T> &GramMat1,
    const MyMatrix<T> &EXT2, const MyMatrix<T> &GramMat2, std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: |EXT1|=" << EXT1.rows() << " |EXT2|=" << EXT2.rows() << "\n";
#endif
  auto f_eval=[&](size_t threshold) -> std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> {
    std::vector<Tidx> CanonicReord1 = LinPolytope_CanonicOrdering_GramMat<T, Tidx>(EXT1, GramMat1, threshold, os);
    std::vector<Tidx> CanonicReord2 = LinPolytope_CanonicOrdering_GramMat<T, Tidx>(EXT2, GramMat2, threshold, os);
    return IsomorphismFromCanonicReord_GramMat<T, Tfield, Tidx>(EXT1, GramMat1, EXT2, GramMat2, CanonicReord1, CanonicReord2, os);
  };
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo = f_eval(THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
#ifdef SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_CANONIC
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo_B = f_eval(THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
  check_iso_info_coherence(IsoInfo, IsoInfo_B, "LinPolytope_Isomorphism_GramMat");
#endif
  if (!IsoInfo)
    return {};
  return IsoInfo->first;
}

template <typename T>
bool is_family_symmmetric(std::vector<MyMatrix<T>> const &ListMat) {
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      return false;
    }
  }
  return true;
}

template <typename T> struct ListMatSymm_Vdiag_WeightMat {
  MyMatrix<T> const &EXT;
  std::vector<MyMatrix<T>> const &ListMat;
  std::vector<T> const &Vdiag;
  int nbRow;
  int nbCol;
  int nMat;
  MyMatrix<T> MatV;
  MyMatrix<T> MatVtr;
  std::vector<T> LScal;
  int i_set;
  int i_set_tr;
  ListMatSymm_Vdiag_WeightMat(MyMatrix<T> const &_EXT,
                              std::vector<MyMatrix<T>> const &_ListMat,
                              std::vector<T> const &_Vdiag)
      : EXT(_EXT), ListMat(_ListMat), Vdiag(_Vdiag), nbRow(EXT.rows()),
        nbCol(EXT.cols()), nMat(ListMat.size()), MatV(nMat, nbCol),
        MatVtr(nMat, nbCol), LScal(nMat + 1) {}
  void f1(int i) {
    for (int iMat = 0; iMat < nMat; iMat++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eSum(0);
        for (int jCol = 0; jCol < nbCol; jCol++) {
          eSum += ListMat[iMat](jCol, iCol) * EXT(i, jCol);
        }
        MatV(iMat, iCol) = eSum;
      }
    }
    i_set = i;
  }
  std::vector<T> f2(int j) {
    for (int iMat = 0; iMat < nMat; iMat++) {
      T eSum(0);
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += MatV(iMat, iCol) * EXT(j, iCol);
      LScal[iMat] = eSum;
    }
    T eVal(0);
    if (i_set == j)
      eVal = Vdiag[j];
    LScal[nMat] = eVal;
    return LScal;
  }
  void f1tr(int i) {
    for (int iMat = 0; iMat < nMat; iMat++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eSum(0);
        for (int jCol = 0; jCol < nbCol; jCol++) {
          eSum += ListMat[iMat](iCol, jCol) * EXT(i, jCol);
        }
        MatVtr(iMat, iCol) = eSum;
      }
    }
    i_set_tr = i;
  }
  std::vector<T> f2tr(int j) {
    for (int iMat = 0; iMat < nMat; iMat++) {
      T eSum(0);
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += MatVtr(iMat, iCol) * EXT(j, iCol);
      LScal[iMat] = eSum;
    }
    T eVal(0);
    if (i_set_tr == j)
      eVal = Vdiag[j];
    LScal[nMat] = eVal;
    return LScal;
  }
};

template <typename T> struct ListMat_Vdiag_WeightMat {
  MyMatrix<T> const &EXT;
  std::vector<MyMatrix<T>> const &ListMat;
  std::vector<T> const &Vdiag;
  int nbRow;
  int nbCol;
  int nMat;
  MyMatrix<T> MatV;
  std::vector<T> LScal;
  int i_set;
  ListMat_Vdiag_WeightMat(MyMatrix<T> const &_EXT,
                          std::vector<MyMatrix<T>> const &_ListMat,
                          std::vector<T> const &_Vdiag)
      : EXT(_EXT), ListMat(_ListMat), Vdiag(_Vdiag), nbRow(EXT.rows()),
        nbCol(EXT.cols()), nMat(ListMat.size()), MatV(nMat, nbCol),
        LScal(nMat + 1) {}
  void f1(int i) {
    if (i < nbRow) {
      // Needed only for the upper part of the matrix
      for (int iMat = 0; iMat < nMat; iMat++) {
        for (int iCol = 0; iCol < nbCol; iCol++) {
          T eSum = 0;
          for (int jCol = 0; jCol < nbCol; jCol++) {
            eSum += ListMat[iMat](jCol, iCol) * EXT(i, jCol);
          }
          MatV(iMat, iCol) = eSum;
        }
      }
    }
    i_set = i;
  }
  std::vector<T> f2(int j) {
    int pos_j = j / nbRow;
    int j_red = j % nbRow;
    int pos_i = i_set / nbRow;
    if (pos_i == pos_j) {
      for (int iMat = 1; iMat < nMat; iMat++)
        LScal[iMat] = 0;
      if (pos_i == 0) {
        LScal[0] = 0;
      } else {
        LScal[0] = 1;
      }
    } else {
      for (int iMat = 0; iMat < nMat; iMat++) {
        T eSum = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eSum += MatV(iMat, iCol) * EXT(j_red, iCol);
        LScal[iMat] = eSum;
      }
    }
    T eVal = 0;
    if (i_set == j)
      eVal = Vdiag[j_red];
    LScal[nMat] = eVal;
    return LScal;
  }
};

template <typename T, typename Tfield, typename Tidx>
DataMapping<Tidx> ExtendPartialAutomorphism(
    MyMatrix<T> const &EXT, const std::vector<Tidx> &Vsubset,
    const std::vector<Tidx> &Vin,
    const std::vector<std::vector<Tidx>> &ListBlocks,
    [[maybe_unused]] const std::vector<MyMatrix<T>> &ListMat,
    [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: Before FindMatrixTransformationTest_Subset\n";
#endif
  std::optional<MyMatrix<Tfield>> test1 =
      FindMatrixTransformationTest_Subset<T, Tfield, Tidx>(EXT, Vsubset, Vin);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: After test1=" << test1.has_value() << "\n";
#endif
  Face block_status(ListBlocks.size());
  if (!test1) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
    os << "PES: f4 exit false 1\n";
#endif
    return {false, block_status, {}};
  }
  const MyMatrix<Tfield> &P = *test1;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  for (auto &eMat : ListMat) {
    MyMatrix<Tfield> eMat_F = UniversalMatrixConversion<Tfield, T>(eMat);
    MyMatrix<Tfield> eProd = P * eMat_F * TransposedMat(P);
    if (!TestEqualityMatrix(eProd, eMat_F)) {
      std::cerr << "PES: The matrix P should preserve the matrices at this point\n";
      throw TerminalException{1};
    }
  }
#endif
  return RepresentVertexPermutationTest_Blocks<T, Tfield, Tidx>(
      EXT, P, Vsubset, Vin, ListBlocks);
}

// This a piece of code used for purposes of
// ---Computation of automorphism group
// ---Computing canonicalization
// ---(Testing isomorphism) unsure because canonicalization
//    does the job.
// ---Computation of the weight matrix.
//
// The code isusing a strategy of subset where subsets are
// used for equivalence and then we see if the conclusion
// reached can be used to get a general conclusion.
//
// The input of the function is:
// ---TheEXT: The polytope in input.
// ---ListMat: The list of matrices
//    Whether they are symmetric or not is an issue.
// ---Vdiag: The list of diagonal values.
// ---f: The funciton doing the actual processing
//    That takes nbRow and 5 functions f1, f2, f3,
//    f4, f5:
//    ---f1: Set up for the computation of the row
//       of a WeightMatrix
//    ---f2: The computation of the row of the
//       WeightMatrix
//    ---f3: Test whether a subset is acceptable
//       or not for testing isomorphism
//    ---f4: for a transformation on a subset
//       returns the following.
//       --correct: bool, whether it is indeed an
//         equivalence of the full set
//       --block_status: Whether the blocks are
//         preserved or not. So, to know which can be
//         used for extension
//       --eGen: The global transformation of the index
//         if correct=true.
//    ---f5: The input is a canonicalization of
//       a subset and returns from this a canonicalization
//       of the full set.
//
// Logic for how to handle the symmetry:
// ---If ListMat is symmetric then the formalism is quite
//    clear.
// ---For non-symmetric matrices, we have to assume that
//    one matrix (assumed to be the first) is non-degenerate.
// ---Shall we assume that one matrix is positive definite?
//    Statement in question: If (v_i), (v'_i) are full
//    dimensional vector family and Q, Q' are non-degenerate
//    matrices. If c_ij = v_i Q v_j and c'_ij = v'_i Q' v'_j
//    if c_ij = c'_ij then there exist a matrix A such that
//    v'_i A = v_i and A Q A^T = Q'.
//    Proof: If P is the matrix formed with rows the vectors
//    v_i and P' the matrix formed with rows the vectors v'_i
//    then we have P Q P^T = P' Q' P'^T.
//    assuming without loss of generality that the first n
//    rows of P are linearly independent, we obtain P_red, P'_red
//    the restricted matrix obtained by selecting those n rows
//    then we get P_red Q P_red^T = P'_red Q' P'_red^T.
//    So, P_red and P'_red are invertible and we define
//    A = P'_red^{-1} P_red we get A Q A^T = Q'.
//    So, we get P'_red A = P_red which results in v'_i A = v_i
//    for the first n vectors.
//    For the other vectors we have if we write
//    delta = v'_k A - v_k the following equality
//    For i <= n we have
//    delta Q v^T_i = v'_k A Q v^T_i - v_k Q v^T_i
//                  = v'_k A Q A^T v'_i^T - c_{k,i}
//                  = v'_k Q' v'_i^T - c_{k,i}
//                  = c'_{k,i} - c_{k,i} = 0
//    So, delta being orthogonal to a basis of v_i and Q
//    non-degenerate, delta = 0 and v'_k A = v_k for all k.
// ---So, the first matrix should be symmetric non-degenerate
// ---How can we handle the non-symmetric case?
//    ---Working with the overgraph with twice as many vertices
//       is the way to go.
//    ---However, with the subset system, we need to see where
//       we are doing that folding.
//    ---After reflection, when we do the subsetting, it has
//       to respect the partitioning. So, the subsetting has
//       to be done before the folding.
//    ---In other words, the folding is done at the same step
//       as when e go from edge colored to vertex colored
//       graphs.
//
// ---ListMat is assumed to be symmetric
// ---Note that EXT does not have to be of full rank.
//    It makes perfect sense to compute some group
//    and get it only as permutation group.
template <typename T, typename Tfield, typename Tidx, typename Treturn,
          typename F>
Treturn FCT_ListMat_Vdiag(MyMatrix<T> const &EXT,
                          std::vector<MyMatrix<T>> const &ListMat,
                          std::vector<T> const &Vdiag, F f,
                          [[maybe_unused]] std::ostream &os) {
  bool is_symm = is_family_symmmetric(ListMat);
  size_t nbRow = EXT.rows();
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: FCT_ListMat_Vdiag: is_symm=" << is_symm << " nbRow=" << nbRow << "\n";
#endif
  size_t max_val = std::numeric_limits<Tidx>::max();
  if (nbRow > max_val) {
    std::cerr
        << "PES: Error in FCT_ListMat_Vdiag due to too small coefficient range\n";
    std::cerr << "PES: nbRow=" << nbRow
              << " std::numeric_limits<Tidx>::max()=" << max_val << "\n";
    throw TerminalException{1};
  }
  // The lambda for the construction of the weight matrix.
  ListMatSymm_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
  auto f1 = [&](size_t iRow) -> void { lms.f1(iRow); };
  auto f2 = [&](size_t jRow) -> std::vector<T> { return lms.f2(jRow); };
  auto f1tr = [&](size_t iRow) -> void { lms.f1tr(iRow); };
  auto f2tr = [&](size_t jRow) -> std::vector<T> { return lms.f2tr(jRow); };
  // Preemptive check that the subset is adequate
  auto f3 = [&](std::vector<Tidx> const &Vsubset) -> bool {
    return IsSubsetFullRank<T, Tfield, Tidx>(EXT, Vsubset);
  };
  // Extension of the partial automorphism
  auto f4 = [&](const std::vector<Tidx> &Vsubset, const std::vector<Tidx> &Vin,
                const std::vector<std::vector<Tidx>> &ListBlocks)
      -> DataMapping<Tidx> {
    return ExtendPartialAutomorphism<T, Tfield, Tidx>(EXT, Vsubset, Vin,
                                                      ListBlocks, ListMat, os);
  };
  // Extension of the partial canonicalization
  auto f5 = [&](std::vector<Tidx> const &Vsubset,
                std::vector<Tidx> const &PartOrd) -> std::vector<Tidx> {
    return ExtendPartialCanonicalization<T, Tfield, Tidx>(EXT, Vsubset,
                                                          PartOrd);
  };
  return f(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, is_symm);
}

template <typename T, typename Tfield, typename Tidx, typename Tidx_value>
WeightMatrix<true, std::vector<T>, Tidx_value>
GetWeightMatrix_ListMat_Vdiag(MyMatrix<T> const &TheEXT,
                              std::vector<MyMatrix<T>> const &ListMat,
                              std::vector<T> const &Vdiag, std::ostream &os) {
  using Treturn = WeightMatrix<true, std::vector<T>, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f1tr,
               [[maybe_unused]] auto f2tr, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5,
               [[maybe_unused]] bool is_symm) -> Treturn {
    return WeightMatrix<true, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
  };
  return FCT_ListMat_Vdiag<T, Tfield, Tidx, Treturn, decltype(f)>(
      TheEXT, ListMat, Vdiag, f, os);
}

template <typename T, typename Tfield, typename Tidx_value>
size_t GetInvariant_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  int nbRow = EXT.rows();
  auto get_wmat = [&]() -> WeightMatrix<true, std::vector<T>, Tidx_value> {
    bool is_symm = is_family_symmmetric(ListMat);
    if (is_symm) {
      ListMatSymm_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
      auto f1 = [&](size_t iRow) -> void { lms.f1(iRow); };
      auto f2 = [&](size_t jRow) -> std::vector<T> { return lms.f2(jRow); };
      return WeightMatrix<true, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
    } else {
      ListMat_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
      auto f1 = [&](size_t iRow) -> void { lms.f1(iRow); };
      auto f2 = [&](size_t jRow) -> std::vector<T> { return lms.f2(jRow); };
      return WeightMatrix<true, std::vector<T>, Tidx_value>(2 * nbRow, f1, f2,
                                                            os);
    }
  };
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat = get_wmat();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: get_wmat|=" << time << "\n";
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: ReorderingSetWeight|=" << time << "\n";
#endif

  size_t e_hash =
      std::hash<WeightMatrix<true, std::vector<T>, Tidx_value>>()(WMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: hash|=" << time << "\n";
#endif
  return e_hash;
}

template <typename T, typename Tfield>
size_t GetInvariant_ListMat_Vdiag(MyMatrix<T> const &EXT,
                                  std::vector<MyMatrix<T>> const &ListMat,
                                  std::vector<T> const &Vdiag,
                                  std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return GetInvariant_ListMat_Vdiag_Tidx_value<T, Tfield, uint8_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return GetInvariant_ListMat_Vdiag_Tidx_value<T, Tfield, uint16_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return GetInvariant_ListMat_Vdiag_Tidx_value<T, Tfield, uint32_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return GetInvariant_ListMat_Vdiag_Tidx_value<T, Tfield, uint64_t>(
        EXT, ListMat, Vdiag, os);
  }
  std::cerr << "PES: Failed to find a matching type\n";
  throw TerminalException{1};
}

template <typename T, typename Tfield, typename Tgroup, typename Tidx_value>
std::vector<std::vector<typename Tgroup::Telt::Tidx>> GetListGenAutomorphism_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Treturn = std::vector<std::vector<Tidx>>;
#ifdef SANITY_CHECK_POLYTOPE_EQUI_STAB
  for (!is_family_symmmetric(ListMat)) {
    std::cerr << "PES: The matrices of ListMat are not symmetric\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f1tr, auto f2tr, auto f3, auto f4,
               [[maybe_unused]] auto f5, bool is_symm) -> Treturn {
    return f_for_stab<std::vector<T>, Tgroup, Tidx_value, decltype(f1),
                      decltype(f2), decltype(f1tr), decltype(f2tr),
                      decltype(f3), decltype(f4)>(nbRow, f1, f2, f1tr, f2tr, f3, f4, is_symm, os);
  };
  Treturn ListGen = FCT_ListMat_Vdiag<T, Tfield, Tidx, Treturn, decltype(f)>(
      EXT, ListMat, Vdiag, f, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetListGenAutomorphism_ListMat_Vdiag|=" << time << "\n";
#endif
  return ListGen;
}

template <typename T, typename Tfield, typename Tgroup>
std::vector<std::vector<typename Tgroup::Telt::Tidx>> GetListGenAutomorphism_ListMat_Vdiag(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_val_poss = nbRow * nbRow / 2 + 1;
  if (max_val_poss < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<T, Tfield, Tgroup,
                                                           uint8_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<T, Tfield, Tgroup,
                                                           uint16_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<T, Tfield, Tgroup,
                                                           uint32_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<T, Tfield, Tgroup,
                                                           uint64_t>(
        EXT, ListMat, Vdiag, os);
  }
  std::cerr << "PES: Failed to find a matching Tidx_value\n";
  throw TerminalException{1};
}

template <typename T, typename Tfield, typename Tidx, typename Tidx_value>
std::vector<Tidx> Canonicalization_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, size_t threshold, std::ostream &os) {
#ifdef SANITY_CHECK_POLYTOPE_EQUI_STAB
  for (!is_family_symmmetric(ListMat)) {
    std::cerr << "PES: The matrices of ListMat are not symmetric\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f1tr, auto f2tr, auto f3, auto f4, auto f5,
               [[maybe_unused]] bool is_symm) -> Treturn {
    return f_for_canonic<std::vector<T>, Tidx, Tidx_value, decltype(f1), decltype(f2),
                         decltype(f1tr), decltype(f2tr), decltype(f3), decltype(f4),
                         decltype(f5)>(nbRow, f1, f2, f1tr, f2tr, f3, f4, f5, is_symm, threshold, os);
  };
  Treturn CanonicReord =
      FCT_ListMat_Vdiag<T, Tfield, Tidx, Treturn, decltype(f)>(EXT, ListMat,
                                                               Vdiag, f, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: Canonicalization_ListMat_Vdiag|=" << time << "\n";
#endif
  return CanonicReord;
}

template <typename T, typename Tfield, typename Tidx>
std::vector<Tidx>
Canonicalization_ListMat_Vdiag(MyMatrix<T> const &EXT,
                               std::vector<MyMatrix<T>> const &ListMat,
                               std::vector<T> const &Vdiag, size_t threshold, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint8_t>(EXT, ListMat, Vdiag, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint16_t>(EXT, ListMat, Vdiag, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint32_t>(EXT, ListMat, Vdiag, threshold, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint64_t>(EXT, ListMat, Vdiag, threshold, os);
  }
  std::cerr << "PES: No matching type for Tidx_value\n";
  throw TerminalException{1};
}

template <typename T, typename Tfield, typename Tidx, typename Tidx_value>
std::optional<std::vector<Tidx>> TestEquivalence_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT1, std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<T> const &Vdiag1, MyMatrix<T> const &EXT2,
    std::vector<MyMatrix<T>> const &ListMat2, std::vector<T> const &Vdiag2,
    std::ostream &os) {
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  size_t nbRow = EXT1.rows();
  auto f_eval=[&](size_t threshold) -> std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> {
    std::vector<Tidx> CanonicReord1 =
    Canonicalization_ListMat_Vdiag<T, Tfield, Tidx>(EXT1, ListMat1, Vdiag1,
                                                    threshold, os);
    std::vector<Tidx> CanonicReord2 =
    Canonicalization_ListMat_Vdiag<T, Tfield, Tidx>(EXT2, ListMat2, Vdiag2,
                                                    threshold, os);
    return IsomorphismFromCanonicReord<T, Tfield, Tidx>(EXT1, EXT2, CanonicReord1, CanonicReord2);
  };
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo = f_eval(THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
#ifdef SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_CANONIC
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo_B = f_eval(THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
  check_iso_info_coherence(IsoInfo, IsoInfo_B, "TestEquivalence_ListMat_Vdiag_Tidx_value");
#endif
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: We have IsoInfo\n";
#endif
  if (!IsoInfo) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
    os << "PES: Not equiv from IsomorphismFromCanonicReord\n";
#endif
    return {};
  }
  const std::vector<Tidx> &eList = IsoInfo->first;
  MyMatrix<Tfield> P = Inverse(IsoInfo->second);
  // The checks below are basically needed:
  // ---If we used a ListMat with ListMat.size() == 1 and it is
  // (sum_i v_i^T v_i)^{-1} then the success of the computation
  // IsomorphismFromCanonicReord is equivalent to the resolution
  // of the problems.
  // ---However, for more than 1 matrix, more is needed.
  // And so in the end it is better to check everything, if only
  // to be safe.
  size_t nMat = ListMat1.size();
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: nMat=" << nMat << "\n";
#endif
  for (size_t iMat = 0; iMat < nMat; iMat++) {
    MyMatrix<Tfield> eMat1 =
        UniversalMatrixConversion<Tfield, T>(ListMat1[iMat]);
    MyMatrix<Tfield> eMat2 =
        UniversalMatrixConversion<Tfield, T>(ListMat2[iMat]);
    MyMatrix<Tfield> eProd = P * eMat1 * TransposedMat(P);
    if (!TestEqualityMatrix(eProd, eMat2)) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
      os << "PES: Not equiv from TestEqualityMatrix from iMat=" << iMat << "/" << nMat << "\n";
#endif
      return {};
    }
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: Pass the ListMat1 / ListMat2 checks\n";
#endif
  Tidx nbRow_tidx = nbRow;
  for (Tidx i1 = 0; i1 < nbRow_tidx; i1++) {
    Tidx i2 = eList[i1];
    if (Vdiag1[i1] != Vdiag2[i2]) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
      os << "PES: Not equiv from Vdiag1 different from Vdiag2\n";
#endif
      return {};
    }
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: Before returning eList\n";
#endif
  return eList;
}

template <typename T, typename Tfield, typename Tidx>
std::optional<std::vector<Tidx>> TestEquivalence_ListMat_Vdiag(
    MyMatrix<T> const &EXT1, std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<T> const &Vdiag1, MyMatrix<T> const &EXT2,
    std::vector<MyMatrix<T>> const &ListMat2, std::vector<T> const &Vdiag2,
    std::ostream &os) {
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: Beginning of TestEquivalence_ListMat_Vdiag nbRow=" << EXT1.rows()
     << "\n";
#endif
  if (EXT1.rows() != EXT2.rows()) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
    os << "PES: Not equiv because |EXT1|=" << EXT1.rows()
       << " |EXT2|=" << EXT2.rows() << "\n";
#endif
    return {};
  }
  size_t nbRow = EXT1.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: max_poss_val=" << max_poss_val << "\n";
#endif
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint8_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint16_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint32_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, uint64_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  std::cerr << "PES: Failed to find a match for Tidx_value\n";
  throw TerminalException{1};
}

//
// The antipodal configurations and the absolute trick
//

template <typename T, typename Tidx_value> struct WeightMatrixAbs {
  Tidx_value positionZero;
  Face ArrSigns;
  WeightMatrix<true, T, Tidx_value> WMat;
};

template <typename T, typename Tidx_value>
WeightMatrixAbs<T, Tidx_value> GetSimpleWeightMatrixAntipodal_AbsTrick(
    MyMatrix<T> const &TheEXT, MyMatrix<T> const &Qmat, std::ostream &os) {
  static_assert(is_totally_ordered<T>::value,
                "Requires T to be a totally ordered field");
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  size_t nbPair = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  size_t n_ent = (nbPair * (nbPair + 1)) / 2;
  std::vector<Tidx_value> INP_TheMat(n_ent);
  Face ArrSigns(n_ent);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  Tidx_value miss_val = std::numeric_limits<Tidx_value>::max();
  Tidx_value positionZero = miss_val;
  //
  auto set_entry = [&](size_t iRow, size_t jRow, Tidx_value pos,
                       bool eChg) -> void {
    size_t idx = weightmatrix_idx<true>(nbPair, iRow, jRow);
    INP_TheMat[idx] = pos;
    ArrSigns[idx] = eChg;
  };
  MyVector<T> V(nbCol);
  for (size_t iPair = 0; iPair < nbPair; iPair++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eSum = 0;
      for (size_t jCol = 0; jCol < nbCol; jCol++)
        eSum += Qmat(iCol, jCol) * TheEXT(iPair, jCol);
      V(iCol) = eSum;
    }
    for (size_t jPair = 0; jPair <= iPair; jPair++) {
      T eScal = 0;
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eScal += V(iCol) * TheEXT(jPair, iCol);
      bool ChgSign = false;
      if (eScal < 0) {
        eScal = -eScal;
        ChgSign = true;
      }
      Tidx_value &value = ValueMap[eScal];
      if (value == 0) {
        // This is a missing value
        if (positionZero == miss_val && eScal == 0)
          positionZero = idxWeight;
        idxWeight++;
        value = idxWeight;
        INP_ListWeight.push_back(eScal);
      }
      Tidx_value pos = value - 1;
      set_entry(iPair, jPair, pos, ChgSign);
    }
  }
  /* This cannot be simplified be a classic constructor
     WeightMatrix(nbRow,f1,f2)
     because we also need to compute the positionZero and the ArrSigns. */
  bool weight_ordered = false;
  WeightMatrix<true, T, Tidx_value> WMat(nbPair, INP_TheMat, INP_ListWeight,
                                         weight_ordered, os);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: Before positionZero=" << positionZero << "\n";
#endif
  positionZero = WMat.ReorderingSetWeight_specificPosition(positionZero);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: After positionZero=" << positionZero << "\n";
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
#endif
  return {positionZero, std::move(ArrSigns), std::move(WMat)};
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetSimpleWeightMatrixAntipodal(MyMatrix<T> const &TheEXT,
                               MyMatrix<T> const &Qmat, std::ostream &os) {
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  size_t nbPair = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  size_t INP_nbRow = 2 * nbPair;
  size_t nb = nbPair * (2 * nbPair + 1);
  std::vector<Tidx_value> INP_TheMat(nb);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  //
  auto set_entry = [&](size_t iRow, size_t jRow, Tidx_value pos) -> void {
    size_t idx = weightmatrix_idx<true>(2 * nbPair, iRow, jRow);
    INP_TheMat[idx] = pos;
  };
  MyVector<T> V(nbCol);
  for (size_t iPair = 0; iPair < nbPair; iPair++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eSum = 0;
      for (size_t jCol = 0; jCol < nbCol; jCol++)
        eSum += Qmat(iCol, jCol) * TheEXT(iPair, jCol);
      V(iCol) = eSum;
    }
    for (size_t jPair = 0; jPair <= iPair; jPair++) {
      T eSum1 = 0;
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eSum1 += V(iCol) * TheEXT(jPair, iCol);
      T eSum2 = -eSum1;
      Tidx_value &value1 = ValueMap[eSum1];
      if (value1 == 0) {
        // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eSum1);
      }
      Tidx_value &value2 = ValueMap[eSum2];
      if (value2 == 0) {
        // This is a missing value
        idxWeight++;
        value2 = idxWeight;
        INP_ListWeight.push_back(eSum2);
      }
      Tidx_value pos1 = value1 - 1;
      Tidx_value pos2 = value2 - 1;
      set_entry(2 * iPair, 2 * jPair, pos1);
      set_entry(2 * iPair + 1, 2 * jPair, pos2);
      set_entry(2 * iPair, 2 * jPair + 1, pos2);
      set_entry(2 * iPair + 1, 2 * jPair + 1, pos1);
    }
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetSimpleWeightMatrixAntipodal|=" << time << "\n";
#endif
  bool weight_ordered = false;
  return WeightMatrix<true, T, Tidx_value>(INP_nbRow, INP_TheMat,
                                           INP_ListWeight, weight_ordered, os);
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixAntipodal(MyMatrix<T> const &TheEXT, std::ostream &os) {
  MyMatrix<T> Qmat = GetQmatrix(TheEXT, os);
  return GetSimpleWeightMatrixAntipodal<T, Tidx_value>(TheEXT, Qmat, os);
}

template <typename T> void SignRenormalizationMatrix(MyMatrix<T> &M) {
  size_t nbRow = M.rows();
  size_t n = M.cols();
  auto get_need_chgsign = [&](int const &iRow) -> bool {
    for (size_t i = 0; i < n; i++) {
      T eVal = M(iRow, i);
      if (eVal != 0) {
        return eVal < 0;
      }
    }
    return false;
  };
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (get_need_chgsign(iRow)) {
      for (size_t i = 0; i < n; i++)
        M(iRow, i) = -M(iRow, i);
    }
  }
}

template <typename T> MyMatrix<T> ExpandReducedMatrix(MyMatrix<T> const &M) {
  size_t nbPair = M.rows();
  size_t n = M.cols();
  MyMatrix<T> Mret(2 * nbPair, n);
  for (size_t iPair = 0; iPair < nbPair; iPair++)
    for (size_t i = 0; i < n; i++) {
      Mret(2 * iPair, i) = M(iPair, i);
      Mret(2 * iPair + 1, i) = -M(iPair, i);
    }
  return Mret;
}

/*
  Consider the case of the A2 root system with vectors
  \pm (1,0), \pm (0,1), \pm (1,1).
  If we consider the automorphisms of this vector configuration what we get is:
  ---Rotation by 2pi / 6 : Define subgroup of size 6
  ---Symmetry by axis.
  All together the group is of size 12.
  ----
  If we consider the absolute graph formed by the 3 vectors: (1,0), (0,1) and
  (1,1) then we get that this system defined a complete graph on 3 elements. So
  the group is of size 6. So, we indeed have the equality G = {\pm Id} x
  G_{abs}.
  ---
  The following holds:
  ---The construction of the weight matrix and so on means that orthogonal
  transformation on the vectors are not a problem.
  ---Since the absolute graph is the complete graph, we obtain that any ordering
  of the vector is possible by the canonicalization code.
  ---Thus if we put the vectors (1,0), (0,1) and (1,1)
  then the absolute canonicalization code may return us
  {(1,0), (0,1), (1,1)} or {(1,0), (1,1), (0,1)}.
  I think the hermite normal form of those are different.
  So, the method does not work.
  ---But we may be able to do something better. We still have the signs
  to be assigned.
*/
template <typename Tint, typename Tidx_value>
std::optional<MyMatrix<Tint>>
LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick_Tidx_value(
    MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &Qmat, std::ostream &os) {
  using Tgr = GraphBitset;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  WeightMatrixAbs<Tint, Tidx_value> WMatAbs =
      GetSimpleWeightMatrixAntipodal_AbsTrick<Tint, Tidx_value>(EXT, Qmat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
#endif

  using Tidx = uint32_t;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair =
      GetGroupCanonicalizationVector_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetGroupCanonicalizationVector_Kernel|=" << time << "\n";
#endif

  // We check if the Generating vector eGen can be mapped from the absolute
  // graph to the original one.
  auto TestExistSignVector =
      [&](std::vector<unsigned int> const &eGen) -> bool {
    /* We map a vector v_i to another v_j with sign +-1
       V[i] = 0 for unassigned
              1 for positive sign
              2 for negative sign
              3 for positive sign and treated
              4 for negative sign and treated
     */
    std::vector<uint8_t> V(nbRow, 0);
    V[0] = 1;
    while (true) {
      bool IsFinished = true;
      for (size_t i = 0; i < nbRow; i++) {
        uint8_t val = V[i];
        if (val < 3 && val != 0) {
          IsFinished = false;
          V[i] = val + 2;
          size_t iImg = eGen[i];
          for (size_t j = 0; j < nbRow; j++) {
            size_t jImg = eGen[j];
            Tidx_value pos = WMatAbs.WMat.GetValue(i, j);
            if (pos != WMatAbs.positionZero) {
              size_t idx1 = weightmatrix_idx<true>(nbRow, i, j);
              size_t idx2 = weightmatrix_idx<true>(nbRow, iImg, jImg);
              bool ChgSign1 = WMatAbs.ArrSigns[idx1];
              bool ChgSign2 = WMatAbs.ArrSigns[idx2];
              // ChgSign is true if ChgSign1 != ChgSign2
              bool ChgSign = ChgSign1 ^ ChgSign2;
              uint8_t valJ;
              if ((ChgSign && val == 1) || (!ChgSign && val == 2))
                valJ = 2;
              else
                valJ = 1;
              if (V[j] == 0) {
                V[j] = valJ;
              } else {
                if ((valJ % 2) != (V[j] % 2)) {
                  return false;
                }
              }
            }
          }
        }
      }
      if (IsFinished)
        break;
    }
    return true;
  };
  auto IsCorrectListGen = [&]() -> bool {
    for (auto &eGen : ePair.second) {
      bool test = TestExistSignVector(eGen);
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen())
    return {};
  //
  std::vector<Tidx> const &CanonicOrd = ePair.first;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif

  size_t n_cols = EXT.cols();
  MyMatrix<Tint> EXTreord(nbRow, n_cols);
  std::vector<int> ListSigns(nbRow, 0);
  ListSigns[0] = 1;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  std::string strAssign;
  os << "PES: positionZero=" << WMatAbs.positionZero << "\n";
#endif
  auto SetSign = [&](size_t const &i_row) -> void {
    int i_row_orig = CanonicOrd[i_row];
    for (size_t k_row = 0; k_row < nbRow; k_row++) {
      if (k_row != i_row && ListSigns[k_row] != 0) {
        int k_row_orig = CanonicOrd[k_row];
        if (WMatAbs.WMat.GetValue(i_row_orig, k_row_orig) !=
            WMatAbs.positionZero) {
          size_t idx = weightmatrix_idx<true>(nbRow, i_row_orig, k_row_orig);
          bool ChgSign = WMatAbs.ArrSigns[idx];
          int ValSign = 1 - 2 * static_cast<int>(ChgSign);
          int RetSign = ValSign * ListSigns[k_row];
          ListSigns[i_row] = RetSign;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
          strAssign += " (" + std::to_string(i_row) + " / " +
                       std::to_string(k_row) + ")";
#endif
          return;
        }
      }
    }
  };
  while (true) {
    int nbUndone = 0;
    for (size_t i_row = 0; i_row < nbRow; i_row++)
      if (ListSigns[i_row] == 0) {
        nbUndone++;
        SetSign(i_row);
      }
    if (nbUndone == 0)
      break;
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB_REMOVED
  // We have some crash due to this with the MD5 so, let us
  // outcomment it now.
  size_t eHash2 = MD5_hash_T<size_t>(strAssign);
  os << "PES: strAssign=" << strAssign << "\n";
  os << "PES: eHash2=" << eHash2 << "\n";
  std::string strWMat;
  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int i_rowC = CanonicOrd[i_row];
    for (size_t j_row = 0; j_row < nbRow; j_row++) {
      int j_rowC = CanonicOrd[j_row];
      Tidx_value pos = WMatAbs.WMat.GetValue(i_rowC, j_rowC);
      strWMat += " " + std::to_string(pos);
    }
  }
  for (auto &eVal : WMatAbs.WMat.GetWeight()) {
    strWMat += " " + std::to_string(eVal);
  }
  os << "PES: strWMat=" << strWMat << "\n";
  size_t eHash3 = MD5_hash_T<size_t>(strWMat);
  os << "PES: eHash3=" << eHash3 << "\n";
#endif
  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int j_row = CanonicOrd[i_row];
    int eSign = ListSigns[i_row];
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      EXTreord(i_row, i_col) = eSign * EXT(j_row, i_col);
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "PES: EXTreord=\n";
  WriteMatrix(os, EXTreord);
  WriteMatrixGAP(os, EXTreord);
  os << "\n";
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: EXTreord|=" << time << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: ComputeColHermiteNormalForm|=" << time << "\n";
#endif
  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|PES: SignRenormalizationMatrix|=" << time << "\n";
#endif
  return RedMat;
}

template <typename Tint>
std::optional<MyMatrix<Tint>> LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(
    MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &Qmat, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick_Tidx_value<
        Tint, uint8_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick_Tidx_value<
        Tint, uint16_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick_Tidx_value<
        Tint, uint32_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick_Tidx_value<
        Tint, uint64_t>(EXT, Qmat, os);
  }
  std::cerr << "PES: Failed to match for Tidx_value\n";
  throw TerminalException{1};
}

//
// Several asymmetric matrices
//

// The matrices in ListMat do not have to be symmetric.
template <typename T, typename Tint, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value>
T_TranslateToMatrix_ListMat_SHV(std::vector<MyMatrix<T>> const &ListMat,
                                MyMatrix<Tint> const &SHV, std::ostream &os) {
  size_t nbRow = SHV.rows();
  size_t n = SHV.cols();
  size_t nbMat = ListMat.size();
  std::vector<MyVector<T>> ListV(nbMat);
  auto f1 = [&](size_t iRow) -> void {
    for (size_t iMat = 0; iMat < nbMat; iMat++) {
      MyVector<T> V(n);
      for (size_t i = 0; i < n; i++) {
        T eVal = 0;
        for (size_t j = 0; j < n; j++)
          eVal += ListMat[iMat](j, i) * SHV(iRow, j);
        V(i) = eVal;
      }
      ListV[iMat] = V;
    }
  };
  std::vector<T> ListScal(nbMat);
  auto f2 = [&](size_t iCol) -> std::vector<T> {
    for (size_t iMat = 0; iMat < nbMat; iMat++) {
      T eScal = 0;
      for (size_t i = 0; i < n; i++)
        eScal += ListV[iMat](i) * SHV(iCol, i);
      ListScal[iMat] = eScal;
    }
    return ListScal;
  };
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
}

template <bool is_symmetric, typename T, typename Tidx_value>
WeightMatrix<is_symmetric, std::vector<T>, Tidx_value>
GetWeightMatrix_ListComm(MyMatrix<T> const &TheEXT, MyMatrix<T> const &GramMat,
                         std::vector<MyMatrix<T>> const &ListComm,
                         std::ostream &os) {
  size_t nbRow = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  size_t nbComm = ListComm.size();
  std::vector<MyMatrix<T>> ListProd;
  ListProd.push_back(GramMat);
  for (size_t iComm = 0; iComm < nbComm; iComm++) {
    MyMatrix<T> eProd = ListComm[iComm] * GramMat;
    ListProd.push_back(eProd);
  }
  MyMatrix<T> M(nbComm + 1, nbCol);
  auto f1 = [&](size_t iRow) -> void {
    for (size_t iMat = 0; iMat <= nbComm; iMat++) {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        T eSum = 0;
        for (size_t jCol = 0; jCol < nbCol; jCol++)
          eSum += ListProd[iMat](iCol, jCol) * TheEXT(iRow, jCol);
        M(iMat, iCol) = eSum;
      }
    }
  };
  std::vector<T> eVectSum(nbComm + 1);
  auto f2 = [&](size_t jRow) -> std::vector<T> {
    for (size_t iMat = 0; iMat <= nbComm; iMat++) {
      T eSum = 0;
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eSum += TheEXT(jRow, iCol) * M(iMat, iCol);
      eVectSum[iMat] = eSum;
    }
    return eVectSum;
  };
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
}

template <typename T, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value>
GetWeightMatrix_ListMatrix(std::vector<MyMatrix<T>> const &ListMatrix,
                           MyMatrix<T> const &TheEXT, std::ostream &os) {
  size_t nbRow = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  size_t nbMat = ListMatrix.size();
  MyMatrix<T> M(nbMat, nbCol);
  auto f1 = [&](size_t iRow) -> void {
    for (size_t iMat = 0; iMat < nbMat; iMat++) {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        T eSum = 0;
        for (size_t jCol = 0; jCol < nbCol; jCol++)
          eSum += ListMatrix[iMat](iCol, jCol) * TheEXT(iRow, jCol);
        M(iMat, iCol) = eSum;
      }
    }
  };
  std::vector<T> eVectScal(nbMat);
  auto f2 = [&](size_t jRow) -> std::vector<T> {
    for (size_t iMat = 0; iMat < nbMat; iMat++) {
      T eSum = 0;
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eSum += TheEXT(jRow, iCol) * M(iMat, iCol);
      eVectScal[iMat] = eSum;
    }
    return eVectScal;
  };
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
}

//
// Various construction of weighted matrix
//

template <typename T, typename Tint, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
T_TranslateToMatrix_QM_SHV(MyMatrix<T> const &qMat, MyMatrix<Tint> const &SHV,
                           std::ostream &os) {
  size_t nbRow = SHV.rows();
  size_t n = qMat.rows();
  size_t INP_nbRow = nbRow;
  size_t nbPair = nbRow / 2;
  size_t nb = nbPair * (2 * nbPair + 1);
  std::vector<Tidx_value> INP_TheMat(nb);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  //
  auto set_entry = [&](size_t iRow, size_t iCol, Tidx_value val) -> void {
    size_t idx = weightmatrix_idx<true>(nbRow, iRow, iCol);
    INP_TheMat[idx] = val;
  };
  for (size_t iPair = 0; iPair < nbPair; iPair++) {
    MyVector<T> V(n);
    for (size_t i = 0; i < n; i++) {
      T eVal = 0;
      for (size_t j = 0; j < n; j++)
        eVal += qMat(j, i) * SHV(2 * iPair, j);
      V(i) = eVal;
    }
    for (size_t jPair = 0; jPair <= iPair; jPair++) {
      T eScal = 0;
      for (size_t i = 0; i < n; i++)
        eScal += V(i) * SHV(2 * jPair, i);
      Tidx_value &value1 = ValueMap[eScal];
      if (value1 == 0) {
        // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eScal);
      }
      Tidx_value &value2 = ValueMap[-eScal];
      if (value2 == 0) {
        // This is a missing value
        idxWeight++;
        value2 = idxWeight;
        INP_ListWeight.push_back(-eScal);
      }
      Tidx_value pos1 = value1 - 1;
      Tidx_value pos2 = value2 - 1;
      set_entry(2 * iPair, 2 * jPair, pos1);
      set_entry(2 * iPair + 1, 2 * jPair, pos2);
      set_entry(2 * iPair, 2 * jPair + 1, pos2);
      set_entry(2 * iPair + 1, 2 * jPair + 1, pos1);
    }
  }
  bool weight_ordered = false;
  return WeightMatrix<true, T, Tidx_value>(INP_nbRow, INP_TheMat,
                                           INP_ListWeight, weight_ordered, os);
}

// clang-format off
#endif  // SRC_GROUP_POLYTOPEEQUISTAB_H_
// clang-format on
