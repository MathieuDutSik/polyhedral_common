// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_TEMP_POLYTOPEEQUISTAB_H_
#define SRC_GROUP_TEMP_POLYTOPEEQUISTAB_H_

// clang-format off
#include "MatrixGroup.h"
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
#endif

#ifdef TIMINGS
#define TIMINGS_POLYTOPE_EQUI_STAB
#define TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
#endif

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
  std::vector<Tidx> eList(n);
  for (size_t i = 0; i < n; i++) {
    Tidx eVal = OnPoints(i, *test);
    eList[i] = eVal;
  }
  return Telt(std::move(eList));
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
  os << "|UniversalMatrixConversion1|=" << time << "\n";
#endif

  MyMatrix<Tfield> Q_F = Kernel_GetQmatrix(TheEXT_F);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|Kernel_GetQmatrix|=" << time << "\n";
#endif

  MyMatrix<Tfield> Q_F_red = RemoveFractionMatrix(Q_F);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|RemoveFractionMatrix|=" << time << "\n";
#endif

  MyMatrix<T> RetMat = UniversalMatrixConversion<T, Tfield>(Q_F_red);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|UniversalMatrixConversion2|=" << time << "\n";
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
    std::cerr << "Error in FCT_EXT_Qinput due to too small coefficient range\n";
    std::cerr << "nbRow=" << nbRow
              << " std::numeric_limits<Tidx>::max()=" << max_val << "\n";
    throw TerminalException{1};
  }
  size_t nbCol = TheEXT.cols();
  using Tfield = typename overlying_field<T>::field_type;
  MyVector<T> V(nbCol);
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
  return f(nbRow, f1, f2, f3, f4, f5);
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
                                    const std::vector<Tidx> &CanonicReord2) {
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
    // We fail the Gram test
    return {};
  }
  std::pair<std::vector<Tidx>, MyMatrix<Tfield>> IsoInfo{ListIdx, P};
  return IsoInfo;
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetSimpleWeightMatrix(MyMatrix<T> const &TheEXT, MyMatrix<T> const &Qinput, std::ostream& os) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
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
  std::cerr << "Failed to find matching numeric in GetSimpleWeightMatrix\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetWeightMatrix(MyMatrix<T> const &TheEXT,
                                                  std::ostream &os) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
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
  std::cerr << "Failed to find matching numeric in GetWeightMatrix\n";
  throw TerminalException{1};
}

template <typename T>
WeightMatrixLimited<true, T> GetWeightMatrixLimited(MyMatrix<T> const &TheEXT,
                                                    size_t max_offdiag,
                                                    std::ostream &os) {
  using Treturn = WeightMatrixLimited<true, T>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
    return WeightMatrixLimited<true, T>(nbRow, f1, f2, max_offdiag, os);
  };
  using Tidx = size_t;
  return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f, os);
}

template <typename T, bool use_scheme, typename Tgroup, typename Tidx_value>
Tgroup LinPolytope_Automorphism_GramMat_Tidx_value(MyMatrix<T> const &EXT,
                                                   MyMatrix<T> const &GramMat,
                                                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphListAdj;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  HumanTime time;
#endif
  using Treturn = std::vector<std::vector<Tidx>>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               [[maybe_unused]] auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetStabilizerWeightMatrix_Heuristic<T, Tidx>(nbRow, f1, f2, f3, f4,
                                                          os);
    } else {
      WeightMatrix<true, T, Tidx_value> WMat(nbRow, f1, f2, os);
      return GetStabilizerWeightMatrix_Kernel<T, Tgr, Tidx, Tidx_value>(WMat,
                                                                        os);
    }
  };
  Treturn ListGen =
      FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(EXT, GramMat, f);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|LinPolytope_Aut : FCT_EXT_Qinput|=" << time << "\n";
#endif
  std::vector<Telt> LGen;
  for (auto &eList : ListGen) {
    LGen.push_back(Telt(eList));
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|LinPolytope_Aut : LGen|=" << time << "\n";
#endif
  return Tgroup(LGen, nbRow);
}

template <typename T, bool use_scheme, typename Tgroup>
Tgroup LinPolytope_Automorphism_GramMat(MyMatrix<T> const &EXT,
                                        MyMatrix<T> const &GramMat,
                                        std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, use_scheme, Tgroup,
                                                       uint8_t>(EXT, GramMat,
                                                                os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, use_scheme, Tgroup,
                                                       uint16_t>(EXT, GramMat,
                                                                 os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, use_scheme, Tgroup,
                                                       uint32_t>(EXT, GramMat,
                                                                 os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytope_Automorphism_GramMat_Tidx_value<T, use_scheme, Tgroup,
                                                       uint64_t>(EXT, GramMat,
                                                                 os);
  }
  std::cerr << "Failed to find a matching type\n";
  throw TerminalException{1};
}

template <typename T, bool use_scheme, typename Tgroup>
Tgroup LinPolytope_Automorphism(MyMatrix<T> const &EXT, std::ostream &os) {
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  return LinPolytope_Automorphism_GramMat<T, use_scheme, Tgroup>(EXTred, Qmat,
                                                                 os);
}

template <typename T, typename Tidx, bool use_scheme, typename Tidx_value>
std::vector<Tidx> LinPolytope_CanonicOrdering_GramMat_Tidx_value(
    MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat, std::ostream &os) {
  using Tgr = GraphBitset;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif

  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetGroupCanonicalizationVector_Heuristic<T, Tidx>(nbRow, f1, f2,
                                                               f3, f4, f5, os)
          .first;
    } else {
      WeightMatrix<true, T, Tidx_value> WMat(nbRow, f1, f2, os);
      WMat.ReorderingSetWeight();
      return GetGroupCanonicalizationVector_Kernel<T, Tgr, Tidx, Tidx_value>(
                 WMat, os)
          .first;
    }
  };
  std::vector<Tidx> CanonicOrd =
      FCT_EXT_Qinput<T, Tidx, Treturn, decltype(f)>(EXT, GramMat, f);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|FCT_EXT_Qinput|=" << time << "\n";
#endif
  return CanonicOrd;
}

template <typename T, typename Tidx, bool use_scheme>
std::vector<Tidx> LinPolytope_CanonicOrdering_GramMat(
    MyMatrix<T> const &EXT, MyMatrix<T> const &GramMat, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, use_scheme,
                                                          uint8_t>(EXT, GramMat,
                                                                   os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, use_scheme,
                                                          uint16_t>(
        EXT, GramMat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, use_scheme,
                                                          uint32_t>(
        EXT, GramMat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytope_CanonicOrdering_GramMat_Tidx_value<T, Tidx, use_scheme,
                                                          uint64_t>(
        EXT, GramMat, os);
  }
  std::cerr << "Failed to find a match for Tidx_value\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx, bool use_scheme>
std::vector<Tidx> LinPolytope_CanonicOrdering(MyMatrix<T> const &EXT,
                                              std::ostream &os) {
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  return LinPolytope_CanonicOrdering_GramMat<T, Tidx, use_scheme>(EXTred, Qmat,
                                                                  os);
}

template <typename T, bool use_scheme, typename Tidx>
MyMatrix<T> LinPolytope_CanonicForm_Tidx(MyMatrix<T> const &EXT,
                                         std::ostream &os) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  std::vector<Tidx> CanonicOrd =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT, os);
  MyMatrix<T> EXTreord(n_rows, n_cols);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    size_t j_row = CanonicOrd[i_row];
    EXTreord.row(i_row) = EXT.row(j_row);
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|CanonicOrdering + EXTreord|=" << time << "\n";
#endif

  MyMatrix<T> RedMat = CanonicalizeOrderedMatrix(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|CanonicalizeOrderedMatrix|=" << time << "\n";
#endif
  return RedMat;
}

template <typename T, bool use_scheme>
MyMatrix<T> LinPolytope_CanonicForm(MyMatrix<T> const &EXT, std::ostream &os) {
  size_t n_rows = EXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint8_t>(EXT, os);
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint16_t>(EXT, os);
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint32_t>(EXT, os);
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint64_t>(EXT, os);
#endif
  std::cerr << "LinPolytope_CanonicForm : Failed to find matching numeric\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx, bool use_scheme>
std::optional<std::vector<Tidx>>
LinPolytope_Isomorphism(const MyMatrix<T> &EXT1, const MyMatrix<T> &EXT2,
                        std::ostream &os) {
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT1, os);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT2, os);
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(EXT1, EXT2, CanonicReord1,
                                                   CanonicReord2);
  if (!IsoInfo)
    return {};
  return IsoInfo->first;
}

template <typename T, typename Tidx, bool use_scheme>
std::optional<std::vector<Tidx>> LinPolytope_Isomorphism_GramMat(
    const MyMatrix<T> &EXT1, const MyMatrix<T> &GramMat1,
    const MyMatrix<T> &EXT2, const MyMatrix<T> &GramMat2, std::ostream &os) {
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering_GramMat<T, Tidx, use_scheme>(EXT1, GramMat1,
                                                               os);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering_GramMat<T, Tidx, use_scheme>(EXT2, GramMat2,
                                                               os);
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord_GramMat<T, Tfield, Tidx>(
          EXT1, GramMat1, EXT2, GramMat2, CanonicReord1, CanonicReord2);
  if (!IsoInfo)
    return {};
  return IsoInfo->first;
}

template <typename Tint, typename Tidx, typename Tgroup, typename Tidx_value,
          typename Tgr, bool use_scheme>
std::optional<MyMatrix<Tint>>
LinPolytopeIntegral_Isomorphism(const MyMatrix<Tint> &EXT1,
                                const MyMatrix<Tint> &EXT2, std::ostream &os) {
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering<Tint, Tidx, use_scheme>(EXT1, os);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering<Tint, Tidx, use_scheme>(EXT2, os);
  //
  using Tfield = typename overlying_field<Tint>::field_type;
  using Telt = typename Tgroup::Telt;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<Tint, Tfield, Tidx>(EXT1, EXT2, CanonicReord1,
                                                      CanonicReord2);
  if (!IsoInfo)
    return {};
  Telt ePerm(IsoInfo->first);

  MyMatrix<Tfield> EXT1_T = UniversalMatrixConversion<Tfield, Tint>(EXT1);
  MyMatrix<Tfield> EXT2_T = UniversalMatrixConversion<Tfield, Tint>(EXT2);
  Tgroup GRP1 =
      LinPolytope_Automorphism<Tfield, use_scheme, Tgroup>(EXT1_T, os);
  std::optional<MyMatrix<Tfield>> eRes =
      LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, ePerm, os);
  if (eRes)
    return UniversalMatrixConversion<Tint, Tfield>(*eRes);
  return {};
}

template <typename Tint, typename Tidx, typename Tgroup, typename Tidx_value,
          typename Tgr, bool use_scheme>
Tgroup LinPolytopeIntegral_Automorphism(const MyMatrix<Tint> &EXT,
                                        std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom =
      LinPolytope_Automorphism<Tfield, use_scheme, Tgroup>(EXT_T, os);
  Tgroup GRP = LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom, os);
  return GRP;
}

template <typename Tint, typename Tidx, typename Tgroup, typename Tidx_value,
          typename Tgr, bool use_scheme>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Automorphism_RightCoset(const MyMatrix<Tint> &EXT,
                                            std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  using Telt = typename Tgroup::Telt;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom =
      LinPolytope_Automorphism<Tfield, use_scheme, Tgroup>(EXT_T, os);
  std::pair<Tgroup, std::vector<Telt>> pair =
      LinPolytopeIntegral_Stabilizer_RightCoset_Method8(EXT_T, GRPisom, os);
  return pair;
}

template<typename T>
bool is_family_symmmetric(std::vector<MyMatrix<T>> const &ListMat) {
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      return false;
    }
  }
  return true;
}

template<typename T>
struct ListMatSymm_Vdiag_WeightMat {
  MyMatrix<T> const& EXT;
  std::vector<MyMatrix<T>> const& ListMat;
  std::vector<T> const& Vdiag;
  int nbRow;
  int nbCol;
  int nMat;
  MyMatrix<T> MatV;
  std::vector<T> LScal;
  int i_set;
  ListMatSymm_Vdiag_WeightMat(MyMatrix<T> const& _EXT, std::vector<MyMatrix<T>> const& _ListMat, std::vector<T> const& _Vdiag) : EXT(_EXT), ListMat(_ListMat), Vdiag(_Vdiag), nbRow(EXT.rows()), nbCol(EXT.cols()), nMat(ListMat.size()), MatV(nMat, nbCol), LScal(nMat + 1) {
  }
  void f1(int i) {
    for (int iMat = 0; iMat < nMat; iMat++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eSum = 0;
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
      T eSum = 0;
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += MatV(iMat, iCol) * EXT(j, iCol);
      LScal[iMat] = eSum;
    }
    T eVal = 0;
    if (i_set == j)
      eVal = Vdiag[j];
    LScal[nMat] = eVal;
    return LScal;
  }
};

template<typename T>
struct ListMat_Vdiag_WeightMat {
  MyMatrix<T> const& EXT;
  std::vector<MyMatrix<T>> const& ListMat;
  std::vector<T> const& Vdiag;
  int nbRow;
  int nbCol;
  int nMat;
  MyMatrix<T> MatV;
  std::vector<T> LScal;
  int i_set;
  ListMat_Vdiag_WeightMat(MyMatrix<T> const& _EXT, std::vector<MyMatrix<T>> const& _ListMat, std::vector<T> const& _Vdiag) : EXT(_EXT), ListMat(_ListMat), Vdiag(_Vdiag), nbRow(EXT.rows()), nbCol(EXT.cols()), nMat(ListMat.size()), MatV(nMat, nbCol), LScal(nMat + 1) {
  }
  void f1(int i) {
    int i_red = i % nbRow;
    for (int iMat = 0; iMat < nMat; iMat++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eSum = 0;
        for (int jCol = 0; jCol < nbCol; jCol++) {
          eSum += ListMat[iMat](jCol, iCol) * EXT(i_red, jCol);
        }
        MatV(iMat, iCol) = eSum;
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

template<typename T, typename Tfield, typename Tidx>
DataMapping<Tidx> ExtendPartialAutomorphism(MyMatrix<T> const& EXT,
                                            const std::vector<Tidx> &Vsubset,
                                            const std::vector<Tidx> &Vin,
                                            const std::vector<std::vector<Tidx>> &ListBlocks,
                                            [[maybe_unused]] const std::vector<MyMatrix<T>> & ListMat,
                                            [[maybe_unused]] std::ostream& os) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "Before FindMatrixTransformationTest_Subset\n";
#endif
  std::optional<MyMatrix<Tfield>> test1 =
    FindMatrixTransformationTest_Subset<T, Tfield, Tidx>(EXT, Vsubset, Vin);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "After test1=" << test1.has_value() << "\n";
#endif
  Face block_status(ListBlocks.size());
  if (!test1) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB
    os << "f4 exit false 1\n";
#endif
    return {false, block_status, {}};
  }
  const MyMatrix<Tfield> &P = *test1;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  for (auto &eMat : ListMat) {
    MyMatrix<Tfield> eMat_F = UniversalMatrxConversion<Tfield,T>(eMat);
    MyMatrix<Tfield> eProd = P * eMat_F * TransposedMat(P);
    if (!TestEqualityMatrix(eProd, eMat_F)) {
      std::cerr << "The matrix P should preserve the matrices at this point\n";
      throw TerminalException{1};
    }
  }
#endif
  return RepresentVertexPermutationTest_Blocks<T, Tfield, Tidx>(EXT, P, Vsubset, Vin, ListBlocks);
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
// Some issues to address:



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
#ifdef SANITY_CHECK_POLYTOPE_EQUI_STAB
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      std::cerr << "The matrix eMat should be symmetric\n";
      throw TerminalException{1};
    }
  }
#endif
  size_t nbRow = EXT.rows();
  size_t max_val = std::numeric_limits<Tidx>::max();
  if (nbRow > max_val) {
    std::cerr
        << "Error in FCT_ListMat_Vdiag due to too small coefficient range\n";
    std::cerr << "nbRow=" << nbRow
              << " std::numeric_limits<Tidx>::max()=" << max_val << "\n";
    throw TerminalException{1};
  }
  // The lambda for the construction of the weight matrix.
  ListMatSymm_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
  auto f1 = [&](size_t iRow) -> void {
    lms.f1(iRow);
  };
  auto f2 = [&](size_t jRow) -> std::vector<T> {
    return lms.f2(jRow);
  };
  // Preemptive check that the subset is adequate
  auto f3 = [&](std::vector<Tidx> const &Vsubset) -> bool {
    return IsSubsetFullRank<T, Tfield, Tidx>(EXT, Vsubset);
  };
  // Extension of the partial automorphism
  auto f4 = [&](const std::vector<Tidx> &Vsubset, const std::vector<Tidx> &Vin,
                const std::vector<std::vector<Tidx>> &ListBlocks)
      -> DataMapping<Tidx> {
    return ExtendPartialAutomorphism<T, Tfield, Tidx>(EXT, Vsubset, Vin, ListBlocks, ListMat, os);
  };
  // Extension of the partial canonicalization
  auto f5 = [&](std::vector<Tidx> const &Vsubset,
                std::vector<Tidx> const &PartOrd) -> std::vector<Tidx> {
    return ExtendPartialCanonicalization<T, Tfield, Tidx>(EXT, Vsubset,
                                                          PartOrd);
  };
  return f(nbRow, f1, f2, f3, f4, f5);
}

template <typename T, typename Tfield, typename Tidx, typename Tidx_value>
WeightMatrix<true, std::vector<T>, Tidx_value>
GetWeightMatrix_ListMat_Vdiag(MyMatrix<T> const &TheEXT,
                              std::vector<MyMatrix<T>> const &ListMat,
                              std::vector<T> const &Vdiag, std::ostream &os) {
  using Treturn = WeightMatrix<true, std::vector<T>, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
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
  auto get_wmat=[&]() -> WeightMatrix<true, std::vector<T>, Tidx_value> {
    bool is_symm = is_family_symmmetric(ListMat);
    if (is_symm) {
      ListMatSymm_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
      auto f1 = [&](size_t iRow) -> void {
        lms.f1(iRow);
      };
      auto f2 = [&](size_t jRow) -> std::vector<T> {
        return lms.f2(jRow);
      };
      return WeightMatrix<true, std::vector<T>, Tidx_value>(nbRow, f1, f2, os);
    } else {
      ListMat_Vdiag_WeightMat lms(EXT, ListMat, Vdiag);
      auto f1 = [&](size_t iRow) -> void {
        lms.f1(iRow);
      };
      auto f2 = [&](size_t jRow) -> std::vector<T> {
        return lms.f2(jRow);
      };
      return WeightMatrix<true, std::vector<T>, Tidx_value>(2 * nbRow, f1, f2, os);
    }
  };
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat = get_wmat();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|get_wmat|=" << time << "\n";
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|ReorderingSetWeight|=" << time << "\n";
#endif

  size_t e_hash =
      std::hash<WeightMatrix<true, std::vector<T>, Tidx_value>>()(WMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|hash|=" << time << "\n";
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
  std::cerr << "Failed to find a matching type\n";
  throw TerminalException{1};
}

// The use_scheme boolean means that we used subsets to get
// hopefully faster computation of subsets.
template <typename T, typename Tfield, typename Tidx, bool use_scheme,
          typename Tidx_value>
std::vector<std::vector<Tidx>> GetListGenAutomorphism_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  using Treturn = std::vector<std::vector<Tidx>>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               [[maybe_unused]] auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetStabilizerWeightMatrix_Heuristic<std::vector<T>, Tidx>(
          nbRow, f1, f2, f3, f4, os);
    } else {
      WeightMatrix<true, std::vector<T>, Tidx_value> WMat(nbRow, f1, f2, os);
      return GetStabilizerWeightMatrix_Kernel<std::vector<T>, Tgr, Tidx,
                                              Tidx_value>(WMat, os);
    }
  };
  Treturn ListGen = FCT_ListMat_Vdiag<T, Tfield, Tidx, Treturn, decltype(f)>(
      EXT, ListMat, Vdiag, f, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetListGenAutomorphism_ListMat_Vdiag|=" << time << "\n";
#endif
  return ListGen;
}

template <typename T, typename Tfield, typename Tidx, bool use_scheme>
std::vector<std::vector<Tidx>> GetListGenAutomorphism_ListMat_Vdiag(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_val_poss = nbRow * nbRow / 2 + 1;
  if (max_val_poss < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx,
                                                           use_scheme, uint8_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<
        T, Tfield, Tidx, use_scheme, uint16_t>(EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<
        T, Tfield, Tidx, use_scheme, uint32_t>(EXT, ListMat, Vdiag, os);
  }
  if (max_val_poss < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return GetListGenAutomorphism_ListMat_Vdiag_Tidx_value<
        T, Tfield, Tidx, use_scheme, uint64_t>(EXT, ListMat, Vdiag, os);
  }
  std::cerr << "Failed to find a matching Tidx_value\n";
  throw TerminalException{1};
}

// The use_scheme means that we use subsets in order to compute
// the canonical form. It may or may not work as expected.
template <typename T, typename Tfield, typename Tidx, bool use_scheme,
          typename Tidx_value>
std::vector<Tidx> Canonicalization_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT, std::vector<MyMatrix<T>> const &ListMat,
    std::vector<T> const &Vdiag, std::ostream &os) {
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetGroupCanonicalizationVector_Heuristic<std::vector<T>, Tidx>(
                 nbRow, f1, f2, f3, f4, f5, os)
          .first;
    } else {
      WeightMatrix<true, std::vector<T>, Tidx_value> WMat(nbRow, f1, f2, os);
      return GetCanonicalizationVector_Kernel<std::vector<T>, Tgr, Tidx,
                                              Tidx_value>(WMat, os);
    }
  };
  Treturn CanonicReord =
      FCT_ListMat_Vdiag<T, Tfield, Tidx, Treturn, decltype(f)>(EXT, ListMat,
                                                               Vdiag, f, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|Canonicalization_ListMat_Vdiag|=" << time << "\n";
#endif
  return CanonicReord;
}

template <typename T, typename Tfield, typename Tidx, bool use_scheme>
std::vector<Tidx>
Canonicalization_ListMat_Vdiag(MyMatrix<T> const &EXT,
                               std::vector<MyMatrix<T>> const &ListMat,
                               std::vector<T> const &Vdiag, std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx,
                                                     use_scheme, uint8_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx,
                                                     use_scheme, uint16_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx,
                                                     use_scheme, uint32_t>(
        EXT, ListMat, Vdiag, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return Canonicalization_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx,
                                                     use_scheme, uint64_t>(
        EXT, ListMat, Vdiag, os);
  }
  std::cerr << "No matching type for Tidx_value\n";
  throw TerminalException{1};
}

template <typename T, typename Tfield, typename Tidx, bool use_scheme,
          typename Tidx_value>
std::optional<std::vector<Tidx>> TestEquivalence_ListMat_Vdiag_Tidx_value(
    MyMatrix<T> const &EXT1, std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<T> const &Vdiag1, MyMatrix<T> const &EXT2,
    std::vector<MyMatrix<T>> const &ListMat2, std::vector<T> const &Vdiag2,
    std::ostream &os) {
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif

  size_t nbRow = EXT1.rows();

  // Different scenario depending on the size
  if (nbRow < 2000) {
    WeightMatrix<true, std::vector<T>, Tidx_value> WMat1 =
        GetWeightMatrix_ListMat_Vdiag<T, Tfield, Tidx, Tidx_value>(
            EXT1, ListMat1, Vdiag1, os);
    WeightMatrix<true, std::vector<T>, Tidx_value> WMat2 =
        GetWeightMatrix_ListMat_Vdiag<T, Tfield, Tidx, Tidx_value>(
            EXT2, ListMat2, Vdiag2, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
    os << "|GetWeightMatrix_ListMatrix_Subset|=" << time << "\n";
#endif

    WMat1.ReorderingSetWeight();
    WMat2.ReorderingSetWeight();
    if (WMat1.GetWeight() != WMat2.GetWeight()) {
      return {};
    }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
    os << "|ReorderingSetWeight|=" << time << "\n";
#endif

    std::optional<std::vector<Tidx>> PairTest =
        TestEquivalenceWeightMatrix_norenorm<std::vector<T>, Tidx, Tidx_value>(
            WMat1, WMat2, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
    os << "|TestEquivalence_ListMat_Vdiag|=" << time << "\n";
#endif
    return PairTest;
  }

  std::vector<Tidx> CanonicReord1 =
      Canonicalization_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(
          EXT1, ListMat1, Vdiag1, os);
  std::vector<Tidx> CanonicReord2 =
      Canonicalization_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(
          EXT2, ListMat2, Vdiag2, os);

  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(EXT1, EXT2, CanonicReord1,
                                                   CanonicReord2);
  if (!IsoInfo)
    return {};
  const MyMatrix<Tfield> &P = IsoInfo->second;
  // Now checking the mapping of matrices
  size_t nMat = ListMat1.size();
  for (size_t iMat = 0; iMat < nMat; iMat++) {
    MyMatrix<Tfield> eMat1 =
        UniversalMatrixConversion<Tfield, T>(ListMat1[iMat]);
    MyMatrix<Tfield> eMat2 =
        UniversalMatrixConversion<Tfield, T>(ListMat2[iMat]);
    MyMatrix<Tfield> eProd = P * eMat1 * TransposedMat(P);
    if (!TestEqualityMatrix(eProd, eMat2)) {
      return {};
    }
  }
  const std::vector<Tidx> &eList = IsoInfo->first;
  Tidx nbRow_tidx = nbRow;
  for (Tidx i1 = 0; i1 < nbRow_tidx; i1++) {
    Tidx i2 = eList[i1];
    if (Vdiag1[i1] != Vdiag2[i2])
      return {};
  }
  return eList;
}

template <typename T, typename Tfield, typename Tidx, bool use_scheme>
std::optional<std::vector<Tidx>> TestEquivalence_ListMat_Vdiag(
    MyMatrix<T> const &EXT1, std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<T> const &Vdiag1, MyMatrix<T> const &EXT2,
    std::vector<MyMatrix<T>> const &ListMat2, std::vector<T> const &Vdiag2,
    std::ostream &os) {
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  size_t nbRow = EXT1.rows();
  size_t max_poss_val = nbRow * nbRow / 2 + 1;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, use_scheme,
                                                    uint8_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, use_scheme,
                                                    uint16_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, use_scheme,
                                                    uint32_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return TestEquivalence_ListMat_Vdiag_Tidx_value<T, Tfield, Tidx, use_scheme,
                                                    uint64_t>(
        EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, os);
  }
  std::cerr << "Failed to find a match for Tidx_value\n";
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
  os << "Before positionZero=" << positionZero << "\n";
#endif
  positionZero = WMat.ReorderingSetWeight_specificPosition(positionZero);
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "After positionZero=" << positionZero << "\n";
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
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
  os << "|GetSimpleWeightMatrixAntipodal|=" << time << "\n";
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
  os << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
#endif

  using Tidx = uint32_t;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair =
      GetGroupCanonicalizationVector_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetGroupCanonicalizationVector_Kernel|=" << time << "\n";
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
  os << "|GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif

  size_t n_cols = EXT.cols();
  MyMatrix<Tint> EXTreord(nbRow, n_cols);
  std::vector<int> ListSigns(nbRow, 0);
  ListSigns[0] = 1;
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  std::string strAssign;
  os << "positionZero=" << WMatAbs.positionZero << "\n";
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
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  Tint eHash2 = MD5_hash_T<Tint>(strAssign);
  os << "strAssign=" << strAssign << "\n";
  os << "eHash2=" << eHash2 << "\n";
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
  Tint eHash3 = MD5_hash_T<Tint>(strWMat);
  os << "eHash3=" << eHash3 << "\n";
#endif
  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int j_row = CanonicOrd[i_row];
    int eSign = ListSigns[i_row];
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      EXTreord(i_row, i_col) = eSign * EXT(j_row, i_col);
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB
  os << "EXTreord=\n";
  WriteMatrix(os, EXTreord);
  WriteMatrixGAP(os, EXTreord);
  os << "\n";
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|EXTreord|=" << time << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|ComputeColHermiteNormalForm|=" << time << "\n";
#endif
  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|SignRenormalizationMatrix|=" << time << "\n";
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
  std::cerr << "Failed to match for Tidx_value\n";
  throw TerminalException{1};
}

template <typename Tint, typename Tidx_value>
MyMatrix<Tint>
LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value(MyMatrix<Tint> const &EXT,
                                                    std::ostream &os) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetQmatrix|=" << time << "\n";
#endif

  std::optional<MyMatrix<Tint>> eEquiv =
      LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(EXT, Qmat, os);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick|=" << time << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetWeightMatrixAntipodal|=" << time << "\n";
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|ReorderingSetWeight|=" << time << "\n";
#endif

  std::vector<int> CanonicOrd =
      GetCanonicalizationVector_Kernel<Tint, GraphBitset, int>(WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif

  MyMatrix<Tint> EXTreord(n_rows, n_cols);
  size_t idx = 0;
  Face IsIncluded(n_rows);
  for (size_t i_row = 0; i_row < 2 * n_rows; i_row++) {
    int j_row = CanonicOrd[i_row];
    int res = j_row % 2;
    int pos = j_row / 2;
    if (res == 0) {
      if (IsIncluded[pos] == 0) {
        IsIncluded[pos] = 1;
        EXTreord.row(idx) = EXT.row(pos);
        idx++;
      }
    } else {
      if (IsIncluded[pos] == 0) {
        IsIncluded[pos] = 1;
        EXTreord.row(idx) = -EXT.row(pos);
        idx++;
      }
    }
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|EXTreord 2|=" << time << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|ComputeColHermiteNormalForm 2|=" << time << "\n";
#endif

  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|SignRenormalizationMatrix|=" << time << "\n";
#endif
  return RedMat;
}

template <typename Tint>
MyMatrix<Tint>
LinPolytopeAntipodalIntegral_CanonicForm(MyMatrix<Tint> const &EXT,
                                         std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value<Tint, uint8_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value<Tint, uint16_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value<Tint, uint32_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value<Tint, uint64_t>(
        EXT, os);
  }
  std::cerr << "Failed to find a matching type for Tidx_value\n";
  throw TerminalException{1};
}

template <typename Tint, typename Tidx_value>
std::optional<std::vector<std::vector<unsigned int>>>
LinPolytopeAntipodalIntegral_Automorphism_AbsTrick_Tidx_value(
    MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &Qmat, std::ostream &os) {
  using Tgr = GraphBitset;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  WeightMatrixAbs<Tint, Tidx_value> WMatAbs =
      GetSimpleWeightMatrixAntipodal_AbsTrick<Tint, Tidx_value>(EXT, Qmat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
#endif

  using Tidx = uint32_t;
  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetStabilizerWeightMatrix_Kernel|=" << time << "\n";
#endif

  // We check if the Generating vector eGen can be mapped from the absolute
  // graph to the original one.
  std::vector<std::vector<unsigned int>> ListGenRet;
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
    std::vector<unsigned int> eGenRet(2 * nbRow, 0);
    auto setSign = [&](int const &idx, uint8_t const &val) -> void {
      if (val == 1) {
        eGenRet[idx] = eGen[idx];
        eGenRet[idx + nbRow] = eGen[idx] + nbRow;
      } else {
        eGenRet[idx] = eGen[idx] + nbRow;
        eGenRet[idx + nbRow] = eGen[idx];
      }
      V[idx] = val;
    };
    setSign(0, 1);
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
                setSign(j, valJ);
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
    ListGenRet.push_back(eGenRet);
    return true;
  };
  auto IsCorrectListGen = [&]() -> bool {
    for (auto &eGen : ListGen) {
      bool test = TestExistSignVector(eGen);
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen())
    return {};
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "Check Generators|=" << time << "\n";
#endif
  //
  std::vector<unsigned int> AntipodalGen(2 * nbRow, 0);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    AntipodalGen[iRow] = iRow + nbRow;
    AntipodalGen[nbRow + iRow] = iRow;
  }
  ListGenRet.push_back(AntipodalGen);
  //
  return ListGenRet;
}

template <typename Tint>
std::optional<std::vector<std::vector<unsigned int>>>
LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(MyMatrix<Tint> const &EXT,
                                                   MyMatrix<Tint> const &Qmat,
                                                   std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_AbsTrick_Tidx_value<
        Tint, uint8_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_AbsTrick_Tidx_value<
        Tint, uint16_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_AbsTrick_Tidx_value<
        Tint, uint32_t>(EXT, Qmat, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_AbsTrick_Tidx_value<
        Tint, uint64_t>(EXT, Qmat, os);
  }
  std::cerr << "Failed to find a matching type for Tidx_value\n";
  throw TerminalException{1};
}

template <typename Tint, typename Tidx_value>
std::vector<std::vector<unsigned int>>
LinPolytopeAntipodalIntegral_Automorphism_Tidx_value(MyMatrix<Tint> const &EXT,
                                                     std::ostream &os) {
  using Tidx = uint32_t;
  using Tgr = GraphBitset;
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  SecondTime time;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetQmatrix|=" << time << "\n";
#endif

  std::optional<std::vector<std::vector<unsigned int>>> eEquiv =
      LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(EXT, Qmat, os);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|LinPolytopeAntipodalIntegral_Automorphism|=" << time << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetWeightMatrixAntipodal|=" << time << "\n";
#endif

  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB
  os << "|GetStabilizerWeightMatrix_Kernel|=" << time << "\n";
#endif
  return ListGen;
}

template <typename Tint>
std::vector<std::vector<unsigned int>>
LinPolytopeAntipodalIntegral_Automorphism(MyMatrix<Tint> const &EXT,
                                          std::ostream &os) {
  size_t nbRow = EXT.rows();
  size_t max_poss_val = nbRow * nbRow;
  if (max_poss_val < size_t(std::numeric_limits<uint8_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_Tidx_value<Tint, uint8_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint16_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_Tidx_value<Tint, uint16_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint32_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_Tidx_value<Tint, uint32_t>(
        EXT, os);
  }
  if (max_poss_val < size_t(std::numeric_limits<uint64_t>::max() - 1)) {
    return LinPolytopeAntipodalIntegral_Automorphism_Tidx_value<Tint, uint64_t>(
        EXT, os);
  }
  std::cerr << "Failed to find a matching type for Tidx_value\n";
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
                         std::vector<MyMatrix<T>> const &ListComm, std::ostream& os) {
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
                           MyMatrix<T> const &TheEXT, std::ostream& os) {
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

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::vector<MyMatrix<T>> LinPolytopeIntegralWMat_Automorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep,
    std::ostream &os) {
  using Tgr = GraphBitset;
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second, os);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|GRP1|=" << GRP1.size() << " RankMat(ep.first)=" << RankMat(ep.first)
     << " |ep.first|=" << ep.first.rows() << " / " << ep.first.cols() << "\n";
  bool test = CheckStabilizerWeightMatrix(ep.second, GRP1);
  os << "test=" << test << "\n";
#endif
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer_Method8(ep.first, GRP1, os);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "We have GRPfull\n";
#endif
  std::vector<MyMatrix<T>> ListGenMat;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    MyMatrix<T> eMat_T = FindTransformation(ep.first, ep.first, eGen);
    ListGenMat.push_back(eMat_T);
  }
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "We have ListGenMat\n";
#endif
  return ListGenMat;
}

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::optional<MyMatrix<T>> LinPolytopeIntegralWMat_Isomorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep,
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &fp,
    [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphBitset;
  if (ep.first.rows() != fp.first.rows() || ep.first.cols() != fp.first.cols())
    return {};
  if (ep.second.GetWeight() != fp.second.GetWeight())
    return {};
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|ep.first|=" << ep.first.rows() << " / " << ep.first.cols()
     << " rnk=" << RankMat(ep.first) << "\n";
  os << "|fp.first|=" << fp.first.rows() << " / " << fp.first.cols()
     << " rnk=" << RankMat(fp.first) << "\n";
  os << "ep.first=\n";
  WriteMatrix(os, ep.first);
  os << "fp.first=\n";
  WriteMatrix(os, fp.first);
  os << "ep.second=\n";
  PrintWeightedMatrix(os, ep.second);
  os << "fp.second=\n";
  PrintWeightedMatrix(os, fp.second);
#endif
#ifdef TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
  SecondTime time;
#endif

  std::vector<Tidx> eCanonicReord =
      GetGroupCanonicalizationVector_Kernel<Tval, Tgr, Tidx, Tidx_value>(
          ep.second, os)
          .first;
  std::vector<Tidx> fCanonicReord =
      GetGroupCanonicalizationVector_Kernel<Tval, Tgr, Tidx, Tidx_value>(
          fp.second, os)
          .first;
#ifdef TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|GetGroupCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(
          ep.first, fp.first, eCanonicReord, fCanonicReord);
#ifdef TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|IsomorphismFromCanonicReord|=" << time << "\n";
#endif
  if (!IsoInfo) {
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
    os << "We failed to find IsoInfo\n";
#endif
    return {};
  }
  Telt ePerm(IsoInfo->first);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "ePerm=" << ePerm << "\n";
  os << "det(eMat)=" << DeterminantMat(IsoInfo->second)
     << "  eMat=" << StringMatrixGAP(IsoInfo->second) << "\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second, os);
#ifdef TIMINGS_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|GetStabilizerWeightMatrix|=" << time << "\n";
#endif
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "|GRP1|=" << GRP1.size() << "\n";
  for (auto &eGen : GRP1.GeneratorsOfGroup()) {
    MyMatrix<T> eGen_T = FindTransformation(ep.first, ep.first, eGen);
    os << "det(eGen_T)=" << DeterminantMat(eGen_T)
       << " eGen_T=" << StringMatrixGAP(eGen_T) << "\n";
  }
#endif
  std::optional<MyMatrix<T>> eRes = LinPolytopeIntegral_Isomorphism_Method8(
      ep.first, fp.first, GRP1, ePerm, os);
#ifdef TIMINGS
  os << "|LinPolytopeIntegral_Isomorphism_Method8|=" << time << "\n";
#endif
  if (eRes) {
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
    os << "Found one isomorphism\n";
#endif
    return *eRes;
  }
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  os << "eRes is unassigned\n";
#endif
  return {};
}

/*
  The integral canonical form is needed for many application as it makes things
  faster. What can we do:
  ---We can find a canonical ordering of the vertices of the polytope.
  ---From the ordering we can get from GetZbasis (or a variant of it).
  ---That basis M is of determinant d, this gives us d classes.
  ---The group Lin(P) acts on the different embeddings.
  ---The group acts Lin(P) acts on all the embeddings.
  ---Therefore, we can find the canonical set from the set.
  ---From that, we can get the M + embedding. We can canonicalize it.
  ---This gets us an integral transformation.
  ---But this integral transformation does not belong to the group of isometries
  ---Can we handle that? Maybe.
 */
template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
MyMatrix<T> LinPolytopeIntegralWMat_Canonic(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep) {
  std::cerr << "The code is incomplete\n";
  throw TerminalException{1};
  return ep.first;
}

// clang-format off
#endif  // SRC_GROUP_TEMP_POLYTOPEEQUISTAB_H_
// clang-format on
