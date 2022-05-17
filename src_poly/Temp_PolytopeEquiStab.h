#ifndef SRC_POLY_TEMP_POLYTOPEEQUISTAB_H_
#define SRC_POLY_TEMP_POLYTOPEEQUISTAB_H_

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

//
// Equivalence of subsets and stabilizer of a WeightMatrix
//

template <typename T, typename Telt, typename Tidx_value>
std::optional<Telt>
TestEquivalenceSubset(WeightMatrix<true, T, Tidx_value> const &WMat,
                      Face const &f1, Face const &f2) {
  using Tidx = typename Telt::Tidx;
  size_t siz = WMat.GetWeightSize();
  size_t n = WMat.rows();
  auto g = [&](Face const &f, size_t iRow, size_t iCol) -> int {
    if (iRow < n && iCol < n)
      return WMat.GetValue(iRow, iCol);
    if (iRow == n && iCol == n)
      return siz + 2;
    if (iRow == n) { // Thus iCol < n.
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
      n + 1,
      [&](size_t iRow, size_t iCol) -> int { return g(f1, iRow, iCol); });
  WeightMatrix<true, int, Tidx_value> WMat2(
      n + 1,
      [&](size_t iRow, size_t iCol) -> int { return g(f2, iRow, iCol); });
  std::optional<Telt> test =
      TestEquivalenceWeightMatrix_norenorm_perm<int, Telt>(WMat1, WMat2);
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
                        Face const &f) {
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
  WeightMatrix<true, int, Tidx_value> WMatW(n + 1, g);
  Tgroup GRP = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMatW);
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
      T eSum = 0;
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
GetQmatrix(MyMatrix<T> const &TheEXT) {
  return Kernel_GetQmatrix(TheEXT);
}

template <typename T>
inline typename std::enable_if<(not is_ring_field<T>::value), MyMatrix<T>>::type
GetQmatrix(MyMatrix<T> const &TheEXT) {
  using Tfield = typename overlying_field<T>::field_type;
#ifdef TIMINGS
  SingletonTime time1;
#endif

  MyMatrix<Tfield> TheEXT_F = UniversalMatrixConversion<Tfield, T>(TheEXT);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|UniversalMatrixConversion1|=" << ms(time1, time2) << "\n";
#endif

  MyMatrix<Tfield> Q_F = Kernel_GetQmatrix(TheEXT_F);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|Kernel_GetQmatrix|=" << ms(time2, time3) << "\n";
#endif

  MyMatrix<Tfield> Q_F_red = RemoveFractionMatrix(Q_F);
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|RemoveFractionMatrix|=" << ms(time3, time4) << "\n";
#endif

  MyMatrix<T> RetMat = UniversalMatrixConversion<T, Tfield>(Q_F_red);
#ifdef TIMINGS
  SingletonTime time5;
  std::cerr << "|UniversalMatrixConversion2|=" << ms(time4, time5) << "\n";
#endif
  return RetMat;
}

template <typename T>
MyMatrix<T> GetQmatrix_NotFullRank(MyMatrix<T> const &TheEXT) {
  std::vector<int> eList = ColumnReductionSet(TheEXT);
  int len = eList.size();
  int n_cols = TheEXT.cols();
  MyMatrix<T> TheEXT_red = SelectColumn(TheEXT, eList);
  MyMatrix<T> Qmat_red = GetQmatrix(TheEXT_red);
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
      T eSum = 0;
      for (size_t jCol = 0; jCol < nbCol; jCol++)
        eSum += Qinput(iCol, jCol) * TheEXT(iRow, jCol);
      V(iCol) = eSum;
    }
  };
  auto f2 = [&](size_t jRow) -> T {
    T eSum = 0;
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
Treturn FCT_EXT_Qinv(MyMatrix<T> const &TheEXT, F f) {
  MyMatrix<T> Qmat = GetQmatrix(TheEXT);
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
  //  std::cerr << "IsomorphismFromCanonicReord : CanonicReord1=" <<
  //  CanonicReord1 << "\n"; std::cerr << "IsomorphismFromCanonicReord :
  //  CanonicReord2=" << CanonicReord2 << "\n"; std::cerr <<
  //  "IsomorphismFromCanonicReord : ListIdx=" << ListIdx << "\n";

  // Building the matrix equivalence
  MyMatrix<Tfield> Basis1 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT1, CanonicReord1);
  MyMatrix<Tfield> Basis2 =
      GetBasisFromOrdering<T, Tfield, Tidx>(EXT2, CanonicReord2);
  MyMatrix<Tfield> P = Inverse(Basis1) * Basis2;
  //  std::cerr << "P=\n";
  //  WriteMatrix(std::cerr, P);
  // Now testing the obtained mappings
  bool test = CheckEquivalence(EXT1, EXT2, ListIdx, P);
  if (!test) // We fail the polytope equivalence
    return {};
  std::pair<std::vector<Tidx>, MyMatrix<Tfield>> IsoInfo{ListIdx, P};
  return IsoInfo;
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetSimpleWeightMatrix(MyMatrix<T> const &TheEXT, MyMatrix<T> const &Qinput) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
    return WeightMatrix<true, T, Tidx_value>(nbRow, f1, f2);
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
WeightMatrix<true, T, Tidx_value> GetWeightMatrix(MyMatrix<T> const &TheEXT) {
  using Treturn = WeightMatrix<true, T, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
    return WeightMatrix<true, T, Tidx_value>(nbRow, f1, f2);
  };
  //
  size_t n_rows = TheEXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max())) {
    using Tidx = uint8_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f);
  }
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max())) {
    using Tidx = uint16_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f);
  }
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max())) {
    using Tidx = uint32_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f);
  }
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max())) {
    using Tidx = uint64_t;
    return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f);
  }
#endif
  std::cerr << "Failed to find matching numeric in GetWeightMatrix\n";
  throw TerminalException{1};
}

template <typename T>
WeightMatrixLimited<true, T> GetWeightMatrixLimited(MyMatrix<T> const &TheEXT,
                                                    size_t max_offdiag) {
  using Treturn = WeightMatrixLimited<true, T>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
    return WeightMatrixLimited<true, T>(nbRow, f1, f2, max_offdiag);
  };
  using Tidx = size_t;
  return FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(TheEXT, f);
}

template <typename T, bool use_scheme, typename Tgroup>
Tgroup LinPolytope_Automorphism(MyMatrix<T> const &EXT) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphListAdj;
  using Tidx_value = uint16_t;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS
  SingletonTime time1;
#endif
  MyMatrix<T> EXTred = ColumnReduction(EXT);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|LinPolytope_Aut : ColumnReduction|=" << ms(time1, time2)
            << "\n";
#endif
  using Treturn = std::vector<std::vector<Tidx>>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               [[maybe_unused]] auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetStabilizerWeightMatrix_Heuristic<T, Tidx>(nbRow, f1, f2, f3,
                                                          f4);
    } else {
      WeightMatrix<true, T, Tidx_value> WMat(nbRow, f1, f2);
      return GetStabilizerWeightMatrix_Kernel<T, Tgr, Tidx, Tidx_value>(WMat);
    }
  };
  Treturn ListGen = FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(EXTred, f);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|LinPolytope_Aut : FCT_EXT_Qinv|=" << ms(time2, time3) << "\n";
#endif
  std::vector<Telt> LGen;
  //  std::cerr << "nbRow=" << nbRow << " |ListGen|=" << ListGen.size() << "\n";
  for (auto &eList : ListGen) {
    //    std::cerr << "|eList|=" << eList.size() << "\n";
    LGen.push_back(Telt(eList));
  }
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|LinPolytope_Aut : LGen|=" << ms(time3, time4) << "\n";
#endif
  return Tgroup(LGen, nbRow);
}

template <typename T, typename Tidx, bool use_scheme>
std::vector<Tidx> LinPolytope_CanonicOrdering(MyMatrix<T> const &EXT) {
  using Tidx_value = uint16_t;
  using Tgr = GraphBitset;
#ifdef TIMINGS
  SingletonTime time1;
#endif

  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetGroupCanonicalizationVector_Heuristic<T, Tidx>(nbRow, f1, f2,
                                                               f3, f4, f5)
          .first;
    } else {
      WeightMatrix<true, T, Tidx_value> WMat(nbRow, f1, f2);
      WMat.ReorderingSetWeight();
      return GetGroupCanonicalizationVector_Kernel<T, Tgr, Tidx, Tidx_value>(
                 WMat)
          .first;
    }
  };
  std::vector<Tidx> CanonicOrd =
      FCT_EXT_Qinv<T, Tidx, Treturn, decltype(f)>(EXT, f);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|FCT_EXT_Qinv|=" << ms(time1, time2) << "\n";
#endif
  return CanonicOrd;
}

template <typename T, bool use_scheme, typename Tidx>
MyMatrix<T> LinPolytope_CanonicForm_Tidx(MyMatrix<T> const &EXT) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS
  SingletonTime time1;
#endif
  std::vector<Tidx> CanonicOrd =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT);
  MyMatrix<T> EXTreord(n_rows, n_cols);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    size_t j_row = CanonicOrd[i_row];
    EXTreord.row(i_row) = EXT.row(j_row);
  }
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|CanonicOrdering + EXTreord|=" << ms(time1, time2) << "\n";
#endif

  //  MyMatrix<T> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
  MyMatrix<T> RedMat = CanonicalizeOrderedMatrix(EXTreord);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|CanonicalizeOrderedMatrix|=" << ms(time2, time3) << "\n";
#endif
  return RedMat;
}

template <typename T, bool use_scheme>
MyMatrix<T> LinPolytope_CanonicForm(MyMatrix<T> const &EXT) {
  size_t n_rows = EXT.rows();
  if (n_rows < size_t(std::numeric_limits<uint8_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint8_t>(EXT);
  if (n_rows < size_t(std::numeric_limits<uint16_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint16_t>(EXT);
  if (n_rows < size_t(std::numeric_limits<uint32_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint32_t>(EXT);
#if !defined __APPLE__
  if (n_rows < size_t(std::numeric_limits<uint64_t>::max()))
    return LinPolytope_CanonicForm_Tidx<T, use_scheme, uint64_t>(EXT);
#endif
  std::cerr << "LinPolytope_CanonicForm : Failed to find matching numeric\n";
  throw TerminalException{1};
}

template <typename T, typename Tidx, bool use_scheme>
std::optional<std::vector<Tidx>>
LinPolytope_Isomorphism(const MyMatrix<T> &EXT1, const MyMatrix<T> &EXT2) {
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT1);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering<T, Tidx, use_scheme>(EXT2);
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(EXT1, EXT2, CanonicReord1,
                                                   CanonicReord2);
  if (!IsoInfo)
    return {};
  return IsoInfo->first;
}

template <typename Tint, typename Tidx, typename Tgroup, typename Tidx_value,
          typename Tgr, bool use_scheme>
std::optional<MyMatrix<Tint>>
LinPolytopeIntegral_Isomorphism(const MyMatrix<Tint> &EXT1,
                                const MyMatrix<Tint> &EXT2) {
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering<Tint, Tidx, use_scheme>(EXT1);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering<Tint, Tidx, use_scheme>(EXT2);
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
  Tgroup GRP1 = LinPolytope_Automorphism<Tfield, use_scheme, Tgroup>(EXT1_T);
  std::optional<MyMatrix<Tfield>> eRes =
      LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, ePerm);
  if (eRes)
    return UniversalMatrixConversion<Tint, Tfield>(*eRes);
  return {};
}

//
// The Lorentzian case that cause us so much trouble.
//

// ---ListMat is assumed to be symmetric
// ---Note that TheEXT does not have to be of full rank.
//    It makes perfect sense to compute some group
//    and get it only as permutation group.
template <typename T, typename Tidx, typename Treturn, typename F>
Treturn FCT_ListMat_Vdiag(MyMatrix<T> const &TheEXT,
                          std::vector<MyMatrix<T>> const &ListMat,
                          std::vector<T> const &Vdiag, F f) {
  using Tfield = typename overlying_field<T>::field_type;
#ifdef SANITY_CHECK
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      std::cerr << "The matrix eMat should be symmetric\n";
      throw TerminalException{1};
    }
  }
#endif
  size_t nbRow = TheEXT.rows();
  size_t max_val = std::numeric_limits<Tidx>::max();
  if (nbRow > max_val) {
    std::cerr
        << "Error in FCT_ListMat_Vdiag due to too small coefficient range\n";
    std::cerr << "nbRow=" << nbRow
              << " std::numeric_limits<Tidx>::max()=" << max_val << "\n";
    throw TerminalException{1};
  }
  size_t nbCol = TheEXT.cols();
  size_t nMat = ListMat.size();
  std::vector<MyMatrix<Tfield>> ListMat_F;
  for (auto &eMat : ListMat) {
    MyMatrix<Tfield> eMat_F = UniversalMatrixConversion<Tfield, T>(eMat);
    ListMat_F.push_back(eMat_F);
  }
  //
  MyMatrix<T> MatV(nMat, nbCol);
  std::vector<T> LScal(nMat + 1);
  size_t iRow_stor = 0;
  auto f1 = [&](size_t iRow) -> void {
    for (size_t iMat = 0; iMat < nMat; iMat++) {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        T eSum = 0;
        for (size_t jCol = 0; jCol < nbCol; jCol++)
          eSum += ListMat[iMat](jCol, iCol) * TheEXT(iRow, jCol);
        MatV(iMat, iCol) = eSum;
      }
    }
    iRow_stor = iRow;
  };
  auto f2 = [&](size_t jRow) -> std::vector<T> {
    for (size_t iMat = 0; iMat < nMat; iMat++) {
      T eSum = 0;
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eSum += MatV(iMat, iCol) * TheEXT(jRow, iCol);
      LScal[iMat] = eSum;
    }
    T eVal = 0;
    if (iRow_stor == jRow)
      eVal = Vdiag[jRow];
    LScal[nMat] = eVal;
    return LScal;
  };
  // Preemptive check that the subset is adequate
  auto f3 = [&](std::vector<Tidx> const &Vsubset) -> bool {
    return IsSubsetFullRank<T, Tfield, Tidx>(TheEXT, Vsubset);
  };
  // Extension of the partial automorphism
  auto f4 = [&](const std::vector<Tidx> &Vsubset, const std::vector<Tidx> &Vin,
                const std::vector<std::vector<Tidx>> &ListBlocks)
      -> DataMapping<Tidx> {
#ifdef DEBUG_TEMP_POLYTOPE_EQUI_STAB
    std::cerr << "Before FindMatrixTransformationTest_Subset\n";
#endif
    std::optional<MyMatrix<Tfield>> test1 =
        FindMatrixTransformationTest_Subset<T, Tfield, Tidx>(TheEXT, Vsubset,
                                                             Vin);
#ifdef DEBUG_TEMP_POLYTOPE_EQUI_STAB
    std::cerr << "After test1=" << test1.has_value() << "\n";
#endif
    Face block_status(ListBlocks.size());
    if (!test1) {
#ifdef DEBUG_TEMP_POLYTOPE_EQUI_STAB
      std::cerr << "f4 exit false 1\n";
#endif
      return {false, block_status, {}};
    }
    const MyMatrix<Tfield> &P = *test1;
    for (auto &eMat_F : ListMat_F) {
      MyMatrix<Tfield> eProd = P * eMat_F * TransposedMat(P);
      if (!TestEqualityMatrix(eProd, eMat_F)) {
#ifdef DEBUG_TEMP_POLYTOPE_EQUI_STAB
        std::cerr << "f4 exit false 2\n";
#endif
        return {false, block_status, {}};
      }
    }
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

template <typename T, typename Tidx, typename Tidx_value>
WeightMatrix<true, std::vector<T>, Tidx_value>
GetWeightMatrix_ListMat_Vdiag(MyMatrix<T> const &TheEXT,
                              std::vector<MyMatrix<T>> const &ListMat,
                              std::vector<T> const &Vdiag) {
  using Treturn = WeightMatrix<true, std::vector<T>, Tidx_value>;
  auto f = [&](size_t nbRow, auto f1, auto f2, [[maybe_unused]] auto f3,
               [[maybe_unused]] auto f4, [[maybe_unused]] auto f5) -> Treturn {
    return WeightMatrix<true, std::vector<T>, Tidx_value>(nbRow, f1, f2);
  };
  return FCT_ListMat_Vdiag<T, Tidx, Treturn, decltype(f)>(TheEXT, ListMat,
                                                          Vdiag, f);
}

template <typename T>
size_t GetInvariant_ListMat_Vdiag(MyMatrix<T> const &EXT,
                                  std::vector<MyMatrix<T>> const &ListMat,
                                  std::vector<T> const &Vdiag) {
  using Tidx_value = uint16_t;
  using Tidx = unsigned int;
#ifdef TIMINGS
  SingletonTime time1;
#endif

  WeightMatrix<true, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListMat_Vdiag<T, Tidx, Tidx_value>(EXT, ListMat, Vdiag);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << ms(time1, time2)
            << "\n";
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|ReorderingSetWeight|=" << ms(time2, time3) << "\n";
#endif

  size_t e_hash =
      std::hash<WeightMatrix<true, std::vector<T>, Tidx_value>>()(WMat);
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|hash|=" << ms(time3, time4) << "\n";
#endif
  return e_hash;
}

template <typename T, typename Tidx, bool use_scheme>
std::vector<std::vector<Tidx>>
GetListGenAutomorphism_ListMat_Vdiag(MyMatrix<T> const &EXT,
                                     std::vector<MyMatrix<T>> const &ListMat,
                                     std::vector<T> const &Vdiag) {
  using Tidx_value = uint16_t;
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
#ifdef TIMINGS
  SingletonTime time1;
#endif
  using Treturn = std::vector<std::vector<Tidx>>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               [[maybe_unused]] auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetStabilizerWeightMatrix_Heuristic<std::vector<T>, Tidx>(
          nbRow, f1, f2, f3, f4);
    } else {
      WeightMatrix<true, std::vector<T>, Tidx_value> WMat(nbRow, f1, f2);
      return GetStabilizerWeightMatrix_Kernel<std::vector<T>, Tgr, Tidx,
                                              Tidx_value>(WMat);
    }
  };
  Treturn ListGen =
      FCT_ListMat_Vdiag<T, Tidx, Treturn, decltype(f)>(EXT, ListMat, Vdiag, f);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetListGenAutomorphism_ListMat_Vdiag|=" << ms(time1, time2)
            << "\n";
#endif
  return ListGen;
}

template <typename T, typename Tidx, bool use_scheme>
std::vector<Tidx>
Canonicalization_ListMat_Vdiag(MyMatrix<T> const &EXT,
                               std::vector<MyMatrix<T>> const &ListMat,
                               std::vector<T> const &Vdiag) {
  using Tidx_value = uint16_t;
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
#ifdef TIMINGS
  SingletonTime time1;
#endif
  using Treturn = std::vector<Tidx>;
  auto f = [&](size_t nbRow, auto f1, auto f2, auto f3, auto f4,
               auto f5) -> Treturn {
    if constexpr (use_scheme) {
      return GetGroupCanonicalizationVector_Heuristic<std::vector<T>, Tidx>(
                 nbRow, f1, f2, f3, f4, f5)
          .first;
    } else {
      WeightMatrix<true, std::vector<T>, Tidx_value> WMat(nbRow, f1, f2);
      return GetCanonicalizationVector_Kernel<std::vector<T>, Tgr, Tidx,
                                              Tidx_value>(WMat);
    }
  };
  Treturn CanonicReord =
      FCT_ListMat_Vdiag<T, Tidx, Treturn, decltype(f)>(EXT, ListMat, Vdiag, f);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|Canonicalization_ListMat_Vdiag|=" << ms(time1, time2) << "\n";
#endif
  return CanonicReord;
}

template <typename T, typename Tidx, bool use_scheme>
std::optional<std::vector<Tidx>> TestEquivalence_ListMat_Vdiag(
    MyMatrix<T> const &EXT1, std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<T> const &Vdiag1, MyMatrix<T> const &EXT2,
    std::vector<MyMatrix<T>> const &ListMat2, std::vector<T> const &Vdiag2) {
  using Tidx_value = uint16_t;
  using Tfield = typename overlying_field<T>::field_type;
#ifdef TIMINGS
  SingletonTime time1;
#endif

  size_t nbRow1 = EXT1.rows();
  size_t nbRow2 = EXT2.rows();
  if (nbRow1 !=
      nbRow2) { // At this point it should be equal, but better to check
    return {};
  }

  // Different scenario depending on the size
  if (nbRow1 < 2000) {
    WeightMatrix<true, std::vector<T>, Tidx_value> WMat1 =
        GetWeightMatrix_ListMat_Vdiag<T, Tidx, Tidx_value>(EXT1, ListMat1,
                                                           Vdiag1);
    WeightMatrix<true, std::vector<T>, Tidx_value> WMat2 =
        GetWeightMatrix_ListMat_Vdiag<T, Tidx, Tidx_value>(EXT2, ListMat2,
                                                           Vdiag2);
#ifdef TIMINGS
    SingletonTime time2;
    std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << ms(time1, time2)
              << "\n";
#endif

    WMat1.ReorderingSetWeight();
    WMat2.ReorderingSetWeight();
#ifdef TIMINGS
    SingletonTime time3;
    std::cerr << "|ReorderingSetWeight|=" << ms(time2, time3) << "\n";
#endif

    std::optional<std::vector<Tidx>> PairTest =
        TestEquivalenceWeightMatrix_norenorm<std::vector<T>, Tidx, Tidx_value>(
            WMat1, WMat2);
#ifdef TIMINGS
    SingletonTime time4;
    std::cerr << "|TestEquivalence_ListMat_Vdiag|=" << ms(time3, time4) << "\n";
#endif
    return PairTest;
  }

  std::vector<Tidx> CanonicReord1 =
      Canonicalization_ListMat_Vdiag<T, Tidx, use_scheme>(EXT1, ListMat1,
                                                          Vdiag1);
  std::vector<Tidx> CanonicReord2 =
      Canonicalization_ListMat_Vdiag<T, Tidx, use_scheme>(EXT2, ListMat2,
                                                          Vdiag2);

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
  Tidx nbRow1_tidx = nbRow1;
  for (Tidx i1 = 0; i1 < nbRow1_tidx; i1++) {
    Tidx i2 = eList[i1];
    if (Vdiag1[i1] != Vdiag2[i2])
      return {};
  }
  return eList;
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
WeightMatrixAbs<T, Tidx_value>
GetSimpleWeightMatrixAntipodal_AbsTrick(MyMatrix<T> const &TheEXT,
                                        MyMatrix<T> const &Qmat) {
  static_assert(is_totally_ordered<T>::value,
                "Requires T to be a totally ordered field");
#ifdef TIMINGS
  SingletonTime time1;
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
      if (value == 0) { // This is a missing value
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
                                         weight_ordered);
#ifdef DEBUG
  std::cerr << "Before positionZero=" << positionZero << "\n";
#endif
  positionZero = WMat.ReorderingSetWeight_specificPosition(positionZero);
#ifdef DEBUG
  std::cerr << "After positionZero=" << positionZero << "\n";
#endif
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << ms(time1, time2)
            << "\n";
#endif
  return {positionZero, std::move(ArrSigns), std::move(WMat)};
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetSimpleWeightMatrixAntipodal(MyMatrix<T> const &TheEXT,
                               MyMatrix<T> const &Qmat) {
#ifdef TIMINGS
  SingletonTime time1;
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
      if (value1 == 0) { // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eSum1);
      }
      Tidx_value &value2 = ValueMap[eSum2];
      if (value2 == 0) { // This is a missing value
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
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetSimpleWeightMatrixAntipodal|=" << ms(time1, time2) << "\n";
#endif
  bool weight_ordered = false;
  return WeightMatrix<true, T, Tidx_value>(INP_nbRow, INP_TheMat,
                                           INP_ListWeight, weight_ordered);
}

template <typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
GetWeightMatrixAntipodal(MyMatrix<T> const &TheEXT) {
  MyMatrix<T> Qmat = GetQmatrix(TheEXT);
  return GetSimpleWeightMatrixAntipodal<T, Tidx_value>(TheEXT, Qmat);
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
template <typename Tint>
std::optional<MyMatrix<Tint>>
LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(MyMatrix<Tint> const &EXT,
                                                  MyMatrix<Tint> const &Qmat) {
  using Tidx_value = uint16_t;
  using Tgr = GraphBitset;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS
  SingletonTime time1;
#endif
  WeightMatrixAbs<Tint, Tidx_value> WMatAbs =
      GetSimpleWeightMatrixAntipodal_AbsTrick<Tint, Tidx_value>(EXT, Qmat);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << ms(time1, time2)
            << "\n";
#endif

  using Tidx = unsigned int;
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair =
      GetGroupCanonicalizationVector_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|GetGroupCanonicalizationVector_Kernel|=" << ms(time2, time3)
            << "\n";
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
              bool ChgSign =
                  ChgSign1 ^ ChgSign2; // true if ChgSign1 != ChgSign2
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
      //      std::cerr << "test=" << test << "\n";
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen())
    return {};
  //
  std::vector<Tidx> const &CanonicOrd = ePair.first;
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|GetCanonicalizationVector_Kernel|=" << ms(time3, time4)
            << "\n";
#endif

  size_t n_cols = EXT.cols();
  MyMatrix<Tint> EXTreord(nbRow, n_cols);
  std::vector<int> ListSigns(nbRow, 0);
  ListSigns[0] = 1;
#ifdef DEBUG
  std::string strAssign;
  std::cerr << "positionZero=" << WMatAbs.positionZero << "\n";
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
          int ValSign = 1 - 2 * int(ChgSign);
          int RetSign = ValSign * ListSigns[k_row];
          ListSigns[i_row] = RetSign;
#ifdef DEBUG
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
  };
#ifdef DEBUG
  mpz_class eHash2 = MD5_hash_mpz(strAssign);
  std::cerr << "strAssign=" << strAssign << "\n";
  std::cerr << "eHash2=" << eHash2 << "\n";
#endif
#ifdef DEBUG
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
  mpz_class eHash3 = MD5_hash_mpz(strWMat);
  std::cerr << "eHash3=" << eHash3 << "\n";
#endif

  for (size_t i_row = 0; i_row < nbRow; i_row++) {
    int j_row = CanonicOrd[i_row];
    int eSign = ListSigns[i_row];
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      EXTreord(i_row, i_col) = eSign * EXT(j_row, i_col);
  }
#ifdef DEBUG
  std::cerr << "EXTreord=\n";
  WriteMatrix(std::cerr, EXTreord);
  WriteMatrixGAP(std::cerr, EXTreord);
  std::cerr << "\n";
#endif
#ifdef TIMINGS
  SingletonTime time5;
  std::cerr << "|EXTreord|=" << ms(time4, time5) << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS
  SingletonTime time6;
  std::cerr << "|ComputeColHermiteNormalForm|=" << ms(time5, time6) << "\n";
#endif
  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS
  SingletonTime time7;
  std::cerr << "|SignRenormalizationMatrix|=" << ms(time6, time7) << "\n";
#endif
  return RedMat;
}

template <typename Tint>
MyMatrix<Tint>
LinPolytopeAntipodalIntegral_CanonicForm(MyMatrix<Tint> const &EXT) {
  using Tidx_value = uint16_t;
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS
  SingletonTime time1;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetQmatrix|=" << ms(time1, time2) << "\n";
#endif

  std::optional<MyMatrix<Tint>> eEquiv =
      LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(EXT, Qmat);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick|="
            << ms(time2, time3) << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT);
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|GetWeightMatrixAntipodal|=" << ms(time3, time4) << "\n";
#endif
  //  std::cerr << "After direct construction WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMat);

  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  SingletonTime time5;
  std::cerr << "|ReorderingSetWeight|=" << ms(time4, time5) << "\n";
#endif

  std::vector<int> CanonicOrd =
      GetCanonicalizationVector_Kernel<Tint, GraphBitset, int>(WMat);
#ifdef TIMINGS
  SingletonTime time6;
  std::cerr << "|GetCanonicalizationVector_Kernel|=" << ms(time5, time6)
            << "\n";
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
#ifdef TIMINGS
  SingletonTime time7;
  std::cerr << "|EXTreord 2|=" << ms(time6, time7) << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS
  SingletonTime time8;
  std::cerr << "|ComputeColHermiteNormalForm 2|=" << ms(time7, time8) << "\n";
#endif

  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS
  SingletonTime time9;
  std::cerr << "|SignRenormalizationMatrix|=" << ms(time8, time9) << "\n";
#endif
  return RedMat;
}

template <typename Tint>
std::optional<std::vector<std::vector<unsigned int>>>
LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(MyMatrix<Tint> const &EXT,
                                                   MyMatrix<Tint> const &Qmat) {
  using Tidx_value = uint16_t;
  using Tgr = GraphBitset;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS
  SingletonTime time1;
#endif
  WeightMatrixAbs<Tint, Tidx_value> WMatAbs =
      GetSimpleWeightMatrixAntipodal_AbsTrick<Tint, Tidx_value>(EXT, Qmat);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << ms(time1, time2)
            << "\n";
#endif
  //  std::cerr << "WMatAbs.positionZero=" << WMatAbs.positionZero << "\n";
  //  std::cerr << "WMatAbs.WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMatAbs.WMat);

  using Tidx = unsigned int;
  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|GetStabilizerWeightMatrix_Kernel|=" << ms(time2, time3)
            << "\n";
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
              bool ChgSign =
                  ChgSign1 ^ ChgSign2; // true if ChgSign1 != ChgSign2
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
      //      std::cerr << "test=" << test << "\n";
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen())
    return {};
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|Check Generators|=" << ms(time3, time4) << "\n";
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
std::vector<std::vector<unsigned int>>
LinPolytopeAntipodalIntegral_Automorphism(MyMatrix<Tint> const &EXT) {
  using Tidx = unsigned int;
  using Tidx_value = uint16_t;
  using Tgr = GraphBitset;
#ifdef TIMINGS
  SingletonTime time1;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetQmatrix|=" << ms(time1, time2) << "\n";
#endif

  std::optional<std::vector<std::vector<unsigned int>>> eEquiv =
      LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(EXT, Qmat);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|LinPolytopeAntipodalIntegral_Automorphism_AbsTrick|="
            << ms(time2, time3) << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT);
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|GetWeightMatrixAntipodal|=" << ms(time3, time4) << "\n";
#endif

  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(WMat);
#ifdef TIMINGS
  SingletonTime time5;
  std::cerr << "|GetStabilizerWeightMatrix_Kernel|=" << ms(time4, time5)
            << "\n";
#endif
  return ListGen;
}

//
// Several asymmetric matrices
//

// The matrices in ListMat do not have to be symmetric.
template <typename T, typename Tint, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value>
T_TranslateToMatrix_ListMat_SHV(std::vector<MyMatrix<T>> const &ListMat,
                                MyMatrix<Tint> const &SHV) {
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
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}

template <bool is_symmetric, typename T, typename Tidx_value>
WeightMatrix<is_symmetric, std::vector<T>, Tidx_value>
GetWeightMatrix_ListComm(MyMatrix<T> const &TheEXT, MyMatrix<T> const &GramMat,
                         std::vector<MyMatrix<T>> const &ListComm) {
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
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}

template <typename T, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value>
GetWeightMatrix_ListMatrix(std::vector<MyMatrix<T>> const &ListMatrix,
                           MyMatrix<T> const &TheEXT) {
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
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}

//
// Various construction of weighted matrix
//

template <typename T, typename Tint, typename Tidx_value>
WeightMatrix<true, T, Tidx_value>
T_TranslateToMatrix_QM_SHV(MyMatrix<T> const &qMat, MyMatrix<Tint> const &SHV) {
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
      if (value1 == 0) { // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eScal);
      }
      Tidx_value &value2 = ValueMap[-eScal];
      if (value2 == 0) { // This is a missing value
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
                                           INP_ListWeight, weight_ordered);
}

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::vector<MyMatrix<T>> LinPolytopeIntegralWMat_Automorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep) {
  using Tgr = GraphBitset;
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "|GRP1|=" << GRP1.size()
            << " RankMat(ep.first)=" << RankMat(ep.first)
            << " |ep.first|=" << ep.first.rows() << " / " << ep.first.cols()
            << "\n";
  bool test = CheckStabilizerWeightMatrix(ep.second, GRP1);
  std::cerr << "test=" << test << "\n";
#endif
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer_Method8(ep.first, GRP1);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "We have GRPfull\n";
#endif
  std::vector<MyMatrix<T>> ListGenMat;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    MyMatrix<T> eMat_T = FindTransformation(ep.first, ep.first, eGen);
    ListGenMat.push_back(eMat_T);
  }
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "We have ListGenMat\n";
#endif
  return ListGenMat;
}

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::optional<MyMatrix<T>> LinPolytopeIntegralWMat_Isomorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep,
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &fp) {
  //#define DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphBitset;
  if (ep.first.rows() != fp.first.rows() || ep.first.cols() != fp.first.cols())
    return {};
  if (ep.second.GetWeight() != fp.second.GetWeight())
    return {};
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "|ep.first|=" << ep.first.rows() << " / " << ep.first.cols()
            << " rnk=" << RankMat(ep.first) << "\n";
  std::cerr << "|fp.first|=" << fp.first.rows() << " / " << fp.first.cols()
            << " rnk=" << RankMat(fp.first) << "\n";
  std::cerr << "ep.first=\n";
  WriteMatrix(std::cerr, ep.first);
  std::cerr << "fp.first=\n";
  WriteMatrix(std::cerr, fp.first);
  std::cerr << "ep.second=\n";
  PrintWeightedMatrix(std::cerr, ep.second);
  std::cerr << "fp.second=\n";
  PrintWeightedMatrix(std::cerr, fp.second);
#endif
#ifdef TIMINGS
  SingletonTime time1;
#endif

  //  std::cerr << "Before eCanonicReord\n";
  std::vector<Tidx> eCanonicReord =
      GetGroupCanonicalizationVector_Kernel<Tval, Tgr, Tidx, Tidx_value>(
          ep.second)
          .first;
  //  std::cerr << "Before fCanonicReord\n";
  std::vector<Tidx> fCanonicReord =
      GetGroupCanonicalizationVector_Kernel<Tval, Tgr, Tidx, Tidx_value>(
          fp.second)
          .first;
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|GetGroupCanonicalizationVector_Kernel|=" << ms(time1, time2)
            << "\n";
#endif
  using Tfield = typename overlying_field<T>::field_type;
  //  std::cerr << "Before IsomorphismFromCanonicReord\n";
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(
          ep.first, fp.first, eCanonicReord, fCanonicReord);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|IsomorphismFromCanonicReord|=" << ms(time2, time3) << "\n";
#endif
  if (!IsoInfo) {
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
    std::cerr << "We failed to find IsoInfo\n";
#endif
    return {};
  }
  Telt ePerm(IsoInfo->first);
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "ePerm=" << ePerm << "\n";
  std::cerr << "det(eMat)=" << DeterminantMat(IsoInfo->second)
            << "  eMat=" << StringMatrixGAP(IsoInfo->second) << "\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second);
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|GetStabilizerWeightMatrix|=" << ms(time3, time4) << "\n";
#endif
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "|GRP1|=" << GRP1.size() << "\n";
  for (auto &eGen : GRP1.GeneratorsOfGroup()) {
    MyMatrix<T> eGen_T = FindTransformation(ep.first, ep.first, eGen);
    std::cerr << "det(eGen_T)=" << DeterminantMat(eGen_T)
              << " eGen_T=" << StringMatrixGAP(eGen_T) << "\n";
  }
#endif
  std::optional<MyMatrix<T>> eRes =
      LinPolytopeIntegral_Isomorphism_Method8(ep.first, fp.first, GRP1, ePerm);
#ifdef TIMINGS
  SingletonTime time5;
  std::cerr << "|LinPolytopeIntegral_Isomorphism_Method8|=" << ms(time4, time5)
            << "\n";
#endif
  if (eRes) {
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
    std::cerr << "Found one isomorphism\n";
#endif
    return *eRes;
  }
#ifdef DEBUG_LIN_POLYTOPE_INTEGRAL_WMAT
  std::cerr << "eRes is unassigned\n";
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
#endif  // SRC_POLY_TEMP_POLYTOPEEQUISTAB_H_
// clang-format on
