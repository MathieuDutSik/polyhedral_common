// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_POLYTOPEEQUISTABINT_H_
#define SRC_GROUP_POLYTOPEEQUISTABINT_H_

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
#include "MatrixGroup.h"
#include "PolytopeEquiStab.h"
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_POLYTOPE_EQUI_STAB_INT
#endif

#ifdef DEBUG
#define DEBUG_POLYTOPE_EQUI_STAB_INT
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_INT_CANONIC
#endif

template <typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LinPolytopeIntegral_Isomorphism(const MyMatrix<Tint> &EXT1,
                                const MyMatrix<Tint> &EXT2, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tfield = typename overlying_field<Tint>::field_type;
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  auto f_eval=[&](size_t threshold) -> std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> {
    std::vector<Tidx> CanonicReord1 = LinPolytope_CanonicOrdering<Tint, Tidx>(EXT1, threshold, os);
    std::vector<Tidx> CanonicReord2 = LinPolytope_CanonicOrdering<Tint, Tidx>(EXT2, threshold, os);
    //
    return IsomorphismFromCanonicReord<Tint, Tfield, Tidx>(EXT1, EXT2, CanonicReord1, CanonicReord2, os);
  };
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo = f_eval(THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
#ifdef SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_INT_CANONIC
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo_B = f_eval(THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
  check_iso_info_coherence(IsoInfo, IsoInfo_B, "LinPolytopeIntegral_Isomorphism");
#endif
  if (!IsoInfo)
    return {};
  Telt ePerm(IsoInfo->first);

  MyMatrix<Tfield> EXT1_T = UniversalMatrixConversion<Tfield, Tint>(EXT1);
  MyMatrix<Tfield> EXT2_T = UniversalMatrixConversion<Tfield, Tint>(EXT2);
  Tgroup GRP1 = LinPolytope_Automorphism<Tfield, Tgroup>(EXT1_T, os);
  std::optional<MyMatrix<Tfield>> eRes =
      LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, ePerm, os);
  if (eRes)
    return UniversalMatrixConversion<Tint, Tfield>(*eRes);
  return {};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> LinPolytopeIntegral_Isomorphism_GramMat(
    const MyMatrix<Tint> &EXT1, const MyMatrix<T> &GramMat1,
    const MyMatrix<Tint> &EXT2, const MyMatrix<T> &GramMat2, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  MyMatrix<T> EXT1_T = UniversalMatrixConversion<T, Tint>(EXT1);
  MyMatrix<T> EXT2_T = UniversalMatrixConversion<T, Tint>(EXT2);
  auto f_eval=[&](size_t threshold) -> std::optional<std::pair<std::vector<Tidx>, MyMatrix<T>>> {
    std::vector<Tidx> CanonicReord1 = LinPolytope_CanonicOrdering_GramMat<T, Tidx>(EXT1_T, GramMat1, threshold, os);
    //    os << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n";
    std::vector<Tidx> CanonicReord2 = LinPolytope_CanonicOrdering_GramMat<T, Tidx>(EXT2_T, GramMat2, threshold, os);
    //    os << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n";
    return IsomorphismFromCanonicReord_GramMat<T, T, Tidx>(EXT1_T, GramMat1, EXT2_T, GramMat2, CanonicReord1, CanonicReord2, os);
  };
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<T>>> IsoInfo = f_eval(THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
#ifdef SANITY_CHECK_THRESHOLD_SUBSET_SCHEME_INT_CANONIC
  //  os << "---------------------------------------------------------------------------\n";
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<T>>> IsoInfo_B = f_eval(THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
  check_iso_info_coherence(IsoInfo, IsoInfo_B, "LinPolytopeIntegral_Isomorphism_GramMat");
#endif
  if (!IsoInfo)
    return {};
  Telt ePerm(IsoInfo->first);

  Tgroup GRP1 = LinPolytope_Automorphism<T, Tgroup>(EXT1_T, os);
  std::optional<MyMatrix<T>> eRes =
      LinPolytopeIntegral_Isomorphism_Method8(EXT1_T, EXT2_T, GRP1, ePerm, os);
  if (eRes)
    return UniversalMatrixConversion<Tint, T>(*eRes);
  return {};
}

template <typename Tint, typename Tgroup>
Tgroup LinPolytopeIntegral_Automorphism(const MyMatrix<Tint> &EXT,
                                        std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom = LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, os);
  Tgroup GRP = LinPolytopeIntegral_Stabilizer_Method8(EXT_T, GRPisom, os);
  return GRP;
}

template <typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Automorphism_RightCoset(const MyMatrix<Tint> &EXT,
                                            std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom = LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, os);
  return LinPolytopeIntegral_Stabilizer_RightCoset_Method8(EXT_T, GRPisom, os);
}

template <typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<typename Tgroup::Telt>>
LinPolytopeIntegral_Automorphism_DoubleCoset(const MyMatrix<Tint> &EXT,
                                             Tgroup const& GrpV,
                                             std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom = LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, os);
  return LinPolytopeIntegral_Stabilizer_DoubleCoset_Method8(EXT_T, GRPisom, GrpV, os);
}

template <typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<PairCosetStabGens<typename Tgroup::Telt>>>
LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer(const MyMatrix<Tint> &EXT,
                                                       Tgroup const& GrpV,
                                                       std::ostream &os) {
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> EXT_T = UniversalMatrixConversion<Tfield, Tint>(EXT);
  Tgroup GRPisom = LinPolytope_Automorphism<Tfield, Tgroup>(EXT_T, os);
  return LinPolytopeIntegral_Stabilizer_DoubleCosetStabilizer_Method8(EXT_T, GRPisom, GrpV, os);
}

template <typename Tint, typename Tidx_value>
MyMatrix<Tint>
LinPolytopeAntipodalIntegral_CanonicForm_Tidx_value(MyMatrix<Tint> const &EXT,
                                                    std::ostream &os) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  SecondTime time;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetQmatrix|=" << time << "\n";
#endif

  std::optional<MyMatrix<Tint>> eEquiv =
      LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(EXT, Qmat, os);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick|=" << time
     << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetWeightMatrixAntipodal|=" << time << "\n";
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: ReorderingSetWeight|=" << time << "\n";
#endif

  std::vector<int> CanonicOrd =
      GetCanonicalizationVector_Kernel<Tint, GraphBitset, int>(WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetCanonicalizationVector_Kernel|=" << time << "\n";
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
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: EXTreord 2|=" << time << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm_second(EXTreord);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: ComputeColHermiteNormalForm 2|=" << time << "\n";
#endif

  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: SignRenormalizationMatrix|=" << time << "\n";
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
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  SecondTime time;
#endif
  WeightMatrixAbs<Tint, Tidx_value> WMatAbs =
      GetSimpleWeightMatrixAntipodal_AbsTrick<Tint, Tidx_value>(EXT, Qmat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetSimpleWeightMatrixAntipodal_AbsTrick|=" << time << "\n";
#endif

  using Tidx = uint32_t;
  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(
          WMatAbs.WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetStabilizerWeightMatrix_Kernel|=" << time << "\n";
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
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: Check Generators|=" << time << "\n";
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
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  SecondTime time;
#endif
  MyMatrix<Tint> Qmat = GetQmatrix(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetQmatrix|=" << time << "\n";
#endif

  std::optional<std::vector<std::vector<unsigned int>>> eEquiv =
      LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(EXT, Qmat, os);
  if (eEquiv) {
    return *eEquiv;
  }
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: LinPolytopeAntipodalIntegral_Automorphism|=" << time << "\n";
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat =
      GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetWeightMatrixAntipodal|=" << time << "\n";
#endif

  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<Tint, Tgr, Tidx, Tidx_value>(WMat, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetStabilizerWeightMatrix_Kernel|=" << time << "\n";
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

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::vector<MyMatrix<T>> LinPolytopeIntegralWMat_Automorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep,
    std::ostream &os) {
  using Tgr = GraphBitset;
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second, os);
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: |GRP1|=" << GRP1.size()
     << " RankMat(ep.first)=" << RankMat(ep.first)
     << " |ep.first|=" << ep.first.rows() << " / " << ep.first.cols() << "\n";
  bool test = CheckStabilizerWeightMatrix(ep.second, GRP1);
  os << "PES: test=" << test << "\n";
#endif
  Tgroup GRPfull = LinPolytopeIntegral_Stabilizer_Method8(ep.first, GRP1, os);
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: We have GRPfull\n";
#endif
  std::vector<MyMatrix<T>> ListGenMat;
  for (auto &eGen : GRPfull.SmallGeneratingSet()) {
    MyMatrix<T> eMat_T = FindTransformation(ep.first, ep.first, eGen);
    ListGenMat.push_back(eMat_T);
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: We have ListGenMat\n";
#endif
  return ListGenMat;
}

template <typename T, typename Tgroup, typename Tval, typename Tidx_value>
std::optional<MyMatrix<T>> LinPolytopeIntegralWMat_Isomorphism(
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &ep,
    std::pair<MyMatrix<T>, WeightMatrix<true, Tval, Tidx_value>> const &fp,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tgr = GraphBitset;
  if (ep.first.rows() != fp.first.rows() || ep.first.cols() != fp.first.cols())
    return {};
  if (ep.second.GetWeight() != fp.second.GetWeight())
    return {};
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: |ep.first|=" << ep.first.rows() << " / " << ep.first.cols()
     << " rnk=" << RankMat(ep.first) << "\n";
  os << "PES: |fp.first|=" << fp.first.rows() << " / " << fp.first.cols()
     << " rnk=" << RankMat(fp.first) << "\n";
  os << "PES: ep.first=\n";
  WriteMatrix(os, ep.first);
  os << "PES: fp.first=\n";
  WriteMatrix(os, fp.first);
  os << "PES: ep.second=\n";
  PrintWeightedMatrix(os, ep.second);
  os << "PES: fp.second=\n";
  PrintWeightedMatrix(os, fp.second);
#endif
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
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
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetGroupCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(ep.first, fp.first, eCanonicReord, fCanonicReord, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: IsomorphismFromCanonicReord|=" << time << "\n";
#endif
  if (!IsoInfo) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
    os << "PES: We failed to find IsoInfo\n";
#endif
    return {};
  }
  Telt ePerm(IsoInfo->first);
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: ePerm=" << ePerm << "\n";
  os << "PES: det(eMat)=" << DeterminantMat(IsoInfo->second)
     << "  eMat=" << StringMatrixGAP(IsoInfo->second) << "\n";
#endif
  Tgroup GRP1 =
      GetStabilizerWeightMatrix<Tval, Tgr, Tgroup, Tidx_value>(ep.second, os);
#ifdef TIMINGS_POLYTOPE_EQUI_STAB_INT
  os << "|PES: GetStabilizerWeightMatrix|=" << time << "\n";
#endif
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: |GRP1|=" << GRP1.size() << "\n";
  for (auto &eGen : GRP1.GeneratorsOfGroup()) {
    MyMatrix<T> eGen_T = FindTransformation(ep.first, ep.first, eGen);
    os << "PES: det(eGen_T)=" << DeterminantMat(eGen_T)
       << " eGen_T=" << StringMatrixGAP(eGen_T) << "\n";
  }
#endif
  std::optional<MyMatrix<T>> eRes = LinPolytopeIntegral_Isomorphism_Method8(
      ep.first, fp.first, GRP1, ePerm, os);
#ifdef TIMINGS
  os << "|PES: LinPolytopeIntegral_Isomorphism_Method8|=" << time << "\n";
#endif
  if (eRes) {
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
    os << "PES: Found one isomorphism\n";
#endif
    return *eRes;
  }
#ifdef DEBUG_POLYTOPE_EQUI_STAB_INT
  os << "PES: eRes is unassigned\n";
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
#endif  // SRC_GROUP_POLYTOPEEQUISTABINT_H_
// clang-format on
