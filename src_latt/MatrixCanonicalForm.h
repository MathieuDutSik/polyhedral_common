// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_MATRIXCANONICALFORM_H_
#define SRC_LATT_MATRIXCANONICALFORM_H_

// clang-format off
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"
#include "ShortestUniversal.h"
#include "Temp_PolytopeEquiStab.h"
#include <utility>
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint> struct Canonic_PosDef {
  MyMatrix<Tint> Basis;
  MyMatrix<Tint> SHV;
  MyMatrix<T> Mat;
};

template <typename T, typename Tint>
Canonic_PosDef<T, Tint> ComputeCanonicalForm(MyMatrix<T> const &inpMat) {
  //
  // Computing the Z-basis on which the computation relies.
  //
  //  std::cerr << "inpMat=\n";
  //  WriteMatrix(std::cerr, inpMat);
#ifdef TIMINGS
  std::cerr << "Begining of ComputeCanonicalForm\n";
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat);
#ifdef TIMINGS
  std::cerr << "|ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  int nbRow = SHV.rows();
  int n = SHV.cols();
#ifdef DEBUG_CANONIC
  std::cerr << "nbRow=" << nbRow << " n=" << n << "\n";
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr << "The set of vector does not respect the central symmetry "
                 "condition\n";
    throw TerminalException{1};
  }
#endif
  //  std::cerr << "SHV=\n";
  //  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  using Tidx_value = int16_t;
  WeightMatrix<true, T, Tidx_value> WMat =
      T_TranslateToMatrix_QM_SHV<T, Tint, Tidx_value>(inpMat, SHV, std::cerr);
#ifdef TIMINGS
  std::cerr << "|WMat|=" << time << "\n";
#endif
#ifdef DEBUG_CANONIC
  std::cerr << "Original WMat=\n";
  PrintWeightedMatrix(std::cerr, WMat);
#endif
  WMat.ReorderingSetWeight();
#ifdef DEBUG_CANONIC
  std::cerr << "Weight reordering WMat=\n";
  PrintWeightedMatrix(std::cerr, WMat);
#endif
#ifdef TIMINGS
  std::cerr << "|ReorderingSetWeight|=" << time << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::vector<int> CanonicOrd =
      GetCanonicalizationVector_Kernel<T, GraphBitset, int>(WMat, std::cerr);
#ifdef DEBUG_CANONIC
  std::cerr << "CanonicOrd=";
  for (auto &eV : CanonicOrd)
    std::cerr << " " << eV;
  std::cerr << "\n";
#endif
#ifdef TIMINGS
  std::cerr << "|GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan_Tint(n, nbRow);
  for (int iRowCan = 0; iRowCan < nbRow; iRowCan++) {
    int iRowNative = CanonicOrd[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixCol(SHVcan_Tint, iRowCan, eRow_Tint);
  }
#ifdef TIMINGS
  std::cerr << "|SHVcan|=" << time << "\n";
#endif
#ifdef DEBUG_CANONIC
  std::cerr << "SHVred=\n";
  std::pair<MyMatrix<Tint>, MyMatrix<Tint>> pair_calc =
      ComputeRowHermiteNormalForm(SHVcan_Tint);
  MyMatrix<Tint> PrtMat = pair_calc.second;
  MyMatrix<Tint> eDiff =
      Inverse(pair_calc.first) * pair_calc.second - SHVcan_Tint;
  WriteMatrix(std::cerr, PrtMat);
  if (!IsZeroMatrix(eDiff)) {
    std::cerr << "The matrix eDiff should be zero\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> BasisCan_Tint_pre =
      ComputeRowHermiteNormalForm(SHVcan_Tint).first;
  MyMatrix<Tint> BasisCan_Tint = TransposedMat(Inverse(BasisCan_Tint_pre));
  //  MyMatrix<Tint> BasisCan_Tint = TransposedMat(Inverse(BasisCan_Tint_pre));
  //  std::cerr << "SHVcan_Tint=\n";
  //  WriteMatrix(std::cerr, SHVcan_Tint);
#ifdef TIMINGS
  std::cerr << "|ReductionMatrix|=" << time << "\n";
#endif
  MyMatrix<T> BasisCan_T = UniversalMatrixConversion<T, Tint>(BasisCan_Tint);
#ifdef DEBUG_CANONIC
  T eDet = DeterminantMat(BasisCan_T);
  T eDet_abs = T_abs(eDet);
  if (eDet_abs != 1) {
    std::cerr << "The matrix should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> RetMat = BasisCan_T * inpMat * BasisCan_T.transpose();
#ifdef TIMINGS
  std::cerr << "|Matrix products|=" << time << "\n";
#endif
#ifdef DEBUG_CANONIC
  WeightMatrix<true, T, Tidx_value> WMat_B =
      T_TranslateToMatrix_QM_SHV<T, Tint, Tidx_value>(
          RetMat, TransposedMat(SHVcan_Tint), std::cerr);
  WMat_B.ReorderingSetWeight();
#endif
  return {BasisCan_Tint, SHVcan_Tint, RetMat};
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint>
ComputeCanonicalFormMultiple(std::vector<MyMatrix<T>> const &ListMat) {
  //
  // Computing the Z-basis on which the computation relies.
  //
  //  std::cerr << "inpMat=\n";
  //  WriteMatrix(std::cerr, inpMat);
#ifdef TIMINGS
  std::cerr << "Begining of ComputeCanonicalForm\n";
  MicrosecondTime time;
#endif
  MyMatrix<T> inpMat = ListMat[0];
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat);
#ifdef TIMINGS
  std::cerr << "|ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  int nbRow = SHV.rows();
  int n = SHV.cols();
#ifdef DEBUG_CANONIC
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr << "The set of vector does not respect the central symmetry "
                 "condition\n";
    throw TerminalException{1};
  }
#endif
  //  std::cerr << "SHV=\n";
  //  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  using Tidx_value = int16_t;
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      T_TranslateToMatrix_ListMat_SHV<T, Tint, Tidx_value>(ListMat, SHV,
                                                           std::cerr);
#ifdef TIMINGS
  std::cerr << "|WMat|=" << time << "\n";
#endif
  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::cerr << "|ReorderingSetWeight|=" << time << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  WeightMatrix<true, std::vector<T>, Tidx_value> WMatSymm =
      WMat.GetSymmetricWeightMatrix();
  std::vector<int> CanonicOrdSymm =
      GetCanonicalizationVector_Kernel<std::vector<T>, GraphBitset, int>(
          WMatSymm, std::cerr);
  std::vector<int> CanonicOrd =
      GetCanonicalizationFromSymmetrized(CanonicOrdSymm);
#ifdef TIMINGS
  std::cerr << "|GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan_Tint(n, nbRow);
  for (int iRowCan = 0; iRowCan < nbRow; iRowCan++) {
    int iRowNative = CanonicOrd[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixCol(SHVcan_Tint, iRowCan, eRow_Tint);
  }
#ifdef TIMINGS
  std::cerr << "|SHVcan|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan_Tint_pre =
      ComputeRowHermiteNormalForm(SHVcan_Tint).first;
  MyMatrix<Tint> BasisCan_Tint = TransposedMat(Inverse(BasisCan_Tint_pre));
  //  std::cerr << "SHVcan_Tint=\n";
  //  WriteMatrix(std::cerr, SHVcan_Tint);
#ifdef TIMINGS
  std::cerr << "|ReductionMatrix|=" << time << "\n";
#endif
  MyMatrix<T> BasisCan_T = UniversalMatrixConversion<T, Tint>(BasisCan_Tint);
#ifdef DEBUG_CANONIC
  T eDet = DeterminantMat(BasisCan_T);
  T eDet_abs = T_abs(eDet);
  if (eDet_abs != 1) {
    std::cerr << "The matrix should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> RetMat = BasisCan_T * inpMat * TransposedMat(BasisCan_T);
#ifdef TIMINGS
  std::cerr << "|Matrix products|=" << time << "\n";
#endif
  return {BasisCan_Tint, SHVcan_Tint, RetMat};
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint>
ComputeCanonicalFormSymplectic(MyMatrix<T> const &inpMat) {
  int n_tot = inpMat.rows();
  if (n_tot % 2 == 1) {
    std::cerr << "The dimension is odd\n";
    throw TerminalException{1};
  }
  int n = n_tot / 2;
  MyMatrix<Tint> SympFormMat = ZeroMatrix<Tint>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SympFormMat(i, n + i) = 1;
    SympFormMat(n + i, i) = -1;
  }
  Canonic_PosDef<T, Tint> CanPosDef =
      ComputeCanonicalFormMultiple<T, Tint>({inpMat, SympFormMat});
  MyMatrix<Tint> BasisSymp_Tint = SYMPL_ComputeSymplecticBasis(CanPosDef.SHV);
  MyMatrix<T> BasisSymp_T = UniversalMatrixConversion<T, Tint>(BasisSymp_Tint);
  MyMatrix<T> RetMat = BasisSymp_T * inpMat * TransposedMat(BasisSymp_T);
  return {BasisSymp_Tint, CanPosDef.SHV, RetMat};
}

// clang-format off
#endif  // SRC_LATT_MATRIXCANONICALFORM_H_
// clang-format on
