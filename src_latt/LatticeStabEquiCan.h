// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_MATRIXCANONICALFORM_H_
#define SRC_LATT_MATRIXCANONICALFORM_H_

// clang-format off
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include <utility>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef DISABLE_DEBUG_LATTICE_STAB_EQUI_CAN
#undef DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_STAB_EQUI_CAN
#endif

template <typename T, typename Tint> struct Canonic_PosDef {
  MyMatrix<Tint> Basis;
  MyMatrix<Tint> SHV;
  MyMatrix<T> Mat;
};

//
// The canonic form
//

template<typename T>
void check_determinant_one([[maybe_unused]] MyMatrix<T> const& M) {
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  T eDet = DeterminantMat(M);
  T eDet_abs = T_abs(eDet);
  if (eDet_abs != 1) {
    std::cerr << "LSEC: The matrix should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint> ComputeCanonicalForm_inner(std::vector<MyMatrix<T>> const &ListMat,
                                                   MyMatrix<Tint> const &SHV,
                                                   std::ostream &os) {
  using Tgr = GraphListAdj;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<T> const& inpMat = ListMat[0];
  int nbRow = SHV.rows();
  int n = SHV.cols();
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: nbRow=" << nbRow << " n=" << n << "\n";
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr
        << "LSEC: The set of vector does not respect the central symmetry "
           "condition\n";
    throw TerminalException{1};
  }
#endif
  //
  // Computing the scalar product matrix
  //
  using Tidx_value = int16_t;
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      T_TranslateToMatrix_ListMat_SHV<T, Tint, Tidx_value>(ListMat, SHV, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: WMat|=" << time << "\n";
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Original WMat=\n";
  PrintWeightedMatrix(os, WMat);
#endif
  WMat.ReorderingSetWeight();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ReorderingSetWeight|=" << time << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::vector<int> CanonicOrd =
    GetCanonicalizationVector_Kernel<std::vector<T>, Tgr, int>(WMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CNA
  os << "|LSEC: GetCanonicalizationVector_Kernel|=" << time << "\n";
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
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: SHVcan|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan_Tint_pre =
      ComputeRowHermiteNormalForm(SHVcan_Tint).first;
  MyMatrix<Tint> BasisCan_Tint = TransposedMat(Inverse(BasisCan_Tint_pre));
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ReductionMatrix|=" << time << "\n";
#endif
  MyMatrix<T> BasisCan_T = UniversalMatrixConversion<T, Tint>(BasisCan_Tint);
  check_determinant_one(BasisCan_T);
  MyMatrix<T> RetMat = BasisCan_T * inpMat * BasisCan_T.transpose();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: Matrix products|=" << time << "\n";
#endif
  return {BasisCan_Tint, SHVcan_Tint, RetMat};
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint> ComputeCanonicalForm(MyMatrix<T> const &inpMat,
                                             std::ostream &os) {
  //
  // Computing the Z-basis on which the computation relies.
  //
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ComputeCanonicalForm_inner(ListMat, SHV, os);
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint>
ComputeCanonicalFormMultiple(std::vector<MyMatrix<T>> const &ListMat,
                             std::ostream &os) {
  //
  // Computing the Z-basis on which the computation relies.
  //
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<T> inpMat = ListMat[0];
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  return ComputeCanonicalForm_inner<T, Tint>(ListMat, SHV, os);
}

//
// Automorphism code
//

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> ArithmeticAutomorphismGroupMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat, MyMatrix<Tint> const &SHV,
    std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, 0);

  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(SHV_T, ListMat, Vdiag,
                                                         os);

  std::vector<MyMatrix<Tint>> ListGenRet;
  for (auto &eList : ListGen) {
    std::optional<MyMatrix<T>> opt =
        FindMatrixTransformationTest(SHV_T, SHV_T, eList);
    if (!opt) {
      std::cerr << "LSEC: Failed to find the matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &M_T = *opt;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
    if (!IsIntegralMatrix(M_T)) {
      std::cerr << "LSEC: Bug: The matrix should be integral\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
    ListGenRet.push_back(M);
  }
  return ListGenRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup_inner(MyMatrix<T> const &inpMat,
                                  MyMatrix<Tint> const &SHV, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(ListMat,
                                                                    SHV, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroupMultiple(std::vector<MyMatrix<T>> const &ListMat,
                                    std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis|=" << time << "\n";
#endif
  return ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(ListMat,
                                                                    SHV, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup(MyMatrix<T> const &inpMat, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
}

//
// Equivalence code
//

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> ArithmeticEquivalenceMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat1, MyMatrix<Tint> const &SHV1,
    std::vector<MyMatrix<T>> const &ListMat2, MyMatrix<Tint> const &SHV2,
    std::ostream &os) {
  using Tidx = uint32_t;
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: |SHV1|=" << SHV1.rows() << " |SHV2|=" << SHV2.rows() << "\n";
#endif

  if (SHV1.rows() != SHV2.rows())
    return {};
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: After the test\n";
#endif

  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(SHV2);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: SHV1_T / SHV2_T|=" << time << "\n";
#endif

  int n_rows = SHV1_T.rows();
  std::vector<T> Vdiag1(n_rows, 0);
  std::vector<T> Vdiag2(n_rows, 0);
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Before the TestEquivalence_ListMat_Vdiag\n";
#endif
  std::optional<std::vector<Tidx>> opt =
      TestEquivalence_ListMat_Vdiag<T, T, Tidx>(SHV2_T, ListMat2, Vdiag2,
                                                SHV1_T, ListMat1, Vdiag1, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: TestEquivalence_ListMat_Vdiag|=" << time << "\n";
#endif

  if (!opt)
    return {};
  std::optional<MyMatrix<T>> optB =
      FindMatrixTransformationTest(SHV2_T, SHV1_T, *opt);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: FindMatrixTransformationTest|=" << time << "\n";
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  if (!optB) {
    std::cerr << "LSEC: We have a matrix bug\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> const &M_T = *optB;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  if (!IsIntegralMatrix(M_T)) {
    std::cerr << "LSEC: The matrix should be integral\n";
    throw TerminalException{1};
  }
  for (size_t i = 0; i < ListMat1.size(); i++) {
    MyMatrix<T> eMat1 = ListMat1[i];
    MyMatrix<T> eMat2 = ListMat2[i];
    MyMatrix<T> eProd = M_T * eMat1 * M_T.transpose();
    if (eProd != eMat2) {
      std::cerr << "LSEC: Inconsistency error in the code\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ArithmeticEquivalenceMultiple_inner, M|=" << time << "\n";
#endif
  return M;
}

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence_inner(
    MyMatrix<T> const &inpMat1, MyMatrix<Tint> const &SHV1,
    MyMatrix<T> const &inpMat2, MyMatrix<Tint> const &SHV2, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple_inner(ListMat1, SHV1, ListMat2, SHV2,
                                             os);
}

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>>
ArithmeticEquivalenceMultiple(std::vector<MyMatrix<T>> const &ListMat1,
                              std::vector<MyMatrix<T>> const &ListMat2,
                              std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHV1 =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat1[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis1|=" << time << "\n";
#endif
  MyMatrix<Tint> SHV2 =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat2[0], os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ExtractInvariantVectorFamilyZbasis2|=" << time << "\n";
#endif
  return ArithmeticEquivalenceMultiple_inner(ListMat1, SHV1, ListMat2, SHV2,
                                             os);
}

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence(MyMatrix<T> const &inpMat1,
                                                    MyMatrix<T> const &inpMat2,
                                                    std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple<T, Tint>(ListMat1, ListMat2, os);
}

template <typename T, typename Tint>
Canonic_PosDef<T, Tint>
ComputeCanonicalFormSymplectic(MyMatrix<T> const &inpMat, std::ostream &os) {
  int n_tot = inpMat.rows();
  if (n_tot % 2 == 1) {
    std::cerr << "LSEC: The dimension is odd\n";
    throw TerminalException{1};
  }
  int n = n_tot / 2;
  MyMatrix<Tint> SympFormMat = ZeroMatrix<Tint>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SympFormMat(i, n + i) = 1;
    SympFormMat(n + i, i) = -1;
  }
  Canonic_PosDef<T, Tint> CanPosDef =
      ComputeCanonicalFormMultiple<T, Tint>({inpMat, SympFormMat}, os);
  MyMatrix<Tint> BasisSymp_Tint = SYMPL_ComputeSymplecticBasis(CanPosDef.SHV);
  MyMatrix<T> BasisSymp_T = UniversalMatrixConversion<T, Tint>(BasisSymp_Tint);
  MyMatrix<T> RetMat = BasisSymp_T * inpMat * TransposedMat(BasisSymp_T);
  return {BasisSymp_Tint, CanPosDef.SHV, RetMat};
}

// clang-format off
#endif  // SRC_LATT_MATRIXCANONICALFORM_H_
// clang-format on
