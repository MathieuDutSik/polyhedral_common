// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SHORT_SHORT_SHORTESTCONFIGSTABEQUICANINV_H_
#define SRC_SHORT_SHORT_SHORTESTCONFIGSTABEQUICANINV_H_

// clang-format off
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
// clang-format on

template <typename T, typename Tint> struct ShortIso {
  MyMatrix<T> GramMat;
  MyMatrix<Tint> SHVdisc;
};

template <typename T, typename Tint>
ShortIso<T, Tint> SHORT_GetInformation(MyMatrix<Tint> const &M,
                                       std::ostream &os) {
  int n = M.cols();
  int nbVect = M.rows();
#ifdef DEBUG_SHORTEST_CONFIG
  std::cerr << "nbVect=" << nbVect << "\n";
#endif
  MyMatrix<T> Amat = ZeroMatrix<T>(n, n);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<Tint> eVect_Tint = GetMatrixRow(M, iVect);
    MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect_Tint);
    MyMatrix<T> pMat = RankOneMatrix(eVect_T);
    Amat += pMat;
  }
  MyMatrix<T> TheGramMat = Inverse(Amat);
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(TheGramMat, os);
  MyMatrix<Tint> Mc1 = 2 * M;
  MyMatrix<Tint> Mc2 = -2 * M;
  MyMatrix<Tint> Mcopy = Concatenate(Mc1, Mc2);
  MyMatrix<Tint> SHVret = Concatenate(SHV, Mcopy);
  return {std::move(TheGramMat), std::move(SHVret)};
}

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> SHORT_TestEquivalence(MyMatrix<Tint> const &M1,
                                                    MyMatrix<Tint> const &M2,
                                                    std::ostream &os) {
  ShortIso<T, Tint> eRec1 = SHORT_GetInformation<T, Tint>(M1, os);
  ShortIso<T, Tint> eRec2 = SHORT_GetInformation<T, Tint>(M2, os);
  return ArithmeticEquivalence_inner<T, Tint>(eRec1.GramMat, eRec1.SHVdisc,
                                              eRec2.GramMat, eRec2.SHVdisc, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> SHORT_GetStabilizer(MyMatrix<Tint> const &M,
                                                std::ostream &os) {
  ShortIso<T, Tint> eRec1 = SHORT_GetInformation<T, Tint>(M, os);
  return ArithmeticAutomorphismGroup_inner<T, Tint, Tgroup>(eRec1.GramMat,
                                                            eRec1.SHVdisc, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> SHORT_Canonicalize(MyMatrix<Tint> const &M, std::ostream &os) {
  ShortIso<T, Tint> eRec1 = SHORT_GetInformation<T, Tint>(M, os);
  Canonic_PosDef<T, Tint> eRec2 =
      ComputeCanonicalForm_inner<T, Tint>(eRec1.GramMat, eRec1.SHVdisc, os);
  MyMatrix<Tint> BasisInv = Inverse(eRec2.Basis);
  MyMatrix<Tint> M_can = M * BasisInv;
  return M_can;
}

// clang-format off
#endif  // SRC_SHORT_SHORT_SHORTESTCONFIGSTABEQUICANINV_H_
// clang-format on
