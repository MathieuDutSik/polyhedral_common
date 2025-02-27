// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SIGN_SIGNATURE_H_
#define SRC_LATT_SIGN_SIGNATURE_H_

// clang-format off
#include "MAT_Matrix.h"
#include "COMB_Vectors.h"
// clang-format on

template <typename T>
MyVector<T> FindNonIsotropeVector(MyMatrix<T> const &SymMat) {
  int n = SymMat.rows();
  MyVector<T> eVect = ZeroVector<T>(n);
  for (int i = 0; i < n; i++)
    if (SymMat(i, i) != 0) {
      eVect(i) = 1;
      return eVect;
    }
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++)
      if (SymMat(i, j) != 0) {
        eVect(i) = 1;
        eVect(j) = 1;
        return eVect;
      }
  std::cerr << "SymMat=\n";
  WriteMatrix(std::cerr, SymMat);
  std::cerr << "Clear error in FindNonIsotropeVector\n";
  throw TerminalException{1};
}

template <typename T> struct DiagSymMat {
  MyMatrix<T> Transform;
  MyMatrix<T> RedMat;
  int nbZero;
  int nbPlus;
  int nbMinus;
};

template <typename T>
DiagSymMat<T>
DiagonalizeNonDegenerateSymmetricMatrix(MyMatrix<T> const &SymMat) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n = SymMat.rows();
#ifdef DEBUG_POSITIVITY
  std::cerr << "Beginning of DiagonalizeNonDegenerateSymmetricMatrix\n";
  WriteMatrix(std::cerr, SymMat);
  int rnk = RankMat(SymMat);
  if (rnk != n) {
    std::cerr << "POS: Error in DiagonalizeNonDegenerateSymmetricMatrix\n";
    std::cerr << "POS: rnk=" << rnk << " n=" << n << "\n";
    throw TerminalException{1};
  }
#endif
  if (n == 0) {
    return {{}, {}, 0, 0, 0};
  }
  std::vector<MyVector<T>> ListVect;
  for (int i = 0; i < n; i++) {
    MyMatrix<T> BasisOrthogonal;
    if (i == 0) {
      BasisOrthogonal = IdentityMat<T>(n);
    } else {
#ifdef DEBUG_POSITIVITY
      std::cerr << "|ListVect|=" << ListVect.size() << "\n";
#endif
      MyMatrix<T> TheBasis = MatrixFromVectorFamily(ListVect);
      MyMatrix<T> eProd = SymMat * TheBasis.transpose();
      BasisOrthogonal = NullspaceMat(eProd);
    }
    MyMatrix<T> RedInBasis =
        BasisOrthogonal * SymMat * BasisOrthogonal.transpose();
    MyVector<T> V = FindNonIsotropeVector(RedInBasis);
    MyVector<T> Vins = (BasisOrthogonal.transpose()) * V;
    ListVect.push_back(Vins);
  }
  MyMatrix<T> TheBasis = MatrixFromVectorFamily(ListVect);
  MyMatrix<T> RedMat = TheBasis * SymMat * TheBasis.transpose();
#ifdef DEBUG_POSITIVITY
  std::cerr << "RedMat=\n";
  WriteMatrix(std::cerr, RedMat);
#endif
  int nbPlus = 0;
  int nbMinus = 0;
  int nbZero = 0;
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0)
      nbPlus++;
    if (RedMat(i, i) < 0)
      nbMinus++;
  }
#ifdef DEBUG_POSITIVITY
  std::cerr << "nbPlus=" << nbPlus << " nbMinus=" << nbMinus << "\n";
#endif
  return {TheBasis, RedMat, nbZero, nbPlus, nbMinus};
}

template <typename T> struct NSPreduction {
  MyMatrix<T> RedMat;
  MyMatrix<T> Transform;
  MyMatrix<T> NonDegenerate;
};

template <typename T>
MyMatrix<T> SymmetricExtractSubMatrix(MyMatrix<T> const &M,
                                      std::vector<int> const &S) {
  int siz = S.size();
  MyMatrix<T> Mret(siz, siz);
  for (int i = 0; i < siz; i++)
    for (int j = 0; j < siz; j++) {
      int i2 = S[i];
      int j2 = S[j];
      Mret(i, j) = M(i2, j2);
    }
  return Mret;
}

template <typename T>
NSPreduction<T> NullspaceReduction(MyMatrix<T> const &SymMat) {
  int n = SymMat.rows();
  MyMatrix<T> PreTransfMat = NullspaceMat(SymMat);
  int DimKern = PreTransfMat.rows();
  MyMatrix<T> TransfMat = ZeroMatrix<T>(n, n);
  for (int iRow = 0; iRow < DimKern; iRow++)
    for (int i = 0; i < n; i++)
      TransfMat(iRow, i) = PreTransfMat(iRow, i);
  std::vector<int> eSet = ColumnReductionSet(PreTransfMat);
  std::vector<int> totSet(n);
  for (int i = 0; i < n; i++)
    totSet[i] = i;
  std::vector<int> diffSet = DifferenceVect(totSet, eSet);
  int idx = DimKern;
  for (auto &eCol : diffSet) {
    TransfMat(idx, eCol) = 1;
    idx++;
  }
  MyMatrix<T> SymMat2 = TransfMat * SymMat * TransfMat.transpose();
  std::vector<int> gSet;
  for (int i = DimKern; i < n; i++)
    gSet.push_back(i);
  MyMatrix<T> MatExt = SymmetricExtractSubMatrix(SymMat2, gSet);
  return {SymMat2, TransfMat, MatExt};
}

template <typename T>
DiagSymMat<T> DiagonalizeSymmetricMatrix(MyMatrix<T> const &SymMat) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n1 = SymMat.rows();
  NSPreduction<T> NSP1 = NullspaceReduction(SymMat);
  MyMatrix<T> RMat1 = NSP1.Transform;
  MyMatrix<T> SymMat2 = NSP1.NonDegenerate;
  DiagSymMat<T> NSP2 = DiagonalizeNonDegenerateSymmetricMatrix(SymMat2);
  int n2 = SymMat2.rows();
  MyMatrix<T> RMat2 = ZeroMatrix<T>(n1, n1);
  for (int i = 0; i < n2; i++)
    for (int j = 0; j < n2; j++)
      RMat2(i + n1 - n2, j + n1 - n2) = NSP2.Transform(i, j);
  for (int i = 0; i < n1 - n2; i++)
    RMat2(i, i) = 1;
  MyMatrix<T> RMat = RMat2 * RMat1;
  MyMatrix<T> RedMat = RMat * SymMat * RMat.transpose();
  int nbPlus = 0;
  int nbMinus = 0;
  int nbZero = 0;
  for (int i = 0; i < n1; i++) {
    if (RedMat(i, i) > 0)
      nbPlus++;
    if (RedMat(i, i) < 0)
      nbMinus++;
    if (RedMat(i, i) == 0)
      nbZero++;
  }
  return {RMat, RedMat, nbZero, nbPlus, nbMinus};
}

// clang-format off
#endif  //  SRC_LATT_SIGN_SIGNATURE_H_
// clang-format on
