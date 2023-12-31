// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_GENERATION_H_
#define SRC_LATT_TSPACE_GENERATION_H_

#ifdef DEBUG
#define DEBUG_TSPACE_GENERATION
#endif

template <typename T> LinSpaceMatrix<T> ComputeCanonicalSpace(int const &n) {
  std::vector<MyMatrix<T>> ListMat;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      MyMatrix<T> eMat(n, n);
      ZeroAssignation(eMat);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListMat.push_back(eMat);
    }
  MyMatrix<T> SuperMat = IdentityMat<T>(n);
  return BuildLinSpace(SuperMat, ListMat, {});
}

template <typename T>
MyMatrix<T> __RealQuadMatSpace(MyMatrix<T> const &eMatB,
                               MyMatrix<T> const &eMatC, int n, T const &eSum,
                               T const &eProd) {
  MyMatrix<T> eMatN = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i2 = 0; i2 < n; i2++)
    for (int j2 = 0; j2 < n; j2++) {
      T bVal = eMatB(i2, j2);
      T cVal = eMatC(i2, j2);
      T eVal1 = 2 * bVal + eSum * cVal;
      T eVal1b = (eSum * eSum - 2 * eProd) * bVal + eSum * eProd * cVal;
      eMatN(i2, j2) = eVal1;
      eMatN(i2 + n, j2 + n) = eVal1b;
      T eVal2 = eSum * bVal + 2 * eProd * cVal;
      eMatN(i2, j2 + n) = eVal2;
      eMatN(i2 + n, j2) = eVal2;
    }
  return eMatN;
}

template <typename T>
LinSpaceMatrix<T> ComputeRealQuadraticSpace(int n, T const &eSum,
                                            T const &eProd) {
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<T> eMatB(n, n);
  MyMatrix<T> eMatC(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      ZeroAssignation(eMatB);
      ZeroAssignation(eMatC);
      eMatB(i, j) = 1;
      eMatB(j, i) = 1;
      //
      MyMatrix<T> eMat = __RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
      ListMat.push_back(eMat);
      //
      ZeroAssignation(eMatB);
      ZeroAssignation(eMatC);
      eMatC(i, j) = 1;
      eMatC(j, i) = 1;
      //
      MyMatrix<T> fMat = __RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
      ListMat.push_back(fMat);
    }
  eMatB = IdentityMat<int>(n);
  ZeroAssignation(eMatC);
  MyMatrix<T> SuperMat = __RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
  // The matrix eComm represent multiplication by the primitive element
  MyMatrix<T> eComm = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    eComm(i, n + i) = 1;
    eComm(n + i, i) = -eProd;
    eComm(n + i, n + i) = eSum;
  }
  return BuildLinSpace(SuperMat, ListMat, {eComm});
}

template <typename T>
LinSpaceMatrix<T> ComputeImagQuadraticSpace(int n, T const &eSum,
                                            T const &eProd) {
  std::vector<MyMatrix<T>> ListMat;
  T Discriminant = eSum * eSum - 4 * eProd;
  if (Discriminant >= 0) {
    std::cerr << "Discriminant is positive\n";
    std::cerr << "Impossible for an imaginary space\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(2 * n, 2 * n);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      T eVal = eSum / 2;
      eMat(n + i, j) = eVal;
      eMat(n + j, i) = eVal;
      eMat(i, n + j) = eVal;
      eMat(j, n + i) = eVal;
      //
      eMat(n + i, n + j) = eProd;
      eMat(n + j, n + i) = eProd;
      ListMat.push_back(eMat);
    }
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++) {
      MyMatrix<T> eMat(2 * n, 2 * n);
      ZeroAssignation(eMat);
      eMat(n + i, j) = 1;
      eMat(j, n + i) = 1;
      eMat(n + j, i) = -1;
      eMat(i, n + j) = -1;
      ListMat.push_back(eMat);
    }
  MyMatrix<T> SuperMat(2 * n, 2 * n);
  ZeroAssignation(SuperMat);
  for (int i = 0; i < n; i++) {
    SuperMat(i, i) = 1;
    T eVal = eSum / 2;
    SuperMat(n + i, i) = eVal;
    SuperMat(i, n + i) = eVal;
    //
    SuperMat(n + i, n + i) = eProd;
  }
  // The matrix eComm represent multiplication by the primitive element
  MyMatrix<T> eComm = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    eComm(i, n + i) = 1;
    eComm(n + i, i) = -eProd;
    eComm(n + i, n + i) = eSum;
  }
  return BuildLinSpace(SuperMat, ListMat, {eComm});
}

//
// We search for the set of matrices satisfying g M g^T = M for all g in ListGen
//
template <typename T>
std::vector<MyMatrix<T>>
BasisInvariantForm(int const &n, std::vector<MyMatrix<T>> const &ListGen) {
  std::vector<std::vector<int>> ListCoeff;
  MyMatrix<int> ListCoeffRev(n,n);
  int pos = 0;
  for (int iLin = 0; iLin < n; iLin++) {
    for (int iCol = 0; iCol <= iLin; iCol++) {
      ListCoeff.push_back({iLin, iCol});
      ListCoeffRev(iLin, iCol) = pos;
      ListCoeffRev(iCol, iLin) = pos;
      pos++;
    }
  }
  int nbCoeff = pos;
  auto FuncPos = [&](int const &i, int const &j) -> int {
    return ListCoeffRev(i, j);
  };
  std::vector<MyVector<T>> ListEquations;
  for (auto &eGen : ListGen)
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++) {
        MyVector<T> TheEquation = ZeroVector<T>(nbCoeff);
        TheEquation(FuncPos(i, j)) += 1;
        for (int k = 0; k < n; k++)
          for (int l = 0; l < n; l++) {
            int pos = FuncPos(k, l);
            TheEquation(pos) += -eGen(i, k) * eGen(j, l);
          }
        ListEquations.push_back(TheEquation);
      }
  MyMatrix<T> MatEquations = MatrixFromVectorFamily(ListEquations);
  std::cerr << "Before call to NullspaceTrMat nbGen=" << ListGen.size()
            << " MatEquations(rows/cols)=" << MatEquations.rows() << " / "
            << MatEquations.cols() << "\n";
  MyMatrix<T> NSP = NullspaceTrMat(MatEquations);
  std::cerr << "After call to NullspaceTrMat NSP.rows=" << NSP.rows()
            << " cols=" << NSP.cols() << "\n";
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    MyMatrix<T> eMat(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        eMat(i, j) = eRow(FuncPos(i, j));
    TheBasis[iDim] = eMat;
  }
#ifdef DEBUG_TSPACE_GENERATION
  for (auto & eBasis : TheBasis) {
    for (auto & eGen : ListGen) {
      MyMatrix<T> eProd = eGen * eBasis * eGen.transpose();
      if (eProd != eBasis) {
        std::cerr << "We have eProd <> eBasis, so a linear algebra bug\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return TheBasis;
}

template <typename T> std::vector<MyMatrix<T>> StandardSymmetricBasis(int n) {
  std::vector<MyMatrix<T>> ListMat;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListMat.push_back(eMat);
    }
  return ListMat;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_GENERATION_H_
// clang-format on
