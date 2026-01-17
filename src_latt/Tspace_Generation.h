// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_GENERATION_H_
#define SRC_LATT_TSPACE_GENERATION_H_

// clang-format off
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_TSPACE_GENERATION
#endif

#ifdef DISABLE_DEBUG_TSPACE_GENERATION
#undef DEBUG_TSPACE_GENERATION
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
  std::vector<MyMatrix<T>> ListComm;
  return BuildLinSpace(SuperMat, ListMat, ListComm);
}

template <typename T>
MyMatrix<T> __RealQuadMatSpace(MyMatrix<T> const &eMatB,
                               MyMatrix<T> const &eMatC, int n, T const &eSum,
                               T const &eProd) {
  MyMatrix<T> eMatN = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i2 = 0; i2 < n; i2++) {
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
  }
  return eMatN;
}

template <typename T>
MyMatrix<T> GetCommRealQuadratic(int n, T const &eSum, T const &eProd) {
  MyMatrix<T> Imultiplication = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    Imultiplication(i, n + i) = 1;
    Imultiplication(n + i, i) = -eProd;
    Imultiplication(n + i, n + i) = eSum;
  }
  // Need to write the code for the real quadratic.
  // Maybe it is the same matrix as the one for the imaginary.
  // We need to investigate.
  std::cerr << "We need to code the GetCommRealQuadratic function\n";
  throw TerminalException{1};

  return Imultiplication;
}

template <typename T>
MyMatrix<T> GetCommImagQuadratic(int n, T const &eSum, T const &eProd) {
  MyMatrix<T> Imultiplication = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    Imultiplication(i, n + i) = 1;
    Imultiplication(n + i, i) = -eProd;
    Imultiplication(n + i, n + i) = eSum;
  }
  return Imultiplication;
}

template <typename T>
LinSpaceMatrix<T> ComputeRealQuadraticSpace(int n, T const &eSum,
                                            T const &eProd, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<T> eMatB(n, n);
  MyMatrix<T> eMatC(n, n);
  for (int i = 0; i < n; i++) {
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
  }
  // That choice of SuperMat matters in the self-duality.
  // It represents the identity x identity in the space S^n x S^n
  // where it lives.
  eMatB = IdentityMat<T>(n);
  ZeroAssignation(eMatC);
  MyMatrix<T> SuperMat = __RealQuadMatSpace(eMatB, eMatC, n, eSum, eProd);
  // The matrix eComm represent multiplication by the primitive element
  MyMatrix<T> eComm = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    eComm(i, n + i) = 1;
    eComm(n + i, i) = -eProd;
    eComm(n + i, n + i) = eSum;
  }
  std::vector<MyMatrix<T>> ListMatRet = IntegralSaturationSpace(ListMat, os);
  LinSpaceMatrix<T> LinSpa = BuildLinSpace(SuperMat, ListMatRet, {eComm});
  LinSpa.l_spanning_elements.push_back(eComm);
  return LinSpa;
}

template <typename T>
LinSpaceMatrix<T> ComputeImagQuadraticSpace(int n, T const &eSum,
                                            T const &eProd, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat;
  T Discriminant = eSum * eSum - 4 * eProd;
  if (Discriminant >= 0) {
    std::cerr << "Discriminant is positive\n";
    std::cerr << "Impossible for an imaginary space\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++) {
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
      MyMatrix<T> fMat = RemoveFractionMatrix(eMat);
      ListMat.push_back(fMat);
    }
  }
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      MyMatrix<T> eMat(2 * n, 2 * n);
      ZeroAssignation(eMat);
      eMat(n + i, j) = 1;
      eMat(j, n + i) = 1;
      eMat(n + j, i) = -1;
      eMat(i, n + j) = -1;
      ListMat.push_back(eMat);
    }
  }
  // That choice of SuperMat matters in the self-duality.
  // It corresponds to the identity matrix in the H^n space
  // in which the SuperMat truly lives.
  MyMatrix<T> SuperMat_pre = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SuperMat_pre(i, i) = 1;
    T eVal = eSum / 2;
    SuperMat_pre(n + i, i) = eVal;
    SuperMat_pre(i, n + i) = eVal;
    //
    SuperMat_pre(n + i, n + i) = eProd;
  }
  MyMatrix<T> SuperMat = RemoveFractionMatrix(SuperMat_pre);
  // The matrix eComm represent multiplication by the primitive element
  MyMatrix<T> eComm = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    eComm(i, n + i) = 1;
    eComm(n + i, i) = -eProd;
    eComm(n + i, n + i) = eSum;
  }
  std::vector<MyMatrix<T>> ListMatRet = IntegralSaturationSpace(ListMat, os);
  LinSpaceMatrix<T> LinSpa = BuildLinSpace(SuperMat, ListMatRet, {eComm});
  LinSpa.l_spanning_elements.push_back(eComm);
  return LinSpa;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_GENERATION_H_
// clang-format on
