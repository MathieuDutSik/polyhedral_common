// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_ENUMERATION_K_SPACE_H_
#define SRC_RANKIN_ENUMERATION_K_SPACE_H_

#include "Shvec_exact.h"
#include "MAT_MatrixInt.h"

// Returns the Hermite constant at the n-th power.
// That is the mmaximum of min(A)^n / det(A)
// --
// So we would have min(A)^n / det(A) <= H(n)
// which gets us min(A)^n <= H(n) det(A)
template<typename T>
T GetUpperBoundHermitePower(int n) {
  if (n == 1) {
    return 1;
  }
  if (n == 2) {
    return 4/3;
  }
  if (n == 3) {
    return 2;
  }
  if (n == 4) {
    return 4;
  }
  if (n == 5) {
    return 8;
  }
  if (n == 6) {
    return 64/3;
  }
  if (n == 7) {
    return 64;
  }
  if (n == 8) {
    return 256;
  }
  int h = n*(n-1)/2;
  return (4/3)^h;
}

/*
  orthogonal projector means self adjoint.
  that is that we have
  <x, p(y)> = < p(x), y>
  We have the scalar product <x,y> = X G Y^t
  and so we get
  <p(x), y> = XP G Y^T
  <x, p(y)> = XG P^T Y
  so we have PG = G P^T
  ---
  See the code of __GetOrthogonalProjector in
  SublatticeEnumeration.g
  Though the single formula needs to be checked.
*/
template<typename T, typename Tint>
MyMatrix<T> GetOrthogonalProjector(MyMatrix<T> const& TheGramMat, MyMatrix<Tint> const& TheSubBasis) {
  int n = TheGramMat.rows();
  int hDim = TheSubBasis.rows();
  if (n != TheSubBasis.cols()) {
    std::cerr << "The matrix size does not match\n";
    throw TerminalException{1};
  }
  MyMatrix<T> TheSubBasis_T = UniversalMatrixConversion<T,Tint>(TheSubBasis);
  MyMatrix<T> eGram = TheSubBasis_T * TheGramMat * TheSubBasis_T.transpose();
  MyMatrix<T> eGramInv = Inverse(eGram);

  MyMatrix<T> TheProj = TheSubBasis_T.transpose() * eGramInv * TheSubBasis_T * TheGramMat;
#ifdef DEBUG_RANKIN
  MyMatrix<T> TheProjSqr = TheProj * TheProj;
  if (TheProj * TheProj != TheProj) {
    std::cerr << "The matrix is not a projector\n";
    throw TerminalException{1};
  }
  if (TheProj * TheGramMat != TheGramMat * TheProj) {
    std::cerr << "The obtained projector is not self adjoint\n";
    throw TerminalException{1};
  }
#endif
  return TheProj;
}

template<typename T, typename Tint>
T UpperBoundRankinMinimalDeterminant(MyMatrix<T> const& TheGramMat, int k) {
  int n = TheGramMat.rows();
  T_shvec_info<T, Tint> SHVmin = computeMinimum_GramMat(TheGramMat);
  T rNorm = SHVmin.minimum;
  if (k == 1) {
    return rNorm;
  }
  MyVector<Tint> eVect = SHVmin.short_vectors[0];
  MyVector<T> eVect_T = UniversalVectorConversion<T,Tint>(eVect);
  MyVector<T> eVect_T_TheGramMat = TheGramMat * eVect_T.transpose();
  MyMatrix<T> eVect_M = MatrixFromVector(eVect);
  MyMatrix<Tint> TheCompl = SubspaceCompletionInt(eVect_M, n);
  MyMatrix<T> TheProj(n-1, n);
  for (int i=0; i<n-1; i++) {
    MyVector<Tint> fVect = GetMatrixRow(TheCompl, i);
    MyVector<T> fVect_T = UniversalVectorConversion<T,Tint>(fVect);
    T scal = (eVect_T_TheGramMat.dot(fVect_T)) / rNorm;
    MyVector<T> rVect = fVect - eVect * scal;
    AssignMatrixRow(TheProj, i, rVect);
  }
  MyMatrix<T> ReducedGramMat = TheProj * TheGramMat * TheProj.transpose();
  T upper = UpperBoundRankinMinimalDeterminant(ReducedGramMat, k-1);
  return rNorm * upper;
}



template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> Randin_k_level(MyMatrix<T> const& A, int const& k, T const& MaxDet);


template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> Rankin_k_level(MyMatrix<T> const& A, int const& k, T const& MaxDet) {
  // We use the HermiteNormalForm
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  if (k == 1) {
    T bound = MaxDet;
    T_shvec_info<T, Tint> SHVmin = computeLevel_GramMat(A, bound);
    std::vector<MyMatrix<Tint>> RetList;
    for (auto & eV : SHVmin.short_vectors) {
      MyMatrix<Tint> M = MatrixFromVector(eV);
      RetList.push_back(M);
    }
    return RetList;
  } else {
    // We are now using the Hermite constant to get a bound on the minimum
    // That is we have min(A)^k <= H(n) * MaxDet
    T upper = GetUpperBoundHermitePower(k) * MaxDet;
    // What could be done to get upper bound on min(A)?
    // 
    
    
  }
  std::vector<MyMatrix<Tint>> vec_subspaces;
  for (auto& mat : set_subspaces) {
    vec_subspaces.push_back(mat);
  }
  return vec_subspaces;
}




// clang-format off
#endif  // SRC_RANKIN_ENUMERATION_K_SPACE_H_
// clang-format on
