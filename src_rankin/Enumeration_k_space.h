// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_ENUMERATION_K_SPACE_H_
#define SRC_RANKIN_ENUMERATION_K_SPACE_H_

#include "Shvec_exact.h"

// Returns the Hermite constant at the n-th power.
// That is the minimum of min(A)^n / det(A)
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
std::vector<MyMatrix<Tint>> Compute_k_minimum(MyMatrix<T> const& A, int const& k);


template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> Compute_k_minimum(MyMatrix<T> const& A, int const& k, T const& MaxDet) {
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  if (k == 1) {
    
  } else {
    

    
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
