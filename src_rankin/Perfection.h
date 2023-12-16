// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_PERFECTION_H_
#define SRC_RANKIN_PERFECTION_H_

template<typename T, typename Tint>
MyMatrix<T> GetProjectorMatrix(MyMatrix<Tint> const& P, MyMatrix<T> const& eG) {
  MyMatrix<T> P_T = UniversalMatrixConversion<T,Tint>(P);
  MyMatrix<T> eG_red = P_T * eG * P_T.transpose();
  MyMatrix<T> eGinv = Inverse(eG_red);
  MyMatrix<T> ProjMat = P_T.transpose() * eGinv * P_T;
  return ProjMat;
}


/*
  we define the function Phi_{Pi}(A) = det PAP^T
  We have the equations Phi_{Pi}(A) = 1 for A = A0 + sum_k B_k
  In truth we want to minimize the function
  f = sum_i (Phi_{Pi}(A) -1)^2
  because we could have a number of conditions on the P higher than the
  number of matrices B_k, in other words the problem could be overdetermined
  due to degenerate facets.
  ---
  We need the first and second differentials.
  det(A + dA) = det(A(I_n + A^{-1} dA))
              = det(A) det(I_n + A^{-1} dA)
              = det(A) (1 + Tr(A^{-1} dA))
  So, the differential is
  det(A) Tr(A^{-1} dA)
  Writing A_P = P A P^T we get
  Phi_P(A+dA) = det(A_P + dA_P)
              = det(A_P) Tr(A_P^{-1} dA_P)
  However, all of this is simply linear transformations on A, so we do not need to
  keep track of the A_P and simply work with A.
  writing Psi(A,dB) = det(A) Tr(A^{-1} dB)
  we have
  Phi(A + dA, dB) = det(A + dA) Tr((A+dA)^{-1} dB)
                  = det(A) (1 + Tr(A^{-1}dA)) Tr((A^{-1} - A^{-1} dA A^{-1}) dB)
                  = C + det(A) { Tr(A^{-1} dA) Tr(A^{-1} dB) - Tr(A^{-1} dA A^{-1} dB) }
 */
template<typename T, typename Tint>
MyMatrix<T> FindApproximateMatrix(MyMatrix<T> const& M, std::vector<MyMatrix<Tint>> const& ListP) {

  
}










// clang-format off
#endif  // SRC_RANKIN_PERFECTION_H_
// clang-format on
