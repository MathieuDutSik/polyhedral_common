// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_PERFECTION_H_
#define SRC_RANKIN_PERFECTION_H_

template <typename T, typename Tint>
MyMatrix<T> GetProjectorMatrix(MyMatrix<Tint> const &P, MyMatrix<T> const &eG) {
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
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
  However, all of this is simply linear transformations on A, so we do not need
  to keep track of the A_P and simply work with A. writing Psi(A,dB) = det(A)
  Tr(A^{-1} dB) we have HessPhi(A + dA, dB) = det(A + dA) Tr((A+dA)^{-1} dB) =
  det(A) (1 + Tr(A^{-1}dA)) Tr((A^{-1} - A^{-1} dA A^{-1}) dB) = C + det(A) {
  Tr(A^{-1} dA) Tr(A^{-1} dB) - Tr(A^{-1} dA A^{-1} dB) }

  ---
  Now for the functional f that will get us
  df = sum_i 2 (Phi_{Pi}(A) - 1) dPhi_{Pi}
  and
  Hess(f) = sum_i 2 dPhi_{Pi}(A,deltaA) dPhi_{Pi}(A,deltaA)
          + sum_i 2 (Phi_{Pi}(A) - 1) HessPhi_{Pi}(deltaA, deltaA)
 */
template <typename T, typename Tint>
MyMatrix<T> FindApproximateMatrix(MyMatrix<T> const &A0,
                                  std::vector<MyMatrix<T>> const &B_k,
                                  std::vector<MyMatrix<Tint>> const &ListP,
                                  tol const &maxDeltaDet) {
  int dim_space = B_k.size();
  int n_iter = 0;
  std::vector<MyMatrix<T>> ListP_T;
  for (auto &eP : ListP) {
    MyMatrix<T> eP_T = UniversalMatrixConversion<T, Tint>(eP);
    ListP_T.push_back(EP_T);
  }
  int n_space = ListP.size();
  std::vector<std::vector<MyMatrix<T>>> ll_B_k_P;
  for (int i_space = 0; i_space < n_space; i_space++) {
    std::vector<MyMatrix<T>> l_B_k_P;
    for (int i_basis = 0; i_basis < dim_space; i_basis++) {
      MyMatrix<T> B_trans =
          ListP[i_space] * B_k[i_basis] * ListP[i_space].transpose();
      l_B_k_P.push_back(B_trans);
    }
    ll_B_k_P.push_back(l_B_k_P);
  }
  MyMatrix<T> A = A0;
  while (true) {
    T f = 0;
    MyVector<T> grad_f = ZeroVector<T>(dim_space);
    MyMatrix<T> Hess_f = ZeroMatrix<T>(dim_space, dim_space);
    //
    for (int i_space = 0; i_space < n_space; i_space++) {
      MyMatrix<T> const &P = ListP[i_space];
      MyMatrix<T> A_P = P * A * P.transpose();
      T det = DeterminantMat(A_P);
      MyMatrix<T> Ainv = Inverse(A_P);
      f += (det - 1) ^ 2;
      for (int i_basis = 0; i_basis < dim_space; i_basis++) {
        MyMatrix<T> const &B1 = ll_B_k_P[i_space][i_basis];
        MyMatrix<T> prod1 = Ainv * B1;
        grad_f(i_basis) += 2 * (det - 1) * det * prod1.trace();
        for (int j_basis = 0; j_basis <= i_basis; j_basis++) {
          MyMatrix<T> const &B2 = ll_B_k_P[i_space][j_basis];
          MyMatrix<T> prod2 = Ainv * B1;
          MyMatrix<T> prod12 = prod1 * prod2;
          T h1 = 2 * det * det * prod1.trace() * prod2.trace();
          T h2 = 2 * det * (det - 1) *
                 (prod1.trace() * prod2.trace() - prod12.trace());
          T h = h1 + h2;
          Hess_F(i_basis, j_basis) += h;
          if (i_basis != j_basis) {
            Hess_F(j_basis, i_basis) += h;
          }
        }
      }
    }
  }
}

// clang-format off
#endif  // SRC_RANKIN_PERFECTION_H_
// clang-format on
