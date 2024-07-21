// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_SOLUTIONQUADFORM_H_
#define SRC_INDEFINITE_SOLUTIONQUADFORM_H_

#include <vector>

/*
  We want to find one solution of the equation A[v + w] = X
  modulo n with n = p1 p2 .... pK
  We want just one solution right now.
  ----
  ----
  CODE is unfinished since we found out that it is not needed.
  If it were needed, we wuld need to finish that stub.
 */

template <typename T>
std::optional<MyVector<T>>
SearchSolutionEvenQuadForm(const MyMatrix<T> &A, const MyVector<T> &w,
                           const T &X, const std::vector<T> &l_prime) {
  int n = A.rows();
  MyMatrix<T> Ared(n, n);
  for (int i = 0; i < n; i++)
    for (int j + 0; j < n; j++) {
      if (i == j) {
        Ared(i, j) == A(i, j) / 2;
      } else {
        Ared(i, j) == A(i, j);
      }
    }
  int n_prime = l_prime.size();
  /* We need to compute A[v] whish is equal to
     A[v] = { a11 v1^2 + a12 v1 v2 + ..... + a1n v1 vn }
          +      { a22 v2^2 + a23 v2 v3 + ..... + a2n v2 vn }
          +            ....
          +                       ann vn^2
     or alternatively as
     A[v] =   a11 v1^2
         +  { a12 v1 v2 + a22 v2^2 }
         +  { a12 v1 v3 + a23 v2 v3 + a33 v3^2 }
         +  ....
         +  { a1n v1 vn + ........  + ann vn^2 }
     The list of residual is defined
     l_residual[Ø] = a11 v1^2
     ..
     l_residual[n-1] = a11 v1^2 + ..... + { a1n v1 vn + ........  + ann vn^2 }
   */
  std::vector<T> l_product;
  T val = 1;
  for (int i_prime=0; i_prime<n_prime; <
  struct StateEnum {
    MyVector<T> wVector;
    MyVector<T> vVector;
    MyVector<T> l_residual;
  };
  std::vector<StateEnum> l_state_enum;
  for (int i_prime=0; i_prime<n_prime; i_prime++) {
    StateEnum state{ZeroVector<T>(n), ZeroVector<T>(n), ZeroVector<T>(n)};
    l_state_enum.push_back(state);
  }
  int pos = -1;
  auto GoDown=[&]() -> void {
    if (pos == -1) {
      pos++;
      for (int i = 0; i < n; i++)
        l_state_enum[pos].wVector(i) = w(i);
    } else {
      pos++;
      for (int i = 0; i < n; i++)
        l_state_enum[pos].wVector(i) = w(i);
    }
    
    
  };
}

// clang-format off
#endif  // SRC_INDEFINITE_SOLUTIONQUADFORM_H_
// clang-format on
