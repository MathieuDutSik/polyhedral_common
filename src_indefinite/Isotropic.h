// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_ISOTROPIC_H_
#define SRC_INDEFINITE_ISOTROPIC_H_

// clang-format off
#include "DeterminantMinimization.h"
#include "Legendre_equation.h"
#include "Indefinite_LLL.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

/*
  The method here is inspired by the following
  papers and manuscripts:
  A) "Quadratic equations in dimensions 4, 5 and more" (P1)
  Denis Simon
  B) "Solving Quadratic Equations using reduced unimodular
  quadratic forms" (P2)
  Denis Simon
  ---
  The existing implementations that inspired us:
  1) Hecke.jl : src/QuadForm/Quad/Spaces.jl
  2) Pari: src/basemath/qfsolve.c
  ---
 */


// Try to find isotropic subspace by using Indefinite LLL
template<typename T>
std::optional<MyVector<T>> GetIsotropIndefiniteLLL(MyMatrix<T> const& M) {
  using Tint = typename underlying_ring<T>::ring_type;
  int n = M.rows();
  // Compute the LLL reduction.
  ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T,Tint>(M);
  if (!res.success) {
    return res.Xisotrop;
  }
  MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(res.B);
  MyMatrix<T> Binv_T = Inverse(B_T);
  // Compute the product
  MyMatrix<T> const& Mred = res.Mred;
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      if (Mred(i,j) == 0 && Mred(i,i) + Mred(j, j) == 0) {
        MyVector<T> eV = ZeroVector<T>(n);
        eV(i) = 1;
        eV(j) = 1;
        MyVector<T> fV = Binv_T.transpose() * eV;
        return fV;
      }
    }
  }
  return {};
}




// For rank 1 this is trivial
template<typename T>
std::optional<MyVector<T>> FindIsotropicRankOne(MyMatrix<T> const& M) {
  if (M(0,0) == 0)
    return {};
  MyVector<T> V(1);
  V(0) = 1;
  return V;
}



// For rank 2 this is trivial
template<typename T>
std::optional<MyVector<T>> FindIsotropicRankTwo(MyMatrix<T> const& M) {
  T a = M(0,0);
  T b = M(0,1);
  T c = M(1,1);
  // quadratic form q(x,y) = a x^2 + 2 b xy + c y^2
  // equation q(x,y) = 0
  MyVector<T> V(2);
  // Special case a = 0
  if (a == 0) {
    V(0) = 1;
    return V;
  }
  // Now set y = 1
  // Solutions are
  // x = (-2b + sqrt(Delta)) / 2a
  // with Delta = (2b)^2 - 4ac
  // which leads to
  // x = (-b + sqrt(b^2 - ac)) / a  and y = 1
  // or
  // x = -b + sqrt(b^2 - ac)  and y = a
  T Delta = b * b - a * c;
  std::optional<T> opt = UniversalSquareRoot(Delta);
  if (opt) {
    T const& sqrt_Delta = *opt;
    V(0) = -b + sqrt_Delta;
    V(1) = a;
    return V;
  }
  return {};
}








/*
  The algorithm is purely 
 */
template<typename T>
MyVector<T> Kernel_FindIsotropic(MyMatrix<T> const& Q) {
  int n = Q.rows();
  ResultDetMin<T> res = DeterminantMinimization(Q);
  MyMatrix<T> Pw = res.P;
  MyMatrix<T> Qw = res.Mred;
  while (true) {
    MyMatrix<T> U = get_random_int_matrix<T>(n);
    Pw = U * Pw;
    Qw = U * Qw * U.transpose();
    std::optional<MyVector<T>> opt = GetIsotropIndefiniteLLL(Qw);
    if (opt) {
      MyVector<T> const& eV = *opt;
      MyMatrix<T> Pw_inv = Inverse(Pw);
      MyVector<T> fV = Pw_inv.transpose() * eV;
      return fV;
    }
  }
}




template<typename T>
std::optional<MyVector<T>> FindIsotropic(MyMatrix<T> const& M) {
  int n = M.rows();
  if (n == 1) {
    return FindIsotropicRankOne(M);
  }
  if (n == 2) {
    return FindIsotropicRankTwo(M);
  }
  if (n == 3) {
    // Apply first the Legendre theorem
    bool test = ternary_has_isotropic_vector(M);
    if (test) {
      return Kernel_FindIsotropic(M);
    } else {
      return {};
    }
  }
  if (n == 4) {
    // This is the case where we have a missing ingredient
    std::cerr << "The following call might fail because there is no isotropic vectors\n";
    return Kernel_FindIsotropic(M);
  }
  // Now n is greater than 5, so there is an isotropic vector
  return Kernel_FindIsotropic(M);
}

// clang-format off
#endif  //  SRC_INDEFINITE_ISOTROPIC_H_
// clang-format on
