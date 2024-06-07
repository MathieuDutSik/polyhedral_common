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

#ifdef DEBUG
#define DEBUG_ISOTROPIC
#endif

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
template <typename T>
std::optional<MyVector<T>> GetIsotropIndefiniteLLL(MyMatrix<T> const &Q, [[maybe_unused]] std::ostream& os) {
  using Tint = typename underlying_ring<T>::ring_type;
  int n = Q.rows();
  auto get_norm = [&](MyMatrix<T> const &mat) -> T {
    T sum = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(mat(i, j));
    return sum;
  };
  MyMatrix<T> Pw = IdentityMat<T>(n);
  MyMatrix<T> Qw = Q;
  T curr_norm = get_norm(Q);
  while(true) {
    // Compute the LLL reduction.
    ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(Qw);
    if (!res.success) {
      MyVector<T> V = Pw.transpose() * res.Xisotrop;
#ifdef DEBUG_ISOTROPIC
      os << "ISOTROP: GetIsotropIndefiniteLLL finding isotrop in the process\n";
#endif
      return V;
    }
    // Compute the product
    MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(res.B);
    Pw = B_T * Pw;
    Qw = res.Mred;
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: GetIsotropIndefiniteLLL Qw=\n";
    WriteMatrix(os, Qw);
#endif
    // Checking first for diagonal zeros.
    for (int i = 0; i < n; i++) {
      if (Qw(i, i) == 0) {
        MyVector<T> eV = ZeroVector<T>(n);
        eV(i) = 1;
        MyVector<T> fV = Pw.transpose() * eV;
#ifdef DEBUG_ISOTROPIC
        os << "ISOTROP: GetIsotropIndefiniteLLL finding isotrop in the diagonal\n";
#endif
        return fV;
      }
    }
    // Checking for opposite signs in the diagonal
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        if (Qw(i, j) == 0 && Qw(i, i) + Qw(j, j) == 0) {
          MyVector<T> eV = ZeroVector<T>(n);
          eV(i) = 1;
          eV(j) = 1;
          MyVector<T> fV = Pw.transpose() * eV;
#ifdef DEBUG_ISOTROPIC
        os << "ISOTROP: GetIsotropIndefiniteLLL finding isotrop as sum of two diagonal terms\n";
#endif
          return fV;
        }
      }
    }
    T norm = get_norm(Qw);
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: GetIsotropIndefiniteLLL norm=" << norm << " curr_norm=" << curr_norm << "\n";
#endif
    if (norm > curr_norm) {
      break;
    }
    curr_norm = norm;
    MyMatrix<T> RandUnit = get_random_int_matrix<T>(n);
    Pw = RandUnit * Pw;
    Qw = RandUnit * Qw * RandUnit.transpose();
  }
  return {};
}



// For rank 1 this is trivial
template <typename T>
std::optional<MyVector<T>> FindIsotropicRankOne(MyMatrix<T> const &M) {
  if (M(0, 0) == 0)
    return {};
  MyVector<T> V(1);
  V(0) = 1;
  return V;
}

// For rank 2 this is trivial
template <typename T>
std::optional<MyVector<T>> FindIsotropicRankTwo(MyMatrix<T> const &M) {
  T a = M(0, 0);
  T b = M(0, 1);
  T c = M(1, 1);
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
    T const &sqrt_Delta = *opt;
    V(0) = -b + sqrt_Delta;
    V(1) = a;
    return V;
  }
  return {};
}

/*
  The algorithm is purely
 */
template <typename T> MyVector<T> Kernel_FindIsotropic(MyMatrix<T> const &Q, std::ostream& os) {
  int n = Q.rows();
#ifdef DEBUG_ISOTROPIC
  os << "ISOTROP: Before determinant minimization\n";
#endif
  ResultDetMin<T> res = DeterminantMinimization(Q, os);
#ifdef DEBUG_ISOTROPIC
  os << "ISOTROP: After determinant minimization\n";
  os << "ISOTROP: res.P=\n";
  WriteMatrix(os, res.P);
  os << "ISOTROP: res.Mred=\n";
  WriteMatrix(os, res.Mred);
  os << "ISOTROP: det(Mred)=" << DeterminantMat(res.Mred) << "\n";
#endif
  MyMatrix<T> Pw = res.P;
  MyMatrix<T> Qw = res.Mred;
  while (true) {
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: Before GetIsotropIndefiniteLLL\n";
#endif
    std::optional<MyVector<T>> opt = GetIsotropIndefiniteLLL(Qw, os);
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: We have opt\n";
#endif
    if (opt) {
      MyVector<T> const &eV = *opt;
#ifdef DEBUG_ISOTROPIC
      os << "ISOTROP: Qw=\n";
      WriteMatrix(os, Qw);
      os << "ISOTROP: eV=" << StringVector(eV) << "\n";
#endif
      MyVector<T> fV = Pw.transpose() * eV;
      return fV;
    }
    MyMatrix<T> U = get_random_int_matrix<T>(n);
    Pw = U * Pw;
    Qw = U * Qw * U.transpose();
  }
}

template <typename T>
std::optional<MyVector<T>> FindIsotropic(MyMatrix<T> const &M, std::ostream& os) {
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
#ifdef DEBUG_ISOTROPIC
    os << "ternary_has_isotropic_vector test=" << test << "\n";
#endif
    if (test) {
      return Kernel_FindIsotropic(M, os);
    } else {
      return {};
    }
  }
  if (n == 4) {
    // This is the case where we have a missing ingredient
    std::cerr << "The following call might fail because there is no isotropic "
                 "vectors\n";
    return Kernel_FindIsotropic(M, os);
  }
  // Now n is greater than 5, so there is an isotropic vector
  return Kernel_FindIsotropic(M, os);
}

// clang-format off
#endif  //  SRC_INDEFINITE_ISOTROPIC_H_
// clang-format on
