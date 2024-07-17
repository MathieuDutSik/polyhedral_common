// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_FIFTHANDHIGHERORDER_H_
#define SRC_INDEFINITE_FIFTHANDHIGHERORDER_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include "factorizations.h"
#include "Positivity.h"
#include "Quaternary.h"
#include "Legendre_equation.h"
#include "NumberTheoryPadic.h"
// clang-format on

/*
  We follow here "Algorithms for solving rational quadratic forms"
  by Josef Schicho and Jana Pilnikova for the 5th order.
  It is at least a rigorous method.
  ---
  For higher dimensions, we actually find a right subspace by the
  usual techniques. This is different from SP which uses a strict
  recursion. The technique is still rigorous.
 */

#ifdef DEBUG
#define DEBUG_FIFTH_AND_HIGHER_ORDER
#endif




/*
  By Meyer theorem, we know that there is one isotropic vector.
 */
template<typename T>
MyVector<T> solve_fifth_equation(MyVector<T> const& a, std::ostream& os) {
  T t(1);
  // Quadratic form -a1 x1^2 - a2 x2^2 + t z^2 = 0
  std::vector<T> v_tern{-a(0), -a(1)};
  std::vector<int> poss_tern = get_poss_signs(v_tern);
  // Quadratic form a3 x3^2 + a4 x4^2 + a5 x5^2 + t z^2 = 0
  std::vector<T> v_quad{a(2), a(3), a(4)};
  std::vector<int> poss_quad = get_poss_signs(v_quad);
  std::vector<int> poss = IntersectionVect(poss_tern, poss_quad);
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  if (poss.size() == 0) {
    std::cerr << "FIFTH: No possibilities\n";
    throw TerminalException{1};
  }
#endif
  int sign = poss[0];
  T t(sign);
  MyVector<T> a_tern(3);
  a_tern(0) = a(0);
  a_tern(1) = a(1);
  MyVector<T> a_quad(4);
  a_tern(0) = a(2);
  a_tern(1) = a(3);
  a_tern(2) = a(4);
  while(true) {
    a_tern(2) = -t;
    a_quad(3) = t;
    std::optional<MyVector<T>> opt_tern = TernaryIsotropicVectorDiagonal(a_tern, os);
    if (opt_tern) {
      MyVector<T> iso_tern = *opt_tern;
      std::optional<MyVector<T>> opt_quad = QuaternaryIsotropicVectorDiagonal(a_quad, os);
      if (opt_quad) {
        MyVector<T> iso_quad = *opt_quad;
        T val_tern = iso_tern(2);
        T val_quad = iso_quad(3);
        MyVector<T> iso_fifth(5);
        for (int u=0; u<2; u++) {
          iso_fifth(u) = iso_tern(u) * val_quad;
        }
        for (int u=2; u<5; u++) {
          iso_fifth(u) = iso_quad(u-2) * val_tern;
        }
        return iso_fifth;
      }
    }
  }
}


/*
  Compute a basis of the space 
 */
template<typename T>
MyMatrix<T> compute_fifth_basis(MyMatrix<T> const& Q) {
  int n = M.rows();
  size_t max_one_sign = 4;
  std::vector<MyVector<T>> basis;
  MyMatrix<T> M;
  auto is_linearly_independent=[&](MyVector<T> const& v) -> bool {
    if (M.rows() == 0) {
      return true;
    }
    std::optional<MyVector<T>> opt = SolutionMat(M, v);
    return !opt.has_value();
  };
  auto insert_vector=[&](MyVector<T> const& v) -> void {
    basis.push_back(v);
    M = MatrixFromVectorFamily(basis);
  };
  auto test_vector=[&](MyVector<T> const& v) -> bool {
    if (basis.size() == 5) {
      return false;
    }
    if (!is_linearly_independent(v)) {
      return false;
    }
    std::vector<T> test_basis = basis;
    test_basis.push_back(v);
    MyMatrix<T> test_M = MatrixFromVectorFamily(test_basis);
    MyMatrix<T> QuadM = test_M * Q * test_M.transpose();
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(QuadM);
    if (DiagInfo.nbPlus > max_one_sign) {
      return false;
    }
    if (DiagInfo.nbMinus > max_one_sign) {
      return false;
    }
    return true;
  };
  auto insert_vector_if_ok=[&](MyVector<T> const& v) {
    if (test_vector(v)) {
      insert_vector(v);
    }
  };
  for (int i=0; i<n; i++) {
    MyVector<T> v = ZeroVector<T>(n);
    v(i) = 1;
    insert_vector_if_ok(v);
  }
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      MyVector<T> v = ZeroVector<T>(n);
      v(i) = 1;
      v(j) = 1;
      insert_vector_if_ok(v);
    }
  }
  auto get_sign_vector=[&](int const& sign) -> Myvector<T> {
    MyMatrix<T> ProdM = M * Q;
    MyMatrix<T> NSP = NullspaceTrMat(ProdM);
    MyMatrix<T> QuadQ = sign * NSP * Q * NSP.transpose();
    T CritNorm(0);
    bool StrictIneq = true;
    bool NeedNonZero = true;
    MyVector<T> v1 = GetIntegralVector_allmeth(QuadM, CritNorm, StrictIneq, NeedNonZero, os);
    MyVector<T> v2 = NSP.transpose() * v1;
    return v2;
  };
  while(true) {
    if (M.rows() == 5) {
      return M;
    }
    MyMatrix<T> QuadM = M * Q * M.transpose();
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(QuadM);
    if (DiagInfo.nbPlus < max_one_sign) {
      MyVector<T> v = get_sign_vector(1);
      insert_vector(v);
    } else {
      if (DiagInfo.nbMinus < max_one_sign) {
        MyVector<T> v = get_sign_vector(-1);
        insert_vector(v);
      }
    }
  }
}








// clang-format off
#endif  // FIFTHANDHIGHERORDER_H_
// clang-format on
