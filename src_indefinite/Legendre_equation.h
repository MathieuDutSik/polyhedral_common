// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_LEGENDRE_EQUATION_H_
#define SRC_MATRIX_LEGENDRE_EQUATION_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "factorizations.h"
// clang-format on

/*
  The ternary equation x M x = 0 for L a 3x3 matrix
  is named the Legendre equation.
  The paper we use is LEGENDRE’S THEOREM, LEGRANGE’S DESCENT
  https://public.csusm.edu/aitken_html/notes/legendre.pdf
  Refered to as P1 below.
  -
  
 */



template<typename T>
MyVector<T> solve_ternary_equation(MyVector<T> const& aV) {
  MyMatrix<T> RedMatrix = IdentityMat<T>(3);
  while(true) {
    

    
  }



}





/*
  We apply Corollary 4 of P1 to case of abc is square-free.
  The conditions are:
  * a,b,c are not all of the same sign
  * −bc is a square modulo |a|
  * −ac is a square modulo |b|
  * −ab is a square modulo |c|.
 */
template<typename T>
bool determine_solvability(MyVector<T> const& aReduced) {
  size_t n_plus = 0;
  size_t n_minus = 0;
  for (int i=0; i<3; i++) {
    if (aReduced(i) > 0) {
      n_plus++;
    }
    if (aReduced(i) < 0) {
      n_minus++;
    }
  }
  if (n_plus == 0 && n_minus == 0) {
    return false;
  }
  T a = aReduced(0);
  T b = aReduced(1);
  T c = aReduced(2);
  T a_abs = T_abs(a);
  T b_abs = T_abs(b);
  T c_abs = T_abs(c);
  T a_cnt = - b * c;
  T b_cnt = - a * c;
  T c_cnt = - a * b;
  if (!is_quadratic_residue(a_cnt, a_abs)) {
    return false;
  }
  if (!is_quadratic_residue(b_cnt, b_abs)) {
    return false;
  }
  if (!is_quadratic_residue(c_cnt, c_abs)) {
    return false;
  }
  return false;
}


/*
  We apply some prime reduction so that no prime factor is shared between the coefficients.
  We return something so that:
  * The vector aRet expressing the matrix
  * The diagonal matrix B such that
  B Diag(aRet) B = u Diag(aV)  for some coefficient u.
 */
template<typename T>
std::pair<MyMatrix<T>, MyVector<T>> reduction_information(MyVector<T> const& aV) {
  T a = aV(0);
  T b = aV(1);
  T c = aV(2);
  std::map<T,size_t> a_map = FactorsIntMap(T_abs(a));
  std::map<T,size_t> b_map = FactorsIntMap(T_abs(b));
  std::map<T,size_t> c_map = FactorsIntMap(T_abs(c));
  MyMatrix<T> TransMat = IdentityMat<T>(3);
  //
  // Eliminating the even prime powers.
  //
  auto get_set=[&](std::map<T,size_t> const& map, int const& idx) -> std::set<T> {
    std::set<T> eset;
    for (auto & kv : a_map) {
      T const& val = kv.first;
      size_t mult = kv.second;
      size_t r = ResInt(mult, 2);
      size_t q = QuoInt(mult, 2);
      if (r > 0)
        eset.insert(val);
      TransMat(idx, idx) *= MyPow(val, q);
    }
    return eset;
  };
  std::set<T> a_set = get_set(a_map, 0);
  std::set<T> b_set = get_set(b_map, 1);
  std::set<T> c_set = get_set(c_map, 2);
  //
  // Now reducing so that no primes is shared
  //
  std::set<T> primes;
  for (auto & kv : a_map) {
    primes.insert(kv.first);
  }
  for (auto & kv : b_map) {
    primes.insert(kv.first);
  }
  for (auto & kv : c_map) {
    primes.insert(kv.first);
  }
  T a_prod(1);
  T b_prod(1);
  T c_prod(1);
  T ab_prod(1);
  T ac_prod(1);
  T bc_prod(1);
  for (auto & p : primes) {
    bool a_in = a_set.count(p) == 1;
    bool b_in = b_set.count(p) == 1;
    bool c_in = c_set.count(p) == 1;
    if (a_in && !b_in && !c_in) {
      a_prod *= p;
    }
    if (!a_in && b_in && c_in) {
      bc_prod *= p;
    }
    if (!a_in && b_in && !c_in) {
      b_prod *= p;
    }
    if (a_in && !b_in && c_in) {
      ac_prod *= p;
    }
    if (!a_in && !b_in && c_in) {
      c_prod *= p;
    }
    if (a_in && b_in && !c_in) {
      ab_prod *= p;
    }
  }
  TransMat(0,0) *= bc_prod;
  TransMat(1,1) *= ac_prod;
  TransMat(2,2) *= ab_prod;
  MyVector<T> aRet(3);
  aRet(0) = T_sign(a) * a_prod * bc_prod;
  aRet(1) = T_sign(b) * b_prod * ac_prod;
  aRet(2) = T_sign(c) * c_prod * ab_prod;
  return {TransMat, aRet};
}



template<typename T>
std::optional<MyVector<T>> get_ternary_isotropic_vector(MyMatrix<T> const& M) {
  

}





// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXFLT_H_
// clang-format on
