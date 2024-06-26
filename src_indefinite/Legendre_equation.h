// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_LEGENDRE_EQUATION_H_
#define SRC_INDEFINITE_LEGENDRE_EQUATION_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include "factorizations.h"
#include "Positivity.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_LEGENDRE
#endif

/*
  The ternary equation x M x = 0 for L a 3x3 matrix
  is named the Legendre equation.
  The paper we use is LEGENDRE’S THEOREM, LEGRANGE’S DESCENT
  https://public.csusm.edu/aitken_html/notes/legendre.pdf
  Refered to as P1 below.
  -
  Here we just put the test of existence.
 */

/*
  We apply Corollary 4 of P1 to case of abc is square-free.
  The conditions are:
  * a,b,c are not all of the same sign
  * −bc is a square modulo |a|
  * −ac is a square modulo |b|
  * −ab is a square modulo |c|.
 */
template <typename T> bool determine_solvability_dim3(MyVector<T> const &aReduced, [[maybe_unused]] std::ostream& os) {
  size_t n_plus = 0;
  size_t n_minus = 0;
  for (int i = 0; i < 3; i++) {
    if (aReduced(i) > 0) {
      n_plus++;
    }
    if (aReduced(i) < 0) {
      n_minus++;
    }
  }
#ifdef DEBUG_LEGENDRE
  os << "LEG: n_plus=" << n_plus << " n_minus=" << n_minus << "\n";
#endif
  if (n_plus == 0 || n_minus == 0) {
#ifdef DEBUG_LEGENDRE
    os << "LEG: Returning false by positivity condition\n";
#endif
    return false;
  }
  T a = aReduced(0);
  T b = aReduced(1);
  T c = aReduced(2);
#ifdef DEBUG_LEGENDRE
  T prod = a * b * c;
  os << "LEG: prod=" << prod << "\n";
#endif
  T a_abs = T_abs(a);
  T b_abs = T_abs(b);
  T c_abs = T_abs(c);
  T a_cnt = -b * c;
  T b_cnt = -a * c;
  T c_cnt = -a * b;
#ifdef DEBUG_LEGENDRE
  os << "LEG: a_abs=" << a_abs << " a_cnt=" << a_cnt << "\n";
  os << "LEG: b_abs=" << b_abs << " b_cnt=" << b_cnt << "\n";
  os << "LEG: c_abs=" << c_abs << " c_cnt=" << c_cnt << "\n";
#endif
  //
  // Test early termination by using Jacobi.
  // This is preferable because Jacobi is super fast to compute.
  //
  bool test_jacobi_a = compute_jacobi_symbol(a_cnt, a_abs);
  bool test_jacobi_b = compute_jacobi_symbol(b_cnt, b_abs);
  bool test_jacobi_c = compute_jacobi_symbol(c_cnt, c_abs);
#ifdef DEBUG_LEGENDRE
  os << "LEG: test_jacobi, a=" << test_jacobi_a << " b=" << test_jacobi_b << " c=" << test_jacobi_c << "\n";
#endif
  if (!test_jacobi_a || !test_jacobi_b || !test_jacobi_c) {
#ifdef DEBUG_LEGENDRE
    os << "LEG: early termination from Jacobi criterion\n";
#endif
    return false;
  }
  //
  bool test_a = is_quadratic_residue(a_cnt, a_abs);
#ifdef DEBUG_LEGENDRE
  os << "QUADTEST: [" << a_cnt << "," << a_abs << "," << test_a << "],\n";
#endif
  if (!test_a) {
#ifdef DEBUG_LEGENDRE
    os << "LEG: Returning false by quadratic residue for a\n";
#endif
    return false;
  }
  //
  //
  bool test_b = is_quadratic_residue(b_cnt, b_abs);
#ifdef DEBUG_LEGENDRE
  os << "QUADTEST: [" << b_cnt << "," << b_abs << "," << test_b << "],\n";
#endif
  if (!test_b) {
#ifdef DEBUG_LEGENDRE
    os << "LEG: Returning false by quadratic residue for b\n";
#endif
    return false;
  }
  //
  //
  bool test_c = is_quadratic_residue(c_cnt, c_abs);
#ifdef DEBUG_LEGENDRE
  os << "QUADTEST: [" << c_cnt << "," << c_abs << "," << test_c << "],\n";
#endif
  if (!test_c) {
#ifdef DEBUG_LEGENDRE
    os << "LEG: Returning false by quadratic residue for c\n";
#endif
    return false;
  }
#ifdef DEBUG_LEGENDRE
  os << "LEG: Returning true as the condition appears to be satisfied\n";
#endif
  return true;
}

/*
  We apply some prime reduction so that no prime factor is shared between the
  coefficients. We return something so that:
  * The vector aRet expressing the matrix
  * The diagonal matrix B such that
  B Diag(aRet) B = u Diag(aV)  for some coefficient u.
 */
template <typename T>
std::pair<MyMatrix<T>, MyVector<T>>
reduction_information(MyVector<T> const &aV) {
#ifdef DEBUG_LEGENDRE
  if (aV.size() != 3) {
    std::cerr << "The length should be exactly 3\n";
    throw TerminalException{1};
  }
#endif
  T a = aV(0);
  T b = aV(1);
  T c = aV(2);
  std::map<T, size_t> a_map = FactorsIntMap(T_abs(a));
  std::map<T, size_t> b_map = FactorsIntMap(T_abs(b));
  std::map<T, size_t> c_map = FactorsIntMap(T_abs(c));
  MyMatrix<T> TransMat = IdentityMat<T>(3);
  //
  // Eliminating the even prime powers.
  //
  auto get_set = [&](std::map<T, size_t> const &map,
                     int const &idx) -> std::set<T> {
    std::set<T> eset;
    for (auto &kv : map) {
      T const &val = kv.first;
      size_t mult = kv.second;
      size_t r = mult % 2;
      size_t q = mult / 2;
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
  for (auto &kv : a_map) {
    primes.insert(kv.first);
  }
  for (auto &kv : b_map) {
    primes.insert(kv.first);
  }
  for (auto &kv : c_map) {
    primes.insert(kv.first);
  }
  T a_prod(1);
  T b_prod(1);
  T c_prod(1);
  T ab_prod(1);
  T ac_prod(1);
  T bc_prod(1);
  for (auto &p : primes) {
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
  TransMat(0, 0) *= bc_prod;
  TransMat(1, 1) *= ac_prod;
  TransMat(2, 2) *= ab_prod;
  MyVector<T> aRet(3);
  aRet(0) = T_sign(a) * a_prod * bc_prod;
  aRet(1) = T_sign(b) * b_prod * ac_prod;
  aRet(2) = T_sign(c) * c_prod * ab_prod;
  return {TransMat, aRet};
}

template<typename T>
MyVector<typename underlying_ring<T>::ring_type> get_reduced_diagonal(MyMatrix<T> const &M, [[maybe_unused]] std::ostream& os) {
  using Tring = typename underlying_ring<T>::ring_type;
  DiagSymMat<T> dsm = DiagonalizeNonDegenerateSymmetricMatrix(M);
  MyVector<T> V1 = GetDiagonal(dsm.RedMat);
#ifdef DEBUG_LEGENDRE
  os << "V1=" << StringVectorGAP(V1) << "\n";
#endif
  MyVector<T> V2 = RemoveFractionVector(V1);
#ifdef DEBUG_LEGENDRE
  os << "V2=" << StringVectorGAP(V2) << "\n";
#endif
  MyVector<Tring> V3 = UniversalVectorConversion<Tring, T>(V2);
#ifdef DEBUG_LEGENDRE
  os << "V3=" << StringVectorGAP(V3) << "\n";
#endif
  return V3;
}



template <typename T> bool ternary_has_isotropic_vector(MyMatrix<T> const &M, std::ostream& os) {
  using Tring = typename underlying_ring<T>::ring_type;
  MyVector<Tring> red_diag_A = get_reduced_diagonal(M, os);
  MyVector<Tring> red_diag_B = reduction_information(red_diag_A).second;
  return determine_solvability_dim3(red_diag_B, os);
}

// clang-format off
#endif  // SRC_INDEFINITE_LEGENDRE_EQUATION_H_
// clang-format on
