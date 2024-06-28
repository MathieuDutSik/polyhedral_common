// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_QUATERNARY_H_
#define SRC_INDEFINITE_QUATERNARY_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include "factorizations.h"
#include "Positivity.h"
#include "Legendre_equation.h"
#include "NumberTheoryPadic.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_QUATERNARY
#endif

/*
  We follow here "Algorithms for solving rational quadratic forms"
  by Josef Schicho and Jana Pilnikova
  ---
  We only provide the determination of whether they are isotropic or not.
 */


/*
  We should have a1 a2 a3 is square free.
  We combine Lemma 7 and Lemma 8.
 */
template<typename T>
bool Padic_isotropy_ternary(MyVector<T> const& a, T const& p) {
  let two(2);
  int miss_val = -1;
  auto get_idx=[&] -> int {
    for (int i=0; i<3; i++) {
      T res = ResInt(a(i), p);
      if (res == 0) {
        return i;
      }
    }
    return miss_val;
  };
  if (p != two) {
    int pos = get_idx();
    if (pos == miss_val) {
      // Lemma 7
      return true;
    }
    int pos1 = ResInt(pos+1, 3);
    int pos2 = ResInt(pos+2, 3);
    T coeff1 = a(pos1);
    T coeff2 = a(pos2);
    T coeff1_inv = mod_inv(coeff1, p);
    T a = - coeff2 * coeff1_inv;
    // Lemma 8, (i)
    return is_quadratic_residue(a, p);
  } else {
    int pos = get_idx();
    if (pos == miss_val) {
      // Lemma 8, (iii)
      T four(4);
      for (int pos1=0; pos1<3; pos1++) {
        int pos2 = ResInt(pos1+1,3);
        T sum = a(pos1) + a(pos2]);
        T res = ResInt(sum, four);
        if (res == 0) {
          return true;
        }
      }
      return false;
    } else {
      // Lemma 8, (ii)
      int pos1 = ResInt(pos+1, 3);
      int pos2 = ResInt(pos+2, 3);
      T a1 = a(pos1);
      T a2 = a(pos2);
      T a3 = a(pos);
      T eight(8);
      for (int s=0; s<2; s++) {
        T sum = a1 + a2 + a3 * s;
        T res = ResInt(sum, eight);
        if (res == 0) {
          return true;
        }
      }
      return false;
    }
  }
}


/*
  This is the main theorem of SP. We do not go over trying to
  find an explicit solution.
 */
template<typename T>
bool determine_solvability_dim4(MyVector<T> const& a) {
  size_t n_plus = 0;
  size_t n_minus = 0;
  for (int i=0; i<4; i++) {
    if (a(i) > 0) {
      n_plus += 1;
    }
    if (a(i) < 0) {
      n_minus += 1;
    }
  }
#ifdef DEBUG_QUATERNARY
  std::cerr << "QUAD: n_plus=" << n_plus << " n_minus=" << n_minus << "\n";
#endif
  if (n_plus == 0 || n_minus == 0) {
#ifdef DEBUG_QUATERNARY
    std::cerr << "QUAD: Returning false by positivity condition\n";
#endif
    return false;
  }
  std::set<T> primes;
  T two(2);
  primes.insert(two);
  for (int i=0; i<4; i++) {
    std::map<T, size_t> map = FactorsIntMap(T_abs(a(i)));
    for (auto & kv : map) {
      if (kv.second > 0) {
        T const& p = kv.first;
        primes.insert(p);
      }
    }
  }
  //
  T a1 = a(0);
  T a2 = a(1);
  T a3 = a(2);
  T a4 = a(3);
  MyVector<T> a12(3);
  MyVector<T> a12(0) = a1;
  MyVector<T> a12(1) = a2;
  MyVector<T> a12(2) = -1;
  MyVector<T> a12_red = reduction_information(a12).second;
  MyVector<T> a34(3);
  MyVector<T> a34(0) = -a3;
  MyVector<T> a34(1) = -a4;
  MyVector<T> a34(2) = -1;
  MyVector<T> a34_red = reduction_information(a34).second;
  /*
    This is Lemma 13 of SP.
   */
  auto Padic_anisotropy_quaternary=[&](p const& val) -> bool {
    T prod12_A = - a1 * a2;
    Padic<T> prod12_B = Padic_from_integer(prod12_A, p);
    Padic<T> prod12_C = Padic_reduce_precision(prod12_B, 3);
    T prod34_A = - a3 * a4;
    Padic<T> prod34_B = Padic_from_integer(prod34_A, p);
    Padic<T> prod34_C = Padic_reduce_precision(prod34_B, 3);
    Padic<T> prod34inv = Padic_inverse(prod34_C, p);
    Padic<T> proc12_34inv = Padic_product(prod12_C, prod34inv, p);
    bool test = Padic_is_square(prod12_34inv, p);
    if (!test) {
      return false;
    }
    auto get_hilbert_a12=[&]() -> int {
      bool test12 = Padic_isotropy_ternary(a12_red, p);
      if (test12) {
        return 1;
      } else {
        return -1;
      }
    };
    auto get_minus_hilbert_a34=[&]() -> int {
      bool test34 = Padic_isotropy_ternary(a34_red, p);
      if (test34) {
        return -1;
      } else {
        return 1;
      }
    };
    int hilbertA = get_hilbert_a12();
    int hilbertB = get_minus_hilbert_a34();
    return hilbertA == hilbertB;
  };
  for (auto & p : primes) {
    bool test = Padic_anisotropy_quaternary(p);
    if (test) {
      // The Lemma is about anisotropy not isotropy
      return false;
    }
  }
  return true;
}

template <typename T> bool quaternary_has_isotropic_vector(MyMatrix<T> const &M) {
  using Tring = typename underlying_ring<T>::ring_type;
  MyVector<Tring> red_diag = get_reduced_diagonal(M);
  return determine_solvability_dim4(red_diag);
}


// clang-format off
#endif  // SRC_INDEFINITE_QUATERNARY_H_
// clang-format on
