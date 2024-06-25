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
bool Padic_isotropy_ternary(std::vector<T> const& a, T const& p) {
  let two(2);
  int miss_val = -1;
  auto get_idx=[&] -> int {
    for (int i=0; i<3; i++) {
      T res = ResInt(a[i], p);
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
    T coeff1 = a[pos1];
    T coeff2 = a[pos2];
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
        T sum = a[pos1] + a[pos2];
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
      T a1 = a[pos1];
      T a2 = a[pos2];
      T a3 = a[pos];
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






template<typename T>
bool determine_solvability_dim4(MyVector<T> const& aReduced) {
  size_t n_plus = 0;
  size_t n_minus = 0;
  for (int i=0; i<4; i++) {
    if (aReduced(i) > 0) {
      n_plus += 1;
    }
    if (aReduced(i) < 0) {
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
    std::map<T, size_t> map = FactorsIntMap(T_abs(aReduced(i)));
    for (auto & kv : map) {
      if (kv.second > 0) {
        T const& p = kv.first;
        primes.insert(p);
      }
    }
  }
  //
  
}







// clang-format off
#endif  // SRC_INDEFINITE_QUATERNARY_H_
// clang-format on
