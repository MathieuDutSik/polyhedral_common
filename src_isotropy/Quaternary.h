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

#ifdef TIMINGS
#define TIMINGS_QUATERNARY
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
template <typename T>
bool Padic_isotropy_ternary(MyVector<T> const &a, T const &p,
                            [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_QUATERNARY
  os << "QUAD: Padic_isotropy_ternary, beginning\n";
#endif
  T two(2);
  int miss_val = -1;
  auto get_idx = [&]() -> int {
    for (int i = 0; i < 3; i++) {
      T res = ResInt(a(i), p);
      if (res == 0) {
        return i;
      }
    }
    return miss_val;
  };
  if (p != two) {
    int pos = get_idx();
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, pos=" << pos << "\n";
#endif
    if (pos == miss_val) {
      // Lemma 7
      return true;
    }
    int pos1 = ResInt(pos + 1, 3);
    int pos2 = ResInt(pos + 2, 3);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, pos1=" << pos1 << " pos2=" << pos2
       << "\n";
#endif
    T coeff1 = ResInt(a(pos1), p);
    T coeff2 = a(pos2);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, coeff1=" << coeff1
       << " coeff2=" << coeff2 << "\n";
#endif
    T coeff1_inv = mod_inv(coeff1, p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, coeff1_inv=" << coeff1_inv << "\n";
#endif
    T a = -coeff2 * coeff1_inv;
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, a=" << a << "\n";
#endif
    // Lemma 8, (i)
    bool test = is_quadratic_residue(a, p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: Padic_isotropy_ternary, test=" << test << "\n";
#endif
    return test;
  } else {
    int pos = get_idx();
    if (pos == miss_val) {
      // Lemma 8, (iii)
      T four(4);
      for (int pos1 = 0; pos1 < 3; pos1++) {
        int pos2 = ResInt(pos1 + 1, 3);
        T sum = a(pos1) + a(pos2);
        T res = ResInt(sum, four);
        if (res == 0) {
          return true;
        }
      }
      return false;
    } else {
      // Lemma 8, (ii)
      int pos1 = ResInt(pos + 1, 3);
      int pos2 = ResInt(pos + 2, 3);
      T a1 = a(pos1);
      T a2 = a(pos2);
      T a3 = a(pos);
      T eight(8);
      for (int s = 0; s < 2; s++) {
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

template <typename T>
std::vector<T> get_local_primes(MyVector<T> const &a,
                                [[maybe_unused]] std::ostream &os) {
  std::set<T> primes;
  T two(2);
  primes.insert(two);
  for (int i = 0; i < 4; i++) {
#ifdef DEBUG_QUATERNARY
    os << "QUAD: get_local_primes i=" << i << " a(i)=" << a(i) << "\n";
#endif
    std::map<T, size_t> map = FactorsIntMap(T_abs(a(i)));
    for (auto &kv : map) {
      if (kv.second > 0) {
        T const &p = kv.first;
        primes.insert(p);
      }
    }
  }
  std::vector<T> v;
  for (auto &p : primes) {
    v.push_back(p);
  }
  return v;
}

/*
  This is the main theorem of SP. We do not go over trying to
  find an explicit solution.
 */
template <typename T>
bool determine_solvability_dim4(MyVector<T> const &a,
                                std::vector<T> const &primes,
                                std::ostream &os) {
#ifdef DEBUG_QUATERNARY
  os << "QUAD: determine_solvability_dim4, beginning |a|=" << a.size() << "\n";
#endif
  size_t n_plus = 0;
  size_t n_minus = 0;
  for (int i = 0; i < 4; i++) {
    if (a(i) > 0) {
      n_plus += 1;
    }
    if (a(i) < 0) {
      n_minus += 1;
    }
  }
#ifdef DEBUG_QUATERNARY
  os << "QUAD: n_plus=" << n_plus << " n_minus=" << n_minus << "\n";
#endif
  if (n_plus == 0 || n_minus == 0) {
#ifdef DEBUG_QUATERNARY
    std::cerr << "QUAD: Returning false by positivity condition\n";
#endif
    return false;
  }
  //
  T a1 = a(0);
  T a2 = a(1);
  T a3 = a(2);
  T a4 = a(3);
#ifdef DEBUG_QUATERNARY
  os << "QUAD: a =" << a1 << " " << a2 << " " << a3 << " " << a4 << "\n";
#endif
  MyVector<T> a12(3);
  a12(0) = a1;
  a12(1) = a2;
  a12(2) = -1;
  MyVector<T> a12_red = reduction_information(a12, os).aReduced;
#ifdef DEBUG_QUATERNARY
  os << "QUAD: a12_red=" << StringVector(a12_red) << "\n";
#endif
  MyVector<T> a34(3);
  a34(0) = -a3;
  a34(1) = -a4;
  a34(2) = -1;
  MyVector<T> a34_red = reduction_information(a34, os).aReduced;
#ifdef DEBUG_QUATERNARY
  os << "QUAD: a34_red=" << StringVector(a34_red) << "\n";
#endif
  /*
    This is Lemma 13 of SP.
   */
  auto Padic_anisotropy_quaternary = [&](T const &p) -> bool {
#ifdef DEBUG_QUATERNARY
    os << "QUAD: ---------------- p=" << p << " -----------------------\n";
#endif
    T prod12_A = -a1 * a2;
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod12_A=" << prod12_A << "\n";
#endif
    Padic<T> prod12_B = Padic_from_integer(prod12_A, p, 3);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod12_B=" << Padic_to_string(prod12_B) << "\n";
#endif
    T prod34_A = -a3 * a4;
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod34_A=" << prod34_A << "\n";
#endif
    Padic<T> prod34_B = Padic_from_integer(prod34_A, p, 3);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod34_B=" << Padic_to_string(prod34_B) << "\n";
#endif
    Padic<T> prod34inv = Padic_inverse(prod34_B, p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod34inv=" << Padic_to_string(prod34inv) << "\n";
#endif
    Padic<T> prod12_34inv = Padic_product(prod12_B, prod34inv, p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: prod12_34inv=" << Padic_to_string(prod12_34inv) << "\n";
#endif
    bool test = Padic_is_square(prod12_34inv, p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: test=" << test << "\n";
#endif
    if (!test) {
      return false;
    }
    auto get_hilbert_a12 = [&]() -> int {
      bool test12 = Padic_isotropy_ternary(a12_red, p, os);
#ifdef DEBUG_QUATERNARY
      os << "QUAD: test12=" << test12 << "\n";
#endif
      if (test12) {
        return 1;
      } else {
        return -1;
      }
    };
    auto get_minus_hilbert_a34 = [&]() -> int {
      bool test34 = Padic_isotropy_ternary(a34_red, p, os);
#ifdef DEBUG_QUATERNARY
      os << "QUAD: test34=" << test34 << "\n";
#endif
      if (test34) {
        return -1;
      } else {
        return 1;
      }
    };
    int hilbertA = get_hilbert_a12();
    int hilbertB = get_minus_hilbert_a34();
#ifdef DEBUG_QUATERNARY
    os << "QUAD: hilbertA=" << hilbertA << " hilbertB=" << hilbertB << "\n";
#endif
    return hilbertA == hilbertB;
  };
  for (auto &p : primes) {
    bool test = Padic_anisotropy_quaternary(p);
#ifdef DEBUG_QUATERNARY
    os << "QUAD: p=" << p << " Padic_anisotropy_quaternary=" << test << "\n";
#endif
    if (test) {
      // The Lemma is about anisotropy not isotropy
      return false;
    }
  }
  return true;
}

template <typename T>
int hilbert_symbol(T const &a, T const &b, T const &p, std::ostream &os) {
  MyVector<T> aV(3);
  aV(0) = a;
  aV(1) = b;
  aV(2) = -1;
  MyVector<T> aV_red = reduction_information(aV, os).aReduced;
  bool test = Padic_isotropy_ternary(aV_red, p, os);
  if (test) {
    return 1;
  } else {
    return -1;
  }
}

template <typename T>
std::vector<T> get_tp_classes(MyVector<T> const &a,
                              std::vector<T> const &primes, std::ostream &os) {
  std::vector<T> classes;
  T a1 = a(0);
  T a2 = a(1);
  T a3 = a(2);
  T a4 = a(3);
  for (auto &p : primes) {
    auto is_matching_a_b = [&](T const &t, T const &a, T const &b) -> bool {
      // The t should satisfy (t, -a a) = (a, a)
      T prod = -a * b;
      int val1 = hilbert_symbol(t, prod, p, os);
      int val2 = hilbert_symbol(a, b, p, os);
      return val1 == val2;
    };
    std::vector<T> candidates = Padic_get_residue_classes(p);
#ifdef DEBUG_QUATERNARY
    size_t n_cond_12 = 0, n_cond_34 = 0;
#endif
    std::vector<T> matches;
    for (auto &t : candidates) {
      bool test12 = is_matching_a_b(t, a1, a2);
      bool test34 = is_matching_a_b(t, -a3, -a4);
#ifdef DEBUG_QUATERNARY
      if (test12) {
        n_cond_12 += 1;
      }
      if (test34) {
        n_cond_34 += 1;
      }
#endif
      if (test12 && test34) {
        matches.push_back(t);
      }
    }
#ifdef DEBUG_QUATERNARY
    os << "QUAD: p=" << p << "\n";
    os << "QUAD: n_cond_12=" << n_cond_12 << " n_cond_34=" << n_cond_34
       << " |matches|=" << matches.size() << "\n";
    if (matches.size() == 0) {
      std::cerr << "QUAD: matches should not be empty\n";
      throw TerminalException{1};
    }
#endif
    T tp = matches[0];
    classes.push_back(tp);
  }
  return classes;
}

template <typename T> std::vector<int> get_poss_signs(std::vector<T> const &a) {
  // Quadratic form is sum a_i x_i^2 + t z^2
  // This returns the possible signs for t
  size_t len = a.size();
  size_t n_plus = 0;
  for (auto &val : a) {
    if (val > 0) {
      n_plus += 1;
    }
  }
  if (n_plus == 0) {
    return {1};
  }
  if (n_plus == len) {
    return {-1};
  }
  return {-1, 1};
}

// The q is defined from section 4.2, item 4.
template <typename T>
T get_q_val(MyVector<T> const &a, std::vector<T> const &classes,
            std::vector<T> const &primes, [[maybe_unused]] std::ostream &os) {
  // We need a1 x1^2 + a2 x2^2 - t u^2 = 0.
  //
  T a1 = a(0);
  T a2 = a(1);
  T a3 = a(2);
  T a4 = a(3);
  std::vector<T> vect12{-a1, -a2};
  std::vector<T> vect34{a3, a4};
  std::vector<int> set12 = get_poss_signs(vect12);
  std::vector<int> set34 = get_poss_signs(vect34);
  std::vector<int> set = IntersectionVect(set12, set34);
#ifdef DEBUG_QUATERNARY
  if (set.size() == 0) {
    std::cerr << "set should not be empty\n";
    throw TerminalException{1};
  }
#endif
  int sign = set[0];
  T q(sign);
  for (size_t u = 0; u < primes.size(); u++) {
    T tp = classes[u];
    T p = primes[u];
    std::pair<size_t, T> pair = Padic_decompose(tp, p);
    T expo = MyPow(p, pair.first);
    q *= expo;
  }
  return q;
}

template <typename T>
MyVector<T> dim4_pair_legendre_iterate_solution(MyVector<T> const &a,
                                                std::vector<T> const &primes,
                                                T const &q, std::ostream &os) {
  size_t len = primes.size();
  T a1 = a(0);
  T a2 = a(1);
  T a3 = a(2);
  T a4 = a(3);
  T ma3 = -a3;
  T ma4 = -a4;
  T prod_12 = -a1 * a2;
  T prod_34 = -a3 * a4;
  std::vector<int> l_symb12, l_symb34;
  for (auto &p : primes) {
    int symb12 = hilbert_symbol(a1, a2, p, os);
    int symb34 = hilbert_symbol(ma3, ma4, p, os);
    l_symb12.push_back(symb12);
    l_symb34.push_back(symb34);
  }
  auto padic_test = [&](T const &t) -> bool {
    for (size_t u = 0; u < len; u++) {
      T const &p = primes[u];
      int symb12 = hilbert_symbol(t, prod_12, p, os);
      int symb34 = hilbert_symbol(t, prod_34, p, os);
      if (symb12 != l_symb12[u] || symb34 != l_symb34[u]) {
        return false;
      }
    }
    return true;
  };
  auto get_solution = [&](T const &t) -> MyVector<T> {
#ifdef DEBUG_QUATERNARY
    os << "QUAD: t=" << t << "\n";
#endif
    MyVector<T> a12(3);
    a12(0) = a1;
    a12(1) = a2;
    a12(2) = -t;
#ifdef DEBUG_QUATERNARY
    os << "QUAD: a12=" << StringVectorGAP(a12) << "\n";
#endif
    std::optional<MyVector<T>> opt12 = TernaryIsotropicVectorDiagonal(a12, os);
    MyVector<T> V12 = unfold_opt(opt12, "opt12 should be solvable");
#ifdef DEBUG_QUATERNARY
    os << "QUAD: V12=" << StringVectorGAP(V12) << "\n";
    T sum12 = a1 * V12(0) * V12(0) + a2 * V12(1) * V12(1) - t * V12(2) * V12(2);
    if (sum12 != 0) {
      std::cerr << "V12 should be an isotropic vector\n";
      throw TerminalException{1};
    }
#endif
    // Special case allows early termination of the computation.
    // If we did not terminate early, we could get a zero vector
    if (V12(2) == 0) {
      MyVector<T> Vret = ZeroVector<T>(4);
      Vret(0) = V12(0);
      Vret(1) = V12(1);
      return Vret;
    }
    MyVector<T> a34(3);
    a34(0) = a3;
    a34(1) = a4;
    a34(2) = t;
#ifdef DEBUG_QUATERNARY
    os << "QUAD: a34=" << StringVectorGAP(a34) << "\n";
#endif
    std::optional<MyVector<T>> opt34 = TernaryIsotropicVectorDiagonal(a34, os);
    MyVector<T> V34 = unfold_opt(opt34, "opt34 should be solvable");
#ifdef DEBUG_QUATERNARY
    os << "QUAD: V34=" << StringVectorGAP(V34) << "\n";
    T sum34 = a3 * V34(0) * V34(0) + a4 * V34(1) * V34(1) + t * V34(2) * V34(2);
    if (sum34 != 0) {
      std::cerr << "V34 should be an isotropic vector\n";
      throw TerminalException{1};
    }
#endif
    //
    T val12 = V12(2);
    T val34 = V34(2);
    MyVector<T> Vret(4);
    Vret(0) = V12(0) * val34;
    Vret(1) = V12(1) * val34;
    Vret(2) = V34(0) * val12;
    Vret(3) = V34(1) * val12;
#ifdef DEBUG_QUATERNARY
    T sum = 0;
    for (int u = 0; u < 4; u++) {
      sum += a(u) * Vret(u) * Vret(u);
    }
    if (sum != 0) {
      std::cerr << "QUAD: Vret should be an isotropic vector\n";
      throw TerminalException{1};
    }
#endif
    return Vret;
  };
  T p0(2);
  while (true) {
    if (IsPrime(p0)) {
      T t = p0 * q;
      if (padic_test(t)) {
        return get_solution(t);
      }
    }
    p0 += 1;
  }
}

template <typename T>
MyVector<T> dim4_compute_solution(MyVector<T> const &a,
                                  std::vector<T> const &primes,
                                  std::ostream &os) {
  std::vector<T> classes = get_tp_classes(a, primes, os);
  T q = get_q_val(a, classes, primes, os);
#ifdef DEBUG_QUATERNARY
  os << "QUAD: q=" << q << "\n";
#endif
  MyVector<T> V = dim4_pair_legendre_iterate_solution(a, primes, q, os);
#ifdef DEBUG_QUATERNARY
  if (IsZeroVector(V)) {
    std::cerr << "QUAD: We found a zero vector\n";
    throw TerminalException{1};
  }
#endif
  return V;
}

template <typename T>
bool quaternary_has_isotropic_vector(MyMatrix<T> const &M, std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
#ifdef DEBUG_QUATERNARY
  os << "QUAD: quaternary_has_isotropic_vector, beginning\n";
#endif
  MyVector<T> red_diag = get_reduced_diagonal(M, os).second;
  MyVector<Tring> red_diagB = UniversalVectorConversion<Tring, T>(red_diag);
#ifdef DEBUG_QUATERNARY
  os << "QUAD: quaternary_has_isotropic_vector, we have red_diag\n";
#endif
  std::vector<Tring> primes = get_local_primes(red_diagB, os);
  return determine_solvability_dim4(red_diagB, primes, os);
}

template <typename T>
std::optional<MyVector<T>>
QuaternaryIsotropicVectorDiagonal(MyVector<T> const &a, std::ostream &os) {
#ifdef DEBUG_QUATERNARY
  os << "QUAD: QuaternaryIsotropicVector, we have diag\n";
#endif
  std::vector<T> primes = get_local_primes(a, os);
#ifdef DEBUG_QUATERNARY
  os << "QUAD: primes =";
  for (auto &p : primes) {
    os << " " << p;
  }
  os << "\n";
#endif
  bool test = determine_solvability_dim4(a, primes, os);
  if (!test) {
    return {};
  }
  return dim4_compute_solution(a, primes, os);
}

template <typename T>
std::optional<MyVector<T>> QuaternaryIsotropicVector(MyMatrix<T> const &M,
                                                     std::ostream &os) {
#ifdef TIMINGS_QUATERNARY
  MicrosecondTime time;
#endif
  using Tring = typename underlying_ring<T>::ring_type;
#ifdef DEBUG_QUATERNARY
  os << "QUAD: QuaternaryIsotropicVector, beginning\n";
#endif
  std::pair<MyMatrix<T>, MyVector<T>> pair1 = get_reduced_diagonal(M, os);
  MyVector<Tring> a = UniversalVectorConversion<Tring, T>(pair1.second);
  std::optional<MyVector<Tring>> opt = QuaternaryIsotropicVectorDiagonal(a, os);
  if (!opt) {
#ifdef TIMINGS_QUATERNARY
    os << "|QUAD: QuaternaryIsotropicVector(None)|=" << time << "\n";
#endif
    return {};
  }
  MyVector<Tring> const &sol1 = *opt;
  MyVector<T> sol2 = UniversalVectorConversion<T, Tring>(sol1);
  MyVector<T> sol3 = pair1.first.transpose() * sol2;
#ifdef DEBUG_QUATERNARY
  os << "QUAD: sol3=" << StringVector(sol3) << "\n";
  T sum3 = EvaluationQuadForm(M, sol3);
  if (sum3 != 0) {
    std::cerr << "QUAD: sol3 is not a solution of the equation\n";
    throw TerminalException{1};
  }
  if (IsZeroVector(sol3)) {
    std::cerr << "QUAD: sol3 should be a non-zero vector\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_QUATERNARY
  os << "|QUAD: QuaternaryIsotropicVector(Some)|=" << time << "\n";
#endif
  return sol3;
}

// clang-format off
#endif  // SRC_INDEFINITE_QUATERNARY_H_
// clang-format on
