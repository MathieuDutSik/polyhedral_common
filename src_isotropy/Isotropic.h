// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_ISOTROPIC_H_
#define SRC_INDEFINITE_ISOTROPIC_H_

// clang-format off
#include "DeterminantMinimization.h"
#include "Legendre_equation.h"
#include "Quaternary.h"
#include "FifthAndHigherOrder.h"
#include "Indefinite_LLL.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ISOTROPIC
#endif

#ifdef TIMINGS
#define TIMINGS_ISOTROPIC
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
    It is known not to be proved correct algorithm.
    * Function is is_isotropic_with_vector
    * It depends on _isotropic_subspace
    * Which uses the Representative(Genus(R,s)) to find a representative
    * The uses the maximal_even_lattice
    * Then uses the maximal_isotropic_subspace_unimodular (which rely on heuristics
  2) Pari: src/basemath/qfsolve.c
    It is a code that is known to have issues as bugs were found.
  3) The algorithm by "Algorithms for solving rational quadratic forms" Josef Schicho and
    Jana Pilnikova (P3) seems relatively more sensible. The induction is dimension is a
    very tall thing. But we could reduce to a 5 dim case by finding the right subspace.
  ---
  The Indefinite LLL with a few tricks are indeed working mostly fine
  except when they do not.
  So, we could have
  * A very sophisticated algorithm by P3 that decides existence of 
  * An Indefinite LLL iteration working.
  Both algorithms could be run in parallel and an early termination of either could
  finish the enumeration.
 */

template<typename T>
MyMatrix<T> GetAnMatrix(int const& k, int const& n) {
  MyMatrix<T> M = IdentityMat<T>(n);
  for (int i=0; i<k; i++) {
    M(i,i) = 2;
  }
  for (int i=0; i<k-1; i++) {
    M(i, i+1) = 1;
    M(i+1, i) = 1;
  }
  return M;
}

template<typename T>
T random_T(T const& q) {
  size_t val1 = random();
  T val2 = UniversalScalarConversion<T,size_t>(val1);
  return ResInt(val2, q);
}

template<typename T>
MyMatrix<T> GetSmithEntry(std::map<T, size_t> const& map, int const& n) {
  MyMatrix<T> M = IdentityMat<T>(n);
  for (auto& kv: map) {
    int pos = random() % n;
    M(pos,pos) *= kv.first;
  }
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      T val = T_min(M(i,i), M(j,j));
      M(i,j) = random_T(val);
    }
  }
  return M;
}

template<typename T>
std::optional<MyVector<T>> get_isotropic_easy_method(MyMatrix<T> const& Q, [[maybe_unused]] std::ostream& os) {
  int n = Q.rows();
  // Checking first for diagonal zeros.
  for (int i = 0; i < n; i++) {
    if (Q(i, i) == 0) {
      MyVector<T> eV = ZeroVector<T>(n);
      eV(i) = 1;
#ifdef DEBUG_ISOTROPIC
      os << "ISOTROP: get_isotropic_easy_method finding isotrop in the diagonal\n";
#endif
      return eV;
    }
  }
  // Checking for opposite signs in the diagonal
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (Q(i, j) == 0 && Q(i, i) + Q(j, j) == 0) {
        MyVector<T> eV = ZeroVector<T>(n);
        eV(i) = 1;
        eV(j) = 1;
#ifdef DEBUG_ISOTROPIC
        os << "ISOTROP: get_isotropic_easy_method finding isotrop as sum of two diagonal terms\n";
#endif
        return eV;
      }
    }
  }
  return {};
}



// Try to find isotropic subspace by using Indefinite LLL
// We iterate by taking random integral matrices and we iterate until
// we cannot reduce anymore the norm of the matrix.
//
// The isotropic vectors can be found using 3 methods:
// --- The Indefinite LLL that found an entry.
// --- 0 entry in the diagonal of the reduction.
// --- Entries (+a, -a) in the diagonal and Q(i,j) = 0 that allow to find an isotropic vector.
template <typename T>
std::optional<MyVector<T>> GetIsotropIndefiniteLLL(MyMatrix<T> const &Q, std::ostream& os) {
  using Tint = typename underlying_ring<T>::ring_type;
  int n = Q.rows();
  auto get_norm = [&](MyMatrix<T> const &mat) -> T {
    T sum = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(mat(i, j));
    return sum;
  };
#ifdef DEBUG_ISOTROPIC
  T det = DeterminantMat(Q);
  os << "ISOTROP: GetIsotropIndefiniteLLL det=" << det << "\n";
  std::map<T, size_t> map = FactorsIntMap(T_abs(det));
  for (auto & kv : map) {
    os << "kv.first=" << kv.first << " kv.second=" << kv.second << "\n";
  }
#endif
  // Direct try at the beginning.
  std::optional<MyVector<T>> optA = get_isotropic_easy_method(Q, os);
  if (optA) {
    MyVector<T> const& eV = *optA;
    return eV;
  }
  MyMatrix<T> Pw = IdentityMat<T>(n);
  MyMatrix<T> Qw = Q;
  T curr_norm = get_norm(Q);
#ifdef DEBUG_ISOTROPIC
  size_t iter = 0;
#endif
  while(true) {
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: GetIsotropIndefiniteLLL while loop, iter=" << iter << "\n";
    iter += 1;
#endif
    // Compute the LLL reduction.
    ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(Qw);
    if (res.Xisotrop) {
      MyVector<T> const& Xisotrop = *res.Xisotrop;
      MyVector<T> V = Pw.transpose() * Xisotrop;
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
    std::optional<MyVector<T>> optB = get_isotropic_easy_method(Qw, os);
    if (optB) {
      MyVector<T> const& eV = *optB;
      MyVector<T> fV = Pw.transpose() * eV;
      return fV;
    }
    T norm = get_norm(Qw);
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: GetIsotropIndefiniteLLL norm=" << norm << " curr_norm=" << curr_norm << "\n";
#endif
    if (norm >= curr_norm) {
      break;
    }
    curr_norm = norm;
    MyMatrix<T> RandM = get_random_int_matrix<T>(n);
    Pw = RandM * Pw;
    Qw = RandM * Qw * RandM.transpose();
#ifdef DEBUG_ISOTROPIC
    os << "ISOTROP: GetIsotropIndefiniteLLL det(Pw)=" << DeterminantMat(Pw) << " det(Qw)=" << DeterminantMat(Qw) << "\n";
#endif
  }
  return {};
}



// For rank 1 this is trivial
template <typename T>
std::optional<MyVector<T>> FindIsotropicRankOne(MyMatrix<T> const &M) {
  if (M(0, 0) != 0)
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
  MyVector<T> V = ZeroVector<T>(2);
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
#ifdef DEBUG_ISOTROPIC
    if (IsZeroVector(V)) {
      std::cerr << "ISOTROP: V should be non-zero\n";
      throw TerminalException{1};
    }
#endif
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
std::optional<MyVector<T>> FindIsotropic_LLL_nfixed(MyMatrix<T> const &Q, std::ostream& os) {
  int n = Q.rows();
  MyMatrix<T> Pw = IdentityMat<T>(n);
  MyMatrix<T> Qw = Q;
  // Iterating 2 times seems right from simulation. Little benefit from going farther.
  int n_iter = 2;
  for (int iter=0; iter<n_iter; iter++) {
    std::optional<MyVector<T>> opt = GetIsotropIndefiniteLLL(Qw, os);
    if (opt) {
#ifdef DEBUG_ISOTROPIC
      os << "ISOTROP: RESULT rec(iter:=" << iter << "),\n";
#endif
      MyVector<T> const &eV = *opt;
      if (!IsZeroVector(eV)) {
        MyVector<T> fV = Pw.transpose() * eV;
        return fV;
      }
    }
    MyMatrix<T> U = get_random_int_matrix<T>(n);
    Pw = U * Pw;
    Qw = U * Qw * U.transpose();
  }
  return {};
}



template <typename T>
std::optional<MyVector<T>> FindIsotropicExact(MyMatrix<T> const &M, std::ostream& os) {
  int n = M.rows();
#ifdef DEBUG_ISOTROPIC
  os << "ISOTROP: STARTEXACT M=" << StringMatrixGAP(M) << "\n";
  int rnk = RankMat(M);
  if (rnk != n) {
    std::cerr << "ISOTROP: FindIsotropicExact failing with rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "ISOTROP: The matrix should be non-degenerate\n";
    throw TerminalException{1};
  }
#endif
  if (n == 1) {
    return FindIsotropicRankOne(M);
  }
  if (n == 2) {
    return FindIsotropicRankTwo(M);
  }
  if (n == 3) {
    return TernaryIsotropicVector(M, os);
  }
  if (n == 4) {
    return QuaternaryIsotropicVector(M, os);
  }
  return FifthAndHigherOrderIsotropicVector(M, os);
}

template <typename T>
std::optional<MyVector<T>> FindIsotropic(MyMatrix<T> const &M, std::ostream& os) {
#ifdef TIMINGS_ISOTROPIC
  MicrosecondTime time;
#endif
  using Tint = typename underlying_ring<T>::ring_type;
  int n = M.rows();
#ifdef DEBUG_ISOTROPIC
  os << "ISOTROP: START M=" << StringMatrixGAP(M) << "\n";
  int rnk = RankMat(M);
  if (rnk != n) {
    std::cerr << "ISOTROP: FindIsotropic failing with rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "ISOTROP: The matrix should be non-degenerate\n";
    throw TerminalException{1};
  }
#endif
  // For low dimension, no point using LLL stuff.
  if (n <= 2) {
    return FindIsotropicExact(M, os);
  }
  // Trying to solve by using the iterated Indefinite-LLL.
  std::optional<MyVector<T>> optA = FindIsotropic_LLL_nfixed(M, os);
#ifdef TIMINGS_LEGENDRE
  os << "|FindIsotropic(optA)|=" << time << "\n";
#endif
  if (optA) {
    return *optA;
  }
  // Computing the Indefinite-LLL reduction. Could get you an isotrop vector
  ResultIndefiniteLLL<T, Tint> res = ComputeReductionIndefinite<T,Tint>(M, os);
#ifdef TIMINGS_LEGENDRE
  os << "|FindIsotropic(ComputeReductionIndefinite)|=" << time << "\n";
#endif
  if (res.Xisotrop) {
    MyVector<T> const& eV = *res.Xisotrop;
    return eV;
  }
  // Now calling the exact algorithm on the reduced matrix
  std::optional<MyVector<T>> optB = FindIsotropicExact(res.Mred, os);
#ifdef TIMINGS_LEGENDRE
  os << "|FindIsotropic(FindIsotropicExact)|=" << time << "\n";
#endif
  if (optB) {
    MyVector<T> const& eV = *optB;
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(res.B);
    MyVector<T> fV = B_T.transpose() * eV;
    return fV;
  } else {
    return {};
  }
}

template <typename T>
bool is_isotropic(MyMatrix<T> const &M, std::ostream& os) {
#ifdef TIMINGS_ISOTROPIC
  MicrosecondTime time;
#endif
  int n = M.rows();
#ifdef DEBUG_ISOTROPIC
  int rnk = RankMat(M);
  if (rnk != n) {
    std::cerr << "ISOTROP: is_isotropic failing with rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "ISOTROP: The matrix should be non-degenerate\n";
    throw TerminalException{1};
  }
#endif
  auto get_val=[&]() -> bool {
    if (n == 1) {
      std::optional<MyVector<T>> opt = FindIsotropicRankOne(M);
      return opt.has_value();
    }
    if (n == 2) {
      std::optional<MyVector<T>> opt = FindIsotropicRankTwo(M);
      return opt.has_value();
    }
    if (n == 3) {
      return ternary_has_isotropic_vector(M, os);
    }
    if (n == 4) {
      return quaternary_has_isotropic_vector(M, os);
    }
    // By Meyer theorem, for n > 4 there is an isotropic vector
    return true;
  };
  bool test = get_val();
#ifdef TIMINGS_ISOTROPIC
  os << "|is_isotropic|=" << time << "\n";
#endif
  return test;
}


// clang-format off
#endif  //  SRC_INDEFINITE_ISOTROPIC_H_
// clang-format on
