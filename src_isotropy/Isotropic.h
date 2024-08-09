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
  * A very sophisticated algorithm by P3 that
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
  T det = DeterminantMat(Q);
  std::map<T, size_t> map = FactorsIntMap(T_abs(det));
#ifdef DEBUG_ISOTROPIC
  os << "ISOTROP: GetIsotropIndefiniteLLL det=" << det << "\n";
  for (auto & kv : map) {
    os << "kv.first=" << kv.first << " kv.second=" << kv.second << "\n";
  }
#endif
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
    auto get_rand_invmat=[&]() -> MyMatrix<T> {
      int choice = random() % 5;
#ifdef DEBUG_ISOTROPIC
      os << "ISOTROP: GetIsotropIndefiniteLLL choice=" << choice << "\n";
#endif
      if (choice == 0) {
        int k = random() % n;
        return GetAnMatrix<T>(k, n);
        //        return GetSmithEntry<T>(map, n);
      } else {
        return get_random_int_matrix<T>(n);
      }
    };
    MyMatrix<T> RandM = get_rand_invmat();
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
  int n_iter = 5;
  for (int iter=0; iter<n_iter; iter++) {
    std::optional<MyVector<T>> opt = GetIsotropIndefiniteLLL(Qw, os);
    if (opt) {
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
  int n = M.rows();
#ifdef DEBUG_ISOTROPIC
  int rnk = RankMat(M);
  if (rnk != n) {
    std::cerr << "ISOTROP: FindIsotropic failing with rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "ISOTROP: The matrix should be non-degenerate\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 2) {
    return FindIsotropicExact(M, os);
  }
  std::optional<MyVector<T>> opt = FindIsotropic_LLL_nfixed(M, os);
  if (opt) {
    return *opt;
  }
  return FindIsotropicExact(M, os);
}

template <typename T>
bool is_isotropic(MyMatrix<T> const &M, std::ostream& os) {
  int n = M.rows();
#ifdef DEBUG_ISOTROPIC
  int rnk = RankMat(M);
  if (rnk != n) {
    std::cerr << "ISOTROP: is_isotropic failing with rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "ISOTROP: The matrix should be non-degenerate\n";
    throw TerminalException{1};
  }
#endif
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
  return true;
}


// clang-format off
#endif  //  SRC_INDEFINITE_ISOTROPIC_H_
// clang-format on
