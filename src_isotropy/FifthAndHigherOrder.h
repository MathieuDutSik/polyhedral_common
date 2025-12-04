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

#ifdef DISABLE_DEBUG_FIFTH_AND_HIGHER_ORDER
#undef DEBUG_FIFTH_AND_HIGHER_ORDER
#endif

#ifdef TIMINGS
#define TIMINGS_FIFTH_AND_HIGHER_ORDER
#endif

/*
  By Meyer theorem, we know that there is one isotropic vector.
 */
template <typename T>
MyVector<T> solve_fifth_equation(MyVector<T> const &a, std::ostream &os) {
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
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: solve_fifth_equation, sign=" << sign << "\n";
#endif
  T t(sign);
  MyVector<T> a_tern(3);
  a_tern(0) = a(0);
  a_tern(1) = a(1);
  MyVector<T> a_quad(4);
  a_quad(0) = a(2);
  a_quad(1) = a(3);
  a_quad(2) = a(4);
  MyVector<T> iso_fifth = ZeroVector<T>(5);
  while (true) {
    a_tern(2) = -t;
    a_quad(3) = t;
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
    os << "FIFTH: solve_fifth_equation, a_tern=" << StringVector(a_tern)
       << "\n";
    os << "FIFTH: solve_fifth_equation, a_quad=" << StringVector(a_quad)
       << "\n";
#endif
    std::optional<MyVector<T>> opt_tern =
        TernaryIsotropicVectorDiagonal(a_tern, os);
    if (opt_tern) {
      MyVector<T> const &iso_tern = *opt_tern;
      // If the last entry is zero there then we can terminate early.
      // Also, continuing could lead us to a zero vector which would be
      // unfortunate
      T val_tern = iso_tern(2);
      if (val_tern == 0) {
        iso_fifth(0) = iso_tern(0);
        iso_fifth(1) = iso_tern(1);
        return iso_fifth;
      }
      std::optional<MyVector<T>> opt_quad =
          QuaternaryIsotropicVectorDiagonal(a_quad, os);
      if (opt_quad) {
        MyVector<T> const &iso_quad = *opt_quad;
        T val_quad = iso_quad(3);
        for (int u = 0; u < 2; u++) {
          iso_fifth(u) = iso_tern(u) * val_quad;
        }
        for (int u = 2; u < 5; u++) {
          iso_fifth(u) = iso_quad(u - 2) * val_tern;
        }
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
        T sum(0);
        for (int u = 0; u < 5; u++) {
          T val = iso_fifth(u);
          sum += a(u) * val * val;
        }
        os << "FIFTH: Before return sum=" << sum << "\n";
        if (sum != 0) {
          std::cerr << "FIFTH: Error, a=" << StringVector(a) << "\n";
          std::cerr << "FIFTH: Error, t=" << t << "\n";
          std::cerr << "FIFTH: Error, iso_fifth=" << StringVector(iso_fifth)
                    << "\n";
          throw TerminalException{1};
        }
#endif
        return iso_fifth;
      }
    }
    t += sign;
  }
}

/*
  Compute a basis of the space
 */
template <typename T>
MyMatrix<T> compute_fifth_basis(MyMatrix<T> const &Q, std::ostream &os) {
  int n = Q.rows();
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: compute_fifth_basis, beginning, n=" << n << "\n";
#endif
  int max_one_sign = 4;
  std::vector<MyVector<T>> basis;
  MyMatrix<T> M;
  auto is_linearly_independent = [&](MyVector<T> const &v) -> bool {
    if (M.rows() == 0) {
      return true;
    }
    std::optional<MyVector<T>> opt = SolutionMat(M, v);
    return !opt.has_value();
  };
  auto insert_vector = [&](MyVector<T> const &v) -> void {
    basis.push_back(v);
    M = MatrixFromVectorFamily(basis);
  };
  auto test_vector = [&](MyVector<T> const &v) -> bool {
    if (basis.size() == 5) {
      return false;
    }
    if (!is_linearly_independent(v)) {
      return false;
    }
    std::vector<MyVector<T>> test_basis = basis;
    test_basis.push_back(v);
    MyMatrix<T> test_M = MatrixFromVectorFamily(test_basis);
    MyMatrix<T> QuadM = test_M * Q * test_M.transpose();
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(QuadM, os);
    if (DiagInfo.nbPlus > max_one_sign) {
      return false;
    }
    if (DiagInfo.nbMinus > max_one_sign) {
      return false;
    }
    return true;
  };
  auto insert_vector_if_ok = [&](MyVector<T> const &v) {
    if (test_vector(v)) {
      insert_vector(v);
    }
  };
  for (int i = 0; i < n; i++) {
    MyVector<T> v = ZeroVector<T>(n);
    v(i) = 1;
    insert_vector_if_ok(v);
  }
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      MyVector<T> v = ZeroVector<T>(n);
      v(i) = 1;
      v(j) = 1;
      insert_vector_if_ok(v);
    }
  }
  auto get_sign_vector = [&](int const &sign) -> MyVector<T> {
    MyMatrix<T> ProdM = M * Q;
    MyMatrix<T> NSP = NullspaceTrMat(ProdM);
    MyMatrix<T> QuadQ = sign * NSP * Q * NSP.transpose();
    T CritNorm(0);
    bool StrictIneq = true;
    MyVector<T> v1 =
        GetIntegralVector_allmeth<T, T>(QuadQ, CritNorm, StrictIneq, os);
    MyVector<T> v2 = NSP.transpose() * v1;
    return v2;
  };
  while (true) {
    if (M.rows() == 5) {
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
      os << "FIFTH: compute_fifth_basis, returning M\n";
#endif
      return M;
    }
    MyMatrix<T> QuadM = M * Q * M.transpose();
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(QuadM, os);
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

template <typename T>
MyVector<T> FifthOrderIsotropicVector(MyMatrix<T> const &M, std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: FifthOrderIsotropicVector, beginning\n";
  if (M.rows() != 5) {
    std::cerr << "FIFTH: M should be a 5x5 matrix\n";
    throw TerminalException{1};
  }
#endif
  std::pair<MyMatrix<T>, MyVector<T>> pair1 = get_reduced_diagonal(M, os);
  MyVector<Tring> a = UniversalVectorConversion<Tring, T>(pair1.second);
  MyVector<Tring> sol1 = solve_fifth_equation(a, os);
  MyVector<T> sol2 = UniversalVectorConversion<T, Tring>(sol1);
  MyVector<T> sol3 = pair1.first.transpose() * sol2;
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: sol3=" << StringVector(sol3) << "\n";
  T sum3 = EvaluationQuadForm(M, sol3);
  if (sum3 != 0) {
    std::cerr << "FIFTH: sol3 is not a solution of the equation\n";
    throw TerminalException{1};
  }
  if (IsZeroVector(sol3)) {
    std::cerr << "FIFTH: sol3 should be non-zero\n";
    throw TerminalException{1};
  }
#endif
  return sol3;
}

template <typename T>
MyVector<T> FifthAndHigherOrderIsotropicVector(MyMatrix<T> const &Q,
                                               std::ostream &os) {
#ifdef TIMINGS_FIFTH_AND_HIGHER_ORDER
  MicrosecondTime time;
#endif
  auto get_basis = [&]() -> MyMatrix<T> {
    if (Q.rows() == 5) {
      return IdentityMat<T>(5);
    } else {
      return compute_fifth_basis(Q, os);
    }
  };
  MyMatrix<T> basis = get_basis();
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: basis=\n";
  WriteMatrix(os, basis);
#endif
  MyMatrix<T> Q2 = basis * Q * basis.transpose();
  if (RankMat(Q2) < 5) {
    MyMatrix<T> NSP = NullspaceMat(Q2);
    MyVector<T> v1 = GetMatrixRow(NSP, 0);
    MyVector<T> v2 = basis.transpose() * v1;
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
    os << "FIFTH: From NSP, v2=" << StringVector(v2) << "\n";
    T sum2 = EvaluationQuadForm(Q, v2);
    if (sum2 != 0) {
      std::cerr << "FIFTH: v2 is not a solution of the equation\n";
      throw TerminalException{1};
    }
    if (IsZeroVector(v2)) {
      std::cerr << "FIFTH: v2 should be non-zero 1\n";
      throw TerminalException{1};
    }
#endif
#ifdef TIMINGS_FIFTH_AND_HIGHER_ORDER
    os << "|FIFTH: FifthAndHigherOrderIsotropicVector(rank)|=" << time << "\n";
#endif
    return v2;
  }
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: Q2=\n";
  WriteMatrix(os, Q2);
#endif
  MyVector<T> v1 = FifthOrderIsotropicVector(Q2, os);
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: v1=" << StringVector(v1) << "\n";
#endif
  MyVector<T> v2 = basis.transpose() * v1;
#ifdef DEBUG_FIFTH_AND_HIGHER_ORDER
  os << "FIFTH: From fifth order v2=" << StringVector(v2) << "\n";
  T sum2 = EvaluationQuadForm(Q, v2);
  if (sum2 != 0) {
    std::cerr << "FIFTH: v2 is not a solution of the equation\n";
    throw TerminalException{1};
  }
  if (IsZeroVector(v2)) {
    std::cerr << "FIFTH: v2 should be non-zero 2\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_FIFTH_AND_HIGHER_ORDER
  os << "|FIFTH: FifthAndHigherOrderIsotropicVector(fifth)|=" << time << "\n";
#endif
  return v2;
}

// clang-format off
#endif  // FIFTHANDHIGHERORDER_H_
// clang-format on
