// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
#define SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_

// clang-format off
#include "factorizations.h"
#include "MAT_MatrixMod.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <map>
// clang-format on

template <typename T> struct ResultDetMin {
  MyMatrix<T> P;
  MyMatrix<T> Mred;
};

/*
  We apply a numbr of ideas from the preprint
  "Quadratic equations in dimensions 4, 5 and more" (P1)
  Denis Simon
  ---
  The following notions are used
  * v_p(Q) to be the multiplicity of p as a prime factor of the determinant of
  Q.
  * d = dim Ker_{F_p}(Q).
  We have the basic result d <= v.
  ---
  The ideas are clear from the manuscript except for Lemma 9, which I do not
  understand yet.
  Things are not done strictly in the same way, but the same ideas
  are implemented.
 */
template <typename T>
ResultDetMin<T> DeterminantMinimization(MyMatrix<T> const &Q) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  using Tring = typename underlying_ring<T>::ring_type;
  if (!IsIntegralMatrix(Q)) {
    std::cerr << "The matrix M should be integral\n";
    throw TerminalException{1};
  }
  int n = Q.rows();
  T det = DeterminantMat(Q);
  T det_abs = T_abs(det);
  Tring det_ai = UniversalScalarConversion<Tring, T>(det_abs);
  std::map<Tring, size_t> map = FactorsIntMap(det_ai);
  MyMatrix<T> Qw = Q;
  MyMatrix<T> Pw = IdentityMat<T>(n);
  while (true) {
    std::vector<Tring> list_P_erase;
    bool DoSomethingGlobal = false;
    for (auto &kv : map) {
      Tring p_ring = kv.first;
      T p = UniversalScalarConversion<T, Tring>(p_ring);
      size_t &v_mult_s = kv.second;
      int v_mult_i = v_mult_s;
      ResultNullspaceMod<T> res = NullspaceMatMod(Qw, p);
      int d_mult_i = res.dimNSP;
      auto change_basis = [&]() -> void {
        Pw = res.BasisTot * Pw;
        Qw = res.BasisTot * Qw * res.BasisTot.transpose();
#ifdef DEBUG_DET_MINIMIZATION
        for (int i = 0; i < d_multi_i; i++) {
          for (int j = 0; j < n; j++) {
            T res = ResInt(Qw(i, j), p);
            if (res != 0) {
              std::cerr << "The coefficient M(i,j) is not reduced modulo p\n";
              std::cerr << "i=" << i << " j=" << j << " p=" << p << "\n";
              throw TerminalException{1};
            }
          }
        }
        int p = n - d_mult_i;
        MyMatrix<T> U(p, p);
        for (int i = 0; i < p; i++) {
          for (int j = 0; j < p; j++) {
            U(i, j) = Qw(i + d_mult_i, j + d_mult_i);
          }
        }
        T det = DeterminantMat(U);
        T res_det = ResInt(det, p);
        if (res_det == 0) {
          std::cerr << "The matrix U is not invertible\n";
          throw TerminalException{1};
        }
#endif
      };
      auto get_qtilde = [&]() -> MyMatrix<T> {
        MyMatrix<T> Qtilde(d_mult_i, d_mult_i);
        for (int i = 0; i < d_mult_i; i++) {
          for (int j = 0; j < d_mult_i; j++) {
            Qtilde(i, j) = Qw(i, j) / p;
          }
        }
        return Qtilde;
      };
      bool DoSomething = false;
      // Apply Lemma 4
      if (d_mult_i == n && !DoSomething) {
        Qw = Qw / p;
        DoSomething = true;
        DoSomethingGlobal = true;
        v_mult_i -= n;
        v_mult_s -= n;
      }
      // Apply Lemma 5
      if (d_mult_i < v_mult_i && !DoSomething) {
        change_basis();
        MyMatrix<T> Qtilde = get_qtilde();
        ResultNullspaceMod<T> resB = NullspaceMatMod(Qtilde, p);
        MyMatrix<T> Hmat = IdentityMat<T>(n);
        for (int i = 0; i < d_mult_i; i++) {
          for (int j = 0; j < d_mult_i; j++) {
            Hmat(i, j) = resB.BasisTot(i, j);
          }
        }
        Pw = Hmat * Pw;
        Qw = Hmat * Qw * Hmat.transpose();
        int dimNSPB = resB.dimNSP;
#ifdef DEBUG_DET_MINIMIZATION
        for (int i = 0; dimNSPB; i++) {
          for (int j = 0; j < dimNSPB; j++) {
            T res = ResInt(Qw(i, j), p * p);
            if (res != 0) {
              std::cerr << "Qw is not divisible by p * p as expected\n";
              throw TerminalException{1};
            }
          }
        }
#endif
        MyMatrix<T> U = IdentityMat<T>(n);
        for (int i = 0; i < dimNSPB; i++)
          U(i, i) = 1 / p;
        Pw = U * Pw;
        Qw = U * Qw * U.transpose();
#ifdef DEBUG_DET_MINIMIZATION
        if (!IsIntegralMatrix(Qw)) {
          std::cerr << "The matrix Qw is not integral\n";
          throw TerminalException{1};
        }
#endif
        DoSomething = true;
        DoSomethingGlobal = true;
        v_mult_i -= dimNSPB;
        v_mult_s -= dimNSPB;
      }
      // Apply Lemma 6
      if (d_mult_i == v_mult_i && n > 2 * d_mult_i && !DoSomething) {
        change_basis();
        MyMatrix<T> U = IdentityMat<T>(n);
        for (int u = d_mult_i; u < n; u++) {
          U(u, u) = p;
        }
        Pw = U * Pw;
        Qw = U * Qw * U.transpose() / p;
        DoSomething = true;
        DoSomethingGlobal = true;
        int dec = n - 2 * d_mult_i;
        v_mult_i -= dec;
        v_mult_s -= dec;
      }
      // Apply Lemma 7 (joined with Lemma 8)
      if (d_mult_i == v_mult_i && !DoSomething) {
        change_basis();
        MyMatrix<T> Qtilde = get_qtilde();
        std::optional<MyVector<T>> opt = FindIsotropicVectorMod(Qtilde, p);
        if (opt) {
          MyVector<T> const &eV = *opt;
          MyMatrix<T> M(1, n);
          for (int i = 0; i < n; i++)
            M(0, i) = eV(i);
          MyMatrix<T> BasisCompl = SubspaceCompletionInt(M, n);
          MyMatrix<T> U = Concatenate(M, BasisCompl);
          Pw = U * Pw;
          Qw = U * Qw * U.transpose();
#ifdef DEBUG_DET_MINIMIZATION
          T res = ResInt(Qw(0, 0), p * p);
          if (res != 0) {
            std::cerr << "We do not have tildeQ(0,0) divisible by p^2\n";
            throw TerminalException{1};
          }
#endif
          MyMatrix<T> V = IdentityMat<T>(n);
          V(0, 0) = 1 / p;
          Pw = V * Pw;
          Qw = V * Qw * V.transpose();
          DoSomething = true;
          DoSomethingGlobal = true;
          v_mult_i -= 2;
          v_mult_s -= 2;
        }
      }
      if (v_mult_s == 0) {
        list_P_erase.push_back(p_ring);
      }
    }
    for (auto &p : list_P_erase) {
      map.erase(p);
    }
    if (!DoSomethingGlobal) {
      break;
    }
  }
  return {Pw, Qw};
}

// clang-format off
#endif  //  SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
// clang-format on
