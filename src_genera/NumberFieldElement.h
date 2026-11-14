// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_NUMBERFIELDELEMENT_H_
#define SRC_GENERA_NUMBERFIELDELEMENT_H_

// clang-format off
#include "PolynomialModQ.h"
#include "GeneraArithmetic.h"
#include "MAT_Matrix.h"
#include <cstdint>
#include <functional>
#include <vector>
// clang-format on

/*
  Elements of a number field E = Q(theta) on the basis of the powers of
  theta (Q(zeta_d) on the powers of zeta, or a field given by the minimal
  polynomial of theta and the complex conjugate of theta): products,
  inverses, conjugation, valuations and residues at a prime.
 */

template <typename T> using ElemE = std::vector<T>;

template <typename T> struct CyclotomicField {
  // d > 0: Q(zeta_d) on the powers of zeta; d = 0: a field Q(theta) on the
  // powers of theta, with Phi the minimal polynomial of theta and conj the
  // complex conjugate of theta on the same basis
  int64_t d;
  int deg;
  std::vector<int64_t> Phi;
  std::vector<T> conj;
  // v_P = v_q(N(a)) / val_div at a prime P alone over q, of residue
  // degree val_div
  int val_div = 1;
  // when set, v_P (any prime P, see LocalPrimeK<T>)
  std::function<int(std::vector<T> const &)> valuation;
};

template <typename T>
inline ElemE<T> EReduce(std::vector<T> a, CyclotomicField<T> const &K) {
  int deg = K.deg;
  for (int k = static_cast<int>(a.size()) - 1; k >= deg; k--) {
    T c = a[k];
    if (c != 0) {
      for (int j = 0; j <= deg; j++) {
        a[k - deg + j] -= c * UniversalScalarConversion<T, int64_t>(K.Phi[j]);
      }
    }
  }
  a.resize(deg, T(0));
  return a;
}

template <typename T>
inline ElemE<T> EMul(ElemE<T> const &a, ElemE<T> const &b, CyclotomicField<T> const &K) {
  std::vector<T> c(a.size() + b.size(), T(0));
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] != 0) {
      for (size_t j = 0; j < b.size(); j++) {
        c[i + j] += a[i] * b[j];
      }
    }
  }
  return EReduce(c, K);
}

template <typename T>
inline ElemE<T> ESub(ElemE<T> a, ElemE<T> const &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] -= b[i];
  }
  return a;
}

template <typename T>
inline ElemE<T> EZetaPower(int64_t t, CyclotomicField<T> const &K) {
  if (K.d == 0) {
    // theta^t, t >= 0
    if (t < 0) {
      std::cerr << "CYCLO: a negative power of theta\n";
      throw TerminalException{1};
    }
    std::vector<T> c(std::max<int64_t>(t + 1, K.deg), T(0));
    c[t] = 1;
    return EReduce(c, K);
  }
  std::vector<T> c(K.d, T(0));
  c[((t % K.d) + K.d) % K.d] = 1;
  return EReduce(c, K);
}

template <typename T>
inline ElemE<T> EConj(ElemE<T> const &a, CyclotomicField<T> const &K) {
  if (K.d == 0) {
    // a(conj(theta)) by Horner
    ElemE<T> res(K.deg, T(0));
    for (int i = K.deg - 1; i >= 0; i--) {
      res = EMul(res, K.conj, K);
      res[0] += a[i];
    }
    return res;
  }
  std::vector<T> c(K.d, T(0));
  for (int t = 0; t < K.deg; t++) {
    c[(K.d - t) % K.d] += a[t];
  }
  return EReduce(c, K);
}

template <typename T>
inline MyMatrix<T> EMulMatrix(ElemE<T> const &a, CyclotomicField<T> const &K) {
  MyMatrix<T> M(K.deg, K.deg);
  for (int j = 0; j < K.deg; j++) {
    ElemE<T> c = EMul(a, EZetaPower(j, K), K);
    for (int i = 0; i < K.deg; i++) {
      M(j, i) = c[i];
    }
  }
  return M;
}

template <typename T>
inline bool EIsZero(ElemE<T> const &a) {
  for (auto &x : a) {
    if (x != 0) {
      return false;
    }
  }
  return true;
}

template <typename T>
inline ElemE<T> EInv(ElemE<T> const &a, CyclotomicField<T> const &K) {
  // x a = 1: x M_a = e_0 with the rows of M_a the a zeta^j
  MyMatrix<T> M = EMulMatrix(a, K);
  MyMatrix<T> Minv = Inverse(M);
  ElemE<T> x(K.deg);
  for (int i = 0; i < K.deg; i++) {
    x[i] = Minv(0, i);
  }
  return x;
}

// v_q(N(a)), large for zero
template <typename T>
inline int EValuation(ElemE<T> const &a, CyclotomicField<T> const &K, int64_t q) {
  if (EIsZero(a)) {
    return 1 << 28;
  }
  if (K.valuation) {
    return K.valuation(a);
  }
  return RationalValuation(DeterminantMat(EMulMatrix(a, K)), q) / K.val_div;
}

// the residue at zeta = 1 modulo q of a P-integral element
template <typename T>
inline typename underlying_z_ring<T>::ring_type EResidue(ElemE<T> const &a,
                                                        int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  Tint qz = UniversalScalarConversion<Tint, int64_t>(q);
  Tint res(0);
  for (auto &x : a) {
    Tint num = ResInt(GetNumerator_z(x), qz);
    Tint den = ResInt(GetDenominator_z(x), qz);
    Tint prod = num * mod_inv(den, qz);
    res += prod;
    res = ResInt(res, qz);
  }
  return res;
}

template <typename T>
inline ElemE<T> EAddE(ElemE<T> a, ElemE<T> const &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] += b[i];
  }
  return a;
}

template <typename T>
inline ElemE<T> EScal(ElemE<T> a, T const &c) {
  for (auto &x : a) {
    x *= c;
  }
  return a;
}

template <typename T>
inline bool EIsOne(ElemE<T> const &a) {
  ElemE<T> b = a;
  b[0] -= 1;
  return EIsZero(b);
}

template <typename T>
inline ElemE<T> EOne(CyclotomicField<T> const &K) { return EZetaPower(0, K); }
// g(theta) in E = Q(theta), g with integer coefficients.
template <typename T>
inline ElemE<T> EPolyOfTheta(PolyQ const &g, CyclotomicField<T> const &K) {
  ElemE<T> res(K.deg, T(0));
  ElemE<T> pw = EZetaPower(0, K);
  ElemE<T> theta = EZetaPower(1, K);
  for (size_t i = 0; i < g.size(); i++) {
    if (g[i] != 0) {
      T g_i = UniversalScalarConversion<T, int64_t>(g[i]);
      for (int t = 0; t < K.deg; t++) {
        res[t] += g_i * pw[t];
      }
    }
    pw = EMul(pw, theta, K);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_NUMBERFIELDELEMENT_H_
// clang-format on
