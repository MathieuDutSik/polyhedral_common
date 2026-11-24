// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SCALARNUMBERFIELD_H_
#define SRC_GENERA_SCALARNUMBERFIELD_H_

// clang-format off
#include "ScalarPolynomialModQ.h"
#include "GeneraArithmetic.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixSolutionMat.h"
#include <initializer_list>
#include <optional>
#include <cstdint>
#include <vector>
// clang-format on

/*
  The scalar type of the number fields: the elements of E = Q(theta) on
  the basis of the powers of theta (Q(zeta_d) on the powers of zeta, or a
  field given by the minimal polynomial of theta and the complex conjugate
  of theta): products, inverses, conjugation, the norm valuation, the
  power traces. The primes of E, their valuations and residues are in
  PrimeDecomposition.h, the abelian fields as subfields of Q(zeta_n) in
  ScalarAbelianField.h.
 */

/*
  An element of E: its deg coordinates on the basis of the powers of theta
  (zeta). The ring operations needing E (EMul, EInv, EConj, ...) take the
  field; the sum, the difference and the product by a rational are
  operators.
 */
template <typename T> struct ElemE {
  std::vector<T> c;

  ElemE() = default;
  ElemE(size_t n, T const &x) : c(n, x) {}
  explicit ElemE(size_t n) : c(n, T(0)) {}
  ElemE(std::initializer_list<T> l) : c(l) {}
  explicit ElemE(std::vector<T> v) : c(std::move(v)) {}

  size_t size() const { return c.size(); }
  T &operator[](size_t i) { return c[i]; }
  T const &operator[](size_t i) const { return c[i]; }
  auto begin() { return c.begin(); }
  auto end() { return c.end(); }
  auto begin() const { return c.begin(); }
  auto end() const { return c.end(); }
  bool is_zero() const {
    for (auto const &x : c) {
      if (x != 0) {
        return false;
      }
    }
    return true;
  }
  ElemE &operator+=(ElemE const &b) {
    for (size_t i = 0; i < c.size(); i++) {
      c[i] += b.c[i];
    }
    return *this;
  }
  ElemE &operator-=(ElemE const &b) {
    for (size_t i = 0; i < c.size(); i++) {
      c[i] -= b.c[i];
    }
    return *this;
  }
  ElemE &operator*=(T const &x) {
    for (auto &y : c) {
      y *= x;
    }
    return *this;
  }
  bool operator==(ElemE const &b) const { return c == b.c; }
};

template <typename T> ElemE<T> operator+(ElemE<T> a, ElemE<T> const &b) {
  return a += b;
}
template <typename T> ElemE<T> operator-(ElemE<T> a, ElemE<T> const &b) {
  return a -= b;
}
template <typename T> ElemE<T> operator-(ElemE<T> a) {
  for (auto &x : a) {
    x = -x;
  }
  return a;
}
template <typename T> ElemE<T> operator*(ElemE<T> a, T const &x) {
  return a *= x;
}
template <typename T> ElemE<T> operator*(T const &x, ElemE<T> a) {
  return a *= x;
}

template <typename T> struct NumberField {
  // d > 0: Q(zeta_d) on the powers of zeta; d = 0: a field Q(theta) on the
  // powers of theta, with Phi the minimal polynomial of theta and conj the
  // complex conjugate of theta on the same basis
  int64_t d;
  int deg;
  std::vector<int64_t> Phi;
  ElemE<T> conj;
};

template <typename T>
inline ElemE<T> EReduce(std::vector<T> a, NumberField<T> const &K) {
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
  return ElemE<T>(std::move(a));
}

template <typename T>
inline ElemE<T> EMul(ElemE<T> const &a, ElemE<T> const &b, NumberField<T> const &K) {
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
inline ElemE<T> EZetaPower(int64_t t, NumberField<T> const &K) {
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
inline ElemE<T> EConj(ElemE<T> const &a, NumberField<T> const &K) {
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
inline MyMatrix<T> EMulMatrix(ElemE<T> const &a, NumberField<T> const &K) {
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
inline ElemE<T> EInv(ElemE<T> const &a, NumberField<T> const &K) {
  // x a = 1: x M_a = e_0 with the rows of M_a the a zeta^j
  MyVector<T> e0 = ZeroVector<T>(K.deg);
  e0(0) = 1;
  std::optional<MyVector<T>> sol = SolutionMat(EMulMatrix(a, K), e0);
  if (!sol) {
    std::cerr << "NUMBERFIELD: the inverse of zero\n";
    throw TerminalException{1};
  }
  ElemE<T> x(K.deg);
  for (int i = 0; i < K.deg; i++) {
    x[i] = (*sol)(i);
  }
  return x;
}

// v_q(N(a)), large for zero (the valuations at the primes of E are in
// PrimeDecomposition.h)
template <typename T>
inline int ENormValuation(ElemE<T> const &a, NumberField<T> const &K,
                          int64_t q) {
  if (a.is_zero()) {
    return 1 << 28;
  }
  return RationalValuation(DeterminantMat(EMulMatrix(a, K)), q);
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
inline bool EIsOne(ElemE<T> const &a) {
  ElemE<T> b = a;
  b[0] -= 1;
  return b.is_zero();
}

template <typename T>
inline ElemE<T> EOne(NumberField<T> const &K) { return EZetaPower(0, K); }
// g(theta) in E = Q(theta), g with integer coefficients.
template <typename T>
inline ElemE<T> EPolyOfTheta(PolyQ const &g, NumberField<T> const &K) {
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

// Tr(theta^k), k = 0, ..., kmax, theta a root of the monic f (by
// increasing degree), from the coefficients of f by Newton's identities.
template <typename T>
inline std::vector<T> PowerTraces(std::vector<int64_t> const &f, int kmax) {
  int deg = static_cast<int>(f.size()) - 1;
  // f = x^deg + c_{deg-1} x^{deg-1} + ... + c_0
  std::vector<T> p(kmax + 1, T(0));
  p[0] = deg;
  for (int k = 1; k <= kmax; k++) {
    T s(0);
    // Newton: p_k + c_{deg-1} p_{k-1} + ... = -k c_{deg-k} for k <= deg
    for (int i = 1; i <= std::min(k - 1, deg); i++) {
      s += UniversalScalarConversion<T, int64_t>(f[deg - i]) * p[k - i];
    }
    if (k <= deg) {
      s += UniversalScalarConversion<T, int64_t>(k) *
           UniversalScalarConversion<T, int64_t>(f[deg - k]);
    }
    p[k] = -s;
  }
  return p;
}

// The trace form on the powers of theta: the matrix (Tr(theta^{a + b}))
// for 0 <= a, b < deg f, of determinant disc(f).
template <typename T>
inline MyMatrix<T> PowerTraceMatrix(std::vector<int64_t> const &f) {
  int deg = static_cast<int>(f.size()) - 1;
  std::vector<T> pt = PowerTraces<T>(f, 2 * deg);
  MyMatrix<T> M(deg, deg);
  for (int a = 0; a < deg; a++) {
    for (int b = 0; b < deg; b++) {
      M(a, b) = pt[a + b];
    }
  }
  return M;
}

// clang-format off
#endif  // SRC_GENERA_SCALARNUMBERFIELD_H_
// clang-format on
