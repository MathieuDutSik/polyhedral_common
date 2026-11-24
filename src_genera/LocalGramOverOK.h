// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALGRAMOVEROK_H_
#define SRC_GENERA_LOCALGRAMOVEROK_H_

// clang-format off
#include "PrimeDecomposition.h"
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_LOCAL_GRAM_OVER_OK
#endif

/*
  A lattice X with an action of O_E written at q as a lattice over O_E (x)
  Z_q: the Gram matrix over E of the form whose trace is the form of X, on
  vectors spanning X (x) Z_q over O_E (x) Z_q (LocalOKGram), for the local
  factors of the masses (quadratic over a totally real field, Hermitian
  over a CM field).
 */

/*
  The Gram matrix over K of phi on the vectors v_i of X: phi(v_i, v_j) on
  the powers of theta, from the traces Tr(theta^t phi(v_i, v_j)) = b(theta^t
  v_i, v_j), theta being any generator of K.
 */
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
LocalOKGramFromVectors(MyMatrix<T> const &B, LocalGenerator<T> const &lg,
                       std::vector<MyVector<T>> const &chosen) {
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = static_cast<int>(chosen.size());
  std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  MyMatrix<T> Tm = PowerTraceMatrix<T>(lg.f);
  MyMatrix<T> Tinv = Inverse(Tm);
  std::vector<std::vector<ElemE<T>>> H(m, std::vector<ElemE<T>>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE<T> val(deg, T(0));
      for (int t = 0; t < deg; t++) {
        MyVector<T> w = Spow[t].transpose() * chosen[i];
        T bt = w.dot(B * chosen[j]);
        for (int s2 = 0; s2 < deg; s2++) {
          val[s2] += bt * Tinv(t, s2);
        }
      }
      H[i][j] = val;
    }
  }
  return H;
}

/*
  The Gram matrix over K of phi on m vectors of X spanning X (x) Z_q over
  O_K (x) Z_q: standard vectors first, then small combinations.
 */
template <typename T> struct LocalOKBasisGram {
  std::vector<MyVector<T>> basis;
  std::vector<std::vector<ElemE<T>>> H;
};

template <typename T>
inline LocalOKBasisGram<T>
LocalOKGramBasis(MyMatrix<T> const &B, LocalGenerator<T> const &lg,
                 int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = r / deg;
  T_small_prime q_u = static_cast<T_small_prime>(q);
  auto red = [&](T const &x) -> T_residue {
    return ResModOfRational<T, Tint>(x, q_u);
  };
  std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  auto orbit_rows = [&](MyVector<T> const &v) -> ModMatrix {
    ModMatrix rows(deg, r);
    for (int t = 0; t < deg; t++) {
      MyVector<T> w = Spow[t].transpose() * v;
      for (int j = 0; j < r; j++) {
        rows(t, j) = red(w(j));
      }
    }
    return rows;
  };
  std::vector<MyVector<T>> chosen;
  ModMatrix acc(0, r);
  auto try_add = [&](MyVector<T> const &v) -> void {
    ModMatrix rows = orbit_rows(v);
    ModMatrix test(acc.rows() + rows.rows(), r);
    test.topRows(acc.rows()) = acc;
    test.bottomRows(rows.rows()) = rows;
    if (RankModPrime(test, q_u) == static_cast<size_t>(test.rows())) {
      chosen.push_back(v);
      acc = test;
    }
  };
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    MyVector<T> v = ZeroVector<T>(r);
    v(k) = 1;
    try_add(v);
  }
  uint64_t seed = 1;
  for (int it = 0; it < 1000 && static_cast<int>(chosen.size()) < m; it++) {
    MyVector<T> v = ZeroVector<T>(r);
    for (int j = 0; j < r; j++) {
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      v(j) = static_cast<long>((seed >> 33) % 5) - 2;
    }
    try_add(v);
  }
#ifdef SANITY_CHECK_LOCAL_GRAM_OVER_OK
  if (static_cast<int>(chosen.size()) != m) {
    std::cerr << "ABQUAD: no local O_K-basis found\n";
    throw TerminalException{1};
  }
#endif
  return {chosen, LocalOKGramFromVectors(B, lg, chosen)};
}

template <typename T>
inline std::vector<std::vector<ElemE<T>>>
LocalOKGram(MyMatrix<T> const &B, LocalGenerator<T> const &lg, int64_t q) {
  return LocalOKGramBasis(B, lg, q).H;
}

// clang-format off
#endif  // SRC_GENERA_LOCALGRAMOVEROK_H_
// clang-format on
