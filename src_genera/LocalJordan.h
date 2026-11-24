// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALJORDAN_H_
#define SRC_GENERA_LOCALJORDAN_H_

// clang-format off
#include "PrimeDecomposition.h"
#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_LOCAL_JORDAN
#endif

/*
  The Jordan splittings of a Gram matrix over a number field at a prime P,
  the valuation v_P given by a PrimeValuation: a symmetric bilinear form
  (quadratic lattices over a totally real field K) or a Hermitian form
  (over a CM field E, for the conjugation of E).

  Both repeat the same step: an entry of minimal valuation gives a pivot
  block of rank one (on the diagonal) or two (off it), every other basis
  vector is made orthogonal to the block by the congruence x_l <- x_l -
  c_1 x_i - c_2 x_j (the coefficients P-integral, the block having the
  smallest valuation), and the block is a piece of the component of its
  scale. The scales do not decrease along the splitting, so the components
  come by increasing scale. The rules of the pivot:
  - Quadratic, as jordan_decomposition of Hecke.jl (Lattices.jl): the
    diagonal first; an off-diagonal entry of smaller valuation is moved to
    the diagonal by x_i <- x_i + x_j, except at a dyadic prime, where it
    gives a block of rank two normalized to an off-diagonal entry pi^v
    (pi the uniformizer of P);
  - Hermitian: a block of rank two when an off-diagonal entry has a
    smaller valuation than the diagonal ones, the diagonal first among
    the entries of equal valuation.
  The splittings over Q itself, on rational matrices, are those of
  SiegelSeries.h and GenusSymbol.h; the Jordan data from the lattice with
  its isometry (orders of kernels, at the unramified primes) are in the
  Hermitian masses.
 */

template <typename T> using SmallGram = std::vector<std::vector<ElemE<T>>>;

template <typename T> struct LocalJordanComponent {
  int s;
  std::vector<SmallGram<T>> pieces;
  int rank() const {
    int r = 0;
    for (auto const &G : pieces) {
      r += static_cast<int>(G.size());
    }
    return r;
  }
};

enum class JordanForm { Quadratic, Hermitian };

// The determinant of a piece of rank one or two.
template <typename T>
inline ElemE<T> JordanPieceDet(SmallGram<T> const &G, NumberField<T> const &K) {
  if (G.size() == 1) {
    return G[0][0];
  }
  return EMul(G[0][0], G[1][1], K) - EMul(G[0][1], G[1][0], K);
}

// The determinant of a component, the product over its pieces.
template <typename T>
inline ElemE<T> LocalComponentDet(LocalJordanComponent<T> const &C,
                                  NumberField<T> const &K) {
  ElemE<T> det = EZetaPower(0, K);
  for (auto const &G : C.pieces) {
    det = EMul(det, JordanPieceDet(G, K), K);
  }
  return det;
}

/*
  The splitting. pi is the uniformizer of P, used by the quadratic form at
  a dyadic prime only (nullptr otherwise).
 */
template <typename T>
inline std::vector<LocalJordanComponent<T>>
LocalJordanSplit(SmallGram<T> G, NumberField<T> const &K,
                 PrimeValuation<T> const &vP, JordanForm form,
                 std::type_identity_t<ElemE<T>> const *pi) {
  bool hermitian = (form == JordanForm::Hermitian);
  bool dyadic = (vP.q == 2);
  auto inv = [&](ElemE<T> const &c) -> ElemE<T> {
    return hermitian ? EConj(c, K) : c;
  };
  std::vector<LocalJordanComponent<T>> res;
  auto add_piece = [&](int s, SmallGram<T> const &B) {
    if (res.empty() || res.back().s != s) {
      res.push_back({s, {}});
    }
    res.back().pieces.push_back(B);
  };
  // x_l <- x_l + c x_k, as a congruence
  auto add_mult = [&](int l, int k, ElemE<T> const &c) {
    int n = static_cast<int>(G.size());
    ElemE<T> cc = inv(c);
    for (int u = 0; u < n; u++) {
      G[l][u] += EMul(c, G[k][u], K);
    }
    for (int u = 0; u < n; u++) {
      G[u][l] += EMul(G[u][k], cc, K);
    }
  };
  auto scale_idx = [&](int k, ElemE<T> const &c) {
    int n = static_cast<int>(G.size());
    ElemE<T> cc = inv(c);
    for (int u = 0; u < n; u++) {
      G[k][u] = EMul(c, G[k][u], K);
    }
    for (int u = 0; u < n; u++) {
      G[u][k] = EMul(G[u][k], cc, K);
    }
  };
  // the vectors done: k of them in front (quadratic, the pivots swapped
  // there as in Hecke), or erased (Hermitian), so that the remaining ones
  // keep the order of the two former splittings
  int k = 0;
  auto swap_idx = [&](int a, int b) {
    if (a == b) {
      return;
    }
    std::swap(G[a], G[b]);
    for (auto &row : G) {
      std::swap(row[a], row[b]);
    }
  };
  auto remove_idx = [&](std::vector<int> idx) {
    std::sort(idx.begin(), idx.end(), std::greater<int>());
    for (int i : idx) {
      G.erase(G.begin() + i);
      for (auto &row : G) {
        row.erase(row.begin() + i);
      }
    }
  };
  // the other vectors made orthogonal to the pivot block on idx
  auto clear_against = [&](std::vector<int> const &idx) {
    int n = static_cast<int>(G.size());
    if (idx.size() == 1) {
      int i0 = idx[0];
      ElemE<T> ninv = EInv(G[i0][i0], K);
      for (int l = k; l < n; l++) {
        if (l != i0) {
          // c = -G[l][i0] / G[i0][i0]
          add_mult(l, i0, -EMul(G[l][i0], ninv, K));
        }
      }
      return;
    }
    int i0 = idx[0], j0 = idx[1];
    ElemE<T> a11 = G[i0][i0], a12 = G[i0][j0], a21 = G[j0][i0],
             a22 = G[j0][j0];
    ElemE<T> dinv = EInv(EMul(a11, a22, K) - EMul(a12, a21, K), K);
    for (int l = k; l < n; l++) {
      if (l == i0 || l == j0) {
        continue;
      }
      // (c_1, c_2) = (G[l][i0], G[l][j0]) times the inverse of the block
      ElemE<T> x1 = G[l][i0], x2 = G[l][j0];
      ElemE<T> c1 = EMul(EMul(x1, a22, K) - EMul(x2, a21, K), dinv, K);
      ElemE<T> c2 = EMul(EMul(x2, a11, K) - EMul(x1, a12, K), dinv, K);
      add_mult(l, i0, -c1);
      add_mult(l, j0, -c2);
    }
  };
  while (k < static_cast<int>(G.size())) {
    int n = static_cast<int>(G.size());
    int vmin = 1 << 28, bi = -1, bj = -1;
    if (hermitian) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          int v = PrimeVal(G[i][j], K, vP);
          if (v < vmin || (v == vmin && i == j && bi != bj)) {
            vmin = v;
            bi = i;
            bj = j;
          }
        }
      }
    } else {
      for (int i = k; i < n; i++) {
        int v = PrimeVal(G[i][i], K, vP);
        if (v < vmin) {
          vmin = v;
          bi = i;
          bj = i;
        }
      }
      for (int i = k; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          int v = PrimeVal(G[i][j], K, vP);
          if (v < vmin) {
            vmin = v;
            bi = i;
            bj = j;
          }
        }
      }
    }
#ifdef SANITY_CHECK_LOCAL_JORDAN
    if (bi < 0) {
      std::cerr << "LOCAL_JORDAN: degenerate Gram matrix\n";
      throw TerminalException{1};
    }
#endif
    if (bi != bj && !hermitian && !dyadic) {
      // the small off-diagonal entry moved to the diagonal
      add_mult(bi, bj, EZetaPower(0, K));
      bj = bi;
    }
    if (bi == bj) {
      if (hermitian) {
        clear_against({bi});
        add_piece(vmin, {{G[bi][bi]}});
        remove_idx({bi});
      } else {
        swap_idx(bi, k);
        clear_against({k});
        add_piece(vmin, {{G[k][k]}});
        k++;
      }
      continue;
    }
    if (!hermitian) {
      // the pair to the front, as in Hecke
      swap_idx(bi, k);
      swap_idx(bj == k ? bi : bj, k + 1);
      bi = k;
      bj = k + 1;
      // dyadic: x_i <- pi^v / G_ij x_i, the off-diagonal entry pi^v
      int v12 = PrimeVal(G[bi][bj], K, vP);
      ElemE<T> c = EInv(G[bi][bj], K);
      for (int i = 0; i < v12; i++) {
        c = EMul(c, *pi, K);
      }
      ElemE<T> pinv = EInv(*pi, K);
      for (int i = 0; i < -v12; i++) {
        c = EMul(c, pinv, K);
      }
      scale_idx(bi, c);
    }
    clear_against({bi, bj});
    add_piece(vmin, {{G[bi][bi], G[bi][bj]}, {G[bj][bi], G[bj][bj]}});
    if (hermitian) {
      remove_idx({bi, bj});
    } else {
      k += 2;
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_LOCALJORDAN_H_
// clang-format on
