// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PROJECTIVE_SYSTEMS_PICARD_REDUCTION_H_
#define SRC_PROJECTIVE_SYSTEMS_PICARD_REDUCTION_H_

// clang-format off
#include "projective_complex.h"
#include <array>
#include <utility>
#include <vector>
// clang-format on

/*
  The reduction of the cusps of the Picard modular group SU(2,1; Z[i])
  ("walk cusp_reduction" in the description), with the form of D. Yasaki
  (Explicit reduction for SU(2,1; Z[i])): Q(u, v) = u^* Q v,
    Q = [[0, 0, i], [0, -1, 0], [-i, 0, 0]],
  the group acting on the column vectors (n, p, q)^T. The description uses
  the row vectors, on which the matrices act on the right: M = A^T for A
  acting on the columns, and the form is h = Q^T.

  A Euclidean algorithm on an isotropic vector v = (n, p, q), q != 0:
  - a translation N(beta, r) = [[1, beta, r + i |beta|^2 / 2], [0, 1, i
    conj(beta)], [0, 0, 1]] of the integral Heisenberg group (beta in
    (1 + i) Z[i], r in Z) maps v to (n + beta p + (r + i |beta|^2 / 2) q,
    p + i conj(beta) q, q); beta is chosen with |p'/q| <= 1 (the covering
    radius of (1 + i) Z[i]) and r with |Re(n'/q)| <= 1/2, so that the
    isotropy |p'|^2 = 2 Im(n' conj(q)) gives |n'/q|^2 <= 1/2,
  - then w = [[0, 0, -1], [0, 1, 0], [1, 0, 0]] maps it to (-q, p', n'),
    whose last entry is smaller.
  It ends at q = 0, the cusp e_1 = (1, 0, 0). With G_k the product of the
  first k steps, the cusps G_k^{-1} e_1 go from e_1 to v, and two
  consecutive ones are A_k {e_1, e_3} with A_k = G_k^{-1} N_{k+1}^{-1}: a
  path of edges equivalent to {e_1, e_3}, as in the continued fractions
  of Manin.
 */

namespace psys {

using GaussElt = std::pair<T, T>;
using GaussMat = std::array<std::array<GaussElt, 3>, 3>;

inline ImagQuadRing<T> gaussian_ring() { return ImagQuadRing<T>{T(0), T(1)}; }

inline GaussMat gauss_identity() {
  GaussMat M;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      M[i][j] = {T(i == j ? 1 : 0), T(0)};
    }
  }
  return M;
}

inline GaussMat gauss_mul(GaussMat const &A, GaussMat const &B) {
  ImagQuadRing<T> ring = gaussian_ring();
  GaussMat C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      GaussElt s{T(0), T(0)};
      for (int k = 0; k < 3; k++) {
        s = ring.add(s, ring.mul(A[i][k], B[k][j]));
      }
      C[i][j] = s;
    }
  }
  return C;
}

// The nearest integer, the ties rounded down.
inline T round_rational(T const &x) {
  T x2 = x + T(1, 2);
  Tint num = x2.get_num();
  Tint den = x2.get_den();
  Tint q;
  mpz_fdiv_q(q.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());
  return T(q);
}

// N(beta, r), acting on the columns.
inline GaussMat heisenberg(GaussElt const &beta, T const &r) {
  ImagQuadRing<T> ring = gaussian_ring();
  GaussMat N = gauss_identity();
  N[0][1] = beta;
  N[0][2] = {r, ring.norm(beta) / T(2)};
  N[1][2] = ring.mul({T(0), T(1)}, ring.conj(beta));
  return N;
}

/*
  The realification of the transpose of A (A acting on the columns): the
  matrix acting on the right on the realified row vectors (a_1, a_2, a_3,
  b_1, b_2, b_3).
 */
inline MyMatrix<Tint> realify_transpose(GaussMat const &A) {
  MyMatrix<Tint> R(6, 6);
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      // The entry (k, l) of A^T is A[l][k].
      GaussElt c = A[l][k];
      Tint a = UniversalScalarConversion<Tint, T>(c.first);
      Tint b = UniversalScalarConversion<Tint, T>(c.second);
      R(k, l) = a;
      R(k, 3 + l) = b;
      R(3 + k, l) = -b;
      R(3 + k, 3 + l) = a;
    }
  }
  return R;
}

// Whether the ring and the form are those of Yasaki.
inline bool is_picard_gaussian(ProjectiveRing const &ring) {
  if (ring.d != 2 || ring.m != 3 || ring.s != 0 || ring.p != 1) {
    return false;
  }
  MyMatrix<Tint> h_a = ZeroMatrix<Tint>(3, 3);
  MyMatrix<Tint> h_b = ZeroMatrix<Tint>(3, 3);
  h_a(1, 1) = -1;
  h_b(0, 2) = -1;
  h_b(2, 0) = 1;
  return ring.B == get_real_form(3, 0, 1, h_a, h_b);
}

/*
  The result of the reduction of v: the matrices M_0, ..., M_{K-1}
  (acting on the realified row vectors) such that the rows of {e_1, e_3}
  M_k are consecutive edges of a path from e_1 to the line of v, and the
  matrix R_G of the product G of the steps, with v R_G on the line of e_1.
 */
struct PicardReduction {
  std::vector<MyMatrix<Tint>> l_M;
  MyMatrix<Tint> R_G;
};

inline PicardReduction picard_cusp_reduction(ProjectiveRing const &ring, MyVector<Tint> const &v_row) {
  if (!is_picard_gaussian(ring)) {
    std::cerr << "PSYS: the walk cusp_reduction needs the Picard group over Z[i] with the form of Yasaki\n";
    throw TerminalException{1};
  }
  ImagQuadRing<T> gring = gaussian_ring();
  std::array<GaussElt, 3> v;
  for (int j = 0; j < 3; j++) {
    v[j] = {T(v_row(j)), T(v_row(3 + j))};
  }
  GaussMat w_inv = gauss_identity();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      w_inv[i][j] = {T(0), T(0)};
    }
  }
  // w = [[0, 0, -1], [0, 1, 0], [1, 0, 0]], w^{-1} = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]].
  w_inv[0][2] = {T(1), T(0)};
  w_inv[1][1] = {T(1), T(0)};
  w_inv[2][0] = {T(-1), T(0)};
  GaussMat w = gauss_identity();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      w[i][j] = {T(0), T(0)};
    }
  }
  w[0][2] = {T(-1), T(0)};
  w[1][1] = {T(1), T(0)};
  w[2][0] = {T(1), T(0)};
  GaussMat G = gauss_identity();
  GaussMat G_inv = gauss_identity();
  std::vector<MyMatrix<Tint>> l_M;
  GaussElt one_plus_i{T(1), T(1)};
  int n_iter = 0;
  while (!gring.is_zero(v[2])) {
    n_iter += 1;
    if (n_iter > 1000) {
      std::cerr << "PSYS: picard_cusp_edge_path, the reduction does not end\n";
      throw TerminalException{1};
    }
    GaussElt n = v[0], p = v[1], q = v[2];
    // lambda in (1 + i) Z[i] nearest to p / q, beta = -i conj(lambda).
    GaussElt z = gring.div(gring.div(p, q), one_plus_i);
    GaussElt t{round_rational(z.first), round_rational(z.second)};
    GaussElt lambda = gring.mul(one_plus_i, t);
    GaussElt beta = gring.mul({T(0), T(-1)}, gring.conj(lambda));
    // n' = n + beta p + i |beta|^2 / 2 q, then r with |Re(n'' / q)| <= 1/2.
    GaussElt n1 = gring.add(gring.add(n, gring.mul(beta, p)), gring.mul({T(0), gring.norm(beta) / T(2)}, q));
    T r = -round_rational(gring.div(n1, q).first);
    GaussMat N = heisenberg(beta, r);
    GaussMat N_inv = heisenberg(gring.neg(beta), -r);
    // The edge A_k {e_1, e_3}, A_k = G_k^{-1} N^{-1}.
    l_M.push_back(realify_transpose(gauss_mul(G_inv, N_inv)));
    // v -> w N v, G^{-1} -> G^{-1} N^{-1} w^{-1}.
    std::array<GaussElt, 3> v_new;
    for (int i = 0; i < 3; i++) {
      GaussElt s{T(0), T(0)};
      for (int k = 0; k < 3; k++) {
        s = gring.add(s, gring.mul(N[i][k], v[k]));
      }
      v_new[i] = s;
    }
    // w (n, p, q) = (-q, p, n).
    v = {gring.neg(v_new[2]), v_new[1], v_new[0]};
    G_inv = gauss_mul(gauss_mul(G_inv, N_inv), w_inv);
    G = gauss_mul(gauss_mul(w, N), G);
  }
#ifdef SANITY_CHECK_PROJECTIVE_COMPLEX
  for (auto &M : l_M) {
    if (!ring.in_group(M)) {
      std::cerr << "PSYS: picard_cusp_edge_path, an edge transformation is not in the group\n";
      throw TerminalException{1};
    }
  }
#endif
  return {l_M, realify_transpose(G)};
}

}  // namespace psys

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_PICARD_REDUCTION_H_
// clang-format on
