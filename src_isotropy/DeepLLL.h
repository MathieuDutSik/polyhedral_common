// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ISOTROPY_DEEPLLL_H_
#define SRC_ISOTROPY_DEEPLLL_H_
// clang-format off
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_DEEP_LLL
#endif

#ifdef DISABLE_DEBUG_DEEP_LLL
#undef DEBUG_DEEP_LLL
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_DEEP_LLL
#endif

#ifdef TIMINGS
#define TIMINGS_DEEP_LLL
#endif

/*
  Schnorr-Euchner deep insertion.

  LLL keeps a basis size reduced and enforces one inequality per index, the
  Lovasz condition relating b_k to its immediate predecessor. When the
  condition fails the two are swapped. The whole of the discovery is thus
  carried out one adjacent transposition at a time, and a vector that
  geometrically belongs near the front of the basis has to travel there by a
  sequence of them.

  Deep insertion replaces the adjacent swap by a move of unbounded length.
  Write pi_i for the orthogonal projection onto the complement of
  span(b_1, ..., b_{i-1}), so that pi_i(b_i) = b_i^* is the Gram-Schmidt
  vector. The deep condition at index k is

      delta |b_i^*|^2 <= |pi_i(b_k)|^2   for every i < k,

  of which the Lovasz condition is the single case i = k-1. When it fails at
  some i, b_k is removed from its position and inserted at position i, the
  vectors b_i, ..., b_{k-1} shifting up by one:

      (b_1, ..., b_{i-1}, b_k, b_i, ..., b_{k-1}, b_{k+1}, ...).

  The test costs nothing extra. Running i upwards from 1 and maintaining

      c_1 = |b_k|^2,   c_{i+1} = c_i - mu_{i,k}^2 |b_i^*|^2,

  gives c_i = |pi_i(b_k)|^2 at each step, so the whole scan is O(n) on data
  the size reduction has already produced.

  TERMINATION IS A DIFFERENT ARGUMENT FROM LLL'S, and the difference is worth
  stating because it is easy to assume otherwise. Let D_j be the leading
  principal minor of order j of the Gram matrix, so D_j = prod_{l<=j}
  |b_l^*|^2. An insertion at position i from position k:

    * leaves D_1, ..., D_{i-1} unchanged, the first i-1 vectors being
      untouched;
    * leaves D_k, ..., D_n unchanged, the first k vectors being permuted among
      themselves and so spanning the same sublattice;
    * strictly decreases D_i, since the new i-th Gram-Schmidt norm is
      c_i < delta |b_i^*|^2, so D_i' < delta D_i;
    * says nothing whatever about D_j for i < j < k, which may increase.

  For an adjacent swap, i = k-1, the third and fourth points leave no gap and
  the product prod_j D_j strictly decreases: that is the classical LLL
  potential argument. For a genuine deep insertion, k - i >= 2, the product
  may rise, and the potential argument fails. What survives is that the vector
  (D_1, ..., D_{n-1}) strictly decreases in the lexicographic order. On an
  integral form these are positive integers, and the lexicographic order on
  tuples of positive integers is well founded: D_1 is non-increasing so it is
  eventually constant, after which every insertion has i >= 2 and D_2 is
  non-increasing, and so on. Hence the algorithm terminates. The bound this
  argument gives is exponential rather than polynomial, and indeed no
  polynomial bound on the number of deep insertions is known; this is the
  price of the larger move, and it is why the depth of the insertion is
  usually restricted in practice.

  Reference: C. P. Schnorr, M. Euchner, Lattice basis reduction: improved
  practical algorithms and solving subset sum problems, Math. Programming 66
  (1994) 181--199.
 */

/*
  The Gram-Schmidt data of a Gram matrix: B(i) = |b_i^*|^2 and mu(i,j) for
  j < i. Recomputed from row i_start onwards, the earlier rows being taken as
  already correct, which is what an insertion at position i_start leaves
  valid.
 */
template <typename Tfield>
void DeepLLL_UpdateGSO(MyMatrix<Tfield> const &gram, MyMatrix<Tfield> &mu,
                       MyVector<Tfield> &B, int const &i_start) {
  int n = gram.rows();
  for (int r = i_start; r < n; r++) {
    for (int j = 0; j < r; j++) {
      Tfield sum = gram(r, j);
      for (int l = 0; l < j; l++) {
        sum -= mu(r, l) * mu(j, l) * B(l);
      }
      mu(r, j) = sum / B(j);
    }
    Tfield sum = gram(r, r);
    for (int l = 0; l < r; l++) {
      sum -= mu(r, l) * mu(r, l) * B(l);
    }
    B(r) = sum;
  }
}

/*
  Is the position i admissible as an insertion target for index k?

  With depth <= 0 every position is, which is deep insertion as Schnorr and
  Euchner state it. With depth > 0 only the first depth positions and the last
  depth positions before k are tried. That restriction is what makes the
  method usable in higher dimension: the unrestricted version has no
  polynomial bound on its number of insertions, and the positions in the
  middle are empirically the ones that pay least.
 */
inline bool DeepLLL_AdmissiblePosition(int const &i, int const &k,
                                       int const &depth) {
  if (depth <= 0) {
    return true;
  }
  return i < depth || i >= k - depth;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> DeepLLLreducedBasisDepth(MyMatrix<T> const &GramMat,
                                               int const &depth,
                                               [[maybe_unused]]
                                               std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n = GramMat.rows();
#ifdef SANITY_CHECK_DEEP_LLL
  if (n != GramMat.cols()) {
    std::cerr << "DEEPLLL: The matrix should be square\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(GramMat)) {
    std::cerr << "DEEPLLL: The Gram matrix should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 1) {
    return {GramMat, IdentityMat<Tint>(n)};
  }
#ifdef SANITY_CHECK_DEEP_LLL
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "DEEPLLL: The deep LLL reduction needs a positive definite "
                 "matrix\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_DEEP_LLL
  MicrosecondTime time;
#endif
  // The same delta as the classic reduction of this package, so that the two
  // are compared on equal terms and the deep condition is a strict
  // strengthening of the Lovasz one rather than a different trade-off.
  Tfield const delta = Tfield(99) / Tfield(100);
  Tfield const half = Tfield(1) / Tfield(2);
  MyMatrix<Tfield> gram = UniversalMatrixConversion<Tfield, T>(GramMat);
  MyMatrix<Tint> H = IdentityMat<Tint>(n);
  MyMatrix<Tfield> mu = ZeroMatrix<Tfield>(n, n);
  MyVector<Tfield> B(n);
  DeepLLL_UpdateGSO(gram, mu, B, 0);
  //
  // b_k <- b_k - q b_j on the Gram matrix. The row operation is applied to
  // the whole matrix first and the column operation second, reading the
  // already updated entries: that is exactly U gram U^T for U = I - q E_kj,
  // the diagonal entry picking up the -2q gram(k,j) + q^2 gram(j,j) it should.
  //
  auto f_reduce = [&](int const &k, int const &j) -> void {
    if (mu(k, j) <= half && mu(k, j) >= -half) {
      return;
    }
    Tint q = UniversalNearestScalarInteger<Tint, Tfield>(mu(k, j));
    Tfield q_F = UniversalScalarConversion<Tfield, Tint>(q);
    RowSubMul(gram, k, q_F, j);
    ColSubMul(gram, k, q_F, j);
    RowSubMul(H, k, q, j);
    for (int l = 0; l < j; l++) {
      mu(k, l) -= q_F * mu(j, l);
    }
    mu(k, j) -= q_F;
  };
  //
  // Moving index k to position i, the entries in between shifting up. Done as
  // an explicit permutation of the Gram matrix and of the transformation:
  // an in-place rotation of a symmetric matrix is easy to get subtly wrong,
  // and this is not the inner loop.
  //
  auto f_insert = [&](int const &k, int const &i) -> void {
    std::vector<int> perm(n);
    for (int a = 0; a < i; a++) {
      perm[a] = a;
    }
    perm[i] = k;
    for (int a = i + 1; a <= k; a++) {
      perm[a] = a - 1;
    }
    for (int a = k + 1; a < n; a++) {
      perm[a] = a;
    }
    MyMatrix<Tfield> gram_new(n, n);
    MyMatrix<Tint> H_new(n, n);
    for (int a = 0; a < n; a++) {
      for (int b = 0; b < n; b++) {
        gram_new(a, b) = gram(perm[a], perm[b]);
      }
      for (int b = 0; b < n; b++) {
        H_new(a, b) = H(perm[a], b);
      }
    }
    gram = std::move(gram_new);
    H = std::move(H_new);
  };
  //
  [[maybe_unused]] size_t n_insert = 0;
  [[maybe_unused]] size_t n_deep = 0;
  int k = 1;
  while (k < n) {
    for (int j = k - 1; j >= 0; j--) {
      f_reduce(k, j);
    }
    DeepLLL_UpdateGSO(gram, mu, B, k);
    // c runs through |pi_i(b_k)|^2 as i increases, at the cost of one
    // multiplication per step.
    Tfield c = gram(k, k);
    int i_found = -1;
    for (int i = 0; i < k; i++) {
      if (DeepLLL_AdmissiblePosition(i, k, depth) && c < delta * B(i)) {
        i_found = i;
        break;
      }
      c -= mu(k, i) * mu(k, i) * B(i);
    }
    if (i_found >= 0) {
      f_insert(k, i_found);
      DeepLLL_UpdateGSO(gram, mu, B, i_found);
      n_insert++;
      if (i_found < k - 1) {
        n_deep++;
      }
      // Restart at the insertion point: the vectors from there on have moved
      // and their deep conditions have to be established again.
      k = i_found > 1 ? i_found : 1;
    } else {
      k++;
    }
  }
#ifdef DEBUG_DEEP_LLL
  os << "DEEPLLL: n=" << n << " depth=" << depth << " insertions=" << n_insert
     << " of which genuinely deep=" << n_deep << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(H);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(H)};
#ifdef SANITY_CHECK_DEEP_LLL
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_DEEP_LLL
  os << "DEEPLLL: DeepLLLreducedBasisDepth took " << time << "\n";
#endif
  return res;
}

/*
  Deep insertion without restriction on the position, which is the algorithm
  as Schnorr and Euchner state it.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint> DeepLLLreducedBasis(MyMatrix<T> const &GramMat,
                                          std::ostream &os) {
  return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 0, os);
}

template <typename T, typename Tint>
LLLreduction<T, Tint> DeepLLLreducedGeneral(MyMatrix<T> const &GramMat,
                                            std::string const &method,
                                            std::ostream &os) {
  if (method == "full") {
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 0, os);
  }
  if (method == "depth5") {
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 5, os);
  }
  if (method == "depth10") {
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 10, os);
  }
  std::cerr << "DEEPLLL: No matching method for " << method
            << ", allowed are full, depth5, depth10\n";
  throw TerminalException{1};
}

/*
  Tests that a Gram matrix is deep LLL reduced for the given delta and depth:
  size reduced, and satisfying the deep condition at every admissible pair.
  Used by the test program, and available to a caller who has obtained a basis
  by other means and wants to know whether it is already reduced.
 */
template <typename T>
bool IsDeepLLLreduced(MyMatrix<T> const &GramMat, int const &depth,
                      T const &delta, std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n = GramMat.rows();
  if (n <= 1) {
    return true;
  }
  Tfield delta_F = UniversalScalarConversion<Tfield, T>(delta);
  Tfield half = Tfield(1) / Tfield(2);
  MyMatrix<Tfield> gram = UniversalMatrixConversion<Tfield, T>(GramMat);
  MyMatrix<Tfield> mu = ZeroMatrix<Tfield>(n, n);
  MyVector<Tfield> B(n);
  DeepLLL_UpdateGSO(gram, mu, B, 0);
  for (int k = 1; k < n; k++) {
    for (int j = 0; j < k; j++) {
      if (mu(k, j) > half || mu(k, j) < -half) {
        os << "DEEPLLL: not size reduced, mu(" << k << "," << j
           << ")=" << mu(k, j) << "\n";
        return false;
      }
    }
    Tfield c = gram(k, k);
    for (int i = 0; i < k; i++) {
      if (DeepLLL_AdmissiblePosition(i, k, depth) && c < delta_F * B(i)) {
        os << "DEEPLLL: the deep condition fails at k=" << k << " i=" << i
           << ", |pi_i(b_k)|^2=" << c << " against delta |b_i^*|^2="
           << delta_F * B(i) << "\n";
        return false;
      }
      c -= mu(k, i) * mu(k, i) * B(i);
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_ISOTROPY_DEEPLLL_H_
// clang-format on
