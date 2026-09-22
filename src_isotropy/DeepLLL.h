// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ISOTROPY_DEEPLLL_H_
#define SRC_ISOTROPY_DEEPLLL_H_
// clang-format off
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "QuoIntFcts.h"
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
  Schnorr-Euchner deep insertion, fraction free.

  LLL keeps a basis size reduced and enforces one inequality per index, the
  Lovasz condition relating b_k to its immediate predecessor; when it fails the
  two are swapped. The whole of the discovery is thus carried out one adjacent
  transposition at a time, and a vector that geometrically belongs near the
  front of the basis has to travel there by a sequence of them.

  Deep insertion replaces the adjacent swap by a move of unbounded length.
  Write pi_i for the orthogonal projection onto the complement of
  span(b_1, ..., b_{i-1}), so that pi_i(b_i) = b_i^* is the Gram-Schmidt
  vector. The deep condition at index k is

      delta |b_i^*|^2 <= |pi_i(b_k)|^2   for every i < k,

  of which the Lovasz condition is the single case i = k-1. When it fails at
  some i, b_k is removed and inserted at position i, the vectors b_i, ...,
  b_{k-1} shifting up by one.

  TERMINATION IS A DIFFERENT ARGUMENT FROM LLL'S, and the difference is worth
  stating because it is easy to assume otherwise. Let D_j be the leading
  principal minor of order j of the Gram matrix, so D_j = prod_{l<=j}
  |b_l^*|^2. An insertion at position i from position k:

    * leaves D_1, ..., D_{i-1} unchanged, the first i-1 vectors being
      untouched;
    * leaves D_k, ..., D_n unchanged, the first k vectors being permuted among
      themselves and so spanning the same sublattice;
    * strictly decreases D_i, the new i-th Gram-Schmidt norm being
      |pi_i(b_k)|^2 < delta |b_i^*|^2;
    * says nothing whatever about D_j for i < j < k, which may increase.

  For an adjacent swap, i = k-1, the fourth point is vacuous and the product
  prod_j D_j strictly decreases: that is the classical LLL potential argument.
  For a genuine deep insertion, k - i >= 2, the product may rise and the
  argument fails. What survives is that (D_1, ..., D_{n-1}) strictly decreases
  in the LEXICOGRAPHIC order, which over the positive integers is well founded:
  D_1 is non-increasing so eventually constant, after which every insertion has
  i >= 2 and D_2 is non-increasing, and so on. Hence termination. The bound is
  exponential rather than polynomial, and indeed no polynomial bound on the
  number of deep insertions is known; this is the price of the larger move and
  the reason the depth is usually restricted.

  THE ARITHMETIC IS INTEGRAL THROUGHOUT. The obvious implementation keeps the
  mu_{j,k} as rationals, and then pays for it: the numerators and denominators
  grow with the minors of the Gram matrix, and the data has to be recomputed
  after every insertion. Instead we keep de Weger's integral form,

      d_i = D_i  (with d_0 = 1),     lambda_{i,j} = d_j mu_{j,i},

  both integral, from which mu_{j,i} = lambda_{i,j}/d_j and
  |b_i^*|^2 = d_i/d_{i-1}. Three things then become integer operations.

    * The Gram-Schmidt data itself, by the Bareiss-type recursion of
      DeepLLL_IntegralGSO below, whose divisions are exact.

    * The projected norms. Setting S_i = d_{i-1} |pi_i(b_k)|^2, one has
      S_1 = |b_k|^2 and

          S_{i+1} = (S_i d_i - lambda_{k,i}^2) / d_{i-1},

      an exact division. S_i is not an artefact of the algebra: it is the Gram
      determinant of (b_1, ..., b_{i-1}, b_k), which is why it is a positive
      integer and why the division comes out even.

    * The deep test. With delta = num/den,

          |pi_i(b_k)|^2 < delta |b_i^*|^2   <=>   den S_i < num d_i,

      the denominators cancelling. So the condition that drives the whole
      algorithm is a comparison of two integers.

  Size reduction is integral for the same reason: |mu_{j,k}| <= 1/2 reads
  2|lambda_{k,j}| <= d_j, and the rounded quotient is
  floor((2 lambda_{k,j} + d_j) / (2 d_j)).

  Reference: C. P. Schnorr, M. Euchner, Lattice basis reduction: improved
  practical algorithms and solving subset sum problems, Math. Programming 66
  (1994) 181--199. The integral Gram-Schmidt is de Weger's, as presented in
  Cohen, A Course in Computational Algebraic Number Theory, Algorithm 2.6.7.
 */

/*
  The integral Gram-Schmidt data of an integral Gram matrix: d(i) for
  i = 0, ..., n with d(0) = 1 and d(i) the leading principal minor of order i,
  and lambda(i,j) = d(j+1) mu(j,i) for j < i.

  DeepLLL_IntegralGSO_Row computes ONE row, assuming the rows before it
  correct. That is the form the descent wants, and the reason is that the
  descent at index k never looks above row k: the deep test reads lambda(k,i)
  and d(i) for i <= k only. Recomputing rows k, ..., n-1 at every step, which
  is the obvious thing to do and what this did at first, therefore spends
  O(n^3) per step on data that will be recomputed again before it is read,
  and makes a sweep O(n^4) where O(n^3) suffices.
 */
template <typename Tring>
void DeepLLL_IntegralGSO_Row(MyMatrix<Tring> const &gram,
                             MyMatrix<Tring> &lambda, std::vector<Tring> &d,
                             int const &i) {
  d[0] = Tring(1);
  for (int j = 0; j <= i; j++) {
    Tring u = gram(i, j);
    for (int k = 0; k < j; k++) {
      Tring num = d[k + 1] * u - lambda(i, k) * lambda(j, k);
      Tring quot = num / d[k];
#ifdef SANITY_CHECK_DEEP_LLL
      if (quot * d[k] != num) {
        std::cerr << "DEEPLLL: non-exact division in the integral "
                     "Gram-Schmidt, the input is not an integral form over "
                     "an integral domain\n";
        throw TerminalException{1};
      }
#endif
      u = quot;
    }
    if (j < i) {
      lambda(i, j) = u;
    } else {
      d[i + 1] = u;
    }
  }
}

template <typename Tring>
void DeepLLL_IntegralGSO(MyMatrix<Tring> const &gram,
                         MyMatrix<Tring> &lambda, std::vector<Tring> &d,
                         int const &i_start) {
  int n = gram.rows();
  for (int i = i_start; i < n; i++) {
    DeepLLL_IntegralGSO_Row(gram, lambda, d, i);
  }
}

/*
  Is the position i admissible as an insertion target for index k?

  With depth <= 0 every position is, which is deep insertion as Schnorr and
  Euchner state it. With depth > 0 only the first depth positions and the last
  depth positions before k are tried. The tail positions include i = k-1, the
  LLL swap, so the restricted algorithm still produces an LLL reduced basis for
  every depth >= 1. That restriction is what makes the method usable in higher
  dimension, the unrestricted version having no polynomial bound on its number
  of insertions, and the middle positions being the ones that discard the most
  work when they fire.
 */
inline bool DeepLLL_AdmissiblePosition(int const &i, int const &k,
                                       int const &depth) {
  if (depth <= 0) {
    return true;
  }
  return i < depth || i >= k - depth;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> DeepLLLreducedBasisDepthDelta(MyMatrix<T> const &GramMat,
                                                    int const &depth,
                                                    int const &delta_num,
                                                    int const &delta_den,
                                                    [[maybe_unused]]
                                                    std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
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
  if (delta_num <= 0 || delta_den <= 0 || delta_num >= delta_den) {
    std::cerr << "DEEPLLL: delta must lie strictly between 0 and 1, got "
              << delta_num << "/" << delta_den << "\n";
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
  // Both conditions of the algorithm are invariant under a positive rescaling
  // of the form, so the descent may be run on an integral rescaling and the
  // transformation it finds is the same. That is what puts the whole
  // computation over the integers.
  MyMatrix<Tring> gram =
      UniversalMatrixConversion<Tring, T>(RemoveFractionMatrix(GramMat));
  MyMatrix<Tint> H = IdentityMat<Tint>(n);
  MyMatrix<Tring> lambda = ZeroMatrix<Tring>(n, n);
  std::vector<Tring> d(n + 1, Tring(0));
  // Rows 0, ..., n_valid-1 of the Gram-Schmidt data are correct, together with
  // d(0), ..., d(n_valid). Rows are brought up only as the descent reaches
  // them, and invalidated only as far down as a move actually reaches.
  int n_valid = 0;
  auto f_ensure = [&](int const &r) -> void {
    for (int rr = n_valid; rr <= r; rr++) {
      DeepLLL_IntegralGSO_Row(gram, lambda, d, rr);
    }
    if (r + 1 > n_valid) {
      n_valid = r + 1;
    }
  };
  Tring const two(2);
  Tring const num(delta_num);
  Tring const den(delta_den);
  //
  // b_k <- b_k - q b_j on the Gram matrix. The row operation is applied to the
  // whole matrix first and the column operation second, reading the already
  // updated entries: that is exactly U gram U^T for U = I - q E_kj, the
  // diagonal entry picking up the -2q gram(k,j) + q^2 gram(j,j) it should.
  // The minors d are untouched, size reduction not changing the flag.
  //
  auto f_reduce = [&](int const &k, int const &j) -> void {
    Tring abs_lam = T_abs(lambda(k, j));
    if (two * abs_lam <= d[j + 1]) {
      return;
    }
    // The nearest integer to mu = lambda/d, with ties resolved DOWNWARDS so
    // as to agree with NearestInteger of the package, which the classic
    // reduction uses: q = ceil(mu - 1/2) = -floor((d - 2 lambda) / (2 d)).
    // Rounding ties the other way is just as correct -- both leave
    // |mu| <= 1/2 -- but would give a different, equally reduced, basis, and
    // the two reductions should not disagree on so small a thing.
    Tring quo_num = d[j + 1] - two * lambda(k, j);
    Tring quo_den = two * d[j + 1];
    Tring q = -QuoInt(quo_num, quo_den);
    if (q == 0) {
      return;
    }
    RowSubMul(gram, k, q, j);
    ColSubMul(gram, k, q, j);
    Tint q_int = UniversalScalarConversion<Tint, Tring>(q);
    RowSubMul(H, k, q_int, j);
    for (int l = 0; l < j; l++) {
      lambda(k, l) -= q * lambda(j, l);
    }
    lambda(k, j) -= q * d[j + 1];
  };
  //
  // Moving index k to position i, the entries in between shifting up. Done as
  // an explicit permutation of the Gram matrix and of the transformation: an
  // in-place rotation of a symmetric matrix is easy to get subtly wrong, and
  // this is not the inner loop.
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
    MyMatrix<Tring> gram_new(n, n);
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
    f_ensure(k);
    for (int j = k - 1; j >= 0; j--) {
      f_reduce(k, j);
    }
    // Size reduction of row k leaves that row's data correct, f_reduce
    // updating lambda(k,.) by the same transvection and d being untouched, the
    // flag being unchanged. Rows above k read gram(.,k), which did change.
    if (n_valid > k + 1) {
      n_valid = k + 1;
    }
    // S runs through d(i) |pi_i(b_k)|^2 as i increases: the Gram determinant
    // of (b_0, ..., b_{i-1}, b_k), an integer, at the cost of one
    // multiplication and one exact division per step.
    Tring S = gram(k, k);
    int i_found = -1;
    for (int i = 0; i < k; i++) {
      if (DeepLLL_AdmissiblePosition(i, k, depth) && den * S < num * d[i + 1]) {
        i_found = i;
        break;
      }
      Tring next = S * d[i + 1] - lambda(k, i) * lambda(k, i);
      Tring quot = next / d[i];
#ifdef SANITY_CHECK_DEEP_LLL
      if (quot * d[i] != next) {
        std::cerr << "DEEPLLL: non-exact division in the projected norm "
                     "recursion\n";
        throw TerminalException{1};
      }
      if (quot <= 0) {
        std::cerr << "DEEPLLL: the projected norm " << quot
                  << " is not positive, which a positive definite form "
                     "forbids\n";
        throw TerminalException{1};
      }
#endif
      S = quot;
    }
    if (i_found >= 0) {
      f_insert(k, i_found);
      if (n_valid > i_found) {
        n_valid = i_found;
      }
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
  os << "DEEPLLL: n=" << n << " depth=" << depth << " delta=" << delta_num
     << "/" << delta_den << " insertions=" << n_insert
     << " of which genuinely deep=" << n_deep << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(H);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(H)};
#ifdef SANITY_CHECK_DEEP_LLL
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_DEEP_LLL
  os << "DEEPLLL: DeepLLLreducedBasisDepthDelta took " << time << "\n";
#endif
  return res;
}

/*
  The same delta as the classic reduction of this package, so that the two are
  compared on equal terms and the deep condition is a strict strengthening of
  the Lovasz one rather than a different trade-off.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint> DeepLLLreducedBasisDepth(MyMatrix<T> const &GramMat,
                                               int const &depth,
                                               std::ostream &os) {
  return DeepLLLreducedBasisDepthDelta<T, Tint>(GramMat, depth, 99, 100, os);
}

/*
  Deep insertion without restriction on the position, which is the algorithm as
  Schnorr and Euchner state it.
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
  Tests that a Gram matrix is LLL reduced for the given delta: size reduced and
  satisfying the Lovasz condition, and nothing more.

  This is NOT IsDeepLLLreduced at depth 1. The admissible set at depth 1 is
  {0, k-1}, which adds the deep condition at i = 0, that is
  delta |b_1^*|^2 <= |b_k|^2 for every k. A deep reduced basis satisfies it but
  an LLL reduced one need not: LLL is perfectly content to leave b_3 shorter
  than b_1. Using the depth-1 test as a test for LLL reducedness therefore
  rejects legitimate bases, and did.
 */
template <typename T>
bool IsLLLreduced(MyMatrix<T> const &GramMat, int const &delta_num,
                  int const &delta_den, std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = GramMat.rows();
  if (n <= 1) {
    return true;
  }
  MyMatrix<Tring> gram =
      UniversalMatrixConversion<Tring, T>(RemoveFractionMatrix(GramMat));
  MyMatrix<Tring> lambda = ZeroMatrix<Tring>(n, n);
  std::vector<Tring> d(n + 1, Tring(0));
  DeepLLL_IntegralGSO(gram, lambda, d, 0);
  Tring const two(2);
  Tring const num(delta_num);
  Tring const den(delta_den);
  for (int k = 1; k < n; k++) {
    for (int j = 0; j < k; j++) {
      if (two * T_abs(lambda(k, j)) > d[j + 1]) {
        os << "DEEPLLL: not size reduced, 2|lambda(" << k << "," << j
           << ")|=" << two * T_abs(lambda(k, j)) << " exceeds d=" << d[j + 1]
           << "\n";
        return false;
      }
    }
    // The Lovasz condition delta |b_{k-1}^*|^2 <= |pi_{k-1}(b_k)|^2 reads
    // num d(k) <= den S_{k-1}, the common factor d(k-1) cancelling.
    Tring S = gram(k, k);
    for (int i = 0; i + 1 < k; i++) {
      S = (S * d[i + 1] - lambda(k, i) * lambda(k, i)) / d[i];
    }
    if (den * S < num * d[k]) {
      os << "DEEPLLL: the Lovasz condition fails at k=" << k << ", " << den
         << "*" << S << " < " << num << "*" << d[k] << "\n";
      return false;
    }
  }
  return true;
}

/*
  Tests that a Gram matrix is deep LLL reduced for the given delta and depth:
  size reduced, and satisfying the deep condition at every admissible pair.
  Recomputes the integral Gram-Schmidt data from scratch, so it is an
  independent check of a descent and not a restatement of it. Also available to
  a caller who has obtained a basis by other means and wants to know whether it
  is already reduced.
 */
template <typename T>
bool IsDeepLLLreduced(MyMatrix<T> const &GramMat, int const &depth,
                      int const &delta_num, int const &delta_den,
                      std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = GramMat.rows();
  if (n <= 1) {
    return true;
  }
  MyMatrix<Tring> gram =
      UniversalMatrixConversion<Tring, T>(RemoveFractionMatrix(GramMat));
  MyMatrix<Tring> lambda = ZeroMatrix<Tring>(n, n);
  std::vector<Tring> d(n + 1, Tring(0));
  DeepLLL_IntegralGSO(gram, lambda, d, 0);
  Tring const two(2);
  Tring const num(delta_num);
  Tring const den(delta_den);
  for (int k = 1; k < n; k++) {
    for (int j = 0; j < k; j++) {
      if (two * T_abs(lambda(k, j)) > d[j + 1]) {
        os << "DEEPLLL: not size reduced, 2|lambda(" << k << "," << j
           << ")|=" << two * T_abs(lambda(k, j)) << " exceeds d=" << d[j + 1]
           << "\n";
        return false;
      }
    }
    Tring S = gram(k, k);
    for (int i = 0; i < k; i++) {
      if (DeepLLL_AdmissiblePosition(i, k, depth) && den * S < num * d[i + 1]) {
        os << "DEEPLLL: the deep condition fails at k=" << k << " i=" << i
           << ", " << den << "*" << S << " < " << num << "*" << d[i + 1]
           << "\n";
        return false;
      }
      S = (S * d[i + 1] - lambda(k, i) * lambda(k, i)) / d[i];
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_ISOTROPY_DEEPLLL_H_
// clang-format on
