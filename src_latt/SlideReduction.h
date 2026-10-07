// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SLIDEREDUCTION_H_
#define SRC_LATT_SLIDEREDUCTION_H_
// clang-format off
#include "BKZ.h"
#include "ClassicLLL.h"
#include "DeepLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixInverse.h"
#include "Shvec_exact.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SLIDE
#endif

#ifdef DISABLE_DEBUG_SLIDE
#undef DEBUG_SLIDE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SLIDE
#endif

#ifdef TIMINGS
#define TIMINGS_SLIDE
#endif

/*
  Slide reduction (Gama-Nguyen).

  BKZ enforces one family of conditions, on overlapping blocks, and pays for it
  with a number of tours that has no known polynomial bound. Slide reduction
  enforces TWO families, on blocks that do not overlap, and buys back the
  polynomial bound. With n = p k the two families are, indices being 0-based:

    PRIMAL, for i = 0, ..., p-1: the projected block on indices
      [i k, i k + k - 1] is HKZ-reduced, that is BKZ-reduced at block size k,
      which for a block of rank k is the same thing.

    DUAL, for i = 0, ..., p-2: the projected block on indices
      [i k + 1, (i+1) k] -- the same blocks shifted right by one -- is
      DUAL-reduced: its LAST Gram-Schmidt norm |b^*_{(i+1)k}| is as large as
      possible over the bases of that projected block.

  The primal blocks tile the basis; the dual blocks are the same tiling shifted
  by one, so each straddles exactly one primal boundary. That is the whole
  design. A primal step rearranges vectors inside one block, so it preserves
  the span of every prefix ending at a block boundary and changes none of the
  determinants D_{ik}. A dual step on the block [ik+1, (i+1)k] increases
  |b^*_{(i+1)k}| while the determinant of that block is fixed, so it decreases
  the product of the earlier norms in the block and therefore D_{(i+1)k}
  strictly; and it changes no other D_{jk}, the indices it touches lying
  strictly between the boundaries ik and (i+2)k.

  Hence the potential

      Pi = prod_{i=1}^{p-1} D_{ik},

  a positive integer on an integral form, is untouched by primal steps and
  strictly decreased by every dual step. THAT IS THE POINT OF THE METHOD: the
  number of dual steps is bounded as in the LLL analysis, polynomially in the
  size of the input, where for BKZ no bound on the number of tours is known.
  The quality bound is also slightly better than BKZ's at the same block size,
  the exponent being (n-k)/(2(k-1)) rather than (n-1)/(2(k-1)).

  Two implementation points. The reduction interleaved after each step is a
  SIZE reduction and not a full LLL: size reduction changes no Gram-Schmidt
  norm, so it cannot disturb a block condition just established, whereas LLL is
  free to swap across a block boundary and undo the work. And the dual step is
  carried out on the reversed dual of the block: if the block has Gram matrix G
  and J is the reversal permutation, the reversed dual basis has Gram matrix
  J G^{-1} J, its first Gram-Schmidt norm is the reciprocal of the block's
  last, and a transformation U of the reversed dual corresponds to
  J U^{-T} J of the block. So maximising the last norm of the block is finding
  a shortest vector of J adj(G) J and putting it first, which is the same
  enumeration the primal step uses.

  Reference: N. Gama, P. Q. Nguyen, Finding short lattice vectors within
  Mordell's inequality, STOC 2008, 207--216.
 */

/*
  Slide reduction is defined for a block size k >= 2 dividing the dimension n.
  A form of dimension at most one is reduced whatever the block size.

  The block size is never adjusted to fit: a basis returned as slide reduced
  at block size k is slide reduced at block size k, and a k that does not
  divide n is an error rather than an invitation to run something else.
 */
inline bool SlideIsApplicable(int const &n, int const &k) {
  return n <= 1 || (k >= 2 && n % k == 0);
}

/*
  The dual step on the block of size m starting at index j: make the last
  Gram-Schmidt norm of that projected block maximal. Returns true when the
  block was changed, that is when the condition was violated by more than the
  slack delta.

  blk is the projected block Gram matrix up to a positive scale, with its
  content divided out (BKZ_ProjectedBlockGram); taking it as an argument lets a
  caller that already has every block from one elimination, as self-dual BKZ
  does, avoid recomputing it.
 */
template <typename T, typename Tring, typename Tint>
bool SlideDualStepOnBlock(MyMatrix<Tring> const &blk, MyMatrix<Tring> &gram,
                          MyMatrix<Tint> &H, int const &j, int const &m,
                          Tring const &num, Tring const &den,
                          std::ostream &os) {
  std::pair<MyMatrix<Tring>, Tring> pair = AdjugateDeterminant(blk);
#ifdef SANITY_CHECK_SLIDE
  if (pair.second <= 0) {
    std::cerr << "SLIDE: the projected block is not positive definite, det="
              << pair.second << "\n";
    throw TerminalException{1};
  }
#endif
  // The reversed dual Gram matrix, scaled by det(blk), which changes no
  // shortest vector.
  MyMatrix<Tring> rdual_raw(m, m);
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      rdual_raw(a, b) = pair.first(m - 1 - a, m - 1 - b);
    }
  }
  // Again the content is divided out; the adjugate has raised whatever scale
  // the block carried to the power m - 1.
  MyMatrix<Tring> rdual = RemoveFractionMatrix(rdual_raw);
  MyMatrix<T> rdual_T = UniversalMatrixConversion<T, Tring>(rdual);
  Tshortest<T, Tint> shv = T_ShortestVector<T, Tint>(rdual_T, os);
  Tring min_r = UniversalScalarConversion<Tring, T>(shv.min);
  // rdual(0,0) is the current first squared norm of the reversed dual, whose
  // reciprocal is the block's last Gram-Schmidt norm. Both carry the same
  // scaling, so the comparison is between two integers.
  if (den * min_r >= num * rdual(0, 0)) {
    return false;
  }
  MyVector<Tint> y = GetMatrixRow(shv.SHV, 0);
  MyMatrix<Tint> U = BKZ_UnimodularWithFirstRow(y);
  // A transformation U of the reversed dual is J U^{-T} J of the block.
  MyMatrix<Tint> Uinv = Inverse(U);
  MyMatrix<Tint> V(m, m);
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      V(a, b) = Uinv(m - 1 - b, m - 1 - a);
    }
  }
#ifdef SANITY_CHECK_SLIDE
  Tint det_V = DeterminantMat(V);
  if (det_V != 1 && det_V != -1) {
    std::cerr << "SLIDE: the dual block transformation is not unimodular, det="
              << det_V << "\n";
    throw TerminalException{1};
  }
#endif
  BKZ_ApplyBlockTransformation(gram, H, V, j);
  return true;
}

template <typename T, typename Tring, typename Tint>
bool SlideDualStep(MyMatrix<Tring> &gram, MyMatrix<Tint> &H, int const &j,
                   int const &m, Tring const &num, Tring const &den,
                   std::ostream &os) {
  MyMatrix<Tring> blk = BKZ_ProjectedBlockGram(gram, j, m);
  return SlideDualStepOnBlock<T, Tring, Tint>(blk, gram, H, j, m, num, den, os);
}

template <typename T, typename Tint>
LLLreduction<T, Tint> SlideReducedBasisDelta(MyMatrix<T> const &GramMat,
                                             int const &k,
                                             int const &delta_num,
                                             int const &delta_den,
                                             std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = GramMat.rows();
#ifdef SANITY_CHECK_SLIDE
  if (n != GramMat.cols()) {
    std::cerr << "SLIDE: The matrix should be square\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(GramMat)) {
    std::cerr << "SLIDE: The Gram matrix should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 1) {
    return {GramMat, IdentityMat<Tint>(n)};
  }
  if (!SlideIsApplicable(n, k)) {
    std::cerr << "SLIDE: the block size must be at least 2 and divide the "
                 "dimension, but k="
              << k << " and n=" << n << ". The admissible block sizes are:";
    for (int k_div = 2; k_div <= n; k_div++) {
      if (n % k_div == 0) {
        std::cerr << " " << k_div;
      }
    }
    std::cerr << "\n";
    throw TerminalException{1};
  }
#ifdef SANITY_CHECK_SLIDE
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "SLIDE: The slide reduction needs a positive definite "
                 "matrix\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_SLIDE
  MicrosecondTime time;
#endif
  int p = n / k;
  LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<Tint> H = lll.Pmat;
  MyMatrix<Tring> gram = UniversalMatrixConversion<Tring, T>(
      RemoveFractionMatrix(lll.GramMatRed));
  Tring const num(delta_num);
  Tring const den(delta_den);
  [[maybe_unused]] size_t n_primal = 0;
  [[maybe_unused]] size_t n_dual = 0;
  while (true) {
    // Primal tour: HKZ-reduce each of the p disjoint blocks. This changes no
    // D_{ik} and so cannot make the potential rise.
    for (int i = 0; i < p; i++) {
      int j = i * k;
      MyMatrix<Tring> blk = BKZ_ProjectedBlockGram(gram, j, k);
      MyMatrix<T> blk_T = UniversalMatrixConversion<T, Tring>(blk);
      LLLreduction<T, Tint> hkz = BKZreducedBasis<T, Tint>(blk_T, k, os);
      MyMatrix<Tint> Id = IdentityMat<Tint>(k);
      if (hkz.Pmat != Id) {
        BKZ_ApplyBlockTransformation(gram, H, hkz.Pmat, j);
        IntegralSizeReduce(gram, H);
        n_primal++;
      }
    }
    // Dual tour: the same tiling shifted by one. Each accepted step strictly
    // decreases the potential, which is what bounds the whole computation.
    bool did_dual = false;
    for (int i = 0; i + 1 < p; i++) {
      int j = i * k + 1;
      if (SlideDualStep<T, Tring, Tint>(gram, H, j, k, num, den, os)) {
        IntegralSizeReduce(gram, H);
        n_dual++;
        did_dual = true;
        break;
      }
    }
#ifdef DEBUG_SLIDE
    {
      // diagnostic: the potential must strictly decrease on every dual step
      MyMatrix<Tring> w = gram;
      Tring pv(1);
      Tring pot(1);
      size_t nbits = 0;
      for (int i2 = 0; i2 + 1 < n; i2++) {
        if (i2 > 0 && i2 % k == 0) {
          pot *= w(i2, i2);
        }
        BKZ_BareissStep(w, i2, pv);
      }
      for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
          Tring v = T_abs(gram(a, b));
          size_t nb = 0;
          Tring two_r(2);
          while (v > 0) { v = QuoInt(v, two_r); nb++; }
          if (nb > nbits) nbits = nb;
        }
      }
      os << "SLIDE: iter primal=" << n_primal << " dual=" << n_dual
         << " max_gram_bits=" << nbits << "\n";
    }
#endif
    if (!did_dual) {
      break;
    }
  }
  IntegralSizeReduce(gram, H);
#ifdef DEBUG_SLIDE
  os << "SLIDE: n=" << n << " k=" << k << " p=" << p
     << " primal steps=" << n_primal << " dual steps=" << n_dual << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(H);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(H)};
#ifdef SANITY_CHECK_SLIDE
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_SLIDE
  os << "SLIDE: SlideReducedBasisDelta k=" << k << " took " << time << "\n";
#endif
  return res;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> SlideReducedBasis(MyMatrix<T> const &GramMat,
                                        int const &k, std::ostream &os) {
  return SlideReducedBasisDelta<T, Tint>(GramMat, k, 99, 100, os);
}

// clang-format off
#endif  // SRC_LATT_SLIDEREDUCTION_H_
// clang-format on
