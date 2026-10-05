// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SELFDUALBKZ_H_
#define SRC_LATT_SELFDUALBKZ_H_
// clang-format off
#include "BKZ.h"
#include "ClassicLLL.h"
#include "DeepLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "SlideReduction.h"
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SDBKZ
#endif

#ifdef DISABLE_DEBUG_SDBKZ
#undef DEBUG_SDBKZ
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SDBKZ
#endif

#ifdef TIMINGS
#define TIMINGS_SDBKZ
#endif

/*
  Self-dual BKZ (Micciancio-Walter).

  BKZ asks, for each block [j, j+beta-1] of the projected basis, that its FIRST
  Gram-Schmidt norm be as small as possible: b_j^* is a shortest vector of the
  projected block. The dual condition asks that its LAST Gram-Schmidt norm be as
  large as possible, which is the same condition on the reversed dual of the
  block. Self-dual BKZ alternates the two:

    FORWARD TOUR, j = 0, ..., n-beta: if the block [j, j+beta-1] has a vector
      shorter than delta |b_j^*|, insert a shortest one at position j.

    BACKWARD TOUR, j = n-beta, ..., 0: if the last Gram-Schmidt norm of the
      block [j, j+beta-1] can be made larger than |b_{j+beta-1}^*| / delta,
      make it maximal.

  Only full blocks of size beta are used. The forward tour is the primal step
  of BKZ, and the backward tour is the dual step of slide reduction, carried
  out through the reversed dual J adj(G) J of the block (see
  SlideReduction.h); what self-dual BKZ adds is to apply the dual step to every
  block, overlapping as in BKZ, rather than to a fixed tiling as slide
  reduction does. Its interest, in the analysis of Micciancio and Walter, is
  that it has proven bounds comparable to slide reduction's and BKZ's
  practical behaviour:
  measured with fplll on the SVP challenge lattices of dimension 60, it reaches
  at block size 20 the profile BKZ reaches at block size 30.

  AFTER EACH INSERTION the block alone is LLL reduced, then the whole basis is
  size reduced; there is no global LLL. Neither pass can undo the insertion:
  in a forward step b_j is a shortest vector of the block, so no Lovasz swap
  can move it, and in a backward step the last Gram-Schmidt norm is maximal,
  so no swap can increase it, while a swap at an earlier pair does not change
  it. Size reduction changes no Gram-Schmidt norm at all. Avoiding the global
  LLL matters: in BKZ.h it is the full LLL after each insertion that dominates
  the cost at small block sizes.

  TERMINATION. The termination argument of BKZ does not carry over. A forward
  insertion at j decreases D_j and fixes D_i for i < j, so it decreases
  (D_1, ..., D_{n-1}) lexicographically; a backward step on [j, k] decreases
  D_{k-1} and fixes D_i for i >= k, but may increase D_j, ..., D_{k-2}, so it
  decreases the same vector in the REVERSE lexicographic order. No single
  order serves both, and Micciancio and Walter analyse the convergence of the
  profile rather than the termination of the exact loop. The computation is
  therefore in two phases:
    * the self-dual rounds, a forward and a backward tour each, run while a
      round decreases strictly the LLL potential prod_{i=1}^{n-1} D_i, a
      positive integer on an integral form; the round that fails to decrease
      it is kept, since the potential is not what the rounds improve;
    * forward tours alone, until one makes no insertion. These terminate by
      the argument of BKZ, the lexicographic one above.
  On return every full block [j, j+beta-1] satisfies the condition of BKZ up to
  delta: b_j^* is a shortest vector of the projected block. The backward
  conditions are those of the last self-dual round and are not certified after
  the final phase. The CI section 16B checks the forward conditions;
  DEBUG_SDBKZ reports the rounds and the insertions of each phase.

  Reference: D. Micciancio, M. Walter, Practical, predictable lattice basis
  reduction, EUROCRYPT 2016, LNCS 9665, 820--849.
 */

/*
  The scaled projected block Gram matrices of every block of size m, for
  j = 0, ..., n-m, from one elimination: the block at j is the trailing m x m
  block of the matrix after j steps of the fraction-free elimination, as in
  BKZ.h. One pass costs O(n^3), where recomputing each block on its own would
  cost O(n^4).
 */
template <typename Tring>
std::vector<MyMatrix<Tring>> SDBKZ_ProjectedBlocks(MyMatrix<Tring> const &gram,
                                                   int const &m) {
  int n = gram.rows();
  MyMatrix<Tring> work = gram;
  Tring prev(1);
  std::vector<MyMatrix<Tring>> l_block;
  for (int j = 0; j + m <= n; j++) {
    MyMatrix<Tring> blk(m, m);
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        blk(a, b) = work(j + a, j + b);
      }
    }
    l_block.emplace_back(std::move(blk));
    if (j + m < n) {
      BKZ_BareissStep(work, j, prev);
    }
  }
  return l_block;
}

/*
  The LLL potential prod_{i=1}^{n-1} D_i, D_i the leading principal minors.
  After i steps of the elimination the diagonal entry i is D_{i+1}.
 */
template <typename Tring> Tring SDBKZ_Potential(MyMatrix<Tring> const &gram) {
  int n = gram.rows();
  MyMatrix<Tring> work = gram;
  Tring prev(1);
  Tring pot(1);
  for (int i = 0; i + 1 < n; i++) {
    pot *= work(i, i);
    BKZ_BareissStep(work, i, prev);
  }
  return pot;
}

/*
  LLL reduce the projected block [j, j+m-1] alone and apply the transformation
  to the whole basis. The vectors before j are unchanged, so the conditions
  established on earlier blocks are not disturbed.
 */
template <typename T, typename Tring, typename Tint>
void SDBKZ_LLLBlock(MyMatrix<Tring> &gram, MyMatrix<Tint> &H, int const &j,
                    int const &m, std::ostream &os) {
  MyMatrix<Tring> blk = BKZ_ProjectedBlockGram(gram, j, m);
  MyMatrix<T> blk_T = UniversalMatrixConversion<T, Tring>(blk);
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(blk_T, os);
  if (red.Pmat != IdentityMat<Tint>(m)) {
    BKZ_ApplyBlockTransformation(gram, H, red.Pmat, j);
  }
}

/*
  Self-dual BKZ at block size beta with the slack delta = delta_num /
  delta_den. max_tour caps the number of rounds (a forward and a backward
  tour); 0 means no cap, the loop then stopping by the rule above.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint>
SelfDualBKZReducedBasisDelta(MyMatrix<T> const &GramMat, int const &beta,
                             int const &delta_num, int const &delta_den,
                             int const &max_tour, std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = GramMat.rows();
#ifdef SANITY_CHECK_SDBKZ
  if (n != GramMat.cols()) {
    std::cerr << "SDBKZ: The matrix should be square\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(GramMat)) {
    std::cerr << "SDBKZ: The Gram matrix should be symmetric\n";
    throw TerminalException{1};
  }
  if (beta < 2) {
    std::cerr << "SDBKZ: the block size must be at least 2, got " << beta
              << "\n";
    throw TerminalException{1};
  }
  if (delta_num <= 0 || delta_den <= 0 || delta_num > delta_den) {
    std::cerr << "SDBKZ: delta must lie in (0,1], got " << delta_num << "/"
              << delta_den << "\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 1) {
    return {GramMat, IdentityMat<Tint>(n)};
  }
#ifdef SANITY_CHECK_SDBKZ
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "SDBKZ: The reduction needs a positive definite matrix\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_SDBKZ
  MicrosecondTime time;
#endif
  // A block larger than the dimension is the whole basis.
  int m = beta < n ? beta : n;
  LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<Tint> H = lll.Pmat;
  // The scale fixed here never changes: a unimodular congruence leaves the
  // content of an integral matrix unchanged, so the potentials of successive
  // rounds are compared on one scale.
  MyMatrix<Tring> gram =
      UniversalMatrixConversion<Tring, T>(RemoveFractionMatrix(lll.GramMatRed));
  Tring const num(delta_num);
  Tring const den(delta_den);
  Tring pot = SDBKZ_Potential(gram);
  [[maybe_unused]] size_t n_forward = 0;
  [[maybe_unused]] size_t n_backward = 0;
  [[maybe_unused]] bool stopped_by_potential = false;
  [[maybe_unused]] size_t n_final = 0;
  int n_tour = 0;
  auto after_insertion = [&](int const &j) -> void {
    SDBKZ_LLLBlock<T, Tring, Tint>(gram, H, j, m, os);
    IntegralSizeReduce(gram, H);
  };
  // One forward tour, with the blocks of the current basis; returns the
  // number of insertions.
  auto forward_tour = [&](std::vector<MyMatrix<Tring>> &l_block) -> size_t {
    size_t n_ins = 0;
    for (int j = 0; j + m <= n; j++) {
      std::optional<MyVector<Tint>> opt_z =
          BKZ_ShorterInBlock<T, Tring, Tint>(l_block[j], num, den, os);
      if (opt_z) {
        MyMatrix<Tint> U_blk = BKZ_UnimodularWithFirstRow(*opt_z);
        BKZ_ApplyBlockTransformation(gram, H, U_blk, j);
        after_insertion(j);
        l_block = SDBKZ_ProjectedBlocks(gram, m);
        n_ins++;
      }
    }
    return n_ins;
  };
  auto backward_tour = [&](std::vector<MyMatrix<Tring>> &l_block) -> size_t {
    size_t n_ins = 0;
    for (int j = n - m; j >= 0; j--) {
      // The content is divided out before the adjugate is taken, for the
      // reason given in SlideDualStep.
      MyMatrix<Tring> blk = RemoveFractionMatrix(l_block[j]);
      if (SlideDualStepOnBlock<T, Tring, Tint>(blk, gram, H, j, m, num, den,
                                               os)) {
        after_insertion(j);
        l_block = SDBKZ_ProjectedBlocks(gram, m);
        n_ins++;
      }
    }
    return n_ins;
  };
  // Phase 1: the self-dual rounds, kept while they decrease the potential.
  while (true) {
    std::vector<MyMatrix<Tring>> l_block = SDBKZ_ProjectedBlocks(gram, m);
    size_t n_fwd = forward_tour(l_block);
    size_t n_bwd = backward_tour(l_block);
    n_forward += n_fwd;
    n_backward += n_bwd;
    n_tour++;
    if (n_fwd + n_bwd == 0) {
      break;
    }
    Tring pot_new = SDBKZ_Potential(gram);
    if (pot_new >= pot) {
      stopped_by_potential = true;
      break;
    }
    pot = pot_new;
    if (max_tour > 0 && n_tour >= max_tour) {
      break;
    }
  }
  // Phase 2: forward tours to a fixed point, which certifies the block
  // condition of every full block. This terminates as BKZ does: an insertion
  // at j fixes D_1, ..., D_{j-1} and decreases D_j strictly, and neither the
  // LLL reduction of the block, which keeps the shortest vector it starts
  // with, nor the size reduction increases any D_i.
  while (true) {
    std::vector<MyMatrix<Tring>> l_block = SDBKZ_ProjectedBlocks(gram, m);
    size_t n_fwd = forward_tour(l_block);
    n_final += n_fwd;
    if (n_fwd == 0) {
      break;
    }
  }
#ifdef DEBUG_SDBKZ
  os << "SDBKZ: n=" << n << " beta=" << beta << " rounds=" << n_tour
     << " forward insertions=" << n_forward
     << " backward insertions=" << n_backward
     << " stopped by the potential=" << stopped_by_potential
     << " final forward insertions=" << n_final << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(H);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(H)};
#ifdef SANITY_CHECK_SDBKZ
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_SDBKZ
  os << "SDBKZ: SelfDualBKZReducedBasisDelta beta=" << beta << " took " << time
     << "\n";
#endif
  return res;
}

/*
  The same delta as BKZ and slide reduction in this package, and no cap on the
  number of rounds.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint> SelfDualBKZReducedBasis(MyMatrix<T> const &GramMat,
                                              int const &beta,
                                              std::ostream &os) {
  return SelfDualBKZReducedBasisDelta<T, Tint>(GramMat, beta, 99, 100, 0, os);
}

// clang-format off
#endif  // SRC_LATT_SELFDUALBKZ_H_
// clang-format on
