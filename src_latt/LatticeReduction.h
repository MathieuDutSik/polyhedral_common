// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEREDUCTION_H_
#define SRC_LATT_LATTICEREDUCTION_H_
// clang-format off
#include "BKZ.h"
#include "ClassicLLL.h"
#include "DeepLLL.h"
#include "MAT_Matrix.h"
#include "SeysenReduction.h"
#include "SlideReduction.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  One dispatch over every reduction of a positive definite Gram matrix that
  this package implements, so that a caller naming a method by string has a
  single list to consult and a new reduction has a single place to be added.

  All of them return the same LLLreduction pair (GramMatRed, Pmat) with
  Pmat GramMat Pmat^T = GramMatRed, so they are interchangeable at any call
  site.

  What distinguishes them is worth keeping in mind when choosing:

    * direct and dual are LLL, the second reducing the form through its
      adjugate. On every measure applied in this package the dual variant is
      the weaker of the two, sometimes by orders of magnitude.

    * seysen and its variants minimise sum_i |b_i|^2 |dual b_i|^2, which
      treats all indices symmetrically and reduces the form and its dual in
      one descent. This is the one to try when what matters is that EVERY
      coefficient be small rather than that one vector be short.

    * deep is Schnorr-Euchner deep insertion: LLL's move set widened from the
      adjacent swap to an insertion at any earlier position.

    * bkz and slide change the ORACLE rather than the move set, asking at each
      index for a shortest vector of a projected block of the given size.
      Slide reduction differs from BKZ in using two families of conditions on
      non-overlapping blocks, which buys a polynomial bound on the number of
      steps that BKZ has no analogue of.

  The block methods cost more than the rest, superexponentially in the block
  size, and the quality they buy improves slowly. Above block size eight the
  gains observed in this package were negligible.
 */
inline std::vector<std::string> LatticeReductionSingleMethods() {
  return {"direct", "dual",   "seysen", "seysen_best", "seysen_lll",
          "deep",   "deep5",  "deep10", "bkz4",        "bkz8",
          "bkz12",  "slide4", "slide8"};
}

template <typename T, typename Tint>
LLLreduction<T, Tint> LatticeReducedGeneral(MyMatrix<T> const &GramMat,
                                            std::string const &method,
                                            std::ostream &os) {
  if (method == "direct") {
    return LLLreducedBasis<T, Tint>(GramMat, os);
  }
  if (method == "dual") {
    return LLLreducedBasisDual<T, Tint>(GramMat, os);
  }
  if (method == "seysen") {
    return SeysenReducedBasis<T, Tint>(GramMat, os);
  }
  if (method == "seysen_best") {
    return SeysenReducedBasisBest<T, Tint>(GramMat, os);
  }
  if (method == "seysen_lll") {
    return SeysenLLLreducedBasis<T, Tint>(GramMat, os);
  }
  if (method == "deep") {
    return DeepLLLreducedBasis<T, Tint>(GramMat, os);
  }
  if (method == "deep5") {
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 5, os);
  }
  if (method == "deep10") {
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, 10, os);
  }
  if (method == "bkz4") {
    return BKZreducedBasis<T, Tint>(GramMat, 4, os);
  }
  if (method == "bkz8") {
    return BKZreducedBasis<T, Tint>(GramMat, 8, os);
  }
  if (method == "bkz12") {
    return BKZreducedBasis<T, Tint>(GramMat, 12, os);
  }
  if (method == "slide4") {
    return SlideReducedBasisAuto<T, Tint>(GramMat, 4, os);
  }
  if (method == "slide8") {
    return SlideReducedBasisAuto<T, Tint>(GramMat, 8, os);
  }
  std::cerr << "LATTICE_REDUCTION: unknown method " << method
            << ". Allowed are:";
  for (auto &meth : LatticeReductionSingleMethods()) {
    std::cerr << " " << meth;
  }
  std::cerr << " best\n";
  throw TerminalException{1};
}

/*
  The quality of a reduced Gram matrix, for the "best" search below.

  There is NO canonical measure here, and the choice made is worth stating
  rather than hiding. The primary criterion is the squared orthogonality
  defect prod_i G_ii / det G, which is scale free, is at least one by Hadamard
  and equals one exactly for an orthogonal basis; the sum of absolute values of
  the entries, which is the size of the integers everything downstream will
  compute with, breaks its ties.

  The order of the two matters, and taking it the other way round was a
  mistake made here first. Coefficient size alone is not a measure of
  reduction: a sparse but badly skewed form can have a smaller L1 norm than
  any reduction of it, so a search on L1 first happily returns the UNREDUCED
  input, which is not what a reduction is for. On a scrambled E8 the input had
  L1 norm 34 against 48 for every reduction of it, and defect 512 against 256.

  A caller who cares about something else -- a short first vector, a balanced
  Gram-Schmidt profile, Seysen's measure -- should name the method rather than
  ask for the best, because the best by this criterion need not be the best by
  that one.
 */
template <typename T> struct LatticeReductionQuality {
  T orth_defect_sq;
  T l1_norm;
  bool operator<(LatticeReductionQuality<T> const &other) const {
    if (orth_defect_sq != other.orth_defect_sq) {
      return orth_defect_sq < other.orth_defect_sq;
    }
    return l1_norm < other.l1_norm;
  }
};

template <typename T>
LatticeReductionQuality<T>
ComputeLatticeReductionQuality(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
  T det = DeterminantMat(GramMat);
  T prod_diag(1);
  for (int i = 0; i < n; i++) {
    prod_diag *= GramMat(i, i);
  }
  return {prod_diag / det, L1_norm_mat(GramMat)};
}

template <typename T, typename Tint> struct LatticeReductionResult {
  MyMatrix<T> GramMatRed;
  MyMatrix<Tint> Pmat;
  std::string method;
  LatticeReductionQuality<T> quality;
};

/*
  Run every method and keep whichever minimises the criterion above. The
  unreduced input is among the candidates, so the result is never worse than
  what was handed in.
 */
template <typename T, typename Tint>
LatticeReductionResult<T, Tint>
LatticeReducedBest(MyMatrix<T> const &GramMat, std::ostream &os) {
  LatticeReductionResult<T, Tint> best{
      GramMat, IdentityMat<Tint>(GramMat.rows()), "none",
      ComputeLatticeReductionQuality(GramMat)};
  for (auto &method : LatticeReductionSingleMethods()) {
    LLLreduction<T, Tint> red =
        LatticeReducedGeneral<T, Tint>(GramMat, method, os);
    LatticeReductionQuality<T> quality =
        ComputeLatticeReductionQuality(red.GramMatRed);
#ifdef DEBUG_LATTICE_REDUCTION
    os << "LATTICE_REDUCTION: method " << method << " L1=" << quality.l1_norm
       << " defect^2=" << quality.orth_defect_sq << "\n";
#endif
    if (quality < best.quality) {
      best = {std::move(red.GramMatRed), std::move(red.Pmat), method,
              std::move(quality)};
    }
  }
  return best;
}

template <typename T, typename Tint>
LatticeReductionResult<T, Tint>
LatticeReducedByName(MyMatrix<T> const &GramMat, std::string const &method,
                     std::ostream &os) {
  if (method == "best") {
    return LatticeReducedBest<T, Tint>(GramMat, os);
  }
  LLLreduction<T, Tint> red =
      LatticeReducedGeneral<T, Tint>(GramMat, method, os);
  LatticeReductionQuality<T> quality =
      ComputeLatticeReductionQuality(red.GramMatRed);
  return {std::move(red.GramMatRed), std::move(red.Pmat), method,
          std::move(quality)};
}

// clang-format off
#endif  // SRC_LATT_LATTICEREDUCTION_H_
// clang-format on
