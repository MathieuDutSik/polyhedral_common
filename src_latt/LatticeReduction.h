// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEREDUCTION_H_
#define SRC_LATT_LATTICEREDUCTION_H_
// clang-format off
#include "BKZ.h"
#include "Basic_string.h"
#include "ClassicLLL.h"
#include "DeepLLL.h"
#include "MAT_Matrix.h"
#include "SeysenReduction.h"
#include "SlideReduction.h"
#include <optional>
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
/*
  The three block-like reductions take a parameter and are named "deep-<d>",
  "bkz-<b>" and "slide-<k>" for any admissible value: the parameter is parsed
  from the name rather than being one of a fixed handful of spellings. Bare
  "deep" is deep insertion without restriction on the depth, which is the
  algorithm as Schnorr and Euchner state it; "deep-<d>" restricts it.

  The list below is not the set of accepted methods, which is infinite. It is
  the set of candidates that "best" tries, chosen to span the useful range
  without costing more than it is worth. A caller wanting a value outside it
  names the method.
 */
inline std::vector<std::string> LatticeReductionSingleMethods() {
  return {"direct", "dual",    "seysen",  "seysen_best", "seysen_lll",
          "deep",   "deep-5",  "deep-10", "bkz-4",       "bkz-8",
          "bkz-12", "slide-4", "slide-8"};
}

/*
  A strictly parsed non-negative integer: digits only, non-empty, and short
  enough that the conversion cannot overflow. Anything else yields nothing,
  and the caller reports the whole name as unknown rather than silently
  accepting a parameter it did not mean; std::stoi on its own would take
  "8junk" for 8.
 */
inline std::optional<int> LatticeReduction_ParseInteger(std::string const &str) {
  if (str.empty() || str.size() > 9) {
    return {};
  }
  for (auto &c : str) {
    if (c < '0' || c > '9') {
      return {};
    }
  }
  return StringToInt(str);
}

/*
  The parameter of a method named "<prefix>-<value>", or nothing if the name
  does not have that shape.
 */
inline std::optional<int>
LatticeReduction_PrefixedParameter(std::string const &method,
                                   std::string const &prefix) {
  std::optional<std::string> tail = get_postfix(method, prefix + "-");
  if (!tail) {
    return {};
  }
  return LatticeReduction_ParseInteger(*tail);
}

/*
  The parameter must be at least min_value. This is a check on USER INPUT and
  so is unconditional: it is not a programming invariant and must not be
  compiled out with the sanity checks.
 */
inline void LatticeReduction_CheckParameter(std::string const &method,
                                            int const &value,
                                            int const &min_value) {
  if (value < min_value) {
    std::cerr << "LATTICE_REDUCTION: the method " << method
              << " asks for a parameter of " << value
              << ", but the smallest meaningful value is " << min_value
              << "\n";
    throw TerminalException{1};
  }
}

/*
  The message for a name that matched nothing. The old spellings without the
  hyphen were in use, so a name of that shape gets told what it should be
  rather than just being rejected.
 */
inline void LatticeReduction_UnknownMethod(std::string const &method) {
  std::cerr << "LATTICE_REDUCTION: unknown method " << method << "\n";
  for (auto &prefix : {"deep", "bkz", "slide"}) {
    std::optional<std::string> tail = get_postfix(method, prefix);
    if (tail && LatticeReduction_ParseInteger(*tail)) {
      std::cerr << "  did you mean " << prefix << "-" << *tail
                << "? The parameter is separated by a hyphen.\n";
    }
  }
  std::cerr << "  the methods without a parameter are:";
  for (auto &meth : {"direct", "dual", "seysen", "seysen_best", "seysen_lll",
                     "deep", "best"}) {
    std::cerr << " " << meth;
  }
  std::cerr << "\n";
  std::cerr << "  the methods with one are: deep-<depth>, bkz-<blocksize>, "
               "slide-<blocksize>\n";
  throw TerminalException{1};
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
  // A depth of zero would mean no restriction, which is what bare "deep"
  // already provides, so the smallest meaningful restriction is one.
  std::optional<int> depth = LatticeReduction_PrefixedParameter(method, "deep");
  if (depth) {
    LatticeReduction_CheckParameter(method, *depth, 1);
    return DeepLLLreducedBasisDepth<T, Tint>(GramMat, *depth, os);
  }
  // A block of size one carries no condition; BKZ at two is LLL.
  std::optional<int> beta = LatticeReduction_PrefixedParameter(method, "bkz");
  if (beta) {
    LatticeReduction_CheckParameter(method, *beta, 2);
    return BKZreducedBasis<T, Tint>(GramMat, *beta, os);
  }
  // The value is an upper bound: slide reduction needs the block size to
  // divide the dimension, and SlideBlockSize takes the largest divisor that
  // does not exceed it.
  std::optional<int> k_max =
      LatticeReduction_PrefixedParameter(method, "slide");
  if (k_max) {
    LatticeReduction_CheckParameter(method, *k_max, 2);
    return SlideReducedBasisAuto<T, Tint>(GramMat, *k_max, os);
  }
  LatticeReduction_UnknownMethod(method);
  // Not reached; LatticeReduction_UnknownMethod always throws.
  return LLLnoreduction<T, Tint>(GramMat);
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
