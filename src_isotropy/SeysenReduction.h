// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ISOTROPY_SEYSENREDUCTION_H_
#define SRC_ISOTROPY_SEYSENREDUCTION_H_
// clang-format off
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixInverse.h"
#include "QuoIntFcts.h"
#include <string>
#include <utility>
// clang-format on

#ifdef DEBUG
#define DEBUG_SEYSEN_REDUCTION
#endif

#ifdef DISABLE_DEBUG_SEYSEN_REDUCTION
#undef DEBUG_SEYSEN_REDUCTION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SEYSEN_REDUCTION
#endif

#ifdef TIMINGS
#define TIMINGS_SEYSEN_REDUCTION
#endif

/*
  Seysen reduction.

  Seysen's measure of a basis is

      S(B) = sum_i |b_i|^2 |b^_i|^2,

  where b^_1, ..., b^_n is the reciprocal basis, characterised by
  <b_i, b^_j> = delta_ij. Beware the notation: in the LLL literature b_i^*
  is the Gram-Schmidt orthogonalisation, whereas here the hatted vector is
  the reciprocal, a different object.

  On the Gram side, with A the Gram matrix of B, the Gram matrix of the
  reciprocal basis is A^-1, so the measure is a function of A alone,

      S(A) = sum_i A_ii (A^-1)_ii.

  This has four properties that together make it the natural global objective
  for reduction:

    * S(A) >= n by Cauchy-Schwarz, with equality exactly for an orthogonal
      basis, so it measures what one wants and not a proxy for it;
    * S(lambda A) = S(A), so it is scale free;
    * it depends on B only through A, so the orthogonal factor of the Iwasawa
      decomposition B = KAN never enters and never has to be recovered;
    * S(A) = S(A^-1), so a basis and its reciprocal are reduced by one and the
      same descent. This is the point of Seysen's paper and the property that
      distinguishes the method from everything in the LLL family, where a
      primal pass and a dual pass are different computations.

  The moves are the transvections b_i <- b_i + lambda b_j. Such a move changes
  exactly two of the 2n vectors involved and changes them in opposite
  directions: b_i in the primal, and, since the reciprocal basis transforms
  contragradiently, b^_j <- b^_j - lambda b^_i in the dual. The measure prices
  that trade, and the optimum of the trade is a rounded quotient: with
  A* = A^-1,

      lambda_opt = (A_jj A*_ij - A_ij A*_ii) / (2 A_jj A*_ii),

  and the step is the nearest integer to it. Note what this depends on. The
  size reduction of LLL reads its coefficient off the primal Gram-Schmidt data
  alone and always accepts the step; here the coefficient involves A*_ij and
  A*_ii as well, that is the geometry of the reciprocal lattice, and the step
  is taken only if it lowers the measure. Two bases with identical mu's and
  identical Gram-Schmidt profiles can call for different Seysen steps.

  Everything below runs over the integers. The rational A^-1 is replaced by
  the adjugate, adj(A) = det(A) A^-1, which is integral for an integral A, and
  det(A) cancels from the formula for lambda. The quantity descended on is

      Shat(A) = sum_i A_ii adj(A)_ii = det(A) S(A),

  a positive integer, bounded below by n det(A), with det(A) fixed along the
  whole descent since the moves are unimodular. Each accepted move strictly
  decreases it, so the descent terminates, by the same argument that makes LLL
  terminate. Note that the analogous statement fails for the logarithmic
  potentials one is otherwise tempted to write down: those are real valued,
  have no integrality to fall back on, and cannot even be evaluated in exact
  arithmetic.

  What this does not do is escape a local minimum. Seysen reduction is a
  greedy descent over transvections; it improves the information used to rank
  the neighbourhood, not the size of the neighbourhood. It is also weaker than
  LLL at producing one short vector, which is the natural consequence of a
  measure treating all n indices symmetrically, and correspondingly better at
  keeping the basis balanced and the coefficients small. The two are
  complements, and the useful combination is to alternate them.

  Reference: M. Seysen, Simultaneous reduction of a lattice basis and its
  reciprocal basis, Combinatorica 13 (1993) 363--376.
 */

/*
  Seysen's measure S(A) = sum_i A_ii (A^-1)_ii of a positive definite form,
  as an exact element of T. It is at least the dimension, with equality
  exactly for an orthogonal basis.
 */
template <typename T> T SeysenMeasure(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
  std::pair<MyMatrix<T>, T> pair = AdjugateDeterminant(GramMat);
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  if (pair.second <= 0) {
    std::cerr << "SEYSEN: SeysenMeasure expects a positive definite form, but "
                 "the determinant is "
              << pair.second << "\n";
    throw TerminalException{1};
  }
#endif
  T num(0);
  for (int i = 0; i < n; i++) {
    num += GramMat(i, i) * pair.first(i, i);
  }
  return num / pair.second;
}

/*
  The integral potential Shat(A) = sum_i A_ii adj(A)_ii = det(A) S(A), from a
  form and its adjugate that the caller already has. This is the quantity the
  descent decreases, and it is what a caller should monitor rather than S: it
  is an integer, so a comparison of two states is exact and cheap.
 */
template <typename T>
T SeysenIntegralPotential(MyMatrix<T> const &A, MyMatrix<T> const &adj) {
  int n = A.rows();
  T pot(0);
  for (int i = 0; i < n; i++) {
    pot += A(i, i) * adj(i, i);
  }
  return pot;
}

/*
  The descent itself.

  best_move selects between the two sweep policies. With best_move false a
  sweep applies every improving pair as it is met, which is what one wants in
  practice: a sweep then does real work throughout rather than spending O(n^2)
  evaluations to perform a single move. With best_move true a sweep evaluates
  all ordered pairs and applies only the most profitable one, which is slower
  but makes the trajectory independent of the order in which pairs are
  visited, and is therefore the variant to use when comparing implementations
  or investigating a descent that stalls.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint> SeysenReducedBasisKernel(MyMatrix<T> const &GramMat,
                                               bool const &best_move,
                                               [[maybe_unused]]
                                               std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = GramMat.rows();
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  if (n != GramMat.cols()) {
    std::cerr << "SEYSEN: The matrix should be square\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(GramMat)) {
    std::cerr << "SEYSEN: The Gram matrix should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 1) {
    return {GramMat, IdentityMat<Tint>(n)};
  }
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "SEYSEN: The Seysen reduction needs a positive definite "
                 "matrix\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_SEYSEN_REDUCTION
  MicrosecondTime time;
#endif
  // The measure is scale free, so the descent may be run on an integral
  // rescaling of the form and the transformation found is the same. This is
  // what puts the whole computation over the integers, and with it the
  // termination argument: the potential is then a strictly decreasing
  // positive integer.
  MyMatrix<Tring> A =
      UniversalMatrixConversion<Tring, T>(RemoveFractionMatrix(GramMat));
  std::pair<MyMatrix<Tring>, Tring> pair = AdjugateDeterminant(A);
  MyMatrix<Tring> adj = std::move(pair.first);
  MyMatrix<Tint> P = IdentityMat<Tint>(n);
  Tring const two(2);
#ifdef DEBUG_SEYSEN_REDUCTION
  os << "SEYSEN: n=" << n << " det=" << pair.second
     << " initial potential=" << SeysenIntegralPotential(A, adj) << "\n";
#endif
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  Tring pot_prev = SeysenIntegralPotential(A, adj);
#endif
  //
  // The gain of the move (i, j, lambda) and the optimal lambda. Writing
  //   num = A_jj adj_ij - A_ij adj_ii,   den = 2 A_jj adj_ii > 0,
  // the real optimum is num/den and the gain at an integer lambda is
  //   lambda (lambda den - 2 num),
  // which is the expansion of Seysen's formula with det(A) cancelled. The
  // rounding is the floor of num/den + 1/2; at a tie the two candidates have
  // the same gain, the parabola being symmetric about its vertex, so the
  // choice between them does not matter.
  //
  auto f_candidate = [&](int const &i, int const &j) -> std::pair<Tring, Tring> {
    Tring den = two * A(j, j) * adj(i, i);
    Tring num = A(j, j) * adj(i, j) - A(i, j) * adj(i, i);
    // The arguments are materialised rather than passed as expressions: the
    // multiprecision types use expression templates, and QuoInt deduces its
    // single parameter from both arguments.
    Tring quo_num = two * num + den;
    Tring quo_den = two * den;
    Tring lambda = QuoInt(quo_num, quo_den);
    if (lambda == 0) {
      return {Tring(0), Tring(0)};
    }
    Tring gain = lambda * (lambda * den - two * num);
    return {lambda, gain};
  };
  //
  // Applying a move. Only row and column i of A change, and only row and
  // column j of adj, by Seysen's contragradient action; both updates are
  // O(n) and the diagonal entries must be formed from the old values, hence
  // the two saved quantities.
  //
  auto f_apply = [&](int const &i, int const &j, Tring const &lambda) -> void {
    Tring new_Aii =
        A(i, i) + two * lambda * A(i, j) + lambda * lambda * A(j, j);
    Tring new_adjjj =
        adj(j, j) - two * lambda * adj(i, j) + lambda * lambda * adj(i, i);
    for (int k = 0; k < n; k++) {
      if (k != i) {
        A(i, k) += lambda * A(j, k);
        A(k, i) += lambda * A(k, j);
      }
      if (k != j) {
        adj(j, k) -= lambda * adj(i, k);
        adj(k, j) -= lambda * adj(k, i);
      }
    }
    A(i, i) = new_Aii;
    adj(j, j) = new_adjjj;
    Tint lambda_int = UniversalScalarConversion<Tint, Tring>(lambda);
    RowAddMul(P, i, lambda_int, j);
  };
  //
  size_t n_move = 0;
  [[maybe_unused]] size_t n_sweep = 0;
  while (true) {
    bool did_move = false;
    if (best_move) {
      Tring best_gain(0);
      int best_i = -1, best_j = -1;
      Tring best_lambda(0);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i != j) {
            std::pair<Tring, Tring> cand = f_candidate(i, j);
            if (cand.second < best_gain) {
              best_gain = cand.second;
              best_lambda = cand.first;
              best_i = i;
              best_j = j;
            }
          }
        }
      }
      if (best_i >= 0) {
        f_apply(best_i, best_j, best_lambda);
        n_move++;
        did_move = true;
      }
    } else {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i != j) {
            std::pair<Tring, Tring> cand = f_candidate(i, j);
            if (cand.second < 0) {
              f_apply(i, j, cand.first);
              n_move++;
              did_move = true;
            }
          }
        }
      }
    }
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
    Tring pot_now = SeysenIntegralPotential(A, adj);
    if (did_move && pot_now >= pot_prev) {
      std::cerr << "SEYSEN: a sweep that moved did not decrease the potential, "
                   "it went from "
                << pot_prev << " to " << pot_now << "\n";
      throw TerminalException{1};
    }
    pot_prev = pot_now;
#endif
    if (!did_move) {
      break;
    }
#ifdef DEBUG_SEYSEN_REDUCTION
    n_sweep++;
    os << "SEYSEN: after sweep " << n_sweep << ", moves=" << n_move
       << " potential=" << SeysenIntegralPotential(A, adj) << "\n";
#endif
  }
#ifdef DEBUG_SEYSEN_REDUCTION
  os << "SEYSEN: terminated after " << n_sweep << " sweeps and " << n_move
     << " moves, final potential=" << SeysenIntegralPotential(A, adj) << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(P)};
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_SEYSEN_REDUCTION
  os << "SEYSEN: SeysenReducedBasisKernel took " << time << "\n";
#endif
  return res;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> SeysenReducedBasis(MyMatrix<T> const &GramMat,
                                         std::ostream &os) {
  return SeysenReducedBasisKernel<T, Tint>(GramMat, false, os);
}

template <typename T, typename Tint>
LLLreduction<T, Tint> SeysenReducedBasisBest(MyMatrix<T> const &GramMat,
                                             std::ostream &os) {
  return SeysenReducedBasisKernel<T, Tint>(GramMat, true, os);
}

/*
  Seysen and LLL are complements rather than competitors: the first keeps the
  basis balanced and the coefficients small, the second is the one that
  produces a short leading vector. Alternating them until neither improves the
  Seysen potential gives a basis that is good in both senses, at a cost which
  is a small multiple of either.
 */
template <typename T, typename Tint>
LLLreduction<T, Tint> SeysenLLLreducedBasis(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  MyMatrix<T> work = GramMat;
  MyMatrix<Tint> P = IdentityMat<Tint>(GramMat.rows());
  // The comparison is on the scale free measure and not on the integral
  // potential of the kernel. The potential is relative to whatever rescaling
  // made the form integral, and a congruence can lower the denominators of
  // the entries, so two rounds need not measure their potentials against the
  // same unit. The measure has no such defect, and it still gives
  // termination: its denominator divides det(G), which is invariant under
  // the unimodular congruences applied here, so the values it takes lie in a
  // fixed (1/det)Z and are bounded below by the dimension. A strictly
  // decreasing sequence in such a set is finite.
  T meas = SeysenMeasure(work);
  while (true) {
    LLLreduction<T, Tint> seysen = SeysenReducedBasis<T, Tint>(work, os);
    LLLreduction<T, Tint> lll =
        LLLreducedBasis<T, Tint>(seysen.GramMatRed, os);
    MyMatrix<Tint> step = lll.Pmat * seysen.Pmat;
    T meas_new = SeysenMeasure(lll.GramMatRed);
    if (meas_new >= meas) {
      // The LLL pass undid more than the Seysen pass gained, so the pair of
      // them is not an improvement and the previous state is returned.
      break;
    }
    meas = meas_new;
    work = lll.GramMatRed;
    P = step * P;
  }
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  LLLreduction<T, Tint> res = {P_T * GramMat * P_T.transpose(), std::move(P)};
#ifdef SANITY_CHECK_SEYSEN_REDUCTION
  CheckLLLreduction(res, GramMat);
#endif
  return res;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> SeysenReducedGeneral(MyMatrix<T> const &GramMat,
                                           std::string const &method,
                                           std::ostream &os) {
  if (method == "first") {
    return SeysenReducedBasis<T, Tint>(GramMat, os);
  }
  if (method == "best") {
    return SeysenReducedBasisBest<T, Tint>(GramMat, os);
  }
  if (method == "seysen_lll") {
    return SeysenLLLreducedBasis<T, Tint>(GramMat, os);
  }
  std::cerr << "SEYSEN: No matching method for " << method
            << ", allowed are first, best, seysen_lll\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_ISOTROPY_SEYSENREDUCTION_H_
// clang-format on
