// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_MINKOWSKIREDUCTION_H_
#define SRC_LATT_MINKOWSKIREDUCTION_H_
// clang-format off
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_MINKOWSKI
#endif

#ifdef DISABLE_DEBUG_MINKOWSKI
#undef DEBUG_MINKOWSKI
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MINKOWSKI
#endif

#ifdef TIMINGS
#define TIMINGS_MINKOWSKI
#endif

/*
  Minkowski reduction.

  A basis b_1, ..., b_n of L is Minkowski reduced when, for every i, b_i is a
  SHORTEST vector among those v in L for which (b_1, ..., b_{i-1}, v) extends
  to a basis of L. It is the strongest of the classical notions: where LLL,
  deep insertion and BKZ each ask for a local or blockwise optimality, this
  asks, at every index in turn, for the true optimum over the whole lattice
  subject only to the previous vectors being kept.

  THE CONDITION IS DECIDABLE IN ANY DIMENSION, and it is worth saying so
  because a restriction to dimension seven is often quoted in connection with
  Minkowski reduction. That restriction is a different statement: the explicit
  finite system of inequalities cutting out the Minkowski fundamental domain in
  the cone of positive definite forms is known only up to dimension seven
  (Tammela). Reducing a GIVEN form asks nothing of that domain. It asks only
  for a shortest admissible vector at each index, which is a finite
  computation, and the shortest-vector enumerator performs it.

  Extendability in coordinates. Writing v = sum_j c_j b_j in the current basis,
  (b_1, ..., b_{i-1}, v) extends to a basis of L exactly when

      gcd(c_i, c_{i+1}, ..., c_n) = 1.

  The quotient of L by the span of b_1, ..., b_{i-1} is free on the images of
  b_i, ..., b_n; the pair extends to a basis exactly when the image of v is
  primitive there, which is the stated condition. Note that it does not involve
  c_1, ..., c_{i-1} at all.

  One forward pass suffices. The set of admissible v at index i depends only on
  L and on b_1, ..., b_{i-1}. A later step changes b_j for j > i and nothing
  before it, so it cannot disturb a condition already established. There is no
  outer loop and no termination argument to make: the algorithm is n
  enumerations.

  THE COST IS EXPONENTIAL and this is inherent, not an artefact. Each step
  enumerates the lattice points of norm at most |b_i|^2, and the last index is
  the expensive one. The method belongs where the dimension is small and
  exactness matters more than speed, which is the regime of the perfect-form,
  Ctype and T-space computations this package exists for; it is not a
  replacement for LLL or BKZ at larger sizes and is deliberately left out of
  the candidate list that the "best" searches try.

  The basis is LLL reduced first. That is not required for correctness but
  decides the cost: the enumeration bound at index i is |b_i|^2, so starting
  from a skewed basis enumerates a ball larger than necessary at every step.
 */

/*
  Is the coefficient vector z admissible at index i, that is does
  gcd(z_i, ..., z_{n-1}) equal one? A vector supported entirely on the first i
  coordinates lies in the span of b_1, ..., b_{i-1} and has a gcd of zero
  there, so it is correctly rejected.
 */
template <typename Tint>
bool Minkowski_IsAdmissible(MyVector<Tint> const &z, int const &i) {
  int n = z.size();
  Tint gcd(0);
  for (int j = i; j < n; j++) {
    gcd = GcdPair(gcd, z(j));
  }
  return T_abs(gcd) == 1;
}

/*
  The unimodular matrix putting the admissible vector z at index i while
  leaving b_1, ..., b_{i-1} exactly as they are. Its rows above i are the unit
  vectors, its row i is z, and its rows below complete the primitive tail
  (z_i, ..., z_{n-1}) to a basis of Z^{n-i}. Being block lower triangular with
  an identity block above, its determinant is that of the tail block, which is
  plus or minus one.
 */
template <typename Tint>
MyMatrix<Tint> Minkowski_ReplacementMatrix(MyVector<Tint> const &z,
                                           int const &i) {
  int n = z.size();
  int m = n - i;
  MyVector<Tint> tail(m);
  for (int j = 0; j < m; j++) {
    tail(j) = z(i + j);
  }
#ifdef SANITY_CHECK_MINKOWSKI
  if (!IsVectorPrimitive(tail)) {
    std::cerr << "MINKOWSKI: the tail of the chosen vector is not primitive, "
                 "so it does not extend to a basis\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> U = ZeroMatrix<Tint>(n, n);
  for (int j = 0; j < i; j++) {
    U(j, j) = Tint(1);
  }
  for (int j = 0; j < n; j++) {
    U(i, j) = z(j);
  }
  if (m > 1) {
    MyMatrix<Tint> compl_mat = ComplementToBasis(tail);
    for (int r = 0; r + 1 < m; r++) {
      for (int c = 0; c < m; c++) {
        U(i + 1 + r, i + c) = compl_mat(r, c);
      }
    }
  }
#ifdef SANITY_CHECK_MINKOWSKI
  Tint det_U = DeterminantMat(U);
  if (det_U != 1 && det_U != -1) {
    std::cerr << "MINKOWSKI: the replacement matrix is not unimodular, det="
              << det_U << "\n";
    throw TerminalException{1};
  }
#endif
  return U;
}

template <typename T, typename Tint>
LLLreduction<T, Tint> MinkowskiReducedBasis(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  int n = GramMat.rows();
#ifdef SANITY_CHECK_MINKOWSKI
  if (n != GramMat.cols()) {
    std::cerr << "MINKOWSKI: The matrix should be square\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(GramMat)) {
    std::cerr << "MINKOWSKI: The Gram matrix should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  if (n <= 1) {
    return {GramMat, IdentityMat<Tint>(n)};
  }
#ifdef SANITY_CHECK_MINKOWSKI
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "MINKOWSKI: The Minkowski reduction needs a positive "
                 "definite matrix\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_MINKOWSKI
  MicrosecondTime time;
#endif
  LLLreduction<T, Tint> lll = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<T> gram = lll.GramMatRed;
  MyMatrix<Tint> H = lll.Pmat;
  [[maybe_unused]] size_t n_change = 0;
  for (int i = 0; i < n; i++) {
    T bound = gram(i, i);
    std::vector<MyVector<Tint>> cands =
        computeLevel_GramMat<T, Tint>(gram, bound, os);
    // The current vector is always a candidate and is admissible, its
    // coefficient vector being the unit vector at i; it is seeded here rather
    // than relied upon to come back from the enumeration.
    MyVector<Tint> best = ZeroVector<Tint>(n);
    best(i) = Tint(1);
    T best_norm = bound;
    for (auto &z : cands) {
      if (!Minkowski_IsAdmissible(z, i)) {
        continue;
      }
      T norm = EvaluationQuadForm<T, Tint>(gram, z);
      if (norm < best_norm) {
        best_norm = norm;
        best = z;
      }
    }
    if (best_norm < bound) {
      MyMatrix<Tint> U = Minkowski_ReplacementMatrix(best, i);
      MyMatrix<T> U_T = UniversalMatrixConversion<T, Tint>(U);
      gram = U_T * gram * U_T.transpose();
      H = U * H;
      n_change++;
#ifdef SANITY_CHECK_MINKOWSKI
      if (gram(i, i) != best_norm) {
        std::cerr << "MINKOWSKI: after the replacement the diagonal entry is "
                  << gram(i, i) << " and not the norm " << best_norm
                  << " that was chosen\n";
        throw TerminalException{1};
      }
#endif
    }
#ifdef DEBUG_MINKOWSKI
    os << "MINKOWSKI: index " << i << " candidates=" << cands.size()
       << " |b_i|^2 " << bound << " -> " << gram(i, i) << "\n";
#endif
  }
#ifdef DEBUG_MINKOWSKI
  os << "MINKOWSKI: n=" << n << " replacements=" << n_change << "\n";
#endif
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(H);
  MyMatrix<T> GramMatRed = P_T * GramMat * P_T.transpose();
  LLLreduction<T, Tint> res = {std::move(GramMatRed), std::move(H)};
#ifdef SANITY_CHECK_MINKOWSKI
  CheckLLLreduction(res, GramMat);
#endif
#ifdef TIMINGS_MINKOWSKI
  os << "MINKOWSKI: MinkowskiReducedBasis took " << time << "\n";
#endif
  return res;
}

/*
  Tests the condition directly: at every index, no admissible vector is
  strictly shorter than the one in place. The enumeration is redone from the
  matrix, so this checks the result rather than the descent that produced it.
 */
template <typename T, typename Tint>
bool IsMinkowskiReduced(MyMatrix<T> const &GramMat, std::ostream &os) {
  int n = GramMat.rows();
  if (n <= 1) {
    return true;
  }
  for (int i = 0; i < n; i++) {
    T bound = GramMat(i, i);
    std::vector<MyVector<Tint>> cands =
        computeLevel_GramMat<T, Tint>(GramMat, bound, os);
    for (auto &z : cands) {
      if (!Minkowski_IsAdmissible(z, i)) {
        continue;
      }
      T norm = EvaluationQuadForm<T, Tint>(GramMat, z);
      if (norm < bound) {
        os << "MINKOWSKI: at index " << i << " the vector in place has norm "
           << bound << " but an admissible vector of norm " << norm
           << " exists\n";
        return false;
      }
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_LATT_MINKOWSKIREDUCTION_H_
// clang-format on
