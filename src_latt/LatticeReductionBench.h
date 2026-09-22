// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEREDUCTIONBENCH_H_
#define SRC_LATT_LATTICEREDUCTIONBENCH_H_
// clang-format off
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixDeterminant.h"
#include "MAT_MatrixInt.h"
#include "QuoIntFcts.h"
#include "Shvec_exact.h"
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_REDUCTION_BENCH
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_REDUCTION_BENCH
#endif

/*
  The hidden-good-basis benchmark.

  A reduction algorithm is judged here by the following experiment. A lattice
  is presented by a Gram matrix G_M of a basis known to be good (a root
  lattice, Z^n, or a random well-rounded form). A random unimodular U is drawn
  and the presentation is destroyed:

      G_B = U G_M U^T.

  The algorithm receives G_B alone and returns P with G_red = P G_B P^T. Two
  different questions are then asked of it, and they must not be confused:

    * recovery: is P U a signed permutation matrix, i.e. was the original
      presentation found back;
    * reduction: is G_red a geometrically good presentation, judged by the
      measures below and compared against G_M itself.

  The second is the real lattice-reduction problem. Exact recovery of U is a
  strictly stronger demand and is not geometrically meaningful: any signed
  permutation of a good basis is just as good, and so is any other basis of
  the same quality. An algorithm can fail recovery completely and still be the
  better reducer.

  Everything is computed on the Gram side and in exact arithmetic. The
  orthogonal factor of the Iwasawa decomposition B = KAN never appears: it is
  killed by G = B B^T, and the remaining factors are read off the Bareiss
  decomposition that FullGramInfo already computes, through
  a_i^2 = d(i)/d(i-1) and mu_{i,j} = Nmat(i,j)/d(i).
 */

/*
  The quality measures.

  orth_defect_sq = prod_i G_ii / det G is the squared orthogonality defect. It
  is >= 1 by Hadamard, with equality exactly for an orthogonal basis, and it is
  the least arbitrary scale-free measure of basis quality.

  lll_potential = prod_{i<n} D_i, the product of the proper leading principal
  minors, is the classical LLL potential. It is a positive integer for an
  integral form and it decreases strictly under every profitable swap or
  insertion, which is what gives LLL its termination proof. It is the measure
  to use for descent; the logarithmic potentials are for reading, not for
  deciding.

  hermite_pow_n = G_11^n / det G is the n-th power of the squared Hermite
  factor, kept in that form so that it stays in the ring.

  max_coeff and total_bits measure the size of the integers in the reduced
  form. In a package doing exact arithmetic downstream this is not a cosmetic
  quantity: it is the cost of everything that reads the matrix afterwards.

  log_E_A and log_E_N are the two halves of the KAN energy, the scale imbalance
  sum (log a_i - (1/n) log det)^2 and the shear sum of mu_{i,j}^2. They are
  reported separately on purpose: a transformation routinely improves one while
  making the other worse, so a scalar objective built from them hides exactly
  the information one wants when a descent stalls.
 */
template <typename T> struct ReductionQuality {
  T orth_defect_sq;
  T lll_potential;
  T hermite_pow_n;
  T max_coeff;
  size_t total_bits;
  double log_orth_defect;
  double log_E_A;
  double log_E_N;
};

/*
  The number of bits of an entry. It is counted in the underlying ring rather
  than in T, since a repeated division by 2 in a field never terminates; the
  forms handled here are integral, so the conversion is exact.
 */
template <typename T> size_t BitSizeOfEntry(T const &val) {
  using Tring = typename underlying_ring<T>::ring_type;
  Tring a = T_abs(UniversalScalarConversion<Tring, T>(val));
  Tring two(2);
  Tring zero(0);
  size_t n_bit = 0;
  while (a > zero) {
    a = QuoInt(a, two);
    n_bit++;
  }
  return n_bit;
}

template <typename T>
ReductionQuality<T> ComputeReductionQuality(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
#ifdef SANITY_CHECK_REDUCTION_BENCH
  if (n != GramMat.cols()) {
    std::cerr << "REDUCTION_BENCH: the Gram matrix must be square\n";
    throw TerminalException{1};
  }
#endif
  // Checked unconditionally: the Bareiss decomposition below divides by the
  // leading principal minors, so a singular form is a division by zero and
  // the failure would otherwise be a signal rather than a message.
  if (DeterminantMat(GramMat) == 0) {
    std::cerr << "REDUCTION_BENCH: the Gram matrix is singular, the quality "
                 "measures are not defined for it\n";
    throw TerminalException{1};
  }
  FullGramInfo<T> info(GramMat);
  // d(i) is the leading principal minor of order i+1, so det G = d(n-1) and
  // a_{i+1}^2 = d(i)/d(i-1) with the convention d(-1) = 1.
  T det = info.d(n - 1);
  T prod_diag(1);
  for (int i = 0; i < n; i++) {
    prod_diag *= GramMat(i, i);
  }
  T lll_potential(1);
  for (int i = 0; i < n - 1; i++) {
    lll_potential *= info.d(i);
  }
  // The shortest diagonal entry is an upper bound for lambda_1^2 that does not
  // depend on the ordering of the basis, which the first entry would.
  T min_diag = GramMat(0, 0);
  for (int i = 1; i < n; i++) {
    if (GramMat(i, i) < min_diag) {
      min_diag = GramMat(i, i);
    }
  }
  T hermite_pow_n(1);
  for (int i = 0; i < n; i++) {
    hermite_pow_n *= min_diag;
  }
  hermite_pow_n /= det;
  T max_coeff(0);
  size_t total_bits = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T abs_val = T_abs(GramMat(i, j));
      if (abs_val > max_coeff) {
        max_coeff = abs_val;
      }
      total_bits += BitSizeOfEntry(GramMat(i, j));
    }
  }
  // The two KAN energies. These need logarithms, so they leave exact
  // arithmetic; they are diagnostics and nothing decides on them.
  double log_det = std::log(UniversalScalarConversion<double, T>(det));
  double log_E_A = 0;
  T prev(1);
  for (int i = 0; i < n; i++) {
    T a_sq = info.d(i) / prev;
    double log_a = 0.5 * std::log(UniversalScalarConversion<double, T>(a_sq));
    double dev = log_a - log_det / (2 * n);
    log_E_A += dev * dev;
    prev = info.d(i);
  }
  double log_E_N = 0;
  for (int i = 0; i < n; i++) {
    double d_i = UniversalScalarConversion<double, T>(info.d(i));
    for (int j = i + 1; j < n; j++) {
      double mu = UniversalScalarConversion<double, T>(info.Nmat(i, j)) / d_i;
      log_E_N += mu * mu;
    }
  }
  double log_orth_defect =
      0.5 * (std::log(UniversalScalarConversion<double, T>(prod_diag)) -
             log_det);
  T orth_defect_sq = prod_diag / det;
  return {std::move(orth_defect_sq),
          std::move(lll_potential),
          std::move(hermite_pow_n),
          std::move(max_coeff),
          total_bits,
          log_orth_defect,
          log_E_A,
          log_E_N};
}

/*
  The known-good presentations. The root lattices are the honest test cases:
  their standard Gram matrix is the reduced one, it is known independently of
  any reduction algorithm, and their automorphism groups are large enough that
  recovery up to signed permutation is a genuinely weak demand.
 */
template <typename T> MyMatrix<T> GoodGramZn(int n) {
  return IdentityMat<T>(n);
}

template <typename T> MyMatrix<T> GoodGramAn(int n) {
  MyMatrix<T> G = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n; i++) {
    G(i, i) = 2;
    if (i + 1 < n) {
      G(i, i + 1) = -1;
      G(i + 1, i) = -1;
    }
  }
  return G;
}

template <typename T> MyMatrix<T> GoodGramDn(int n) {
#ifdef SANITY_CHECK_REDUCTION_BENCH
  if (n < 4) {
    std::cerr << "REDUCTION_BENCH: D_n is defined here for n >= 4, not " << n
              << "\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> G = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n; i++) {
    G(i, i) = 2;
  }
  // The chain 1 - 2 - ... - (n-1) with the fork n attached to n-2.
  for (int i = 0; i + 1 < n - 1; i++) {
    G(i, i + 1) = -1;
    G(i + 1, i) = -1;
  }
  G(n - 3, n - 1) = -1;
  G(n - 1, n - 3) = -1;
  return G;
}

template <typename T> MyMatrix<T> GoodGramE8() {
  MyMatrix<T> G = ZeroMatrix<T>(8, 8);
  for (int i = 0; i < 8; i++) {
    G(i, i) = 2;
  }
  // Bourbaki E8: chain 1-3-4-5-6-7-8 with node 2 attached to node 4.
  int chain[7] = {0, 2, 3, 4, 5, 6, 7};
  for (int i = 0; i + 1 < 7; i++) {
    G(chain[i], chain[i + 1]) = -1;
    G(chain[i + 1], chain[i]) = -1;
  }
  G(1, 3) = -1;
  G(3, 1) = -1;
  return G;
}

/*
  A random form of low symmetry, built as G = R R^T with R integral and close
  to the identity, then reduced. It is there so that the benchmark is not run
  exclusively on lattices with enormous automorphism groups, where a reducer
  can look good by accident.

  Two caveats, since this family is weaker evidence than the root lattices.
  First, a perturbation of the identity is easily singular, so a draw with zero
  determinant is rejected and redrawn. Second, the "good" presentation here is
  only the LLL-reduced one, not an independently known optimum: on this family
  the benchmark asks whether an algorithm gets back to where LLL would be, not
  whether it reaches the best basis. Conclusions about beating LLL must come
  from the root lattices, whose reduced Gram matrix is known beforehand.
 */
template <typename T>
MyMatrix<T> GoodGramRandomWellRounded(int n, int spread,
                                      std::mt19937_64 &rng) {
  using Tint = typename underlying_ring<T>::ring_type;
  std::uniform_int_distribution<int> distr(-spread, spread);
  int const max_draw = 1000;
  for (int i_draw = 0; i_draw < max_draw; i_draw++) {
    MyMatrix<T> R = IdentityMat<T>(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j) {
          R(i, j) = distr(rng);
        }
      }
    }
    if (DeterminantMat(R) == 0) {
      continue;
    }
    MyMatrix<T> G = R * R.transpose();
    return LLLreducedBasis<T, Tint>(G, std::cerr).GramMatRed;
  }
  std::cerr << "REDUCTION_BENCH: no nonsingular perturbation of the identity "
               "found in "
            << max_draw << " draws at n=" << n << ", spread=" << spread
            << "\n";
  throw TerminalException{1};
}

template <typename T>
MyMatrix<T> GoodGramByName(std::string const &name, int n,
                           std::mt19937_64 &rng) {
  if (name == "Zn") {
    return GoodGramZn<T>(n);
  }
  if (name == "An") {
    return GoodGramAn<T>(n);
  }
  if (name == "Dn") {
    return GoodGramDn<T>(n);
  }
  if (name == "E8") {
    return GoodGramE8<T>();
  }
  if (name == "random") {
    return GoodGramRandomWellRounded<T>(n, 1, rng);
  }
  std::cerr << "REDUCTION_BENCH: unknown lattice name " << name
            << ", allowed are Zn, An, Dn, E8, random\n";
  throw TerminalException{1};
}

/*
  A random element of GL_n(Z), as a product of elementary generators. The
  entries grow roughly geometrically in the number of transvections, so n_ops
  is the knob that says how far the hidden presentation has been pushed away,
  and it is the parameter the benchmark sweeps: a reducer that copes at
  n_ops = 20 and collapses at n_ops = 200 has been found out.
 */
template <typename Tint>
MyMatrix<Tint> RandomUnimodular(int n, int n_ops, std::mt19937_64 &rng) {
  MyMatrix<Tint> U = IdentityMat<Tint>(n);
  std::uniform_int_distribution<int> distr_idx(0, n - 1);
  std::uniform_int_distribution<int> distr_coef(-2, 2);
  std::uniform_int_distribution<int> distr_kind(0, 9);
  for (int i_op = 0; i_op < n_ops; i_op++) {
    int kind = distr_kind(rng);
    if (kind < 8) {
      // A transvection, the generator that actually makes the basis bad.
      int i = distr_idx(rng);
      int j = distr_idx(rng);
      if (i == j) {
        continue;
      }
      Tint coef = distr_coef(rng);
      if (coef == 0) {
        coef = 1;
      }
      for (int k = 0; k < n; k++) {
        U(i, k) += coef * U(j, k);
      }
    } else if (kind == 8) {
      int i = distr_idx(rng);
      int j = distr_idx(rng);
      if (i == j) {
        continue;
      }
      for (int k = 0; k < n; k++) {
        std::swap(U(i, k), U(j, k));
      }
    } else {
      int i = distr_idx(rng);
      for (int k = 0; k < n; k++) {
        U(i, k) = -U(i, k);
      }
    }
  }
  return U;
}

/*
  A signed permutation matrix is the trivial change of presentation: it leaves
  every quantity this benchmark measures unchanged. Recovery is therefore
  tested modulo it rather than against U itself.
 */
template <typename Tint> bool IsSignedPermutation(MyMatrix<Tint> const &M) {
  int n = M.rows();
  if (n != M.cols()) {
    return false;
  }
  std::vector<int> n_row_nz(n, 0), n_col_nz(n, 0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Tint const &val = M(i, j);
      if (val != 0) {
        if (val != 1 && val != -1) {
          return false;
        }
        n_row_nz[i]++;
        n_col_nz[j]++;
      }
    }
  }
  for (int i = 0; i < n; i++) {
    if (n_row_nz[i] != 1 || n_col_nz[i] != 1) {
      return false;
    }
  }
  return true;
}

/*
  One instance of the experiment: the hidden good presentation, the unimodular
  matrix that destroyed it, and the presentation actually handed to the
  algorithm.
 */
template <typename T, typename Tint> struct HiddenBasisInstance {
  std::string name;
  int n;
  int n_ops;
  MyMatrix<T> GramGood;
  MyMatrix<Tint> U;
  MyMatrix<T> GramBad;
};

template <typename T, typename Tint>
HiddenBasisInstance<T, Tint> MakeHiddenBasisInstance(std::string const &name,
                                                     int n, int n_ops,
                                                     std::mt19937_64 &rng) {
  MyMatrix<T> GramGood = GoodGramByName<T>(name, n, rng);
  int dim = GramGood.rows();
  MyMatrix<Tint> U = RandomUnimodular<Tint>(dim, n_ops, rng);
  MyMatrix<T> U_T = UniversalMatrixConversion<T, Tint>(U);
  MyMatrix<T> GramBad = U_T * GramGood * U_T.transpose();
#ifdef SANITY_CHECK_REDUCTION_BENCH
  T det_good = DeterminantMat(GramGood);
  T det_bad = DeterminantMat(GramBad);
  if (det_good != det_bad) {
    std::cerr << "REDUCTION_BENCH: the destruction step changed the lattice, "
                 "det went from "
              << det_good << " to " << det_bad << "\n";
    throw TerminalException{1};
  }
#endif
  return {name, dim, n_ops, std::move(GramGood), std::move(U),
          std::move(GramBad)};
}

/*
  The outcome of running one algorithm on one instance. recovered answers the
  question of Section 29.1 of the overview, quality that of Section 29.2, and
  the point of keeping both is that they come apart.
 */
template <typename T> struct ReductionOutcome {
  std::string algo;
  bool recovered;
  bool matches_good_quality;
  ReductionQuality<T> quality;
  double runtime_ms;
};

/*
  Runs one reducer on one instance. The reducer is handed the bad Gram matrix
  and must return an LLLreduction, that is a pair (G_red, P) with
  P G_bad P^T = G_red; the check that this actually holds is not optional here,
  since a reducer that quietly leaves the lattice is otherwise indistinguishable
  from a very good one.
 */
template <typename T, typename Tint, typename Freduce>
ReductionOutcome<T> RunOneReduction(HiddenBasisInstance<T, Tint> const &inst,
                                    std::string const &algo, Freduce f_reduce,
                                    std::ostream &os) {
  MicrosecondTime time;
  LLLreduction<T, Tint> res = f_reduce(inst.GramBad, os);
  double runtime_ms = static_cast<double>(time.const_eval_int64()) / 1000.0;
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
  MyMatrix<T> check = P_T * inst.GramBad * P_T.transpose();
  if (check != res.GramMatRed) {
    std::cerr << "REDUCTION_BENCH: algorithm " << algo
              << " returned a transformation inconsistent with its own "
                 "reduced form\n";
    throw TerminalException{1};
  }
  T det_P = DeterminantMat(P_T);
  if (det_P != 1 && det_P != -1) {
    std::cerr << "REDUCTION_BENCH: algorithm " << algo
              << " returned a non-unimodular transformation, det = " << det_P
              << "\n";
    throw TerminalException{1};
  }
  // P U carries the bad presentation back to the reduced one through the good
  // one, so recovery is exactly the statement that it is trivial.
  MyMatrix<Tint> PU = res.Pmat * inst.U;
  bool recovered = IsSignedPermutation(PU);
  ReductionQuality<T> quality = ComputeReductionQuality(res.GramMatRed);
  ReductionQuality<T> quality_good = ComputeReductionQuality(inst.GramGood);
  // Not equality: the demand is that the algorithm did at least as well as the
  // hidden presentation, which is what "reduction succeeded" should mean.
  bool matches_good_quality =
      (quality.orth_defect_sq <= quality_good.orth_defect_sq);
  return {algo, recovered, matches_good_quality, std::move(quality),
          runtime_ms};
}

template <typename T>
void PrintReductionQuality(std::ostream &os, std::string const &label,
                           ReductionQuality<T> const &q) {
  os << "  " << label << ": defect^2=" << q.orth_defect_sq
     << " log(defect)=" << q.log_orth_defect << " pot=" << q.lll_potential
     << " maxcoeff=" << q.max_coeff << " bits=" << q.total_bits
     << " E_A=" << q.log_E_A << " E_N=" << q.log_E_N << "\n";
}

// clang-format off
#endif  // SRC_LATT_LATTICEREDUCTIONBENCH_H_
// clang-format on
