// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_ZERO_ONE_LATTICE_H_
#define SRC_MILP_ZERO_ONE_LATTICE_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include <cmath>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ZERO_ONE_LATTICE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ZERO_ONE_LATTICE
#endif

#ifdef TIMINGS
#define TIMINGS_ZERO_ONE_LATTICE
#endif

/*
  The lattice formulation of the enumeration of the x in {0,1}^n with
  A x = b, in the spirit of the solvediophant of A. Wassermann.

  The right hand side is put into the lattice by homogenizing, which is
  what avoids a coset and an integral particular solution:

      L = ker_Z([A | -b])  in  Z^{n+1}

  a vector (z, s) of L meaning A z = s b. On Z^{n+1} take

      Q(z, s) = || 2z - s 1_n ||^2 + s^2

  which is positive definite, of matrix

      M = [ [ 4 I_n , -2 1_n ] , [ -2 1_n^T , n+1 ] ]

  A 0/1 solution gives Q = n + 1. Conversely a vector of L with
  Q <= n+1 and s = +-1 has all its 2z_j - s odd, hence equal to +-1,
  hence is a solution up to the global sign. So the solutions are
  exactly the vectors of L of Q-norm at most n+1 with s = +-1, and the
  ones with s = 0 are the short vectors of ker_Z(A), which are
  discarded.

  Compared with taking a particular solution x_0 and enumerating the
  coset x_0 + ker_Z(A), this keeps every quantity a small integer: a
  particular solution coming out of a Hermite normal form has entries
  of hundreds of digits, and the centre of the coset enumeration is
  then a rational with a denominator of the size of the determinant.
  Compared with the (n+1+m)-dimensional embedding with an N-scaled
  constraint block, it works in dimension n+1-rank rather than n+1 and
  has no N^2 in the entries.

  Whether the enumeration is affordable is another matter, and
  EstimateEnumerationCost answers it before the enumeration is
  launched: the cost of a Fincke-Pohst is governed by the
  Gram-Schmidt profile of the basis against the radius, and for these
  systems the solutions are much longer than the shortest vectors of
  the lattice, which is the bad case.
*/

template <typename T, typename Tint> struct ZeroOneLattice {
  int n_row;
  int n_col;
  // Dimension of the lattice, that is n_col + 1 - rank([A | -b])
  int dim;
  // The reduced basis, dim x (n_col + 1). A lattice vector is
  // t . Kred, its n_col first entries being z and the last one s.
  MyMatrix<Tint> Kred;
  // The Gram matrix of Kred for the form Q, dim x dim
  MyMatrix<T> Gred;
  // The norm that a solution has, n_col + 1
  T bound;
};

// The matrix of Q on Z^{n+1}
template <typename T> MyMatrix<T> ZeroOneLatticeForm(int n) {
  MyMatrix<T> M = ZeroMatrix<T>(n + 1, n + 1);
  for (int j = 0; j < n; j++) {
    M(j, j) = 4;
    M(j, n) = -2;
    M(n, j) = -2;
  }
  M(n, n) = n + 1;
  return M;
}

template <typename T, typename Tint>
ZeroOneLattice<T, Tint> BuildZeroOneLattice(MyMatrix<T> const &A,
                                            MyVector<T> const &b,
                                            std::ostream &os) {
  int n_row = A.rows();
  int n_col = A.cols();
#ifdef SANITY_CHECK_ZERO_ONE_LATTICE
  if (b.size() != n_row) {
    std::cerr << "ZERO_ONE_LATTICE: A has " << n_row << " rows but b has "
              << b.size() << " entries\n";
    throw TerminalException{1};
  }
#endif
  // Mh = [A | -b]
  MyMatrix<Tint> Mh(n_row, n_col + 1);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++)
      Mh(i, j) = UniversalScalarConversion<Tint, T>(A(i, j));
    Mh(i, n_col) = -UniversalScalarConversion<Tint, T>(b(i));
  }
#ifdef TIMINGS_ZERO_ONE_LATTICE
  MicrosecondTime time;
#endif
  MyMatrix<Tint> K = NullspaceIntTrMat(Mh);
  int dim = K.rows();
#ifdef TIMINGS_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: the kernel of dimension " << dim << " took " << time
     << "\n";
#endif
#ifdef DEBUG_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: n_row=" << n_row << " n_col=" << n_col
     << " dim=" << dim << "\n";
#endif
  T bound = UniversalScalarConversion<T, int>(n_col + 1);
  if (dim == 0) {
    // The homogenized system has no nonzero solution, so A x = b has
    // none either
    MyMatrix<Tint> Kempty(0, n_col + 1);
    MyMatrix<T> Gempty(0, 0);
    return {n_row, n_col, 0, std::move(Kempty), std::move(Gempty), bound};
  }
#ifdef SANITY_CHECK_ZERO_ONE_LATTICE
  MyMatrix<Tint> prod = Mh * K.transpose();
  if (prod.cwiseAbs().maxCoeff() != 0) {
    std::cerr << "ZERO_ONE_LATTICE: the computed kernel does not annihilate "
              << "[A | -b]\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> Mq = ZeroOneLatticeForm<T>(n_col);
  MyMatrix<T> K_T = UniversalMatrixConversion<T, Tint>(K);
  MyMatrix<T> G = K_T * Mq * K_T.transpose();
#ifdef TIMINGS_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: the Gram matrix took " << time << "\n";
#endif
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(G, os);
  MyMatrix<Tint> Kred = rec.Pmat * K;
#ifdef TIMINGS_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: the LLL reduction took " << time << "\n";
#endif
  return {n_row, n_col, dim, std::move(Kred), std::move(rec.GramMatRed), bound};
}

// What a Fincke-Pohst over the lattice would cost, before running it.
struct ZeroOneLatticeEstimate {
  // Squared norms of the Gram-Schmidt vectors of the reduced basis
  double gs_min;
  double gs_max;
  // Decimal logarithm of the number of nodes of the widest level of
  // the enumeration tree, by the Gaussian heuristic, and that level
  double log10_max_node;
  int level;
  // The same for a perfectly flat profile of the same determinant,
  // which no basis reduction can beat
  double log10_max_node_flat;
  int level_flat;
  double log10_det;
};

template <typename T, typename Tint>
ZeroOneLatticeEstimate
EstimateEnumerationCost(ZeroOneLattice<T, Tint> const &lattice) {
  int d = lattice.dim;
  if (d == 0)
    return {0, 0, 0, 0, 0, 0, 0};
  // Gram-Schmidt in double: only an estimate is wanted here, every
  // conclusion drawn from it being about orders of magnitude
  std::vector<std::vector<double>> mu(d, std::vector<double>(d, 0.0));
  std::vector<double> B(d, 0.0);
  std::vector<std::vector<double>> Gd(d, std::vector<double>(d, 0.0));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      Gd[i][j] = UniversalScalarConversion<double, T>(lattice.Gred(i, j));
  for (int i = 0; i < d; i++) {
    B[i] = Gd[i][i];
    for (int j = 0; j < i; j++) {
      double scal = Gd[i][j];
      for (int k = 0; k < j; k++)
        scal -= mu[i][k] * mu[j][k] * B[k];
      mu[i][j] = scal / B[j];
      B[i] -= mu[i][j] * mu[i][j] * B[j];
    }
  }
  double Rsq = UniversalScalarConversion<double, T>(lattice.bound);
  ZeroOneLatticeEstimate est;
  est.gs_min = B[0];
  est.gs_max = B[0];
  est.log10_det = 0.0;
  for (int i = 0; i < d; i++) {
    est.gs_min = std::min(est.gs_min, B[i]);
    est.gs_max = std::max(est.gs_max, B[i]);
    est.log10_det += std::log10(B[i]);
  }
  // Number of lattice points of the projected sublattice of rank k in
  // the ball of squared radius Rsq, by the Gaussian heuristic
  auto log10_ball = [&](int k) -> double {
    return 0.5 * k * std::log10(M_PI) + 0.5 * k * std::log10(Rsq) -
           std::lgamma(1.0 + k / 2.0) / std::log(10.0);
  };
  est.log10_max_node = -1e300;
  est.level = 0;
  est.log10_max_node_flat = -1e300;
  est.level_flat = 0;
  double flat = est.log10_det / d;
  for (int k = 1; k <= d; k++) {
    double val = log10_ball(k);
    for (int i = d - k; i < d; i++)
      val -= 0.5 * std::log10(B[i]);
    if (val > est.log10_max_node) {
      est.log10_max_node = val;
      est.level = k;
    }
    double val_flat = log10_ball(k) - 0.5 * k * flat;
    if (val_flat > est.log10_max_node_flat) {
      est.log10_max_node_flat = val_flat;
      est.level_flat = k;
    }
  }
  return est;
}

// The solutions, by the enumeration of the lattice vectors of Q-norm
// at most n+1. Only affordable when EstimateEnumerationCost says so.
template <typename T, typename Tint>
std::vector<Face> EnumerateZeroOneByLattice(ZeroOneLattice<T, Tint> const &lat,
                                            std::ostream &os) {
  std::vector<Face> ListSol;
  if (lat.dim == 0)
    return ListSol;
  int n = lat.n_col;
#ifdef TIMINGS_ZERO_ONE_LATTICE
  MicrosecondTime time;
#endif
  std::vector<MyVector<Tint>> ListV =
      computeLevel_GramMat<T, Tint>(lat.Gred, lat.bound, os);
#ifdef TIMINGS_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: the enumeration took " << time << "\n";
#endif
#ifdef DEBUG_ZERO_ONE_LATTICE
  size_t n_zero = 0;
  os << "ZERO_ONE_LATTICE: the enumeration returned " << ListV.size()
     << " vectors up to sign\n";
#endif
  for (auto &t : ListV) {
    MyVector<Tint> v = lat.Kred.transpose() * t;
    Tint s = v(n);
    if (s == 0) {
#ifdef DEBUG_ZERO_ONE_LATTICE
      n_zero++;
#endif
      continue;
    }
    if (s != 1 && s != -1)
      continue;
    // The enumeration gives one vector per antipodal pair, so the
    // representative with s = -1 is the negative of a solution
    Face sol(n);
    bool is01 = true;
    for (int j = 0; j < n; j++) {
      Tint val = (s == 1) ? v(j) : -v(j);
      if (val == 1) {
        sol[j] = 1;
      } else {
        if (val != 0)
          is01 = false;
      }
    }
    if (is01)
      ListSol.push_back(std::move(sol));
  }
#ifdef DEBUG_ZERO_ONE_LATTICE
  os << "ZERO_ONE_LATTICE: " << ListSol.size() << " solutions, " << n_zero
     << " vectors of the kernel discarded\n";
#endif
  return ListSol;
}

// clang-format off
#endif  // SRC_MILP_ZERO_ONE_LATTICE_H_
// clang-format on
