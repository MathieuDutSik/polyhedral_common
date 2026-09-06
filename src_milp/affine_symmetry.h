// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_AFFINE_SYMMETRY_H_
#define SRC_MILP_AFFINE_SYMMETRY_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "PolytopeEquiStab.h"
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_AFFINE_SYMMETRY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_AFFINE_SYMMETRY
#endif

#ifdef TIMINGS
#define TIMINGS_AFFINE_SYMMETRY
#endif

/*
  G_aff, the group of the permutations of the n coordinates that
  preserve the affine subspace S = {x : A x = b}.

  It depends on S alone and not on the rows that were written down, so
  it contains the G_mat of system_symmetry.h. Any of its elements maps
  the solutions of the system to solutions, in particular the 0/1 ones.

  The right hand side is handled by homogenizing. The affine subspace
  S of R^n corresponds to the linear subspace

      C(S) = {(x, t) in R^{n+1} : A x = t b}

  and a permutation sigma of the n coordinates, extended by fixing the
  coordinate t, preserves C(S) if and only if it preserves S. Indeed
  the points of C(S) with t != 0 are the (t x, t) for x in S, and those
  with t = 0 are the kernel of A, which sigma preserves as soon as it
  preserves S = x_0 + ker A. So b stops being a constant to deal with
  and becomes the column n+1 of the matrix [A | -b]; the only thing
  left to require is that the group does not move that column, which is
  a vertex colour.

  Then, the coordinate permutations being orthogonal, they preserve a
  subspace exactly when they preserve its orthogonal complement, and
  the complement of C(S) is the row space of [A | -b]. Working there
  rather than in the kernel lowers the dimension from n - rank + 1 to
  rank, which is a large gain when the system has few rows.

  So, with R a basis of the row space of [A | -b], of size r x (n+1),
  the n+1 columns of R are a family of n+1 vectors of R^r and sigma
  preserves the row space if and only if R P_sigma^T = N R for some
  invertible N, that is if and only if sigma preserves that family up
  to a linear transformation. That is the automorphism group of a
  vector configuration, computed by the machinery of PolytopeEquiStab,
  with the last vector coloured apart so that it stays fixed.
*/

/*
  The group of the permutations of the n coordinates preserving
  {x : A x = b}, given as a group on n points.
  ---
  T has to be a field. The system must be consistent; if it is not,
  the returned group is the one of the empty affine subspace, which is
  meaningless, so the caller checks the consistency first.
*/
template <typename T, typename Tgroup>
Tgroup ComputeAffineSymmetry(MyMatrix<T> const &A, MyVector<T> const &b,
                             std::ostream &os) {
  static_assert(is_ring_field<T>::value,
                "ComputeAffineSymmetry requires a field");
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_row = A.rows();
  int n_col = A.cols();
#ifdef SANITY_CHECK_AFFINE_SYMMETRY
  if (b.size() != n_row) {
    std::cerr << "AFFINE_SYMMETRY: A has " << n_row << " rows but b has "
              << b.size() << " entries\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_AFFINE_SYMMETRY
  MicrosecondTime time;
#endif
  // M = [A | -b], of size n_row x (n_col + 1)
  MyMatrix<T> M(n_row, n_col + 1);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++)
      M(i, j) = A(i, j);
    M(i, n_col) = -b(i);
  }
  // A basis of the row space of M
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(M);
  MyMatrix<T> R = SelectRow(M, eSelect.ListRowSelect);
#ifdef DEBUG_AFFINE_SYMMETRY
  os << "AFFINE_SYMMETRY: n_row=" << n_row << " n_col=" << n_col
     << " rank=" << eSelect.TheRank << "\n";
#endif
  // The configuration of the n_col + 1 columns of R, as the rows of a
  // (n_col + 1) x rank matrix
  MyMatrix<T> EXT = R.transpose();
#ifdef TIMINGS_AFFINE_SYMMETRY
  os << "AFFINE_SYMMETRY: building the configuration took " << time << "\n";
#endif
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  std::vector<MyMatrix<T>> ListMat{Qmat};
  // The homogenizing column is the last one and must stay fixed
  std::vector<T> Vdiag(n_col + 1, T(0));
  Vdiag[n_col] = T(1);
  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(EXTred, ListMat, Vdiag,
                                                         os);
#ifdef TIMINGS_AFFINE_SYMMETRY
  os << "AFFINE_SYMMETRY: GetListGenAutomorphism_ListMat_Vdiag took " << time
     << "\n";
#endif
  std::vector<Telt> LGen;
  for (auto &eList : ListGen) {
#ifdef SANITY_CHECK_AFFINE_SYMMETRY
    if (static_cast<int>(eList[n_col]) != n_col) {
      std::cerr << "AFFINE_SYMMETRY: a generator moves the homogenizing "
                << "coordinate to " << eList[n_col] << "\n";
      throw TerminalException{1};
    }
#endif
    std::vector<Tidx> eListRed(n_col);
    for (int j = 0; j < n_col; j++)
      eListRed[j] = eList[j];
    LGen.push_back(Telt(eListRed));
  }
  return Tgroup(LGen, n_col);
}

// clang-format off
#endif  // SRC_MILP_AFFINE_SYMMETRY_H_
// clang-format on
