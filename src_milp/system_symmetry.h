// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_SYSTEM_SYMMETRY_H_
#define SRC_MILP_SYSTEM_SYMMETRY_H_

// clang-format off
#include "MAT_Matrix.h"
#include "PolytopeEquiStab.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SYSTEM_SYMMETRY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SYSTEM_SYMMETRY
#endif

#ifdef TIMINGS
#define TIMINGS_SYSTEM_SYMMETRY
#endif

/*
  G_mat, the automorphism group of the system A x = b as it is written.

  It is the set of permutations sigma of the n columns for which there
  is a permutation tau of the m rows with
     A[tau(i)][sigma(j)] = A[i][j]   and   b[tau(i)] = b[i]
  Such a sigma maps solutions to solutions, whatever the constraint put
  on the entries of x, since it only relabels the system.

  Beware that G_mat depends on the rows that were written down, not
  only on the set of solutions: another generating set of the same
  affine subspace can have a different G_mat. It is a subgroup of the
  stabilizer of the affine subspace, which affine_symmetry.h computes.

  The computation is the automorphism group of the bipartite structure
  of rows and columns, edges coloured by the entries of A, rows
  coloured by their right hand side. It is folded into a single
  symmetric weight matrix on the m + n points

      W(i, i)       = (ROW_DIAG, b_i)          for a row i
      W(j, j)       = (COL_DIAG, colour_j)     for a column j
      W(i, i')      = (ROW_ROW, 0)             two rows
      W(j, j')      = (COL_COL, 0)             two columns
      W(i, j)       = (ENTRY, A[i][j])         a row and a column

  the first component of the weight being what keeps the five blocks
  apart: no weight of one block is ever equal to a weight of another,
  so no automorphism can send a row to a column. The diagonal carries
  the vertex colours, as in the Vdiag convention of the rest of the
  code.

  The column colours are an input so that the group of a partial
  assignment can be asked for: giving the columns already set to 0, to
  1 and still free three distinct colours returns the subgroup of
  G_mat stabilizing that partial assignment, which is what an orbital
  branching needs.
*/

// The blocks of the weight matrix. The value only has to separate
// them, its order is irrelevant.
static const int SYMSYS_ROW_DIAG = 0;
static const int SYMSYS_COL_DIAG = 1;
static const int SYMSYS_ROW_ROW = 2;
static const int SYMSYS_COL_COL = 3;
static const int SYMSYS_ENTRY = 4;

/*
  The group of the system, as permutations of the m + n points. The
  columns are the points m, ..., m + n - 1.
  ---
  ColColour is either empty, meaning that all the columns are
  interchangeable a priori, or of length n.
*/
template <typename T, typename Tgroup>
Tgroup ComputeSystemSymmetryFull(MyMatrix<T> const &A, MyVector<T> const &b,
                                 std::vector<T> const &ColColour,
                                 std::ostream &os) {
  using Tgr = GraphListAdj;
  using Tweight = std::pair<int, T>;
  int n_row = A.rows();
  int n_col = A.cols();
#ifdef SANITY_CHECK_SYSTEM_SYMMETRY
  if (b.size() != n_row) {
    std::cerr << "SYSTEM_SYMMETRY: A has " << n_row << " rows but b has "
              << b.size() << " entries\n";
    throw TerminalException{1};
  }
  if (ColColour.size() != 0 && static_cast<int>(ColColour.size()) != n_col) {
    std::cerr << "SYSTEM_SYMMETRY: ColColour has " << ColColour.size()
              << " entries instead of " << n_col << "\n";
    throw TerminalException{1};
  }
#endif
  bool has_colour = ColColour.size() > 0;
  size_t n_point = static_cast<size_t>(n_row) + static_cast<size_t>(n_col);
  auto f = [&](size_t i, size_t j) -> Tweight {
    bool i_row = i < static_cast<size_t>(n_row);
    bool j_row = j < static_cast<size_t>(n_row);
    if (i == j) {
      if (i_row)
        return {SYMSYS_ROW_DIAG, b(i)};
      T colour = has_colour ? ColColour[i - n_row] : T(0);
      return {SYMSYS_COL_DIAG, colour};
    }
    if (i_row && j_row)
      return {SYMSYS_ROW_ROW, T(0)};
    if (!i_row && !j_row)
      return {SYMSYS_COL_COL, T(0)};
    size_t idx_row = i_row ? i : j;
    size_t idx_col = (i_row ? j : i) - n_row;
    return {SYMSYS_ENTRY, A(idx_row, idx_col)};
  };
#ifdef TIMINGS_SYSTEM_SYMMETRY
  MicrosecondTime time;
#endif
  WeightMatrix<true, Tweight, uint32_t> WMat(n_point, f, os);
#ifdef TIMINGS_SYSTEM_SYMMETRY
  os << "SYSTEM_SYMMETRY: WeightMatrix took " << time << "\n";
#endif
  Tgroup GRP =
      GetStabilizerWeightMatrix<Tweight, Tgr, Tgroup, uint32_t>(WMat, os);
#ifdef TIMINGS_SYSTEM_SYMMETRY
  os << "SYSTEM_SYMMETRY: GetStabilizerWeightMatrix took " << time << "\n";
#endif
#ifdef SANITY_CHECK_SYSTEM_SYMMETRY
  // No generator may send a row to a column, that being what the block
  // separation of the weights is for.
  using Tidx = typename Tgroup::Telt::Tidx;
  for (auto &eGen : GRP.GeneratorsOfGroup()) {
    for (int i = 0; i < n_row; i++) {
      Tidx img = OnPoints(static_cast<Tidx>(i), eGen);
      if (static_cast<int>(img) >= n_row) {
        std::cerr << "SYSTEM_SYMMETRY: a generator sends the row " << i
                  << " to the point " << img << ", which is a column\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return GRP;
}

// The same group, restricted to its action on the n columns. This is
// the group acting on the solutions.
template <typename T, typename Tgroup>
Tgroup ComputeSystemSymmetry(MyMatrix<T> const &A, MyVector<T> const &b,
                             std::vector<T> const &ColColour,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_row = A.rows();
  int n_col = A.cols();
  Tgroup GRPfull = ComputeSystemSymmetryFull<T, Tgroup>(A, b, ColColour, os);
  std::vector<Telt> ListGen;
  for (auto &eGen : GRPfull.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(n_col);
    for (int j = 0; j < n_col; j++) {
      Tidx img = OnPoints(static_cast<Tidx>(n_row + j), eGen);
      eList[j] = static_cast<Tidx>(img - n_row);
    }
    ListGen.push_back(Telt(eList));
  }
  return Tgroup(ListGen, n_col);
}

// clang-format off
#endif  // SRC_MILP_SYSTEM_SYMMETRY_H_
// clang-format on
