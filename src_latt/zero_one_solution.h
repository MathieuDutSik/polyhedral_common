// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ZERO_ONE_SOLUTION_H_
#define SRC_LATT_ZERO_ONE_SOLUTION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "Boost_bitset_kernel.h"
#include "POLY_SimplexClarkson.h"
#include <limits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ZERO_ONE_SOLUTION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ZERO_ONE_SOLUTION
#endif

#ifdef TIMINGS
#define TIMINGS_ZERO_ONE_SOLUTION
#endif

/*
  Enumeration of all the x in {0,1}^n satisfying A x = b, for A an
  integer matrix with n_row rows and n columns and b an integer vector.

  Algorithm: a depth first search over the coordinates, with the two
  prunings below. The state of a node is the assignment of a subset of
  the coordinates; for a row i we keep
  --- resid_i = b_i - sum_{x_j = 1} a_ij, what the still unassigned
      coordinates have to contribute,
  --- lo_i = sum_{j unset} min(a_ij, 0) and hi_i = sum_{j unset}
      max(a_ij, 0), the range that they can contribute.

  --- Bound propagation. The row i is infeasible as soon as resid_i is
      outside of [lo_i, hi_i]. Removing an unset j of coefficient a from
      the row leaves the range [lo_i - min(a,0), hi_i - max(a,0)], so
      x_j = 1 is possible only when resid_i - a lies in it and x_j = 0
      only when resid_i does. When exactly one of the two values passes
      the test for some row, the coordinate is forced, and forcing it
      changes the rows it belongs to, so the propagation is run to a
      fixpoint over a queue of the rows that changed. For a nonnegative
      A this already covers the usual reasoning: a row of residual 0
      forces every remaining coordinate of positive coefficient to 0,
      and a row whose residual equals its maximum forces them to 1.
  --- Linear programming. Any solution of the subtree satisfies the
      relaxation A_free y = resid, 0 <= y <= 1 over the unset
      coordinates, so an infeasible relaxation prunes the subtree. This
      is exact but costly, hence it is only run down to the depth
      lp_max_depth, the shallow nodes being those where a refutation
      pays for itself. It is the same use of exact linear programming as
      in the branch and bound of strongly_semi_eutactic.h.

  Every reported solution is a genuine one, the assignment being
  complete and every row having residual 0; the enumeration is complete
  when the search tree was exhausted, that is when the node budget
  max_node was not hit.
*/

struct ZeroOneOptions {
  // Bound on the number of nodes of the search tree. When it is hit the
  // enumeration stops and the result is reported as not resolved, the
  // solutions found so far having been reported.
  size_t max_node = std::numeric_limits<size_t>::max();
  // Depth down to which the linear programming relaxation is used for
  // pruning. A negative value disables it.
  int lp_max_depth = 0;
  // Branch on the first unset coordinate instead of the most
  // constrained one.
  bool natural_order = false;
};

struct ZeroOneResult {
  // false when the node budget was exhausted, in which case the list of
  // solutions is possibly incomplete
  bool resolved;
  size_t n_node;
  size_t n_solution;
};

// The search itself. F is called on the Face of each solution found, in
// the order of the depth first search.
template <typename T, typename F> struct ZeroOneSolutionSearch {
private:
  static constexpr uint8_t STATE_UNSET = 0;
  static constexpr uint8_t STATE_ZERO = 1;
  static constexpr uint8_t STATE_ONE = 2;
  int n_row;
  int n_col;
  // The nonzero entries of A, by row as (column, coefficient) and by
  // column as (row, coefficient)
  std::vector<std::vector<std::pair<int, T>>> ListRow;
  std::vector<std::vector<std::pair<int, T>>> ListCol;
  std::vector<uint8_t> state;
  std::vector<T> resid;
  std::vector<T> lo;
  std::vector<T> hi;
  std::vector<int> n_free_row;
  int n_free_total;
  // The coordinates assigned since the root, for the undo on backtrack
  std::vector<int> trail;
  // The rows whose (resid, lo, hi) changed and that the propagation has
  // still to examine
  std::vector<int> queue_row;
  std::vector<uint8_t> in_queue;
  ZeroOneOptions options;
  F f;
  std::ostream &os;
  size_t n_node;
  size_t n_solution;
  bool budget_hit;

  void PushQueue(int i) {
    if (in_queue[i] == 0) {
      in_queue[i] = 1;
      queue_row.push_back(i);
    }
  }

  void ClearQueue() {
    for (auto &i : queue_row)
      in_queue[i] = 0;
    queue_row.clear();
  }

  void Assign(int j, uint8_t val) {
    state[j] = val;
    n_free_total--;
    trail.push_back(j);
    for (auto &pair : ListCol[j]) {
      int i = pair.first;
      T const &a = pair.second;
      n_free_row[i]--;
      if (a < 0) {
        lo[i] -= a;
      } else {
        hi[i] -= a;
      }
      if (val == STATE_ONE)
        resid[i] -= a;
      PushQueue(i);
    }
  }

  void UndoTo(size_t trail_size) {
    while (trail.size() > trail_size) {
      int j = trail.back();
      trail.pop_back();
      uint8_t val = state[j];
      state[j] = STATE_UNSET;
      n_free_total++;
      for (auto &pair : ListCol[j]) {
        int i = pair.first;
        T const &a = pair.second;
        n_free_row[i]++;
        if (a < 0) {
          lo[i] += a;
        } else {
          hi[i] += a;
        }
        if (val == STATE_ONE)
          resid[i] += a;
      }
    }
  }

  // Bound propagation to a fixpoint. Returns false when a row is shown
  // to be infeasible.
  bool Propagate() {
    while (!queue_row.empty()) {
      int i = queue_row.back();
      queue_row.pop_back();
      in_queue[i] = 0;
      if (resid[i] < lo[i] || resid[i] > hi[i])
        return false;
      if (n_free_row[i] == 0)
        continue;
      for (auto &pair : ListRow[i]) {
        int j = pair.first;
        if (state[j] != STATE_UNSET)
          continue;
        T const &a = pair.second;
        // The range of the row once j is taken out of it
        T lo_rest = lo[i];
        T hi_rest = hi[i];
        if (a < 0) {
          lo_rest -= a;
        } else {
          hi_rest -= a;
        }
        T res_one = resid[i] - a;
        bool ok_one = lo_rest <= res_one && res_one <= hi_rest;
        bool ok_zero = lo_rest <= resid[i] && resid[i] <= hi_rest;
        if (!ok_one && !ok_zero)
          return false;
        if (ok_one && ok_zero)
          continue;
        Assign(j, ok_one ? STATE_ONE : STATE_ZERO);
      }
    }
    return true;
  }

  // The coordinate to branch on: one of a row with the fewest unset
  // coordinates left, the one of largest coefficient in it. A
  // coordinate belonging to no row is only reached when every row is
  // fully assigned, and then both of its values extend to solutions.
  int SelectBranchVariable() const {
    auto first_unset = [&]() -> int {
      for (int j = 0; j < n_col; j++) {
        if (state[j] == STATE_UNSET)
          return j;
      }
      return -1;
    };
    if (options.natural_order)
      return first_unset();
    int best_row = -1;
    int best_count = 0;
    for (int i = 0; i < n_row; i++) {
      if (n_free_row[i] == 0)
        continue;
      if (best_row == -1 || n_free_row[i] < best_count) {
        best_row = i;
        best_count = n_free_row[i];
      }
    }
    if (best_row == -1)
      return first_unset();
    int best_col = -1;
    T best_abs(0);
    for (auto &pair : ListRow[best_row]) {
      int j = pair.first;
      if (state[j] != STATE_UNSET)
        continue;
      T abs_a = pair.second < 0 ? -pair.second : pair.second;
      if (best_col == -1 || abs_a > best_abs) {
        best_col = j;
        best_abs = abs_a;
      }
    }
    return best_col;
  }

  // Feasibility of A_free y = resid, 0 <= y <= 1 over the unset
  // coordinates. Returns false only when the relaxation is provably
  // empty, which refutes the whole subtree.
  bool IsRelaxationFeasible() {
    std::vector<int> lp_index(n_col, -1);
    int r = 0;
    for (int j = 0; j < n_col; j++) {
      if (state[j] == STATE_UNSET) {
        lp_index[j] = r;
        r++;
      }
    }
    if (r == 0)
      return true;
    int n_cons = 2 * n_row + 2 * r;
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_cons, r + 1);
    int pos = 0;
    for (int i = 0; i < n_row; i++) {
      // resid_i - sum_j a_ij y_j >= 0 and its opposite, for the equality
      ListIneq(pos, 0) = resid[i];
      ListIneq(pos + 1, 0) = -resid[i];
      for (auto &pair : ListRow[i]) {
        int k = lp_index[pair.first];
        if (k >= 0) {
          ListIneq(pos, 1 + k) = -pair.second;
          ListIneq(pos + 1, 1 + k) = pair.second;
        }
      }
      pos += 2;
    }
    for (int k = 0; k < r; k++) {
      // y_k >= 0 and 1 - y_k >= 0
      ListIneq(pos, 1 + k) = 1;
      ListIneq(pos + 1, 0) = 1;
      ListIneq(pos + 1, 1 + k) = -1;
      pos += 2;
    }
    MyVector<T> ToBeMinimized = ZeroVector<T>(r + 1);
    if constexpr (is_ring_field<T>::value) {
      LpSolution<T> eSol =
          SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
      return eSol.DirectSolution.has_value();
    } else {
      LpSolutionScaled<T> eSol =
          SIMPLEX_LinearProgramming_scaled(ListIneq, ToBeMinimized, os);
      return eSol.DirectSolution.has_value();
    }
  }

  void RecordSolution() {
#ifdef SANITY_CHECK_ZERO_ONE_SOLUTION
    for (int i = 0; i < n_row; i++) {
      if (resid[i] != 0) {
        std::cerr << "ZERO_ONE_SOLUTION: the row " << i
                  << " of a reported solution has residual " << resid[i]
                  << " instead of 0\n";
        throw TerminalException{1};
      }
    }
#endif
    Face sol(n_col);
    for (int j = 0; j < n_col; j++) {
      if (state[j] == STATE_ONE)
        sol[j] = 1;
    }
    n_solution++;
    f(sol);
  }

  void Explore(int depth) {
    size_t trail_start = trail.size();
    n_node++;
    if (n_node > options.max_node) {
      budget_hit = true;
    } else {
      if (Propagate()) {
        if (n_free_total == 0) {
          RecordSolution();
        } else if (depth > options.lp_max_depth || IsRelaxationFeasible()) {
          int j = SelectBranchVariable();
#ifdef SANITY_CHECK_ZERO_ONE_SOLUTION
          if (j == -1) {
            std::cerr << "ZERO_ONE_SOLUTION: no coordinate to branch on while "
                      << n_free_total << " are unset\n";
            throw TerminalException{1};
          }
#endif
          for (uint8_t val : {STATE_ONE, STATE_ZERO}) {
            if (budget_hit)
              break;
            size_t sub_start = trail.size();
            Assign(j, val);
            Explore(depth + 1);
            UndoTo(sub_start);
            ClearQueue();
          }
        }
      }
    }
    UndoTo(trail_start);
    ClearQueue();
  }

public:
  ZeroOneSolutionSearch(MyMatrix<T> const &A, MyVector<T> const &b,
                        ZeroOneOptions const &_options, F _f, std::ostream &_os)
      : n_row(A.rows()), n_col(A.cols()), n_free_total(A.cols()),
        options(_options), f(_f), os(_os), n_node(0), n_solution(0),
        budget_hit(false) {
#ifdef SANITY_CHECK_ZERO_ONE_SOLUTION
    if (b.size() != n_row) {
      std::cerr << "ZERO_ONE_SOLUTION: A has " << n_row << " rows but b has "
                << b.size() << " entries\n";
      throw TerminalException{1};
    }
#endif
    ListRow.resize(n_row);
    ListCol.resize(n_col);
    state.assign(n_col, STATE_UNSET);
    resid.resize(n_row);
    lo.assign(n_row, T(0));
    hi.assign(n_row, T(0));
    n_free_row.assign(n_row, 0);
    in_queue.assign(n_row, 0);
    for (int i = 0; i < n_row; i++) {
      resid[i] = b(i);
      for (int j = 0; j < n_col; j++) {
        T const &a = A(i, j);
        if (a == 0)
          continue;
        ListRow[i].push_back({j, a});
        ListCol[j].push_back({i, a});
        n_free_row[i]++;
        if (a < 0) {
          lo[i] += a;
        } else {
          hi[i] += a;
        }
      }
    }
  }

  ZeroOneResult Run() {
#ifdef DEBUG_ZERO_ONE_SOLUTION
    size_t n_nonzero = 0;
    for (int i = 0; i < n_row; i++)
      n_nonzero += ListRow[i].size();
    os << "ZERO_ONE_SOLUTION: n_row=" << n_row << " n_col=" << n_col
       << " n_nonzero=" << n_nonzero << " lp_max_depth=" << options.lp_max_depth
       << "\n";
#endif
    for (int i = 0; i < n_row; i++)
      PushQueue(i);
    Explore(0);
#ifdef DEBUG_ZERO_ONE_SOLUTION
    os << "ZERO_ONE_SOLUTION: n_node=" << n_node
       << " n_solution=" << n_solution << " resolved=" << (!budget_hit) << "\n";
#endif
    return {!budget_hit, n_node, n_solution};
  }
};

// Enumerate the x in {0,1}^n with A x = b, calling f on the Face of
// each of them.
template <typename T, typename F>
ZeroOneResult EnumerateZeroOneSolutions(MyMatrix<T> const &A,
                                        MyVector<T> const &b,
                                        ZeroOneOptions const &options, F f,
                                        std::ostream &os) {
#ifdef TIMINGS_ZERO_ONE_SOLUTION
  MicrosecondTime time;
#endif
  ZeroOneSolutionSearch<T, F> search(A, b, options, f, os);
  ZeroOneResult result = search.Run();
#ifdef TIMINGS_ZERO_ONE_SOLUTION
  os << "ZERO_ONE_SOLUTION: EnumerateZeroOneSolutions took " << time << "\n";
#endif
  return result;
}

// The same, collecting the solutions instead of streaming them.
template <typename T>
std::pair<ZeroOneResult, std::vector<Face>>
GetAllZeroOneSolutions(MyMatrix<T> const &A, MyVector<T> const &b,
                       ZeroOneOptions const &options, std::ostream &os) {
  std::vector<Face> ListSol;
  auto f = [&](Face const &sol) -> void { ListSol.push_back(sol); };
  ZeroOneResult result = EnumerateZeroOneSolutions(A, b, options, f, os);
  return {result, std::move(ListSol)};
}

// The solutions as the rows of a matrix, one solution per row.
template <typename T>
MyMatrix<T> ZeroOneSolutionsAsMatrix(std::vector<Face> const &ListSol,
                                     int n_col) {
  int n_sol = ListSol.size();
  MyMatrix<T> M(n_sol, n_col);
  for (int i_sol = 0; i_sol < n_sol; i_sol++) {
    for (int j = 0; j < n_col; j++)
      M(i_sol, j) = ListSol[i_sol][j];
  }
  return M;
}

// clang-format off
#endif  // SRC_LATT_ZERO_ONE_SOLUTION_H_
// clang-format on
