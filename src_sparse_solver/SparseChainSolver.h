// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SPARSE_SOLVER_SPARSECHAINSOLVER_H_
#define SRC_SPARSE_SOLVER_SPARSECHAINSOLVER_H_

// clang-format off
#include "GampMatlab.h"
#include "SparseExactSolver.h"
#include "Timings.h"
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SPARSE_CHAIN_SOLVER
#endif

#ifdef TIMINGS
#define TIMINGS_SPARSE_CHAIN_SOLVER
#endif

/*
  The exact sparse solutions of x A = b for the large sparse systems of
  the chain complexes (the rows of A are the cells, the columns their
  boundary cells, b a chain): the solution is wanted exact and sparse.

  - SparseSolutionMat_Flow: when the rows with two entries are the edges
    of a graph on the columns (the facets of the 1-dimensional cells), b
    is a 0-chain and the solution is built as a flow.
  - SparseSolutionMat_L1Support: given an approximate solution in floating
    point (the L1 solution, which is sparse), the exact solution on its
    support.
  - SparseSolutionMat_L1Exact: the combination, with a test of the
    consistency modulo a prime and the exact elimination as the last
    resort, every step bounded by the budget.
 */

// The bounds of SparseSolutionMat_L1Exact.
struct SparseSolveBudget {
  size_t max_entries_modular = 6000000; // fill-in of the modular elimination
  size_t max_entries_exact = 3000000;   // fill-in of the rational elimination
  size_t max_nnz_modular = 400000;      // nonzeros for trying the modular test
  int max_yall1_iter = 20000;           // iterations of the L1 solver
  int n_reweight = 1;                   // reweighting rounds of the L1 solver
};

/*
  The flow solution of x A = b when the rows of A with two entries are the
  edges of a graph on the columns: shortest paths between the nodes of
  nonzero residual, each path being solved exactly. Inconsistent if a
  node of nonzero residual is alone in its component.
 */
template <typename T>
SparseSolveResult<T> SparseSolutionMat_Flow(std::vector<Eigen::Triplet<T>> const &tripletList, int const &nbRow,
                                            int const &nbCol, MyVector<T> const &Bvect,
                                            [[maybe_unused]] std::ostream &os) {
  MyVector<T> empty;
  std::vector<std::vector<std::pair<int, int>>> adj(nbCol);
  std::vector<std::vector<int>> row_cols(nbRow);
  for (auto &t : tripletList) {
    row_cols[t.row()].push_back(t.col());
  }
  for (int i1 = 0; i1 < nbRow; i1++) {
    std::vector<int> const &lv = row_cols[i1];
    if (lv.size() == 2 && lv[0] != lv[1]) {
      adj[lv[0]].push_back({lv[1], i1});
      adj[lv[1]].push_back({lv[0], i1});
    }
  }
  MyVector<T> resid = Bvect;
  MyVector<T> sol = ZeroVector<T>(nbRow);
  size_t n_path = 0;
  while (true) {
    if (n_path > static_cast<size_t>(nbCol)) {
      // Each path zeroes one node, so this cannot happen.
      return {SparseSolveStatus::Aborted, empty};
    }
    int v_start = -1;
    for (int v = 0; v < nbCol; v++) {
      if (resid(v) != 0) {
        if (v_start == -1 || T_abs(resid(v)) > T_abs(resid(v_start))) {
          v_start = v;
        }
      }
    }
    if (v_start == -1) {
      break;
    }
    std::vector<int> prev_edge(nbCol, -1), prev_vert(nbCol, -1);
    std::vector<uint8_t> visited(nbCol, 0);
    std::vector<int> queue{v_start};
    visited[v_start] = 1;
    size_t pos = 0;
    int v_end = -1;
    while (pos < queue.size() && v_end == -1) {
      int v = queue[pos];
      pos += 1;
      for (auto &e : adj[v]) {
        if (!visited[e.first]) {
          visited[e.first] = 1;
          prev_edge[e.first] = e.second;
          prev_vert[e.first] = v;
          queue.push_back(e.first);
          if (resid(e.first) != 0) {
            v_end = e.first;
            break;
          }
        }
      }
    }
    if (v_end == -1) {
      return {SparseSolveStatus::Inconsistent, empty};
    }
    std::vector<int> l_edge;
    std::vector<int> l_node;
    int v = v_end;
    while (v != v_start) {
      l_edge.push_back(prev_edge[v]);
      v = prev_vert[v];
      l_node.push_back(v);
    }
    std::unordered_map<int, int> map_edge, map_node;
    for (size_t u = 0; u < l_edge.size(); u++) {
      map_edge[l_edge[u]] = u;
      map_node[l_node[u]] = u;
    }
    MyMatrix<T> Mpath = ZeroMatrix<T>(l_edge.size(), l_node.size());
    for (auto &t : tripletList) {
      auto iter = map_edge.find(t.row());
      if (iter != map_edge.end()) {
        auto iter2 = map_node.find(t.col());
        if (iter2 != map_node.end()) {
          Mpath(iter->second, iter2->second) += t.value();
        }
      }
    }
    MyVector<T> target = ZeroVector<T>(l_node.size());
    target(map_node[v_start]) = resid(v_start);
    std::optional<MyVector<T>> opt_path = SolutionMat(Mpath, target);
    if (!opt_path) {
      return {SparseSolveStatus::Aborted, empty};
    }
    for (size_t u = 0; u < l_edge.size(); u++) {
      sol(l_edge[u]) += (*opt_path)(u);
    }
    for (auto &t : tripletList) {
      auto iter = map_edge.find(t.row());
      if (iter != map_edge.end()) {
        resid(t.col()) -= (*opt_path)(iter->second) * t.value();
      }
    }
    n_path += 1;
  }
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
  os << "SPARSE_CHAIN: SparseSolutionMat_Flow, " << n_path << " paths\n";
#endif
  return {SparseSolveStatus::Solved, sol};
}

/*
  The exact solution of x A = b supported on the rows where the
  approximate solution x_approx is above a threshold, the thresholds of
  l_thr being tried in turn: the system restricted to those rows and to the
  columns they touch (the other columns must have b = 0) is solved by the
  budgeted exact elimination. Aborted if no threshold gives a solution.
 */
template <typename T>
SparseSolveResult<T> SparseSolutionMat_L1Support(MySparseMatrix<T> const &SpMat, MyVector<T> const &Bvect,
                                                 MyVector<double> const &x_approx,
                                                 std::vector<double> const &l_thr,
                                                 size_t const &max_entries_exact,
                                                 [[maybe_unused]] std::ostream &os) {
  int nbRow = SpMat.rows();
  int nbCol = SpMat.cols();
  std::vector<std::vector<std::pair<int, T>>> row_entries(nbRow);
  for (int k = 0; k < SpMat.outerSize(); ++k) {
    for (typename MySparseMatrix<T>::InnerIterator it(SpMat, k); it; ++it) {
      row_entries[it.row()].emplace_back(it.col(), it.value());
    }
  }
  for (auto &thr : l_thr) {
    std::vector<int> l_rows;
    for (int i_row = 0; i_row < nbRow; i_row++) {
      if (std::abs(x_approx(i_row)) > thr) {
        l_rows.push_back(i_row);
      }
    }
    std::vector<int> col_index(nbCol, -1);
    std::vector<int> l_cols;
    for (auto &i_row : l_rows) {
      for (auto &pair : row_entries[i_row]) {
        if (col_index[pair.first] == -1) {
          col_index[pair.first] = l_cols.size();
          l_cols.push_back(pair.first);
        }
      }
    }
    bool is_ok = true;
    for (int j = 0; j < nbCol; j++) {
      if (col_index[j] == -1 && Bvect(j) != 0) {
        is_ok = false;
      }
    }
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
    os << "SPARSE_CHAIN: SparseSolutionMat_L1Support thr=" << thr << " support=" << l_rows.size()
       << " touched columns=" << l_cols.size() << " is_ok=" << is_ok << "\n";
#endif
    if (is_ok) {
      MyVector<T> Bred(l_cols.size());
      for (size_t v = 0; v < l_cols.size(); v++) {
        Bred(v) = Bvect(l_cols[v]);
      }
      std::vector<Eigen::Triplet<T>> trip_red;
      for (size_t u = 0; u < l_rows.size(); u++) {
        for (auto &pair : row_entries[l_rows[u]]) {
          trip_red.emplace_back(u, col_index[pair.first], pair.second);
        }
      }
      MySparseMatrix<T> SpRed(l_rows.size(), l_cols.size());
      SpRed.setFromTriplets(trip_red.begin(), trip_red.end());
#ifdef TIMINGS_SPARSE_CHAIN_SOLVER
      MicrosecondTime time_exact;
#endif
      SparseSolveResult<T> res_red = SparseSolutionMat_Exact_Budget(SpRed, Bred, max_entries_exact, os);
#ifdef TIMINGS_SPARSE_CHAIN_SOLVER
      os << "|SPARSE_CHAIN: exact solve on support " << l_rows.size() << " x " << l_cols.size()
         << " status=" << static_cast<int>(res_red.status) << "|=" << time_exact << "\n";
#endif
      if (res_red.status == SparseSolveStatus::Solved) {
        MyVector<T> sol = ZeroVector<T>(nbRow);
        for (size_t u = 0; u < l_rows.size(); u++) {
          sol(l_rows[u]) = res_red.sol(u);
        }
        return {SparseSolveStatus::Solved, sol};
      }
    }
    if (static_cast<int>(l_rows.size()) == nbRow) {
      break;
    }
  }
  return {SparseSolveStatus::Aborted, MyVector<T>()};
}

/*
  Solution of x A = b for the sparse matrix given by its triplets.
  When is_graph is set (the facet level), the rows with two entries are
  the edges of a graph on the columns and b is a 0-chain: the solution
  is built as a flow (SparseSolutionMat_Flow). Otherwise the consistency
  is tested modulo a prime when the system is small enough, the L1 solver
  in floating point gives a sparse approximate solution which is solved
  exactly on its support and the columns it touches, and the exact
  rational elimination is the last resort for the consistent systems.
  Every step is bounded by the budget; the status says whether the system
  was solved, is inconsistent, or could not be decided within the budget.
 */
template <typename T>
SparseSolveResult<T> SparseSolutionMat_L1Exact(std::vector<Eigen::Triplet<T>> const &tripletList,
                                               int const &nbRow, int const &nbCol, MyVector<T> const &Bvect,
                                               bool const &is_graph, SparseSolveBudget const &budget,
                                               std::ostream &os) {
  MyVector<T> empty;
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
  os << "SPARSE_CHAIN: SparseSolutionMat_L1Exact nbRow=" << nbRow << " nbCol=" << nbCol
     << " nnz=" << tripletList.size() << " is_graph=" << is_graph << "\n";
#endif
#ifdef TIMINGS_SPARSE_CHAIN_SOLVER
  MicrosecondTime time_solve;
#endif
  if (is_graph) {
    return SparseSolutionMat_Flow(tripletList, nbRow, nbCol, Bvect, os);
  }
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  bool small = (tripletList.size() <= budget.max_nnz_modular);
  bool known_consistent = false;
  if (small) {
    SparseSolveStatus status = SparseSystemConsistent_Mod_Budget(SpMat, Bvect, budget.max_entries_modular, os);
    if (status == SparseSolveStatus::Inconsistent) {
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
      os << "SPARSE_CHAIN: SparseSolutionMat_L1Exact inconsistent (mod p)\n";
#endif
      return {SparseSolveStatus::Inconsistent, empty};
    }
    known_consistent = (status == SparseSolveStatus::Solved);
  }
  // The L1 solution and the exact solution on its support
  MySparseMatrix<double> SpMat_d = UniversalSparseMatrixConversion<double, T>(SpMat);
  MyVector<double> Bvect_d = UniversalVectorConversion<double, T>(Bvect);
  MyVector<double> x_out =
      AMP_SolutionSparseSystem_reweighted_d(SpMat_d, Bvect_d, budget.max_yall1_iter, budget.n_reweight, os);
  SparseSolveResult<T> res =
      SparseSolutionMat_L1Support(SpMat, Bvect, x_out, {1e-4, 1e-6, 1e-8}, budget.max_entries_exact, os);
  if (res.status == SparseSolveStatus::Solved) {
#ifdef TIMINGS_SPARSE_CHAIN_SOLVER
    os << "|SPARSE_CHAIN: SparseSolutionMat_L1Exact total|=" << time_solve << "\n";
#endif
    return res;
  }
  if (known_consistent) {
    // Consistent but the L1 support does not carry a solution: the
    // rational elimination, within its budget.
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
    os << "SPARSE_CHAIN: SparseSolutionMat_L1Exact, exact elimination fallback\n";
#endif
    return SparseSolutionMat_Exact_Budget(SpMat, Bvect, budget.max_entries_exact, os);
  }
#ifdef DEBUG_SPARSE_CHAIN_SOLVER
  os << "SPARSE_CHAIN: SparseSolutionMat_L1Exact undecided within the budget\n";
#endif
  return {SparseSolveStatus::Aborted, empty};
}

// clang-format off
#endif  // SRC_SPARSE_SOLVER_SPARSECHAINSOLVER_H_
// clang-format on
