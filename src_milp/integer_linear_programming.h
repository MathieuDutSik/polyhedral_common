// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_INTEGER_LINEAR_PROGRAMMING_H_
#define SRC_MILP_INTEGER_LINEAR_PROGRAMMING_H_

// clang-format off
#include "ilp_fundamental.h"
#include "POLY_LinearProgramming.h"
#include <queue>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#ifdef ENABLE_SCIP_SUPPORT
#include "scip_support.h"
#endif
// clang-format on

#ifdef DEBUG
#define DEBUG_INTEGER_LINEAR_PROGRAMMING
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INTEGER_LINEAR_PROGRAMMING
#endif

#ifdef TIMINGS
#define TIMINGS_INTEGER_LINEAR_PROGRAMMING
#endif

/*
  Integer linear programming.

  The problem is expressed in the conventions of SIMPLEX_LinearProgramming:
  --- ListIneq is a m x (d+1) matrix whose rows encode the inequalities
      ListIneq(i,0) + sum_j ListIneq(i,j) x_j >= 0.
  --- ListEqua is a p x (d+1) matrix whose rows encode the equalities
      ListEqua(i,0) + sum_j ListEqua(i,j) x_j = 0.
  --- ToBeMinimized is a vector of length d+1 encoding the objective
      ToBeMinimized(0) + sum_j ToBeMinimized(j) x_j.
  and the variables x_1, ..., x_d are integers without bounds.

  Two methods are available:
  --- "scip": the SCIP solver of scip_support.h. Only available when
      compiled with ENABLE_SCIP_SUPPORT.
  --- "exact_bb": a best-first branch and bound over the exact simplex
      method of POLY_SimplexClarkson.h. It is slow but exact and needs no
      external library. It is a fallback and a check of the SCIP results.
  --- "default": scip when compiled in, exact_bb otherwise.

  The linear relaxation is expected to be bounded: the callers ensure it
  by adding cutting planes before calling. An unbounded relaxation is
  reported as such without trying to decide the integer problem.
 */

/*
  The exact branch and bound. The nodes are the linear relaxations with
  some bound rows x_j <= k or x_j >= k added. The node of smallest
  relaxation value is processed first, so the first integral relaxation
  solution met is optimal.
 */
template <typename T>
IlpSolution<T> ILP_ExactBranchAndBound(MyMatrix<T> const &ListIneq,
                                       MyMatrix<T> const &ListEqua,
                                       MyVector<T> const &ToBeMinimized,
                                       std::ostream &os) {
  int n_col = ToBeMinimized.size();
  int d = n_col - 1;
  int n_ineq = ListIneq.rows();
  int n_equa = ListEqua.rows();
  MyMatrix<T> BaseIneq(n_ineq + 2 * n_equa, n_col);
  for (int i = 0; i < n_ineq; i++) {
    for (int j = 0; j < n_col; j++) {
      BaseIneq(i, j) = ListIneq(i, j);
    }
  }
  for (int i = 0; i < n_equa; i++) {
    for (int j = 0; j < n_col; j++) {
      BaseIneq(n_ineq + 2 * i, j) = ListEqua(i, j);
      BaseIneq(n_ineq + 2 * i + 1, j) = -ListEqua(i, j);
    }
  }
  // A bound row: (j, value, is_upper) meaning x_j <= value or x_j >= value
  using Tbound = std::tuple<int, T, bool>;
  struct Node {
    T value;
    MyVector<T> x;
    std::vector<Tbound> bounds;
  };
  auto get_ineq = [&](std::vector<Tbound> const &bounds) -> MyMatrix<T> {
    int n_base = BaseIneq.rows();
    int n_bound = bounds.size();
    MyMatrix<T> M = ZeroMatrix<T>(n_base + n_bound, n_col);
    for (int i = 0; i < n_base; i++) {
      for (int j = 0; j < n_col; j++) {
        M(i, j) = BaseIneq(i, j);
      }
    }
    for (int i = 0; i < n_bound; i++) {
      auto const &[jvar, val, is_upper] = bounds[i];
      if (is_upper) {
        M(n_base + i, 0) = val;
        M(n_base + i, jvar + 1) = T(-1);
      } else {
        M(n_base + i, 0) = -val;
        M(n_base + i, jvar + 1) = T(1);
      }
    }
    return M;
  };
  auto cmp = [](Node const &a, Node const &b) -> bool {
    return a.value > b.value;
  };
  std::priority_queue<Node, std::vector<Node>, decltype(cmp)> queue(cmp);
  // Returns true if the relaxation is unbounded.
  auto insert_node = [&](std::vector<Tbound> const &bounds) -> bool {
    MyMatrix<T> M = get_ineq(bounds);
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(M, ToBeMinimized, os);
    if (eSol.DirectSolution && !eSol.DualSolution) {
      return true;
    }
    if (!eSol.DirectSolution) {
      return false;
    }
    Node node{eSol.OptimalValue, *eSol.DirectSolution, bounds};
    queue.push(node);
    return false;
  };
  if (insert_node({})) {
    return {IlpStatus::Unbounded, {}, T(0)};
  }
  size_t n_node = 0;
  size_t max_node = 1000000;
  while (!queue.empty()) {
    Node node = queue.top();
    queue.pop();
    n_node++;
    if (n_node > max_node) {
      std::cerr << "ILP: ILP_ExactBranchAndBound exceeded " << max_node
                << " nodes\n";
      std::cerr << "ILP: Use the scip method for this problem\n";
      throw TerminalException{1};
    }
    int j_frac = -1;
    for (int j = 0; j < d; j++) {
      if (!IsInteger(node.x(j))) {
        j_frac = j;
        break;
      }
    }
    if (j_frac == -1) {
#ifdef DEBUG_INTEGER_LINEAR_PROGRAMMING
      os << "ILP: ILP_ExactBranchAndBound, optimal after n_node=" << n_node
         << "\n";
#endif
      return {IlpStatus::Optimal, node.x, node.value};
    }
    T val_floor = UniversalFloorScalarInteger<T, T>(node.x(j_frac));
    std::vector<Tbound> bounds_down = node.bounds;
    bounds_down.push_back({j_frac, val_floor, true});
    std::vector<Tbound> bounds_up = node.bounds;
    bounds_up.push_back({j_frac, val_floor + 1, false});
    // Bounded at the root means bounded at every node.
    (void)insert_node(bounds_down);
    (void)insert_node(bounds_up);
  }
  return {IlpStatus::Infeasible, {}, T(0)};
}

inline std::string ILP_DefaultMethod() {
#ifdef ENABLE_SCIP_SUPPORT
  return "scip";
#else
  return "exact_bb";
#endif
}

template <typename T>
IlpSolution<T> ILP_IntegerLinearProgramming(MyMatrix<T> const &ListIneq,
                                            MyMatrix<T> const &ListEqua,
                                            MyVector<T> const &ToBeMinimized,
                                            std::string const &method,
                                            std::ostream &os) {
#ifdef TIMINGS_INTEGER_LINEAR_PROGRAMMING
  MicrosecondTime time;
#endif
  std::string method_eff = method;
  if (method_eff == "default") {
    method_eff = ILP_DefaultMethod();
  }
  auto get_solution = [&]() -> IlpSolution<T> {
    if (method_eff == "exact_bb") {
      return ILP_ExactBranchAndBound(ListIneq, ListEqua, ToBeMinimized, os);
    }
    if (method_eff == "scip") {
#ifdef ENABLE_SCIP_SUPPORT
      return ILP_Scip(ListIneq, ListEqua, ToBeMinimized, os);
#else
      std::cerr << "ILP: The scip method requires compiling with "
                   "ENABLE_SCIP_SUPPORT=1\n";
      throw TerminalException{1};
#endif
    }
    std::cerr << "ILP: Unknown method=" << method << "\n";
    std::cerr << "ILP: Allowed methods: default, scip, exact_bb\n";
    throw TerminalException{1};
  };
  IlpSolution<T> sol = get_solution();
#ifdef TIMINGS_INTEGER_LINEAR_PROGRAMMING
  os << "|ILP: ILP_IntegerLinearProgramming(" << method_eff << ")|=" << time
     << "\n";
#endif
#ifdef SANITY_CHECK_INTEGER_LINEAR_PROGRAMMING
  if (sol.status == IlpStatus::Optimal) {
    if (!IsIntegralVector(sol.solution) ||
        !ILP_IsFeasible(ListIneq, ListEqua, sol.solution)) {
      std::cerr << "ILP: The solution is not an integral feasible point\n";
      throw TerminalException{1};
    }
  }
  if (method_eff == "scip") {
    IlpSolution<T> sol_exact =
        ILP_ExactBranchAndBound(ListIneq, ListEqua, ToBeMinimized, os);
    bool is_coherent = sol_exact.status == sol.status;
    if (is_coherent && sol.status == IlpStatus::Optimal) {
      is_coherent = sol_exact.OptimalValue == sol.OptimalValue;
    }
    if (!is_coherent) {
      std::cerr << "ILP: SCIP and the exact branch and bound disagree\n";
      throw TerminalException{1};
    }
  }
#endif
  return sol;
}

// clang-format off
#endif  // SRC_MILP_INTEGER_LINEAR_PROGRAMMING_H_
// clang-format on
