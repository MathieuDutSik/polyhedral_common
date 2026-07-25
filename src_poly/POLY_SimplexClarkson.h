// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_SIMPLEXCLARKSON_H_
#define SRC_POLY_POLY_SIMPLEXCLARKSON_H_

// clang-format off
#include "POLY_LinearProgrammingFund.h"
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SIMPLEX_CLARKSON
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SIMPLEX_CLARKSON
#endif

#ifdef TIMINGS
#define TIMINGS_SIMPLEX_CLARKSON
#endif

/*
  Fraction-free simplex algorithm for exact linear programming.

  This is a from-scratch replacement for the linear programming core of
  POLY_cddlib.h. The eventual goal of this file is:
  --- A fraction-free simplex LP solver (this code).
  --- A floating-point first phase whose terminal basis is certified (and
      if needed repaired) by the exact code.
  --- A Clarkson redundancy elimination built on top of it with warm starts.

  The LP being solved is expressed in the conventions used throughout
  polyhedral_common:
  --- ListIneq is a m x n matrix whose rows encode the inequalities
      f_i(x) = ListIneq(i,0) + sum_j ListIneq(i,j) x_j >= 0.
  --- ToBeMinimized is a vector of length n encoding the objective
      g(x) = ToBeMinimized(0) + sum_j ToBeMinimized(j) x_j to be minimized.
  The returned LpSolution<T> follows the conventions of CDD_LinearProgramming:
  --- Optimal: DirectSolution = optimal point x (length n-1),
      DualSolution = multipliers lambda <= 0 (length m) with
      c_red + sum_i lambda_i a_i = 0 and
      OptimalValue = c_0 + sum_i lambda_i b_i,
      OptimalValue set.
  --- Infeasible: only DualSolution is set. It is a Farkas certificate:
      lambda <= 0, sum_i lambda_i a_i = 0 and sum_i lambda_i b_i > 0,
      that is sum_i (-lambda_i) f_i(x) is identically negative.
  --- Unbounded: only DirectSolution is set. It is a ray d (length n-1)
      with a_i . d >= 0 for all rows and c_red . d < 0.
  Here b_i = ListIneq(i,0), a_i = ListIneq(i,1..) and c_0, c_red similarly
  for the objective.

  The solver is a dictionary simplex method using integer pivoting
  (Edmonds' Q-pivoting, also known as fraction-free or Bareiss pivoting,
  the same scheme as used in lrslib). The dictionary is stored as a matrix
  Tab together with a common denominator det > 0. Row i < m of Tab encodes
  the basic variable of that row as
     basic_i = (1/det) * (Tab(i,0) + sum_j Tab(i,j) u_j)
  over the cobasic variables u_j owning the columns j >= 1. Row m encodes
  the objective in the same way. The pivot on (r,s) updates
     Tab(i,j) <- (Tab(i,j) * Tab(r,s) - Tab(i,s) * Tab(r,j)) / det
  for i != r, j != s, then Tab(r,j) <- -Tab(r,j) for j != s and
  Tab(r,s) <- det, det <- Tab(r,s)_old. By Sylvester's determinant identity
  all divisions are exact: every entry stays equal to a minor of the input
  matrix. For integral input the whole computation is division-free in
  spirit (only exact divisions occur), avoiding the gcd cascades of naive
  rational arithmetic, and the entry sizes are bounded by minors (Hadamard
  bound) instead of blowing up.

  Phases:
  --- Phase 0 pivots the free variables x_j out of the cobasis, choosing for
      each x-column a row whose slack is still basic. Columns that no
      inequality constrains are lineality directions: if the objective is
      nonzero on such a column the problem is unbounded (reported only after
      feasibility has been established), otherwise the column is dropped.
  --- Phase 1 achieves primal feasibility with a single auxiliary variable
      x0 added to the violated rows, minimized by the simplex method. A
      positive minimum yields the Farkas certificate of infeasibility.
  --- Phase 2 is the primal simplex on the objective row.
  Pivoting uses the Dantzig rule and switches permanently to Bland's rule
  after a run of degenerate pivots, which guarantees termination.
 */

template <typename T> struct FractionFreeSimplex {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  // Variable identifiers: 0 is the phase-1 auxiliary variable, 1..m are the
  // slack variables of the rows, m+1..m+n_x are the x variables. Bland's
  // rule uses this numbering, with the auxiliary variable smallest so that
  // it is always preferred to leave the basis on ties.
  int m;
  int n_x;
  int ncol;
  MyMatrix<T> Tab;
  T det;
  std::vector<int> row_owner;
  // col_owner: -2 for the constant column 0, -1 for a dead column,
  // otherwise the variable id owning the column.
  std::vector<int> col_owner;
  bool use_bland;
  int degen_count;
  int degen_threshold;
  MyMatrix<T> const &ListIneq_;
  MyVector<T> const &obj_;
  [[maybe_unused]] std::ostream &os_;

  FractionFreeSimplex(MyMatrix<T> const &ListIneq,
                      MyVector<T> const &ToBeMinimized, std::ostream &os)
      : m(ListIneq.rows()), n_x(ListIneq.cols() - 1), ncol(n_x + 2),
        Tab(m + 1, ncol), det(1), row_owner(m), col_owner(ncol),
        use_bland(false), degen_count(0), degen_threshold(20 + 4 * (m + n_x)),
        ListIneq_(ListIneq), obj_(ToBeMinimized), os_(os) {
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    if (ToBeMinimized.size() != ListIneq.cols()) {
      std::cerr << "SIMPLEX_CLARKSON: |ToBeMinimized|=" << ToBeMinimized.size()
                << " but ListIneq has " << ListIneq.cols() << " columns\n";
      throw TerminalException{1};
    }
#endif
    for (int i = 0; i < m; i++) {
      for (int j = 0; j <= n_x; j++) {
        Tab(i, j) = ListIneq(i, j);
      }
      Tab(i, ncol - 1) = T(0);
      row_owner[i] = 1 + i;
    }
    for (int j = 0; j <= n_x; j++) {
      Tab(m, j) = ToBeMinimized(j);
    }
    Tab(m, ncol - 1) = T(0);
    col_owner[0] = -2;
    for (int k = 0; k < n_x; k++) {
      col_owner[1 + k] = 1 + m + k;
    }
    col_owner[ncol - 1] = -1;
  }

  bool is_slack(int v) const { return v >= 1 && v <= m; }
  bool is_x(int v) const { return v > m; }

  void DoPivot(int r, int s) {
    T piv = Tab(r, s);
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    if (piv == 0 || r < 0 || r >= m || s < 1 || s >= ncol ||
        col_owner[s] == -1) {
      std::cerr << "SIMPLEX_CLARKSON: invalid pivot r=" << r << " s=" << s
                << "\n";
      throw TerminalException{1};
    }
#endif
    for (int i = 0; i <= m; i++) {
      if (i == r) {
        continue;
      }
      T const &c_is = Tab(i, s);
      for (int j = 0; j < ncol; j++) {
        if (j == s || (j > 0 && col_owner[j] == -1)) {
          continue;
        }
        Tab(i, j) = (Tab(i, j) * piv - c_is * Tab(r, j)) / det;
      }
    }
    for (int j = 0; j < ncol; j++) {
      if (j == s || (j > 0 && col_owner[j] == -1)) {
        continue;
      }
      Tab(r, j) = -Tab(r, j);
    }
    Tab(r, s) = det;
    det = piv;
    std::swap(row_owner[r], col_owner[s]);
    if (det < 0) {
      det = -det;
      for (int i = 0; i <= m; i++) {
        for (int j = 0; j < ncol; j++) {
          if (j > 0 && col_owner[j] == -1) {
            continue;
          }
          Tab(i, j) = -Tab(i, j);
        }
      }
    }
  }

  // Entering column for minimizing the basic variable of row crow: a live
  // cobasic column with negative coefficient. Dantzig rule (most negative)
  // by default, Bland rule (smallest variable id) once use_bland is set.
  int SelectEntering(int crow) const {
    int best = -1;
    for (int j = 1; j < ncol; j++) {
      if (col_owner[j] < 1) {
        continue;
      }
      if (!(Tab(crow, j) < 0)) {
        continue;
      }
      if (best == -1) {
        best = j;
        continue;
      }
      if (use_bland) {
        if (col_owner[j] < col_owner[best]) {
          best = j;
        }
      } else {
        if (Tab(crow, j) < Tab(crow, best) ||
            (Tab(crow, j) == Tab(crow, best) &&
             col_owner[j] < col_owner[best])) {
          best = j;
        }
      }
    }
    return best;
  }

  // Leaving row for entering column s: among rows whose basic variable is
  // sign constrained (slack or auxiliary) and whose coefficient is negative,
  // the one of minimum ratio Tab(i,0) / (-Tab(i,s)), ties broken by the
  // smallest basic variable id (Bland).
  int RatioTest(int s) const {
    int best = -1;
    for (int i = 0; i < m; i++) {
      if (row_owner[i] > m) {
        continue;
      }
      if (!(Tab(i, s) < 0)) {
        continue;
      }
      if (best == -1) {
        best = i;
        continue;
      }
      T lhs = Tab(i, 0) * (-Tab(best, s));
      T rhs = Tab(best, 0) * (-Tab(i, s));
      if (lhs < rhs || (lhs == rhs && row_owner[i] < row_owner[best])) {
        best = i;
      }
    }
    return best;
  }

  void TrackDegeneracy(int r) {
    if (Tab(r, 0) == 0) {
      degen_count++;
      if (!use_bland && degen_count > degen_threshold) {
        use_bland = true;
#ifdef DEBUG_SIMPLEX_CLARKSON
        os_ << "SIMPLEX_CLARKSON: switching to Bland rule after "
            << degen_count << " degenerate pivots\n";
#endif
      }
    } else {
      degen_count = 0;
    }
  }

  // A ray of the feasible region obtained by moving along the cobasic
  // variable of column s, expressed in the original x coordinates and
  // scaled by det > 0. sigma is the direction of movement.
  MyVector<T> GetRayOfColumn(int s, int sigma) const {
    MyVector<T> d = ZeroVector<T>(n_x);
    T sigma_T(sigma);
    for (int i = 0; i < m; i++) {
      if (is_x(row_owner[i])) {
        d(row_owner[i] - 1 - m) = sigma_T * Tab(i, s);
      }
    }
    if (is_x(col_owner[s])) {
      d(col_owner[s] - 1 - m) = sigma_T * det;
    }
    return d;
  }

  bool NeedPhase1() const {
    for (int i = 0; i < m; i++) {
      if (is_slack(row_owner[i]) && Tab(i, 0) < 0) {
        return true;
      }
    }
    return false;
  }

  // Insert the auxiliary variable in the violated rows and perform the
  // special pivot that makes the dictionary feasible for the auxiliary
  // problem. Returns the row where the auxiliary variable is basic.
  int ActivatePhase1() {
    int aux_col = ncol - 1;
    int r_min = -1;
    for (int i = 0; i < m; i++) {
      bool viol = is_slack(row_owner[i]) && Tab(i, 0) < 0;
      Tab(i, aux_col) = viol ? det : T(0);
      if (viol && (r_min == -1 || Tab(i, 0) < Tab(r_min, 0))) {
        r_min = i;
      }
    }
    Tab(m, aux_col) = T(0);
    col_owner[aux_col] = 0;
    DoPivot(r_min, aux_col);
    return r_min;
  }

  // Phase 1: minimize the auxiliary variable. Returns the infeasibility
  // solution if the minimum is positive, and an empty optional when primal
  // feasibility has been reached.
  std::optional<LpSolution<T>> RunPhase1() {
    if (!NeedPhase1()) {
      return {};
    }
    int raux = ActivatePhase1();
#ifdef DEBUG_SIMPLEX_CLARKSON
    os_ << "SIMPLEX_CLARKSON: phase 1 started, raux=" << raux << "\n";
#endif
    while (true) {
      int s = SelectEntering(raux);
      if (s == -1) {
        if (Tab(raux, 0) == 0) {
          // The auxiliary variable is basic at value zero: pivot it out on
          // any live column and drop it, the dictionary is feasible.
          int s2 = -1;
          for (int j = 1; j < ncol; j++) {
            if (col_owner[j] >= 1 && Tab(raux, j) != 0) {
              s2 = j;
              break;
            }
          }
          if (s2 == -1) {
            std::cerr << "SIMPLEX_CLARKSON: the auxiliary variable row is "
                         "identically zero which is impossible\n";
            throw TerminalException{1};
          }
          DoPivot(raux, s2);
          col_owner[s2] = -1;
          return {};
        }
        return BuildInfeasible(raux);
      }
      int r = RatioTest(s);
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
      if (r == -1) {
        std::cerr << "SIMPLEX_CLARKSON: phase 1 ratio test failed, the "
                     "auxiliary row should always block\n";
        throw TerminalException{1};
      }
#endif
      TrackDegeneracy(r);
      DoPivot(r, s);
      if (r == raux) {
        // The auxiliary variable left the basis and owns column s now.
        col_owner[s] = -1;
        return {};
      }
    }
  }

  // Phase 2: minimize the objective row over the feasible dictionary.
  LpSolution<T> RunPhase2() {
#ifdef DEBUG_SIMPLEX_CLARKSON
    size_t n_pivot = 0;
#endif
    while (true) {
      int s = SelectEntering(m);
      if (s == -1) {
#ifdef DEBUG_SIMPLEX_CLARKSON
        os_ << "SIMPLEX_CLARKSON: phase 2 finished optimal after " << n_pivot
            << " pivots\n";
#endif
        return BuildOptimal();
      }
      int r = RatioTest(s);
      if (r == -1) {
#ifdef DEBUG_SIMPLEX_CLARKSON
        os_ << "SIMPLEX_CLARKSON: phase 2 finished unbounded after " << n_pivot
            << " pivots\n";
#endif
        return BuildUnbounded(GetRayOfColumn(s, 1));
      }
      TrackDegeneracy(r);
      DoPivot(r, s);
#ifdef DEBUG_SIMPLEX_CLARKSON
      n_pivot++;
#endif
    }
  }

  LpSolution<T> BuildOptimal() const {
    LpSolution<T> sol;
    sol.OptimalValue = Tab(m, 0) / det;
    MyVector<T> DirectSolution = ZeroVector<T>(n_x);
    for (int i = 0; i < m; i++) {
      if (is_x(row_owner[i])) {
        DirectSolution(row_owner[i] - 1 - m) = Tab(i, 0) / det;
      }
    }
    MyVector<T> DualSolution = ZeroVector<T>(m);
    for (int j = 1; j < ncol; j++) {
      if (is_slack(col_owner[j])) {
        DualSolution(col_owner[j] - 1) = -Tab(m, j) / det;
      }
    }
    sol.DirectSolution = DirectSolution;
    sol.DualSolution = DualSolution;
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    for (int i = 0; i < m; i++) {
      T eSum = ListIneq_(i, 0);
      for (int k = 0; k < n_x; k++) {
        eSum += ListIneq_(i, k + 1) * DirectSolution(k);
      }
      if (eSum < 0) {
        std::cerr << "SIMPLEX_CLARKSON: optimal point violates row " << i
                  << "\n";
        throw TerminalException{1};
      }
    }
    if (!CheckDualSolutionGetOptimal(ListIneq_, obj_, sol, os_)) {
      std::cerr << "SIMPLEX_CLARKSON: the dual solution certificate is "
                   "invalid\n";
      throw TerminalException{1};
    }
#endif
    return sol;
  }

  LpSolution<T> BuildInfeasible(int raux) const {
    LpSolution<T> sol;
    MyVector<T> DualSolution = ZeroVector<T>(m);
    for (int j = 1; j < ncol; j++) {
      if (is_slack(col_owner[j])) {
        DualSolution(col_owner[j] - 1) = -Tab(raux, j) / det;
      }
    }
    sol.DualSolution = DualSolution;
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    T eConst(0);
    MyVector<T> eSum = ZeroVector<T>(n_x);
    for (int i = 0; i < m; i++) {
      if (DualSolution(i) > 0) {
        std::cerr << "SIMPLEX_CLARKSON: Farkas certificate has a positive "
                     "entry at row "
                  << i << "\n";
        throw TerminalException{1};
      }
      eConst += DualSolution(i) * ListIneq_(i, 0);
      for (int k = 0; k < n_x; k++) {
        eSum(k) += DualSolution(i) * ListIneq_(i, k + 1);
      }
    }
    if (!(eConst > 0)) {
      std::cerr << "SIMPLEX_CLARKSON: Farkas certificate has nonpositive "
                   "constant\n";
      throw TerminalException{1};
    }
    for (int k = 0; k < n_x; k++) {
      if (eSum(k) != 0) {
        std::cerr << "SIMPLEX_CLARKSON: Farkas certificate does not sum to "
                     "zero at coordinate "
                  << k << "\n";
        throw TerminalException{1};
      }
    }
#endif
    return sol;
  }

  LpSolution<T> BuildUnbounded(MyVector<T> const &d) const {
    LpSolution<T> sol;
    sol.DirectSolution = d;
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    for (int i = 0; i < m; i++) {
      T eSum(0);
      for (int k = 0; k < n_x; k++) {
        eSum += ListIneq_(i, k + 1) * d(k);
      }
      if (eSum < 0) {
        std::cerr << "SIMPLEX_CLARKSON: unbounded ray violates row " << i
                  << "\n";
        throw TerminalException{1};
      }
    }
    T eObj(0);
    for (int k = 0; k < n_x; k++) {
      eObj += obj_(k + 1) * d(k);
    }
    if (!(eObj < 0)) {
      std::cerr << "SIMPLEX_CLARKSON: unbounded ray does not decrease the "
                   "objective\n";
      throw TerminalException{1};
    }
#endif
    return sol;
  }

  LpSolution<T> solve() {
    // Phase 0: pivot the free variables out of the cobasis.
    std::optional<MyVector<T>> pending_ray;
    for (int j = 1; j <= n_x; j++) {
      int r_found = -1;
      for (int i = 0; i < m; i++) {
        if (is_slack(row_owner[i]) && Tab(i, j) != 0) {
          r_found = i;
          break;
        }
      }
      if (r_found >= 0) {
        DoPivot(r_found, j);
        continue;
      }
      // No inequality constrains this direction: a lineality direction.
      if (Tab(m, j) != 0 && !pending_ray) {
        int sigma = (Tab(m, j) > 0) ? -1 : 1;
        pending_ray = GetRayOfColumn(j, sigma);
      }
      col_owner[j] = -1;
    }
#ifdef DEBUG_SIMPLEX_CLARKSON
    os_ << "SIMPLEX_CLARKSON: phase 0 done, det=" << det << "\n";
#endif
    // Phase 1: primal feasibility. Infeasibility takes precedence over the
    // unboundedness detected in phase 0.
    std::optional<LpSolution<T>> infeas = RunPhase1();
    if (infeas) {
#ifdef DEBUG_SIMPLEX_CLARKSON
      os_ << "SIMPLEX_CLARKSON: the problem is infeasible\n";
#endif
      return *infeas;
    }
    if (pending_ray) {
#ifdef DEBUG_SIMPLEX_CLARKSON
      os_ << "SIMPLEX_CLARKSON: unbounded via a lineality direction\n";
#endif
      return BuildUnbounded(*pending_ray);
    }
    return RunPhase2();
  }
};

// Solve the linear program of minimizing ToBeMinimized . (1,x) over the
// polyhedron of x with ListIneq . (1,x) >= 0, exactly, with the fraction
// free simplex method. See the top of the file for the conventions of the
// returned LpSolution<T>.
template <typename T>
LpSolution<T> SIMPLEX_LinearProgramming(MyMatrix<T> const &ListIneq,
                                        MyVector<T> const &ToBeMinimized,
                                        std::ostream &os) {
#ifdef TIMINGS_SIMPLEX_CLARKSON
  MicrosecondTime time;
#endif
  FractionFreeSimplex<T> solver(ListIneq, ToBeMinimized, os);
  LpSolution<T> sol = solver.solve();
#ifdef TIMINGS_SIMPLEX_CLARKSON
  os << "|SIMPLEX_CLARKSON: SIMPLEX_LinearProgramming|=" << time << "\n";
#endif
  return sol;
}

// clang-format off
#endif  // SRC_POLY_POLY_SIMPLEXCLARKSON_H_
// clang-format on
