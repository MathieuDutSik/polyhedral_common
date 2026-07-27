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

  // Perform phase-0 style pivots so that the slack variables of the listed
  // rows become cobasic. This is used to warm start the solver from a basis
  // obtained by other means (typically the floating-point solver). The rows
  // are only a hint: rows that cannot be pivoted are silently skipped and
  // the normal phases repair whatever is left.
  void ApplyBasisHint(std::vector<int> const &rows) {
    for (int r : rows) {
      if (r < 0 || r >= m || !is_slack(row_owner[r])) {
        continue;
      }
      for (int j = 1; j <= n_x; j++) {
        if (is_x(col_owner[j]) && Tab(r, j) != 0) {
          DoPivot(r, j);
          break;
        }
      }
    }
  }

  LpSolution<T> solve() {
    // Phase 0: pivot the free variables out of the cobasis.
    std::optional<MyVector<T>> pending_ray;
    for (int j = 1; j <= n_x; j++) {
      if (!is_x(col_owner[j])) {
        continue;
      }
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
LpSolution<T> SIMPLEX_LinearProgramming_exact(MyMatrix<T> const &ListIneq,
                                              MyVector<T> const &ToBeMinimized,
                                              std::ostream &os) {
#ifdef TIMINGS_SIMPLEX_CLARKSON
  MicrosecondTime time;
#endif
  FractionFreeSimplex<T> solver(ListIneq, ToBeMinimized, os);
  LpSolution<T> sol = solver.solve();
#ifdef TIMINGS_SIMPLEX_CLARKSON
  os << "|SIMPLEX_CLARKSON: SIMPLEX_LinearProgramming_exact|=" << time << "\n";
#endif
  return sol;
}

enum class FloatLpStatus { Optimal, Infeasible, Unbounded, Unreliable };

// The outcome of the floating-point solve: a status and the cobasis at
// termination (the original indices of the rows whose slack variables are
// cobasic). Both are only hints; every conclusion drawn from them is
// certified in exact arithmetic before being returned to the caller.
struct FloatBasisHint {
  FloatLpStatus status;
  std::vector<int> cobasis_rows;
};

/*
  Floating-point version of the dictionary simplex above. It follows the
  same three phases but with the classical normalized pivot update and
  epsilon tolerances instead of exact sign tests. Its only purpose is to
  produce a candidate optimal basis fast; correctness is never assumed.
  Failure modes (cycling, tiny pivots, wrong rank decisions) are handled by
  an iteration cap leading to the Unreliable status, and by the exact
  certification and repair done by the caller.
 */
template <typename Tfloat> struct FloatSimplex {
  int m;
  int n_x;
  int ncol;
  MyMatrix<Tfloat> Tab;
  std::vector<int> row_owner;
  std::vector<int> col_owner;
  bool use_bland;
  int degen_count;
  int degen_threshold;
  size_t iter_count;
  size_t max_iter;
  Tfloat eps_pivot;
  Tfloat eps_cost;
  Tfloat eps_feas;

  // The matrix rows are expected to be prescaled to entries of order 1.
  FloatSimplex(MyMatrix<Tfloat> const &ListIneq,
               MyVector<Tfloat> const &ToBeMinimized)
      : m(ListIneq.rows()), n_x(ListIneq.cols() - 1), ncol(n_x + 2),
        Tab(m + 1, ncol), row_owner(m), col_owner(ncol), use_bland(false),
        degen_count(0), degen_threshold(20 + 4 * (m + n_x)), iter_count(0),
        max_iter(100 * static_cast<size_t>(m + n_x) + 1000),
        eps_pivot(1e-11), eps_cost(1e-9), eps_feas(1e-9) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j <= n_x; j++) {
        Tab(i, j) = ListIneq(i, j);
      }
      Tab(i, ncol - 1) = 0;
      row_owner[i] = 1 + i;
    }
    for (int j = 0; j <= n_x; j++) {
      Tab(m, j) = ToBeMinimized(j);
    }
    Tab(m, ncol - 1) = 0;
    col_owner[0] = -2;
    for (int k = 0; k < n_x; k++) {
      col_owner[1 + k] = 1 + m + k;
    }
    col_owner[ncol - 1] = -1;
  }

  bool is_slack(int v) const { return v >= 1 && v <= m; }
  bool is_x(int v) const { return v > m; }

  std::vector<int> GetCobasisRows() const {
    std::vector<int> rows;
    for (int j = 1; j < ncol; j++) {
      if (is_slack(col_owner[j])) {
        rows.push_back(col_owner[j] - 1);
      }
    }
    return rows;
  }

  FloatBasisHint Result(FloatLpStatus status) const {
    return {status, GetCobasisRows()};
  }

  void DoPivot(int r, int s) {
    iter_count++;
    Tfloat inv = 1 / Tab(r, s);
    for (int j = 0; j < ncol; j++) {
      if (j == s || (j > 0 && col_owner[j] == -1)) {
        continue;
      }
      Tab(r, j) = -Tab(r, j) * inv;
    }
    Tab(r, s) = inv;
    for (int i = 0; i <= m; i++) {
      if (i == r) {
        continue;
      }
      Tfloat c_is = Tab(i, s);
      if (c_is != 0) {
        for (int j = 0; j < ncol; j++) {
          if (j == s || (j > 0 && col_owner[j] == -1)) {
            continue;
          }
          Tab(i, j) += c_is * Tab(r, j);
        }
        Tab(i, s) = c_is * inv;
      }
    }
    std::swap(row_owner[r], col_owner[s]);
  }

  int SelectEntering(int crow) const {
    int best = -1;
    for (int j = 1; j < ncol; j++) {
      if (col_owner[j] < 1) {
        continue;
      }
      if (!(Tab(crow, j) < -eps_cost)) {
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
        if (Tab(crow, j) < Tab(crow, best)) {
          best = j;
        }
      }
    }
    return best;
  }

  // Ratio test with a stability tie-break: among rows whose ratio is within
  // a small window of the minimum, prefer the largest pivot magnitude.
  int RatioTest(int s) const {
    int best = -1;
    Tfloat best_ratio = 0;
    for (int i = 0; i < m; i++) {
      if (row_owner[i] > m) {
        continue;
      }
      Tfloat den = -Tab(i, s);
      if (!(den > eps_pivot)) {
        continue;
      }
      Tfloat num = Tab(i, 0);
      if (num < 0) {
        num = 0;
      }
      Tfloat ratio = num / den;
      if (best == -1) {
        best = i;
        best_ratio = ratio;
        continue;
      }
      Tfloat window = eps_feas * (1 + best_ratio);
      if (ratio < best_ratio - window) {
        best = i;
        best_ratio = ratio;
      } else if (ratio < best_ratio + window) {
        if (use_bland) {
          if (row_owner[i] < row_owner[best]) {
            best = i;
            if (ratio < best_ratio) {
              best_ratio = ratio;
            }
          }
        } else {
          if (-Tab(i, s) > -Tab(best, s)) {
            best = i;
            if (ratio < best_ratio) {
              best_ratio = ratio;
            }
          }
        }
      }
    }
    return best;
  }

  void TrackDegeneracy(int r) {
    if (Tab(r, 0) < eps_feas) {
      degen_count++;
      if (!use_bland && degen_count > degen_threshold) {
        use_bland = true;
      }
    } else {
      degen_count = 0;
    }
  }

  FloatBasisHint solve() {
    // Phase 0 with partial pivoting on the free variable columns.
    bool has_lineality_obj = false;
    for (int j = 1; j <= n_x; j++) {
      if (!is_x(col_owner[j])) {
        continue;
      }
      int r_found = -1;
      Tfloat best_mag = eps_pivot;
      for (int i = 0; i < m; i++) {
        if (!is_slack(row_owner[i])) {
          continue;
        }
        Tfloat mag = Tab(i, j) < 0 ? -Tab(i, j) : Tab(i, j);
        if (mag > best_mag) {
          best_mag = mag;
          r_found = i;
        }
      }
      if (r_found >= 0) {
        DoPivot(r_found, j);
        continue;
      }
      Tfloat obj_mag = Tab(m, j) < 0 ? -Tab(m, j) : Tab(m, j);
      if (obj_mag > eps_cost) {
        has_lineality_obj = true;
      }
      col_owner[j] = -1;
    }
    // Phase 1
    bool need_phase1 = false;
    for (int i = 0; i < m; i++) {
      if (is_slack(row_owner[i]) && Tab(i, 0) < -eps_feas) {
        need_phase1 = true;
        break;
      }
    }
    if (need_phase1) {
      int aux_col = ncol - 1;
      int r_min = -1;
      for (int i = 0; i < m; i++) {
        bool viol = is_slack(row_owner[i]) && Tab(i, 0) < -eps_feas;
        Tab(i, aux_col) = viol ? 1 : 0;
        if (viol && (r_min == -1 || Tab(i, 0) < Tab(r_min, 0))) {
          r_min = i;
        }
      }
      Tab(m, aux_col) = 0;
      col_owner[aux_col] = 0;
      DoPivot(r_min, aux_col);
      int raux = r_min;
      bool feasible = false;
      while (iter_count < max_iter) {
        int s = SelectEntering(raux);
        if (s == -1) {
          if (Tab(raux, 0) < eps_feas * 10) {
            int s2 = -1;
            Tfloat best_mag = eps_pivot;
            for (int j = 1; j < ncol; j++) {
              if (col_owner[j] < 1) {
                continue;
              }
              Tfloat mag = Tab(raux, j) < 0 ? -Tab(raux, j) : Tab(raux, j);
              if (mag > best_mag) {
                best_mag = mag;
                s2 = j;
              }
            }
            if (s2 == -1) {
              return Result(FloatLpStatus::Unreliable);
            }
            DoPivot(raux, s2);
            col_owner[s2] = -1;
            feasible = true;
            break;
          }
          return Result(FloatLpStatus::Infeasible);
        }
        int r = RatioTest(s);
        if (r == -1) {
          return Result(FloatLpStatus::Unreliable);
        }
        TrackDegeneracy(r);
        DoPivot(r, s);
        if (r == raux) {
          col_owner[s] = -1;
          feasible = true;
          break;
        }
      }
      if (!feasible) {
        return Result(FloatLpStatus::Unreliable);
      }
    }
    if (has_lineality_obj) {
      return Result(FloatLpStatus::Unbounded);
    }
    // Phase 2
    while (iter_count < max_iter) {
      int s = SelectEntering(m);
      if (s == -1) {
        return Result(FloatLpStatus::Optimal);
      }
      int r = RatioTest(s);
      if (r == -1) {
        return Result(FloatLpStatus::Unbounded);
      }
      TrackDegeneracy(r);
      DoPivot(r, s);
    }
    return Result(FloatLpStatus::Unreliable);
  }
};

// Convert the LP data to floating point with row scaling and solve it,
// returning the basis hint.
template <typename T, typename Tfloat>
FloatBasisHint ComputeFloatBasisHint(MyMatrix<T> const &ListIneq,
                                     MyVector<T> const &ToBeMinimized) {
  int m = ListIneq.rows();
  int n = ListIneq.cols();
  MyMatrix<Tfloat> A(m, n);
  for (int i = 0; i < m; i++) {
    Tfloat max_mag = 0;
    for (int j = 0; j < n; j++) {
      A(i, j) = UniversalScalarConversion<Tfloat, T>(ListIneq(i, j));
      Tfloat mag = A(i, j) < 0 ? -A(i, j) : A(i, j);
      if (mag > max_mag) {
        max_mag = mag;
      }
    }
    if (max_mag > 0) {
      for (int j = 0; j < n; j++) {
        A(i, j) /= max_mag;
      }
    }
  }
  MyVector<Tfloat> obj(n);
  Tfloat max_mag = 0;
  for (int j = 0; j < n; j++) {
    obj(j) = UniversalScalarConversion<Tfloat, T>(ToBeMinimized(j));
    Tfloat mag = obj(j) < 0 ? -obj(j) : obj(j);
    if (mag > max_mag) {
      max_mag = mag;
    }
  }
  if (max_mag > 0) {
    for (int j = 0; j < n; j++) {
      obj(j) /= max_mag;
    }
  }
  FloatSimplex<Tfloat> solver(A, obj);
  return solver.solve();
}

// Solve the square k x k system M x = rhs by fraction-free Bareiss
// elimination with row pivoting. Returns (xnum, den) with x = xnum / den
// and den = +-det(M), or an empty optional when M is singular. All the
// elimination arithmetic consists of exact divisions, so for integral
// input everything stays integral (no gcd cascades); fractions only appear
// in the O(k^2) back substitution.
template <typename T>
std::optional<std::pair<MyVector<T>, T>>
FractionFreeSolveSquare(MyMatrix<T> W, MyVector<T> rhs) {
  int k = W.rows();
  std::vector<int> pivrow(k);
  std::vector<uint8_t> used(k, 0);
  T prev(1);
  for (int t = 0; t < k; t++) {
    int r_sel = -1;
    for (int r = 0; r < k; r++) {
      if (!used[r] && W(r, t) != 0) {
        r_sel = r;
        break;
      }
    }
    if (r_sel == -1) {
      return {};
    }
    used[r_sel] = 1;
    pivrow[t] = r_sel;
    T const &piv = W(r_sel, t);
    for (int r = 0; r < k; r++) {
      if (used[r]) {
        continue;
      }
      for (int j = t + 1; j < k; j++) {
        W(r, j) = (W(r, j) * piv - W(r, t) * W(r_sel, j)) / prev;
      }
      rhs(r) = (rhs(r) * piv - W(r, t) * rhs(r_sel)) / prev;
      W(r, t) = T(0);
    }
    prev = piv;
  }
  // Rows pivrow[t] form a triangular system: back substitution.
  MyVector<T> x(k);
  for (int t = k - 1; t >= 0; t--) {
    int R = pivrow[t];
    T val = rhs(R);
    for (int j = t + 1; j < k; j++) {
      val -= W(R, j) * x(j);
    }
    x(t) = val / W(R, t);
  }
  // The last pivot is +-det(M), and by Cramer's rule det(M) * x is
  // integral for integral input.
  MyVector<T> xnum(k);
  for (int t = 0; t < k; t++) {
    xnum(t) = x(t) * prev;
  }
  return std::pair<MyVector<T>, T>(std::move(xnum), prev);
}

// Given a candidate optimal basis (a set of n_x row indices), reconstruct
// the exact vertex and the exact dual multipliers by two fraction-free
// linear solves and check the full optimality certificate: the vertex is
// feasible for all rows, the multipliers are nonpositive and supported on
// the basis, and the primal and dual objective values agree. By weak
// duality a success is a complete proof of optimality. Any failure returns
// an empty optional. The feasibility scan is done on denominator-cleared
// values so that for integral input it is pure integer arithmetic.
template <typename T>
std::optional<LpSolution<T>>
SIMPLEX_CertifyBasis(MyMatrix<T> const &ListIneq,
                     MyVector<T> const &ToBeMinimized,
                     std::vector<int> const &rows,
                     [[maybe_unused]] std::ostream &os) {
  int m = ListIneq.rows();
  int n_x = ListIneq.cols() - 1;
  if (static_cast<int>(rows.size()) != n_x) {
    return {};
  }
  // The vertex: x with b_r + a_r . x = 0 for every basis row r, that is
  // M x = -b_B with M(j, i) = a_{r_j, i}.
  MyMatrix<T> M(n_x, n_x);
  MyVector<T> rhs(n_x);
  std::vector<uint8_t> in_basis(m, 0);
  for (int j = 0; j < n_x; j++) {
    int r = rows[j];
    if (r < 0 || r >= m) {
      return {};
    }
    in_basis[r] = 1;
    rhs(j) = -ListIneq(r, 0);
    for (int i = 0; i < n_x; i++) {
      M(j, i) = ListIneq(r, i + 1);
    }
  }
  std::optional<std::pair<MyVector<T>, T>> optX =
      FractionFreeSolveSquare(M, rhs);
  if (!optX) {
    return {};
  }
  MyVector<T> const &xnum = optX->first;
  T const &den = optX->second;
  bool den_pos = den > 0;
  // Feasibility of the vertex: den * f_i(x) = b_i * den + a_i . xnum must
  // have the sign of den. The basis rows are tight by construction.
  for (int i = 0; i < m; i++) {
    if (in_basis[i]) {
      continue;
    }
    T eSum = ListIneq(i, 0) * den;
    for (int k = 0; k < n_x; k++) {
      eSum += ListIneq(i, k + 1) * xnum(k);
    }
    if (den_pos ? (eSum < 0) : (eSum > 0)) {
      return {};
    }
  }
  T primalNum = ToBeMinimized(0) * den;
  for (int k = 0; k < n_x; k++) {
    primalNum += ToBeMinimized(k + 1) * xnum(k);
  }
  T primalValue = primalNum / den;
  // The dual multipliers: pd with sum_j pd_j a_{r_j} = c_red, that is
  // M^T pd = c_red. Then lambda_{r_j} = -pd_j must be nonpositive.
  MyMatrix<T> MT(n_x, n_x);
  MyVector<T> rhs2(n_x);
  for (int j = 0; j < n_x; j++) {
    int r = rows[j];
    for (int i = 0; i < n_x; i++) {
      MT(i, j) = ListIneq(r, i + 1);
    }
  }
  for (int i = 0; i < n_x; i++) {
    rhs2(i) = ToBeMinimized(i + 1);
  }
  std::optional<std::pair<MyVector<T>, T>> optP =
      FractionFreeSolveSquare(MT, rhs2);
  if (!optP) {
    return {};
  }
  MyVector<T> const &pdnum = optP->first;
  T const &dden = optP->second;
  bool dden_pos = dden > 0;
  T dualNum(0);
  for (int j = 0; j < n_x; j++) {
    if (dden_pos ? (pdnum(j) < 0) : (pdnum(j) > 0)) {
      return {};
    }
    dualNum += pdnum(j) * ListIneq(rows[j], 0);
  }
  T dualValue = ToBeMinimized(0) - dualNum / dden;
  if (dualValue != primalValue) {
    return {};
  }
  MyVector<T> x(n_x);
  for (int k = 0; k < n_x; k++) {
    x(k) = xnum(k) / den;
  }
  MyVector<T> DualSolution = ZeroVector<T>(m);
  for (int j = 0; j < n_x; j++) {
    DualSolution(rows[j]) = -pdnum(j) / dden;
  }
  LpSolution<T> sol;
  sol.OptimalValue = primalValue;
  sol.DirectSolution = x;
  sol.DualSolution = DualSolution;
  return sol;
}

// The main entry point: solve the LP in floating point first, certify the
// resulting basis exactly, and fall back to the exact fraction free solver
// warm started from the floating-point basis when certification fails or
// when the floating-point solver reports anything else than optimality.
// The output is always exact and certificate-backed.
template <typename T, typename Tfloat = double>
LpSolution<T> SIMPLEX_LinearProgramming(MyMatrix<T> const &ListIneq,
                                        MyVector<T> const &ToBeMinimized,
                                        std::ostream &os) {
#ifdef TIMINGS_SIMPLEX_CLARKSON
  MicrosecondTime time;
#endif
  int n_x = ListIneq.cols() - 1;
  if (n_x < 3) {
    return SIMPLEX_LinearProgramming_exact(ListIneq, ToBeMinimized, os);
  }
  FloatBasisHint hint =
      ComputeFloatBasisHint<T, Tfloat>(ListIneq, ToBeMinimized);
  if (hint.status == FloatLpStatus::Optimal) {
    std::optional<LpSolution<T>> opt =
        SIMPLEX_CertifyBasis(ListIneq, ToBeMinimized, hint.cobasis_rows, os);
    if (opt) {
#ifdef DEBUG_SIMPLEX_CLARKSON
      os << "SIMPLEX_CLARKSON: floating-point basis certified\n";
#endif
#ifdef TIMINGS_SIMPLEX_CLARKSON
      os << "|SIMPLEX_CLARKSON: SIMPLEX_LinearProgramming(lift)|=" << time
         << "\n";
#endif
      return *opt;
    }
  }
#ifdef DEBUG_SIMPLEX_CLARKSON
  os << "SIMPLEX_CLARKSON: falling back to the exact solver, float status="
     << static_cast<int>(hint.status) << "\n";
#endif
  FractionFreeSimplex<T> solver(ListIneq, ToBeMinimized, os);
  solver.ApplyBasisHint(hint.cobasis_rows);
  LpSolution<T> sol = solver.solve();
#ifdef TIMINGS_SIMPLEX_CLARKSON
  os << "|SIMPLEX_CLARKSON: SIMPLEX_LinearProgramming(fallback)|=" << time
     << "\n";
#endif
  return sol;
}

/*
  Clarkson redundancy elimination.

  Given inequalities f_i(x) = b_i + a_i . x >= 0, the row i is redundant
  when removing it does not change the feasible set, that is when
  min f_i over the other rows is >= 0. The direct algorithm solves one such
  LP per row over the full system. Clarkson's output-sensitive algorithm
  (Clarkson, "More output-sensitive geometric algorithms"; Fukuda,
  "Polyhedral computation", chapter 7) instead maintains a set S of rows
  certified nonredundant and tests each row i only against S:
  --- If min f_i over S (with the relaxation f_i + 1 >= 0 to keep the LP
      bounded) is >= 0 then i is redundant for S, hence for the full system.
  --- Otherwise the LP witness x* satisfies f_S(x*) >= 0, f_i(x*) < 0, and
      shooting a ray from an interior point z towards x* certifies the
      first-hit row as nonredundant (a facet); that row joins S.
  Each LP thus either decides a row or grows S, so at most m + |S| LPs are
  solved, each over roughly |S| rows only, and each solved by the float-first
  certified solver above. The ray shooting uses the same division-free
  minimum-ratio selection with lexicographic tie-breaking as cddlib's
  dd_RayShooting, which guarantees the selected row is a facet (with a
  consistent smallest-index representative for proportional duplicate rows).

  When no interior point exists (the feasible set is not full dimensional)
  the code falls back to the direct algorithm.
 */

// The direct redundancy elimination: one LP per row over the currently
// kept rows. Used as fallback and as a test oracle. The rows are processed
// in decreasing index order, as in dd_RedundantRows, so that among a group
// of mutually redundant rows (duplicates, positive multiples) the smallest
// index is kept -- the same representative the ray shooting tie-break of
// the Clarkson method selects.
template <typename T>
std::vector<int>
SIMPLEX_RedundancyReductionDirect(MyMatrix<T> const &ListIneq,
                                  std::ostream &os) {
  int m = ListIneq.rows();
  int n = ListIneq.cols();
  std::vector<uint8_t> keep(m, 1);
  for (int i = m - 1; i >= 0; i--) {
    int n_kept = 0;
    for (int j = 0; j < m; j++) {
      if (j != i && keep[j]) {
        n_kept++;
      }
    }
    MyMatrix<T> LPmat(n_kept + 1, n);
    int pos = 0;
    for (int j = 0; j < m; j++) {
      if (j != i && keep[j]) {
        for (int k = 0; k < n; k++) {
          LPmat(pos, k) = ListIneq(j, k);
        }
        pos++;
      }
    }
    for (int k = 0; k < n; k++) {
      LPmat(pos, k) = ListIneq(i, k);
    }
    LPmat(pos, 0) += T(1);
    MyVector<T> obj(n);
    for (int k = 0; k < n; k++) {
      obj(k) = ListIneq(i, k);
    }
    LpSolution<T> sol = SIMPLEX_LinearProgramming(LPmat, obj, os);
    if (!sol.DirectSolution || !sol.DualSolution) {
      std::cerr << "SIMPLEX_CLARKSON: the redundancy LP should always be "
                   "feasible and bounded\n";
      throw TerminalException{1};
    }
    if (sol.OptimalValue >= 0) {
      keep[i] = 0;
    }
  }
  std::vector<int> ListIdx;
  for (int i = 0; i < m; i++) {
    if (keep[i]) {
      ListIdx.push_back(i);
    }
  }
  return ListIdx;
}

template <typename T> struct ClarksonRedundancyReduction {
  MyMatrix<T> const &ListIneq_;
  std::ostream &os_;
  int m;
  int n_x;
  std::vector<uint8_t> decided;
  std::vector<uint8_t> redundant;
  std::vector<int> S;
  MyVector<T> z;
  MyVector<T> Fz;
#ifdef DEBUG_SIMPLEX_CLARKSON
  size_t n_lp = 0;
  size_t n_shoot = 0;
#endif

  ClarksonRedundancyReduction(MyMatrix<T> const &ListIneq, std::ostream &os)
      : ListIneq_(ListIneq), os_(os), m(ListIneq.rows()),
        n_x(ListIneq.cols() - 1), decided(m, 0), redundant(m, 0) {}

  // Compute an interior point z with f_i(z) > 0 for all i by maximizing
  // the margin t subject to f_i(x) - t >= 0 and t <= 1. Returns false when
  // the feasible set is not full dimensional.
  bool ComputeInteriorPoint() {
    MyMatrix<T> LPmat(m + 1, n_x + 2);
    for (int i = 0; i < m; i++) {
      for (int k = 0; k <= n_x; k++) {
        LPmat(i, k) = ListIneq_(i, k);
      }
      LPmat(i, n_x + 1) = T(-1);
    }
    for (int k = 0; k <= n_x; k++) {
      LPmat(m, k) = T(0);
    }
    LPmat(m, 0) = T(1);
    LPmat(m, n_x + 1) = T(-1);
    MyVector<T> obj = ZeroVector<T>(n_x + 2);
    obj(n_x + 1) = T(-1);
    LpSolution<T> sol = SIMPLEX_LinearProgramming(LPmat, obj, os_);
    if (!sol.DirectSolution || !sol.DualSolution) {
      std::cerr << "SIMPLEX_CLARKSON: the interior point LP should always "
                   "be feasible and bounded\n";
      throw TerminalException{1};
    }
    if (!(sol.OptimalValue < 0)) {
      return false;
    }
    MyVector<T> const &xt = *sol.DirectSolution;
    z = MyVector<T>(n_x);
    for (int k = 0; k < n_x; k++) {
      z(k) = xt(k);
    }
    Fz = MyVector<T>(m);
    for (int i = 0; i < m; i++) {
      T eSum = ListIneq_(i, 0);
      for (int k = 0; k < n_x; k++) {
        eSum += ListIneq_(i, k + 1) * z(k);
      }
      Fz(i) = eSum;
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
      if (!(Fz(i) > 0)) {
        std::cerr << "SIMPLEX_CLARKSON: the interior point is not interior "
                     "at row "
                  << i << "\n";
        throw TerminalException{1};
      }
#endif
    }
    return true;
  }

  // Shoot the ray z + t d for t > 0 and return the first row hit, that is
  // the row minimizing the ratio (a_i . d) / f_i(z) provided this minimum
  // is negative, with the lexicographic tie-break of dd_RayShooting.
  // Returns -1 when the ray never exits the feasible set.
  int RayShoot(MyVector<T> const &d) {
#ifdef DEBUG_SIMPLEX_CLARKSON
    n_shoot++;
#endif
    int imin = -1;
    T min_num(0);
    T min_den(1);
    for (int i = 0; i < m; i++) {
      T t2(0);
      for (int k = 0; k < n_x; k++) {
        t2 += ListIneq_(i, k + 1) * d(k);
      }
      if (imin == -1) {
        imin = i;
        min_num = t2;
        min_den = Fz(i);
        continue;
      }
      T lhs = t2 * min_den;
      T rhs = min_num * Fz(i);
      if (lhs < rhs) {
        imin = i;
        min_num = t2;
        min_den = Fz(i);
      } else {
        if (lhs == rhs) {
          // Tie: keep the lexicographically smaller normalized row.
          for (int k = 0; k <= n_x; k++) {
            T val_i = ListIneq_(i, k) * min_den;
            T val_min = ListIneq_(imin, k) * Fz(i);
            if (val_i < val_min) {
              imin = i;
              min_num = t2;
              min_den = Fz(i);
              break;
            }
            if (val_i > val_min) {
              break;
            }
          }
        }
      }
    }
    if (imin == -1 || !(min_num < 0)) {
      return -1;
    }
    return imin;
  }

  void CertifyNonRedundant(int r) {
    decided[r] = 1;
    S.push_back(r);
  }

  // Certify some facets cheaply by shooting along the coordinate
  // directions before starting the main loop.
  void Seed() {
    MyVector<T> d = ZeroVector<T>(n_x);
    for (int k = 0; k < n_x; k++) {
      for (int sigma = 0; sigma < 2; sigma++) {
        d(k) = sigma == 0 ? T(1) : T(-1);
        int ired = RayShoot(d);
        if (ired >= 0 && !decided[ired]) {
          CertifyNonRedundant(ired);
        }
      }
      d(k) = T(0);
    }
  }

  // The redundancy test of row i against the current set S, with the
  // relaxation of row i keeping the LP bounded. Returns the solution of
  // the LP minimizing f_i.
  LpSolution<T> SolveRedundancyLP(int i) {
#ifdef DEBUG_SIMPLEX_CLARKSON
    n_lp++;
#endif
    int n = n_x + 1;
    int n_S = S.size();
    MyMatrix<T> LPmat(n_S + 1, n);
    for (int pos = 0; pos < n_S; pos++) {
      int j = S[pos];
      for (int k = 0; k < n; k++) {
        LPmat(pos, k) = ListIneq_(j, k);
      }
    }
    for (int k = 0; k < n; k++) {
      LPmat(n_S, k) = ListIneq_(i, k);
    }
    LPmat(n_S, 0) += T(1);
    MyVector<T> obj(n);
    for (int k = 0; k < n; k++) {
      obj(k) = ListIneq_(i, k);
    }
    LpSolution<T> sol = SIMPLEX_LinearProgramming(LPmat, obj, os_);
    if (!sol.DirectSolution || !sol.DualSolution) {
      std::cerr << "SIMPLEX_CLARKSON: the redundancy LP should always be "
                   "feasible and bounded\n";
      throw TerminalException{1};
    }
    return sol;
  }

  std::vector<int> run() {
    if (m == 0) {
      return {};
    }
    if (!ComputeInteriorPoint()) {
#ifdef DEBUG_SIMPLEX_CLARKSON
      os_ << "SIMPLEX_CLARKSON: no interior point, falling back to the "
             "direct algorithm\n";
#endif
      return SIMPLEX_RedundancyReductionDirect(ListIneq_, os_);
    }
    Seed();
    for (int i = 0; i < m; i++) {
      while (!decided[i]) {
        LpSolution<T> sol = SolveRedundancyLP(i);
        if (sol.OptimalValue >= 0) {
          decided[i] = 1;
          redundant[i] = 1;
        } else {
          MyVector<T> const &xstar = *sol.DirectSolution;
          MyVector<T> d(n_x);
          for (int k = 0; k < n_x; k++) {
            d(k) = xstar(k) - z(k);
          }
          int ired = RayShoot(d);
          if (ired < 0 || decided[ired]) {
            std::cerr << "SIMPLEX_CLARKSON: the ray shooting returned "
                      << ired << " which is impossible since the witness "
                      << "violates row " << i << "\n";
            throw TerminalException{1};
          }
          CertifyNonRedundant(ired);
        }
      }
    }
#ifdef DEBUG_SIMPLEX_CLARKSON
    os_ << "SIMPLEX_CLARKSON: Clarkson done, m=" << m << " |S|=" << S.size()
        << " n_lp=" << n_lp << " n_shoot=" << n_shoot << "\n";
#endif
    std::vector<int> ListIdx;
    for (int i = 0; i < m; i++) {
      if (!redundant[i]) {
        ListIdx.push_back(i);
      }
    }
    return ListIdx;
  }
};

/*
  The block (orbit) accelerated Clarkson redundancy elimination. When a
  group of symmetries of the inequality system permutes the rows, the
  facet-ness of a row is constant on each orbit. Since for a full
  dimensional system the Clarkson notions "redundant" (valid on the
  polyhedron) and "nonredundant" (a facet) are complementary, both
  conclusions propagate to the whole orbit:
  --- a row certified nonredundant by ray shooting certifies its whole
      orbit, which then enters the certified set S in one step;
  --- a row proven redundant against S is redundant for the full system
      (S only relaxes the polyhedron), hence so is its whole orbit.
  Each LP or ray shot thus decides an entire orbit and the number of LPs
  scales with the number of orbits instead of the number of rows.

  The blocks are given by BlockBelong: BlockBelong[i] is the block index
  of row i, with values covering 0..n_block-1. Correctness requires the
  blocks to be orbits (or unions of orbits) of a genuine symmetry group
  of the inequality system and the rows to be pairwise non-proportional;
  otherwise the status propagation is unsound.

  The class reuses the machinery of ClarksonRedundancyReduction (interior
  point, ray shooting, redundancy LP) through a contained core object and
  only the driving loop differs by the orbit-wise propagation.
 */
template <typename T> struct ClarksonRedundancyReductionBlock {
  ClarksonRedundancyReduction<T> core;
  std::vector<std::vector<int>> Blocks;
  std::vector<int> block_id;

  ClarksonRedundancyReductionBlock(MyMatrix<T> const &ListIneq,
                                   std::vector<int> const &BlockBelong,
                                   std::ostream &os)
      : core(ListIneq, os) {
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
    if (static_cast<int>(BlockBelong.size()) != core.m) {
      std::cerr << "SIMPLEX_CLARKSON: |BlockBelong|=" << BlockBelong.size()
                << " but the matrix has " << core.m << " rows\n";
      throw TerminalException{1};
    }
#endif
    int n_block = 0;
    for (int i = 0; i < core.m; i++) {
#ifdef SANITY_CHECK_SIMPLEX_CLARKSON
      if (BlockBelong[i] < 0) {
        std::cerr << "SIMPLEX_CLARKSON: negative block index at row " << i
                  << "\n";
        throw TerminalException{1};
      }
#endif
      if (BlockBelong[i] >= n_block) {
        n_block = BlockBelong[i] + 1;
      }
    }
    block_id = BlockBelong;
    Blocks.resize(n_block);
    for (int i = 0; i < core.m; i++) {
      Blocks[BlockBelong[i]].push_back(i);
    }
  }

  // A certified facet certifies its whole orbit, which enters S.
  void CertifyNonRedundantBlock(int r) {
    for (int j : Blocks[block_id[r]]) {
      if (!core.decided[j]) {
        core.decided[j] = 1;
        core.S.push_back(j);
      }
    }
  }

  // A row proven redundant is redundant for the full system and so is
  // its whole orbit.
  void MarkRedundantBlock(int i) {
    for (int j : Blocks[block_id[i]]) {
      if (!core.decided[j]) {
        core.decided[j] = 1;
        core.redundant[j] = 1;
      }
    }
  }

  std::vector<int> run() {
    if (core.m == 0) {
      return {};
    }
    if (!core.ComputeInteriorPoint()) {
#ifdef DEBUG_SIMPLEX_CLARKSON
      core.os_ << "SIMPLEX_CLARKSON: no interior point, falling back to the "
                  "direct algorithm\n";
#endif
      return SIMPLEX_RedundancyReductionDirect(core.ListIneq_, core.os_);
    }
    // Seeding by coordinate direction shots, with orbit propagation.
    MyVector<T> d = ZeroVector<T>(core.n_x);
    for (int k = 0; k < core.n_x; k++) {
      for (int sigma = 0; sigma < 2; sigma++) {
        d(k) = sigma == 0 ? T(1) : T(-1);
        int ired = core.RayShoot(d);
        if (ired >= 0 && !core.decided[ired]) {
          CertifyNonRedundantBlock(ired);
        }
      }
      d(k) = T(0);
    }
    for (int i = 0; i < core.m; i++) {
      while (!core.decided[i]) {
        LpSolution<T> sol = core.SolveRedundancyLP(i);
        if (sol.OptimalValue >= 0) {
          MarkRedundantBlock(i);
        } else {
          MyVector<T> const &xstar = *sol.DirectSolution;
          MyVector<T> dir(core.n_x);
          for (int k = 0; k < core.n_x; k++) {
            dir(k) = xstar(k) - core.z(k);
          }
          int ired = core.RayShoot(dir);
          if (ired < 0 || core.decided[ired]) {
            std::cerr << "SIMPLEX_CLARKSON: the ray shooting returned "
                      << ired << " which is impossible since the witness "
                      << "violates row " << i << "\n";
            throw TerminalException{1};
          }
          CertifyNonRedundantBlock(ired);
        }
      }
    }
#ifdef DEBUG_SIMPLEX_CLARKSON
    core.os_ << "SIMPLEX_CLARKSON: ClarksonBlock done, m=" << core.m
             << " n_block=" << Blocks.size() << " |S|=" << core.S.size()
             << " n_lp=" << core.n_lp << " n_shoot=" << core.n_shoot << "\n";
#endif
    std::vector<int> ListIdx;
    for (int i = 0; i < core.m; i++) {
      if (!core.redundant[i]) {
        ListIdx.push_back(i);
      }
    }
    return ListIdx;
  }
};

// Clarkson redundancy elimination for rows in inhomogeneous form: row i is
// the inequality ListIneq(i,0) + sum_j ListIneq(i,j) x_j >= 0. Returns the
// ascending list of indices of the nonredundant rows.
template <typename T>
std::vector<int>
SIMPLEX_RedundancyReductionClarkson(MyMatrix<T> const &ListIneq,
                                    std::ostream &os) {
#ifdef TIMINGS_SIMPLEX_CLARKSON
  MicrosecondTime time;
#endif
  ClarksonRedundancyReduction<T> crr(ListIneq, os);
  std::vector<int> ListIdx = crr.run();
#ifdef TIMINGS_SIMPLEX_CLARKSON
  os << "|SIMPLEX_CLARKSON: SIMPLEX_RedundancyReductionClarkson|=" << time
     << "\n";
#endif
  return ListIdx;
}

// The block accelerated variant: BlockBelong[i] is the orbit index of
// row i under a symmetry group of the inequality system. See the comment
// of ClarksonRedundancyReductionBlock for the correctness requirements.
template <typename T>
std::vector<int>
SIMPLEX_RedundancyReductionClarksonBlocks(MyMatrix<T> const &ListIneq,
                                          std::vector<int> const &BlockBelong,
                                          std::ostream &os) {
#ifdef TIMINGS_SIMPLEX_CLARKSON
  MicrosecondTime time;
#endif
  ClarksonRedundancyReductionBlock<T> crrb(ListIneq, BlockBelong, os);
  std::vector<int> ListIdx = crrb.run();
#ifdef TIMINGS_SIMPLEX_CLARKSON
  os << "|SIMPLEX_CLARKSON: SIMPLEX_RedundancyReductionClarksonBlocks|="
     << time << "\n";
#endif
  return ListIdx;
}

// Variant for the homogeneous cone setting: rows of FAC are inequalities
// sum_j FAC(i,j) x_j >= 0. This matches cdd::RedundancyReductionClarksonExt.
template <typename T>
std::vector<int>
SIMPLEX_RedundancyReductionClarksonExt(MyMatrix<T> const &FAC,
                                       std::ostream &os) {
  int m = FAC.rows();
  int n = FAC.cols();
  MyMatrix<T> ListIneq(m, n + 1);
  for (int i = 0; i < m; i++) {
    ListIneq(i, 0) = T(0);
    for (int j = 0; j < n; j++) {
      ListIneq(i, j + 1) = FAC(i, j);
    }
  }
  return SIMPLEX_RedundancyReductionClarkson(ListIneq, os);
}

// clang-format off
#endif  // SRC_POLY_POLY_SIMPLEXCLARKSON_H_
// clang-format on
