// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COPOS_SPN_DECOMPOSITION_H_
#define SRC_COPOS_SPN_DECOMPOSITION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "POLY_LinearProgramming.h"
#include "Positivity.h"
#include "SignatureSymmetric.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SPN_DECOMPOSITION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SPN_DECOMPOSITION
#endif

/*
  Splitting a symmetric matrix as Q = P + N with P positive semidefinite and N
  entrywise non-negative. Such a Q is copositive, since for x >= 0 one has
  Q[x] = P[x] + N[x] >= 0. The converse holds up to n = 4 and fails from n = 5
  on, so the split is a sufficient criterion, not a characterization.

  The set of admissible N is the intersection of the polyhedral cone of the
  entrywise non-negative symmetric matrices with the preimage of the positive
  semidefinite cone. The second one is not polyhedral, so the search is the
  cutting plane loop already used by GetOnePositiveSemiDefiniteMatrix of
  Tspace_SearchPositiveMatrices.h:

  --- the current relaxation is a linear program over N: N_ij >= 0 together
      with the accumulated cuts V^T N V <= V^T Q V, each of which is a
      consequence of Q - N being positive semidefinite;
  --- the linear program is solved exactly, so the candidate N is a rational
      matrix that satisfies the constraints exactly, with no rounding;
  --- if Q - N is positive semidefinite the search is over and the pair is a
      certificate that needs no further checking;
  --- otherwise an integral V with (Q - N)[V] < 0 is produced and the
      corresponding cut is added.

  The three outcomes are genuinely different and the caller has to distinguish
  them:

  --- decomposition_found is a proof. P and N are exact and satisfy
      P + N = Q, P positive semidefinite, N >= 0.
  --- no_decomposition is also a proof. The linear program is infeasible, and
      since every cut is a consequence of the positive semidefiniteness, no
      admissible N exists at all. The simplex being exact, its infeasibility
      is a real Farkas certificate and not a numerical accident.
  --- undetermined is the honest answer of a search that did not settle. The
      positive semidefinite cone is not polyhedral, so the loop can refine the
      relaxation forever without either finding a point or emptying it. This
      is not a defect of the implementation, it is the same intrinsic
      limitation noted in Tspace_SearchPositiveMatrices.h.
 */
enum class SPNstatus {
  decomposition_found,
  no_decomposition,
  undetermined,
};

inline std::string SPNstatus_to_string(SPNstatus const &status) {
  if (status == SPNstatus::decomposition_found) {
    return "decomposition_found";
  }
  if (status == SPNstatus::no_decomposition) {
    return "no_decomposition";
  }
  return "undetermined";
}

template <typename T> struct ResultSPN {
  SPNstatus status;
  // Meaningful only when status is decomposition_found.
  MyMatrix<T> P;
  MyMatrix<T> N;
};

/*
  The entries of N above and on the diagonal are the variables of the linear
  program. The map from a pair of indices to the variable index.
 */
inline int spn_variable_index(int const &n, int const &i, int const &j) {
  // i <= j is assumed.
  return i * n - (i * (i - 1)) / 2 + (j - i);
}

/*
  The coefficients of V^T N V as a linear form in the variables N_ij: the
  diagonal terms appear once and the off diagonal ones twice.
 */
template <typename T>
MyVector<T> spn_quadratic_coefficients(int const &n, MyVector<T> const &V) {
  int n_var = (n * (n + 1)) / 2;
  MyVector<T> coeff = ZeroVector<T>(n_var);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T val = V(i) * V(j);
      if (i != j) {
        val *= 2;
      }
      coeff(spn_variable_index(n, i, j)) = val;
    }
  }
  return coeff;
}

template <typename T>
MyMatrix<T> spn_matrix_from_solution(int const &n,
                                     MyVector<T> const &solution) {
  MyMatrix<T> N(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T val = solution(spn_variable_index(n, i, j));
      N(i, j) = val;
      N(j, i) = val;
    }
  }
  return N;
}

template <typename T> bool spn_is_entrywise_nonnegative(MyMatrix<T> const &M) {
  int n = M.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (M(i, j) < 0) {
        return false;
      }
    }
  }
  return true;
}

/*
  An integral V with M[V] < 0, for a symmetric M that is not positive
  semidefinite. The form may be degenerate, so the search is done on a
  complement of its kernel where it is not, exactly as in
  GetOnePositiveSemiDefiniteMatrix of Tspace_SearchPositiveMatrices.h.

  The vector has to be short, not merely negative. Harvesting instead every
  negative direction of the diagonalization gives several cuts per round, but
  those directions are rational with large denominators and, once scaled to
  integers, produce cuts with huge coefficients that make the exact linear
  program collapse. One short vector per round is worth much more than many
  long ones.
 */
template <typename T, typename Tint>
MyVector<Tint> spn_get_negative_vector(MyMatrix<T> const &M,
                                       std::ostream &os) {
  int n = M.rows();
  MyMatrix<T> NSP_T = NullspaceIntMat(M);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint, T>(NSP_T);
  MyMatrix<Tint> Compl = SubspaceCompletionInt<Tint>(NSP, n);
  MyMatrix<T> Compl_T = UniversalMatrixConversion<T, Tint>(Compl);
  MyMatrix<T> Mred = Compl_T * M * Compl_T.transpose();
  T CritNorm(0);
  bool StrictIneq = true;
  MyVector<Tint> V1 =
      GetShortIntegralVector<T, Tint>(Mred, CritNorm, StrictIneq, os);
  MyVector<Tint> V = Compl.transpose() * V1;
#ifdef SANITY_CHECK_SPN_DECOMPOSITION
  T norm = EvaluationQuadForm<T, Tint>(M, V);
  if (norm >= 0) {
    std::cerr << "SPN: the vector should be a counter example but the norm is "
              << norm << "\n";
    throw TerminalException{1};
  }
#endif
  return V;
}

template <typename T, typename Tint>
ResultSPN<T> SearchSPNdecomposition(MyMatrix<T> const &Q, size_t const &max_iter,
                                    std::ostream &os) {
  int n = Q.rows();
  int n_var = (n * (n + 1)) / 2;
  //
  // The two ansatz that need no linear program.
  //
  if (IsPositiveSemiDefinite(Q, os)) {
    MyMatrix<T> N = ZeroMatrix<T>(n, n);
    return {SPNstatus::decomposition_found, Q, std::move(N)};
  }
  if (spn_is_entrywise_nonnegative(Q)) {
    MyMatrix<T> P = ZeroMatrix<T>(n, n);
    return {SPNstatus::decomposition_found, std::move(P), Q};
  }
  //
  // The cuts. The vectors e_i and e_i +- e_j are put in from the start: they
  // are the cheap consequences of the positive semidefiniteness and they give
  // the program its shape before any separation happens.
  //
  std::vector<MyVector<T>> l_cut;
  for (int i = 0; i < n; i++) {
    MyVector<T> V = ZeroVector<T>(n);
    V(i) = 1;
    l_cut.push_back(std::move(V));
  }
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int sign = -1; sign <= 1; sign += 2) {
        MyVector<T> V = ZeroVector<T>(n);
        V(i) = 1;
        V(j) = sign;
        l_cut.push_back(std::move(V));
      }
    }
  }
  /*
    The margin t is the last variable. A cut asks for V^T (Q - N) V >= t rather
    than >= 0, and t is maximized. Asking only for >= 0 puts the optimum of
    every round exactly on the boundary of the cuts already collected, so the
    candidate is never strictly inside and the separation keeps producing new
    cuts forever: that is what makes the naive loop fail even in dimension 3.
    Maximizing the margin picks the candidate that is as far inside as the
    current relaxation allows, which is what makes the loop converge.

    The bound t <= 1 keeps the objective bounded. The scale of t is arbitrary,
    only its sign matters, and a negative optimum is conclusive: it says that
    for every non-negative N one of the current cuts is violated, so Q - N is
    never positive semidefinite and no decomposition exists.
   */
  int idx_t = n_var;
  int n_col = n_var + 1;
  MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_col);
  ToBeMinimized(1 + idx_t) = -1;
  for (size_t iter = 0; iter < max_iter; iter++) {
    int n_cut = l_cut.size();
    int n_row = n_var + 1 + n_cut;
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_row, 1 + n_col);
    // The non-negativity of the entries of N.
    for (int i_var = 0; i_var < n_var; i_var++) {
      ListIneq(i_var, 1 + i_var) = 1;
    }
    // The bound t <= 1.
    ListIneq(n_var, 0) = 1;
    ListIneq(n_var, 1 + idx_t) = -1;
    // The cuts V^T Q V - V^T N V - t >= 0.
    for (int i_cut = 0; i_cut < n_cut; i_cut++) {
      MyVector<T> const &V = l_cut[i_cut];
      MyVector<T> coeff = spn_quadratic_coefficients(n, V);
      ListIneq(n_var + 1 + i_cut, 0) = EvaluationQuadForm<T, T>(Q, V);
      for (int i_var = 0; i_var < n_var; i_var++) {
        ListIneq(n_var + 1 + i_cut, 1 + i_var) = -coeff(i_var);
      }
      ListIneq(n_var + 1 + i_cut, 1 + idx_t) = -1;
    }
#ifdef DEBUG_SPN_DECOMPOSITION
    os << "SPN: iter=" << iter << " n_cut=" << n_cut << " n_var=" << n_var
       << "\n";
#endif
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.DirectSolution) {
      // Infeasible. N = 0 with a small enough t always satisfies the program,
      // so this should not happen. Reported rather than assumed away.
#ifdef DEBUG_SPN_DECOMPOSITION
      os << "SPN: the linear program came back infeasible, giving up\n";
#endif
      return {SPNstatus::undetermined, MyMatrix<T>(), MyMatrix<T>()};
    }
    if (!eSol.DualSolution) {
      // Unbounded. The objective is -t and t <= 1, so this should not happen.
#ifdef DEBUG_SPN_DECOMPOSITION
      os << "SPN: the linear program came back unbounded, giving up\n";
#endif
      return {SPNstatus::undetermined, MyMatrix<T>(), MyMatrix<T>()};
    }
    MyVector<T> const &solution = *eSol.DirectSolution;
    if (solution(idx_t) < 0) {
      // Every non-negative N violates one of the accumulated cuts, and those
      // are consequences of the positive semidefiniteness. The simplex being
      // exact, this is a proof that no decomposition exists.
#ifdef DEBUG_SPN_DECOMPOSITION
      os << "SPN: the best margin is " << solution(idx_t)
         << " < 0, no decomposition exists\n";
#endif
      return {SPNstatus::no_decomposition, MyMatrix<T>(), MyMatrix<T>()};
    }
    MyMatrix<T> N = spn_matrix_from_solution(n, solution);
    MyMatrix<T> P = Q - N;
#ifdef SANITY_CHECK_SPN_DECOMPOSITION
    if (!spn_is_entrywise_nonnegative(N)) {
      std::cerr << "SPN: the linear program returned a N with a negative "
                   "entry\n";
      throw TerminalException{1};
    }
#endif
    if (IsPositiveSemiDefinite(P, os)) {
#ifdef DEBUG_SPN_DECOMPOSITION
      os << "SPN: found a decomposition after " << (iter + 1) << " rounds\n";
#endif
      return {SPNstatus::decomposition_found, std::move(P), std::move(N)};
    }
    MyVector<Tint> V = spn_get_negative_vector<T, Tint>(P, os);
    l_cut.push_back(UniversalVectorConversion<T, Tint>(V));
  }
#ifdef DEBUG_SPN_DECOMPOSITION
  os << "SPN: the limit of " << max_iter << " rounds was reached\n";
#endif
  return {SPNstatus::undetermined, MyMatrix<T>(), MyMatrix<T>()};
}

// Checks a returned decomposition from scratch. Used by the tests and by the
// sanity checks of the callers.
template <typename T>
bool CheckSPNdecomposition(MyMatrix<T> const &Q, ResultSPN<T> const &result,
                           std::ostream &os) {
  if (result.status != SPNstatus::decomposition_found) {
    return true;
  }
  MyMatrix<T> sum = result.P + result.N;
  if (sum != Q) {
    std::cerr << "SPN: P + N differs from Q\n";
    return false;
  }
  if (!spn_is_entrywise_nonnegative(result.N)) {
    std::cerr << "SPN: N has a negative entry\n";
    return false;
  }
  if (!IsPositiveSemiDefinite(result.P, os)) {
    std::cerr << "SPN: P is not positive semidefinite\n";
    return false;
  }
  return true;
}

// clang-format off
#endif  // SRC_COPOS_SPN_DECOMPOSITION_H_
// clang-format on
