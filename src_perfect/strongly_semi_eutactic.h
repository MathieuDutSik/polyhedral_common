// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_STRONGLY_SEMI_EUTACTIC_H_
#define SRC_PERFECT_STRONGLY_SEMI_EUTACTIC_H_

// clang-format off
#include "MAT_Matrix.h"
#include "Boost_bitset_kernel.h"
#include "POLY_SimplexClarkson.h"
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_STRONGLY_SEMI_EUTACTIC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_STRONGLY_SEMI_EUTACTIC
#endif

/*
  Determination of whether a positive definite quadratic form A is
  strongly semi-eutactic: whether there exist mu > 0 and a subset S of
  Min_h(A) (half of the minimal vectors, one per antipodal pair {v,-v})
  such that
     A^{-1} = mu * sum_{v in S} v v^T
  that is a eutaxy relation with all coefficients lambda_v in {0, mu}.
  The strongly eutactic forms are the special case S = Min_h(A).

  Algorithm: write the eutaxy relations as lambda . B = c with B the
  (m x n(n+1)/2) matrix of rows vec(v v^T) and c = vec(A^{-1}). The
  solution set is lambda(x) = lambda_0 + x . K with K the left kernel
  of B, of dimension d (= m - n(n+1)/2 for a perfect form). Over the
  unknowns y = (x, mu) each condition "lambda_v(x) = 0 or mu" is a
  disjunction of two hyperplane conditions, and a branch and bound
  search is run over these choices:
  --- The current solution subspace is maintained exactly in a
      parametrized form y = q + P t and shrinks by one dimension at
      each imposed condition.
  --- A vector whose lambda_v is constant on the subspace is forced:
      value 0 keeps it out of S, a value cst > 0 forces mu = cst, a
      negative value prunes the branch.
  --- Any solution satisfies the relaxation 0 <= lambda_v <= mu
      together with n/(m * min) <= mu <= 1/min (which follows from
      pairing the relation with A: mu * |S| * min = n and n <= |S| <= m),
      so the branches are pruned by exact linear programming feasibility
      of that relaxation.
  Every returned certificate (S, mu) is directly a proof; refutation is
  the exhaustion of the search tree (the result is then "resolved" with
  no certificate, unless the node budget max_node was hit first).
*/

template <typename T> struct StronglySemiEutacticCert {
  // The rows of the input SHV matrix receiving lambda_v = mu. When the
  // input contains both v and -v of a pair, the first occurrence is the
  // one possibly selected in the subset.
  Face subset;
  T mu;
};

template <typename T> struct StronglySemiEutacticResult {
  // false when the node budget was exhausted before completing the search
  bool resolved;
  std::optional<StronglySemiEutacticCert<T>> cert;
  size_t n_node;
};

// Internal state of the branch and bound search.
template <typename T> struct StronglySemiEutacticSearch {
private:
  static constexpr uint8_t STATE_UNFIXED = 0;
  static constexpr uint8_t STATE_ZERO = 1;
  static constexpr uint8_t STATE_MU = 2;
  // Minimal remaining dimension for running the LP pruning: for smaller
  // dimensions the plain enumeration is cheaper than the LP calls.
  static constexpr int LP_MIN_DIM = 3;
  int n;
  int mh;
  int D;
  // Affine forms over y = (x_1, ..., x_d, mu): the coefficient vectors
  // (length D) and constants of the lambda_v, and the norm data
  std::vector<MyVector<T>> lambda_a;
  std::vector<T> lambda_b;
  T min_norm;
  size_t max_node;
  size_t n_node;
  bool budget_hit;
  std::optional<std::pair<Face, T>> found;
  std::vector<int> map_input_row;
  int n_input_row;
  std::ostream &os;

  // Value of an affine form (a, b) over y on the subspace y = q + P t:
  // returns (alpha, beta) with value = beta + alpha . t
  std::pair<MyVector<T>, T> EvalForm(MyVector<T> const &a, T const &b,
                                     MyMatrix<T> const &P,
                                     MyVector<T> const &q) const {
    int r = P.cols();
    MyVector<T> alpha = ZeroVector<T>(r);
    T beta = b;
    for (int i = 0; i < D; i++) {
      if (a(i) != 0) {
        beta += a(i) * q(i);
        for (int j = 0; j < r; j++)
          alpha(j) += a(i) * P(i, j);
      }
    }
    return {std::move(alpha), std::move(beta)};
  }

  // Impose beta + alpha . t = 0 on the subspace, eliminating one
  // parameter. Returns false when the constraint is contradictory.
  // When the constraint is redundant, P and q are left unchanged.
  bool ImposeConstraint(MyVector<T> const &alpha, T const &beta,
                        MyMatrix<T> &P, MyVector<T> &q) const {
    int r = P.cols();
    int pivot = -1;
    for (int j = 0; j < r; j++) {
      if (alpha(j) != 0) {
        pivot = j;
        break;
      }
    }
    if (pivot == -1)
      return beta == 0;
    T fact = beta / alpha(pivot);
    for (int i = 0; i < D; i++)
      q(i) -= P(i, pivot) * fact;
    MyMatrix<T> Pnew(D, r - 1);
    int pos = 0;
    for (int j = 0; j < r; j++) {
      if (j == pivot)
        continue;
      T coef = alpha(j) / alpha(pivot);
      for (int i = 0; i < D; i++)
        Pnew(i, pos) = P(i, j) - P(i, pivot) * coef;
      pos++;
    }
    P = std::move(Pnew);
    return true;
  }

  // The exact LP feasibility pruning: existence of t with
  // 0 <= lambda_v(t) <= mu(t) for the unfixed v and
  // n/(mh * min) <= mu(t) <= 1/min. Returns false when provably empty.
  bool IsRelaxationFeasible(MyMatrix<T> const &P, MyVector<T> const &q,
                            std::vector<uint8_t> const &state) const {
    int r = P.cols();
    std::vector<std::pair<MyVector<T>, T>> ListConstraint;
    std::pair<MyVector<T>, T> mu_form = GetMuForm(P, q);
    for (int v = 0; v < mh; v++) {
      if (state[v] != STATE_UNFIXED)
        continue;
      std::pair<MyVector<T>, T> lam = EvalForm(lambda_a[v], lambda_b[v], P, q);
      // mu - lambda_v >= 0
      MyVector<T> diff_a = mu_form.first - lam.first;
      T diff_b = mu_form.second - lam.second;
      ListConstraint.push_back({std::move(diff_a), std::move(diff_b)});
      // lambda_v >= 0
      ListConstraint.push_back(std::move(lam));
    }
    // mh * min * mu - n >= 0
    T n_T = UniversalScalarConversion<T, int>(n);
    T mh_T = UniversalScalarConversion<T, int>(mh);
    MyVector<T> low_a = mu_form.first * (mh_T * min_norm);
    T low_b = mu_form.second * (mh_T * min_norm) - n_T;
    ListConstraint.push_back({std::move(low_a), std::move(low_b)});
    // 1 - min * mu >= 0
    MyVector<T> upp_a = -mu_form.first * min_norm;
    T upp_b = 1 - mu_form.second * min_norm;
    ListConstraint.push_back({std::move(upp_a), std::move(upp_b)});
    //
    int n_cons = ListConstraint.size();
    MyMatrix<T> ListIneq(n_cons, r + 1);
    for (int i_cons = 0; i_cons < n_cons; i_cons++) {
      ListIneq(i_cons, 0) = ListConstraint[i_cons].second;
      for (int j = 0; j < r; j++)
        ListIneq(i_cons, 1 + j) = ListConstraint[i_cons].first(j);
    }
    MyVector<T> ToBeMinimized = ZeroVector<T>(r + 1);
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
    return eSol.DirectSolution.has_value();
  }

  std::pair<MyVector<T>, T> GetMuForm(MyMatrix<T> const &P,
                                      MyVector<T> const &q) const {
    int r = P.cols();
    MyVector<T> alpha(r);
    for (int j = 0; j < r; j++)
      alpha(j) = P(D - 1, j);
    return {std::move(alpha), q(D - 1)};
  }

  // Fix the vectors whose lambda_v is constant on the subspace. Returns
  // false when a contradiction is reached.
  bool Propagate(MyMatrix<T> &P, MyVector<T> &q,
                 std::vector<uint8_t> &state) const {
    bool changed = true;
    while (changed) {
      changed = false;
      for (int v = 0; v < mh; v++) {
        if (state[v] != STATE_UNFIXED)
          continue;
        std::pair<MyVector<T>, T> lam =
            EvalForm(lambda_a[v], lambda_b[v], P, q);
        if (!IsZeroVector(lam.first))
          continue;
        if (lam.second == 0) {
          state[v] = STATE_ZERO;
          continue;
        }
        if (lam.second < 0)
          return false;
        // A constant positive lambda_v can only be equal to mu
        std::pair<MyVector<T>, T> mu_form = GetMuForm(P, q);
        MyVector<T> diff_a = -mu_form.first;
        T diff_b = lam.second - mu_form.second;
        if (!ImposeConstraint(diff_a, diff_b, P, q))
          return false;
        state[v] = STATE_MU;
        changed = true;
      }
    }
    // mu must remain possibly positive
    std::pair<MyVector<T>, T> mu_form = GetMuForm(P, q);
    if (IsZeroVector(mu_form.first) && mu_form.second <= 0)
      return false;
    return true;
  }

  // Depth first exploration. Returns true when a certificate was found.
  bool Explore(MyMatrix<T> P, MyVector<T> q, std::vector<uint8_t> state) {
    n_node++;
    if (n_node > max_node) {
      budget_hit = true;
      return false;
    }
    if (!Propagate(P, q, state))
      return false;
    int v_branch = -1;
    for (int v = 0; v < mh; v++) {
      if (state[v] == STATE_UNFIXED) {
        v_branch = v;
        break;
      }
    }
    if (v_branch == -1)
      return TreatLeaf(P, q, state);
    if (P.cols() >= LP_MIN_DIM) {
      if (!IsRelaxationFeasible(P, q, state))
        return false;
    }
    std::pair<MyVector<T>, T> lam =
        EvalForm(lambda_a[v_branch], lambda_b[v_branch], P, q);
    std::pair<MyVector<T>, T> mu_form = GetMuForm(P, q);
    // First branch: lambda_v = mu (finds the strongly eutactic
    // certificates without any descent into the zero branches)
    {
      MyMatrix<T> Pnew = P;
      MyVector<T> qnew = q;
      std::vector<uint8_t> state_new = state;
      state_new[v_branch] = STATE_MU;
      MyVector<T> diff_a = lam.first - mu_form.first;
      T diff_b = lam.second - mu_form.second;
      if (ImposeConstraint(diff_a, diff_b, Pnew, qnew)) {
        if (Explore(std::move(Pnew), std::move(qnew), std::move(state_new)))
          return true;
      }
    }
    if (budget_hit)
      return false;
    // Second branch: lambda_v = 0
    {
      state[v_branch] = STATE_ZERO;
      if (ImposeConstraint(lam.first, lam.second, P, q)) {
        if (Explore(std::move(P), std::move(q), std::move(state)))
          return true;
      }
    }
    return false;
  }

  bool TreatLeaf(MyMatrix<T> const &P, MyVector<T> const &q,
                 std::vector<uint8_t> const &state) {
    // All the lambda_v are fixed: the ones at STATE_MU share the value
    // mu which must be constant and positive
    std::pair<MyVector<T>, T> mu_form = GetMuForm(P, q);
    Face subset(n_input_row);
    bool has_mu = false;
    for (int v = 0; v < mh; v++) {
      if (state[v] == STATE_MU) {
        subset[map_input_row[v]] = 1;
        has_mu = true;
      }
    }
    if (!has_mu)
      return false;
    if (!IsZeroVector(mu_form.first))
      return false;
    T mu = mu_form.second;
    if (mu <= 0)
      return false;
    found = {std::move(subset), std::move(mu)};
    return true;
  }

public:
  StronglySemiEutacticSearch(std::vector<MyVector<T>> const &ListVectHalf,
                             std::vector<int> const &_map_input_row,
                             MyVector<T> const &lambda_part,
                             MyMatrix<T> const &Kern, T const &_min_norm,
                             int _n, int _n_input_row, size_t _max_node,
                             std::ostream &_os)
      : n(_n), mh(ListVectHalf.size()), D(Kern.rows() + 1),
        min_norm(_min_norm), max_node(_max_node), n_node(0),
        budget_hit(false), map_input_row(_map_input_row),
        n_input_row(_n_input_row), os(_os) {
    // The affine form of lambda_v over y = (x, mu) is
    // lambda_part(v) + sum_i x_i Kern(i, v), without mu contribution
    lambda_a.reserve(mh);
    lambda_b.reserve(mh);
    int d = D - 1;
    for (int v = 0; v < mh; v++) {
      MyVector<T> a = ZeroVector<T>(D);
      for (int i = 0; i < d; i++)
        a(i) = Kern(i, v);
      lambda_a.push_back(std::move(a));
      lambda_b.push_back(lambda_part(v));
    }
#ifdef SANITY_CHECK_STRONGLY_SEMI_EUTACTIC
    if (n_input_row < mh) {
      std::cerr << "SSE: incoherent number of input rows\n";
      throw TerminalException{1};
    }
#endif
  }

  StronglySemiEutacticResult<T> Run() {
    MyMatrix<T> P = IdentityMat<T>(D);
    MyVector<T> q = ZeroVector<T>(D);
    std::vector<uint8_t> state(mh, STATE_UNFIXED);
    bool has_found = Explore(std::move(P), std::move(q), std::move(state));
    if (has_found) {
      StronglySemiEutacticCert<T> cert{found->first, found->second};
      return {true, std::move(cert), n_node};
    }
    return {!budget_hit, {}, n_node};
  }
};

// Test whether eGram is strongly semi-eutactic with respect to the
// vector family SHV_T (the minimal vectors, with or without both
// members of the antipodal pairs). max_node bounds the number of nodes
// of the search tree; on exhaustion the result has resolved = false.
template <typename T>
StronglySemiEutacticResult<T>
TestStronglySemiEutactic(MyMatrix<T> const &eGram, MyMatrix<T> const &SHV_T,
                         size_t max_node, std::ostream &os) {
  int n = eGram.rows();
  int n_row = SHV_T.rows();
#ifdef DEBUG_STRONGLY_SEMI_EUTACTIC
  os << "SSE: TestStronglySemiEutactic n=" << n << " n_row=" << n_row << "\n";
#endif
  if (n_row == 0) {
    return {true, {}, 0};
  }
  // Reduction to one vector per antipodal pair
  std::vector<MyVector<T>> ListVectHalf;
  std::vector<int> map_input_row;
  std::map<std::vector<T>, int> map_seen;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV_T, i_row);
    int i_first = 0;
    while (i_first < n && V(i_first) == 0)
      i_first++;
    if (i_first == n)
      continue;
    if (V(i_first) < 0)
      V = -V;
    std::vector<T> V_std(n);
    for (int i = 0; i < n; i++)
      V_std[i] = V(i);
    if (map_seen.try_emplace(std::move(V_std), i_row).second) {
      ListVectHalf.push_back(V);
      map_input_row.push_back(i_row);
    }
  }
  int mh = ListVectHalf.size();
  int dimSymm = (n * (n + 1)) / 2;
  // The linear system lambda . B = c
  MyMatrix<T> B(mh, dimSymm);
  for (int v = 0; v < mh; v++) {
    MyVector<T> const &V = ListVectHalf[v];
    MyMatrix<T> OuterProd(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        OuterProd(i, j) = V(i) * V(j);
    MyVector<T> OuterProdVec = SymmetricMatrixToVector(OuterProd);
    AssignMatrixRow(B, v, OuterProdVec);
  }
  MyMatrix<T> eGramInv = Inverse(eGram);
  MyVector<T> c = SymmetricMatrixToVector(eGramInv);
  std::optional<MyVector<T>> opt_part = SolutionMat(B, c);
  if (!opt_part) {
#ifdef DEBUG_STRONGLY_SEMI_EUTACTIC
    os << "SSE: no eutaxy relation at all\n";
#endif
    return {true, {}, 0};
  }
  MyMatrix<T> Kern = NullspaceMat(B);
#ifdef DEBUG_STRONGLY_SEMI_EUTACTIC
  os << "SSE: mh=" << mh << " kernel dimension=" << Kern.rows() << "\n";
#endif
  T min_norm = EvaluationQuadForm<T, T>(eGram, ListVectHalf[0]);
#ifdef SANITY_CHECK_STRONGLY_SEMI_EUTACTIC
  for (int v = 1; v < mh; v++) {
    T norm = EvaluationQuadForm<T, T>(eGram, ListVectHalf[v]);
    if (norm != min_norm) {
      std::cerr << "SSE: the vector family is not of constant norm\n";
      throw TerminalException{1};
    }
  }
#endif
  StronglySemiEutacticSearch<T> search(ListVectHalf, map_input_row, *opt_part,
                                       Kern, min_norm, n, n_row, max_node, os);
  StronglySemiEutacticResult<T> result = search.Run();
#ifdef SANITY_CHECK_STRONGLY_SEMI_EUTACTIC
  if (result.cert) {
    MyMatrix<T> SumMat = ZeroMatrix<T>(n, n);
    for (int i_row = 0; i_row < n_row; i_row++) {
      if (result.cert->subset[i_row] == 1) {
        MyVector<T> V = GetMatrixRow(SHV_T, i_row);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            SumMat(i, j) += V(i) * V(j);
      }
    }
    SumMat *= result.cert->mu;
    if (SumMat != eGramInv) {
      std::cerr << "SSE: the certificate does not check out\n";
      throw TerminalException{1};
    }
    if (result.cert->mu <= 0) {
      std::cerr << "SSE: the certificate has non positive mu\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef DEBUG_STRONGLY_SEMI_EUTACTIC
  os << "SSE: resolved=" << result.resolved
     << " strongly_semi_eutactic=" << result.cert.has_value()
     << " n_node=" << result.n_node << "\n";
#endif
  return result;
}

// clang-format off
#endif  // SRC_PERFECT_STRONGLY_SEMI_EUTACTIC_H_
// clang-format on
