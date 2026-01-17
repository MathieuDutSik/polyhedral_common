// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_FINITEMATRIXGROUPTEST_H_
#define SRC_GROUP_FINITEMATRIXGROUPTEST_H_

// clang-format off
#include "COMB_Combinatorics.h"
#include "MatrixGroupSimplification.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_FINITE_MATRIX_GROUP
#endif

/*
  For a dimension N, we want to find all the possible integers k such that there
  exist an integer matrix A of order k. The solution is given in
  https://en.wikipedia.org/wiki/Crystallographic_restriction_theorem
  and involves the Psi function
 */
template <typename T>
std::vector<T> GetIntegralMatricesPossibleOrders(T const &N) {
  auto is_prime = [](T const &x) -> bool {
    if (x == 1)
      return false;
    if (x == 2)
      return true;
    T div = 2;
    while (true) {
      T res = ResInt(x, div);
      if (res == 0)
        return false;
      div += 1;
      if (div * div > x)
        break;
    }
    return true;
  };
  std::vector<T> ListPrime;
  for (T val = 2; val <= N + 1; val++) {
    bool test = is_prime(val);
#ifdef DEBUG_FINITE_MATRIX_GROUP
    std::cerr << "FMGT: val=" << val << " test=" << test << "\n";
#endif
    if (test)
      ListPrime.push_back(val);
  }
  //
  struct pair {
    T fact;
    T dim_cost;
  };
  struct Desc {
    T prime;
    std::vector<pair> l_pair;
  };
  auto get_pair = [&](T const &eprime, int const &k) -> pair {
    if (eprime == 2 && k == 1)
      return {2, 0};
    T pow1 = MyPow(eprime, k - 1);
    T fact = pow1 * eprime;
    T dim_cost = pow1 * (eprime - 1);
    return {fact, dim_cost};
  };
  auto get_l_pair = [&](T const &eprime) -> std::vector<pair> {
    std::vector<pair> l_pair;
    l_pair.push_back({1, 0});
    int k = 1;
    while (true) {
      pair epair = get_pair(eprime, k);
      if (epair.dim_cost > N)
        break;
      l_pair.push_back(epair);
      k++;
    }
    return l_pair;
  };
  std::vector<Desc> l_desc;
  for (auto &ePrime : ListPrime)
    l_desc.push_back({ePrime, get_l_pair(ePrime)});
#ifdef DEBUG_FINITE_MATRIX_GROUP
  size_t n_case = 1;
#endif
  std::vector<int> VectSiz;
  for (auto eDesc : l_desc) {
    size_t len = eDesc.l_pair.size();
#ifdef DEBUG_FINITE_MATRIX_GROUP
    n_case *= len;
#endif
    VectSiz.push_back(len);
  }
#ifdef DEBUG_FINITE_MATRIX_GROUP
  std::cerr << "FMGT: n_case=" << n_case << "\n";
#endif
  std::vector<T> l_order;
  for (auto &V : BlockIterationMultiple(VectSiz)) {
    T tot_dim = 0;
    T order = 1;
    for (size_t iPrime = 0; iPrime < ListPrime.size(); iPrime++) {
      int pos = V[iPrime];
      pair epair = l_desc[iPrime].l_pair[pos];
      tot_dim += epair.dim_cost;
      order *= epair.fact;
    }
    if (tot_dim <= N)
      l_order.push_back(order);
  }
  std::sort(l_order.begin(), l_order.end());
  return l_order;
}

template <typename T>
T max_torsion_order_gl_nz(T const &N) {
  std::vector<T> l_order = GetIntegralMatricesPossibleOrders(N);
  size_t len = l_order.size();
  return l_order[len-1];
}

template<typename T>
struct FiniteMatrixGroupTest {
private:
  std::vector<MyMatrix<T>> l_elts;
  std::unordered_set<MyMatrix<T>> set_elts;
  std::vector<MyMatrix<T>> l_gens;
  MyMatrix<T> Id;
  size_t n_done;
  std::optional<bool> is_finite;
  int max_elt_order;
  std::ostream& os;
  bool insert_one_elt(MyMatrix<T> const& x) {
    bool was_inserted = set_elts.insert(x).second;
    if (!was_inserted) { // already present
      return true;
    }
    l_elts.push_back(x);
    return false;
  }
  bool insert_cycle_one_elt(MyMatrix<T> const& x) {
    bool test = insert_one_elt(x);
    if (test) { // Already present, no change
      return false;
    }
    MyMatrix<T> x_work = x;
    for (int u=0; u<max_elt_order; u++) {
      x_work *= x;
      if (x_work == Id) {
        return false;
      }
      insert_one_elt(x_work);
    }
#ifdef DEBUG_FINITE_MATRIX_GROUP
    os << "FMGT: FiniteMatrixGroupTest, setting is_finite=false by bad element\n";
#endif
    is_finite = false; // element of infinite order. Cannot be finite
    return true;
  }
  bool one_level_increase() {
    if (l_gens.size() == 0) { // no generator, no increase to do.
      return true;
    }
    if (is_finite.has_value()) { // Computation has
      return true;
    }
    MyMatrix<T> x_treat = l_elts[n_done]; // First untreated. Copy due to reallocations.
    for (auto & egen: l_gens) {
      MyMatrix<T> prod = x_treat * egen;
      bool test = insert_cycle_one_elt(prod);
      if (test) { // insertion concluded about non-finiteness. Terminating.
        return true;
      }
    }
    n_done += 1;
    if (n_done == l_elts.size()) { // Enumeration has finished.
#ifdef DEBUG_FINITE_MATRIX_GROUP
      os << "FMGT: FiniteMatrixGroupTest, setting is_finite=true by enumeration n_done=" << n_done << "\n";
#endif
      is_finite = true;
    }
    return false;
  }
public:
  FiniteMatrixGroupTest(std::vector<MyMatrix<T>> const& _l_gens, std::ostream& _os) : l_gens(_l_gens), os(_os) {
    if (l_gens.size() == 0) { // no generator, group is trivial and so empty
      is_finite = true;
      return;
    }
    int N = l_gens[0].rows();
    max_elt_order = max_torsion_order_gl_nz(N);
#ifdef DEBUG_FINITE_MATRIX_GROUP
    os << "FMGT: FiniteMatrixGroupTest, N=" << N << " max_elt_order=" << max_elt_order << " |l_gens|=" << l_gens.size() << "\n";
#endif
    Id = IdentityMat<T>(N);
    n_done = 0;
    for (auto & egen: l_gens) {
      bool test = insert_cycle_one_elt(egen);
      if (test) {
        return;
      }
    }
  }
  void timed_increase(int64_t n_nanoseconds) {
    NanosecondTime time;
    while(true) {
      bool test = one_level_increase();
      if (test) { // Nothing more to be done, returning
        return;
      }
      int64_t time_spent = time.const_eval_int64();
      if (time_spent > n_nanoseconds) { // Time budget spent, returning
        return;
      }
    }
  }
  std::optional<bool> get_status() const {
    return is_finite;
  }
  bool full_resolution() {
    while(true) {
      bool test = one_level_increase();
      if (test) { // Nothing more to be done, returning
        break;
      }
    }
    if (!is_finite) {
      std::cerr << "We should have concluded by that point\n";
      throw TerminalException{1};
    }
    return *is_finite;
  }
};

// clang-format off
#endif  // SRC_GROUP_FINITEMATRIXGROUPTEST_H_
// clang-format on
