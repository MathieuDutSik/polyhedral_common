// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_COMMON_H_
#define SRC_PERFECT_HECKE_COMMON_H_

// clang-format off
#include "perfect_complex.h"
#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HECKE_OPERATORS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_HECKE_OPERATORS
#endif

#ifdef TIMINGS
#define TIMINGS_HECKE_OPERATORS
#endif

/*
  The common part of the Hecke operators (see hecke_operators.h for the
  overview): the clean failures and the budgets, the ambient group G of
  the T-space, and the cosets of a Hecke double coset.
 */

/*
  A clean failure: a method exhausted its budget (memory, iterations,
  cells). It is caught by the driver, which reports it and stops,
  instead of letting the memory or the time grow without bound.
 */
struct HeckeFailure {
  std::string message;
};

// The budgets of the methods. They bound the memory and the running time
// of every step; exceeding one of them is a clean failure.
struct HeckeBudget {
  size_t max_cosets = 20000000;         // cosets of G / Gamma
  size_t max_hecke_cosets = 1000000;    // cosets of a Hecke double coset
  size_t max_closure = 1000000;         // elements of a finite closure
  size_t max_walk_steps = 1000000;      // steps of a walk through the cones
  size_t max_peel_moves = 20000;        // moves of the peeling solver
  size_t max_chain = 50000;             // cells of a chain in the peeling
  size_t max_unknowns = 60000;          // unknowns of a linear system
  size_t max_rings = 12;                // rings of unknowns
  SparseSolveBudget solve;              // the sparse exact solutions
  size_t max_sharbly_symbols = 200000;  // symbols of an Ash-Rudolph reduction
};

//
// The ambient group G of the T-space and its generators
//

template <typename T>
bool is_in_ambient_group(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &g) {
  if (!IsIntegralMatrix(g)) {
    return false;
  }
  T det = DeterminantMat(g);
  if (det != T(1) && det != T(-1)) {
    return false;
  }
  return is_stab_space(g, LinSpa);
}

// The stabilizers of the perfect cones together with the adjacency
// elements generate the group of a connected complex.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
get_ambient_group_generators(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                             std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  std::unordered_set<MyMatrix<Tint>> set_gen;
  std::vector<MyMatrix<Tint>> l_gens;
  auto f_insert = [&](MyMatrix<Tint> const &M) -> void {
    if (M == id) {
      return;
    }
    if (set_gen.insert(M).second) {
      l_gens.push_back(M);
    }
  };
  for (auto &cone : fce.pctdi.l_perfect) {
    for (auto &ePerm : cone.GRP_ext.GeneratorsOfGroup()) {
      f_insert(cone.find_matrix(ePerm, os));
    }
    for (auto &adj : cone.l_sing_adj) {
      f_insert(adj.eMat);
    }
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: get_ambient_group_generators |l_gens|=" << l_gens.size() << "\n";
#endif
  return l_gens;
}

//
// The cosets of a Hecke double coset
//

/*
  H x H = union x_k H for H a group given by a membership test and a
  generating set. decompose(y) for y in H x H returns (k, h) with
  y = x_k h.
 */
template <typename T> struct HeckeCosets {
  MyMatrix<T> x;
  std::vector<MyMatrix<T>> l_coset;
  std::vector<MyMatrix<T>> l_coset_inv;
  std::function<bool(MyMatrix<T> const &)> f_member;
  size_t size() const { return l_coset.size(); }
  std::optional<std::pair<int, MyMatrix<T>>>
  find_decompose(MyMatrix<T> const &y) const {
    size_t n_coset = l_coset.size();
    for (size_t k = 0; k < n_coset; k++) {
      MyMatrix<T> h = l_coset_inv[k] * y;
      if (f_member(h)) {
        std::pair<int, MyMatrix<T>> pair{static_cast<int>(k), h};
        return pair;
      }
    }
    return {};
  }
  std::pair<int, MyMatrix<T>> decompose(MyMatrix<T> const &y) const {
    std::optional<std::pair<int, MyMatrix<T>>> opt = find_decompose(y);
    if (!opt) {
      std::cerr << "HECKE: The matrix y does not belong to the double coset\n";
      throw TerminalException{1};
    }
    return *opt;
  }
  /*
    The decomposition for the cosets of G x G with G the ambient group (the
    integral matrices of determinant +-1 preserving the T-space): for y in
    G x G, x_k^{-1} y preserves the T-space and has determinant +-1, so it
    is in G iff it is integral, which is much cheaper than f_member.
   */
  std::pair<int, MyMatrix<T>> decompose_ambient(MyMatrix<T> const &y) const {
    size_t n_coset = l_coset.size();
    for (size_t k = 0; k < n_coset; k++) {
      MyMatrix<T> h = l_coset_inv[k] * y;
      if (IsIntegralMatrix(h)) {
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (!f_member(h)) {
          std::cerr << "HECKE: decompose_ambient, an integral h is not in the ambient group\n";
          throw TerminalException{1};
        }
#endif
        return {static_cast<int>(k), h};
      }
    }
    std::cerr << "HECKE: The matrix y does not belong to the double coset (decompose_ambient)\n";
    throw TerminalException{1};
  }
};

template <typename T, typename Tint>
HeckeCosets<T>
enumerate_hecke_cosets(MyMatrix<T> const &x,
                       std::vector<MyMatrix<Tint>> const &l_gens,
                       std::function<bool(MyMatrix<T> const &)> f_member,
                       [[maybe_unused]] std::ostream &os,
                       size_t const &n_expected = 0) {
  HeckeCosets<T> hc;
  hc.x = x;
  hc.f_member = f_member;
  auto f_insert = [&](MyMatrix<T> const &y) -> void {
    if (hc.find_decompose(y)) {
      return;
    }
    hc.l_coset_inv.push_back(Inverse(y));
    hc.l_coset.push_back(y);
  };
  f_insert(x);
  size_t pos = 0;
  HeckeBudget budget;
  while (pos < hc.l_coset.size()) {
    if (hc.l_coset.size() > budget.max_hecke_cosets) {
      throw HeckeFailure{"enumerate_hecke_cosets: more than " + std::to_string(budget.max_hecke_cosets) + " cosets"};
    }
    MyMatrix<T> y = hc.l_coset[pos];
    for (auto &s : l_gens) {
      // When the number of cosets is known, the enumeration stops as soon
      // as they are all found.
      if (n_expected > 0 && hc.l_coset.size() == n_expected) {
        break;
      }
      MyMatrix<T> s_T = UniversalMatrixConversion<T, Tint>(s);
      MyMatrix<T> prod = s_T * y;
      f_insert(prod);
    }
    pos += 1;
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: enumerate_hecke_cosets |l_coset|=" << hc.l_coset.size() << "\n";
#endif
  return hc;
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_COMMON_H_
// clang-format on
