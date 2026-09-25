// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SPARSE_SOLVER_SPARSEEXACTSOLVER_H_
#define SRC_SPARSE_SOLVER_SPARSEEXACTSOLVER_H_

// clang-format off
#include "MAT_SparseMatrix.h"
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SPARSE_EXACT_SOLVER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SPARSE_EXACT_SOLVER
#endif

/*
  Exact solution of the sparse system x A = b, that is
    sum_i x_i A(i, .) = b,
  by Gaussian elimination in exact arithmetic. The equations are indexed by
  the columns j of A: sum_i A(i,j) x_i = b_j, and the unknowns by the
  rows i. The pivot is chosen in the sparsest active equation, preferring
  unit entries and unknowns occurring in few equations, which keeps the
  fill-in small for the boundary matrices of polyhedral complexes. The
  free unknowns are set to zero. An empty optional is returned when the
  system has no solution.
 */
template <typename T>
std::optional<MyVector<T>>
SparseSolutionMat_Exact(MySparseMatrix<T> const &A, MyVector<T> const &b,
                        [[maybe_unused]] std::ostream &os) {
  int n_unknown = A.rows();
  int n_eq = A.cols();
  std::vector<std::unordered_map<int, T>> eq(n_eq);
  std::vector<T> rhs(n_eq);
  for (int j = 0; j < n_eq; j++) {
    rhs[j] = b(j);
  }
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename MySparseMatrix<T>::InnerIterator it(A, k); it; ++it) {
      if (it.value() != 0) {
        eq[it.col()][it.row()] += it.value();
      }
    }
  }
  // The active equations containing an unknown.
  std::vector<std::unordered_set<int>> occ(n_unknown);
  for (int j = 0; j < n_eq; j++) {
    for (auto it = eq[j].begin(); it != eq[j].end();) {
      if (it->second == 0) {
        it = eq[j].erase(it);
      } else {
        occ[it->first].insert(j);
        ++it;
      }
    }
  }
  std::vector<uint8_t> active(n_eq, 1);
  int n_active = n_eq;
  std::vector<std::pair<int, int>> l_pivot;
  while (n_active > 0) {
    // The sparsest active equation.
    int e_best = -1;
    size_t siz_best = 0;
    for (int j = 0; j < n_eq; j++) {
      if (active[j] == 1) {
        size_t siz = eq[j].size();
        if (e_best == -1 || siz < siz_best) {
          e_best = j;
          siz_best = siz;
          if (siz == 0) {
            break;
          }
        }
      }
    }
    int e = e_best;
    if (eq[e].empty()) {
      if (rhs[e] != 0) {
#ifdef DEBUG_SPARSE_EXACT_SOLVER
        os << "SPARSEEXACT: inconsistent system, n_pivot=" << l_pivot.size() << "\n";
#endif
        return {};
      }
      active[e] = 0;
      n_active -= 1;
      continue;
    }
    // The pivot unknown: unit entries first, then few occurrences.
    int u_best = -1;
    bool unit_best = false;
    size_t occ_best = 0;
    for (auto &kv : eq[e]) {
      int u = kv.first;
      bool is_unit = (kv.second == 1 || kv.second == -1);
      size_t n_occ = occ[u].size();
      bool better = false;
      if (u_best == -1) {
        better = true;
      } else if (is_unit && !unit_best) {
        better = true;
      } else if (is_unit == unit_best && n_occ < occ_best) {
        better = true;
      }
      if (better) {
        u_best = u;
        unit_best = is_unit;
        occ_best = n_occ;
      }
    }
    int u = u_best;
    T pv = eq[e][u];
    // Elimination of u from the other active equations.
    std::vector<int> l_f(occ[u].begin(), occ[u].end());
    for (auto &f : l_f) {
      if (f == e) {
        continue;
      }
      T factor = eq[f][u] / pv;
      for (auto &kv : eq[e]) {
        int i = kv.first;
        if (i == u) {
          continue;
        }
        auto iter = eq[f].find(i);
        if (iter == eq[f].end()) {
          T val = -factor * kv.second;
          eq[f][i] = val;
          occ[i].insert(f);
        } else {
          iter->second -= factor * kv.second;
          if (iter->second == 0) {
            eq[f].erase(iter);
            occ[i].erase(f);
          }
        }
      }
      eq[f].erase(u);
      occ[u].erase(f);
      rhs[f] -= factor * rhs[e];
    }
    // The pivot equation leaves the active set.
    for (auto &kv : eq[e]) {
      occ[kv.first].erase(e);
    }
    active[e] = 0;
    n_active -= 1;
    l_pivot.emplace_back(e, u);
  }
  // Back substitution in reverse order of the pivots.
  MyVector<T> x = ZeroVector<T>(n_unknown);
  for (auto iter = l_pivot.rbegin(); iter != l_pivot.rend(); ++iter) {
    int e = iter->first;
    int u = iter->second;
    T val = rhs[e];
    for (auto &kv : eq[e]) {
      if (kv.first != u) {
        val -= kv.second * x(kv.first);
      }
    }
    x(u) = val / eq[e][u];
  }
#ifdef SANITY_CHECK_SPARSE_EXACT_SOLVER
  MyVector<T> prod = ZeroVector<T>(n_eq);
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename MySparseMatrix<T>::InnerIterator it(A, k); it; ++it) {
      prod(it.col()) += x(it.row()) * it.value();
    }
  }
  for (int j = 0; j < n_eq; j++) {
    if (prod(j) != b(j)) {
      std::cerr << "SPARSEEXACT: The solution does not satisfy the system\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef DEBUG_SPARSE_EXACT_SOLVER
  int n_nz = 0;
  for (int i = 0; i < n_unknown; i++) {
    if (x(i) != 0) {
      n_nz += 1;
    }
  }
  os << "SPARSEEXACT: n_unknown=" << n_unknown << " n_eq=" << n_eq
     << " n_pivot=" << l_pivot.size() << " n_nz(x)=" << n_nz << "\n";
#endif
  return x;
}

/*
  Consistency of x A = b tested by elimination modulo a prime. An
  inconsistency modulo p implies the inconsistency over the rationals;
  a consistency modulo p implies the rational one except with
  probability about 1/p. The arithmetic is on machine words, so the
  fill-in costs nothing like the rational elimination.
 */
inline int64_t residue_mod_p(mpq_class const &x, int64_t const &p) {
  int64_t num = mpz_fdiv_ui(x.get_num_mpz_t(), p);
  int64_t den = mpz_fdiv_ui(x.get_den_mpz_t(), p);
  // The inverse of den by Fermat
  int64_t res = 1, base = den, e = p - 2;
  while (e > 0) {
    if (e & 1) {
      res = (res * base) % p;
    }
    base = (base * base) % p;
    e >>= 1;
  }
  return (num * res) % p;
}

template <typename T>
int64_t residue_mod_p(T const &x, int64_t const &p) {
  mpq_class x_q = UniversalScalarConversion<mpq_class, T>(x);
  return residue_mod_p(x_q, p);
}

template <typename T>
bool SparseSystemConsistent_Mod(MySparseMatrix<T> const &A, MyVector<T> const &b,
                                [[maybe_unused]] std::ostream &os) {
  int64_t const p = 2147483647;
  auto inv_mod = [&](int64_t a) -> int64_t {
    int64_t res = 1, base = ((a % p) + p) % p, e = p - 2;
    while (e > 0) {
      if (e & 1) {
        res = (res * base) % p;
      }
      base = (base * base) % p;
      e >>= 1;
    }
    return res;
  };
  int n_unknown = A.rows();
  int n_eq = A.cols();
  std::vector<std::unordered_map<int, int64_t>> eq(n_eq);
  std::vector<int64_t> rhs(n_eq);
  for (int j = 0; j < n_eq; j++) {
    rhs[j] = residue_mod_p(b(j), p);
  }
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename MySparseMatrix<T>::InnerIterator it(A, k); it; ++it) {
      int64_t val = residue_mod_p(it.value(), p);
      int64_t &ref = eq[it.col()][it.row()];
      ref = (ref + val) % p;
    }
  }
  std::vector<std::unordered_set<int>> occ(n_unknown);
  for (int j = 0; j < n_eq; j++) {
    for (auto it = eq[j].begin(); it != eq[j].end();) {
      if (it->second == 0) {
        it = eq[j].erase(it);
      } else {
        occ[it->first].insert(j);
        ++it;
      }
    }
  }
  std::vector<uint8_t> active(n_eq, 1);
  int n_active = n_eq;
  while (n_active > 0) {
    int e = -1;
    size_t siz_best = 0;
    for (int j = 0; j < n_eq; j++) {
      if (active[j] == 1) {
        size_t siz = eq[j].size();
        if (e == -1 || siz < siz_best) {
          e = j;
          siz_best = siz;
          if (siz == 0) {
            break;
          }
        }
      }
    }
    if (eq[e].empty()) {
      if (rhs[e] != 0) {
        return false;
      }
      active[e] = 0;
      n_active -= 1;
      continue;
    }
    int u = -1;
    size_t occ_best = 0;
    for (auto &kv : eq[e]) {
      size_t n_occ = occ[kv.first].size();
      if (u == -1 || n_occ < occ_best) {
        u = kv.first;
        occ_best = n_occ;
      }
    }
    int64_t pv_inv = inv_mod(eq[e][u]);
    std::vector<int> l_f(occ[u].begin(), occ[u].end());
    for (auto &f : l_f) {
      if (f == e) {
        continue;
      }
      int64_t factor = (eq[f][u] * pv_inv) % p;
      for (auto &kv : eq[e]) {
        int i = kv.first;
        if (i == u) {
          continue;
        }
        int64_t sub = (factor * kv.second) % p;
        auto iter = eq[f].find(i);
        if (iter == eq[f].end()) {
          eq[f][i] = (p - sub) % p;
          occ[i].insert(f);
        } else {
          iter->second = ((iter->second - sub) % p + p) % p;
          if (iter->second == 0) {
            eq[f].erase(iter);
            occ[i].erase(f);
          }
        }
      }
      eq[f].erase(u);
      occ[u].erase(f);
      rhs[f] = ((rhs[f] - (factor * rhs[e]) % p) % p + p) % p;
    }
    for (auto &kv : eq[e]) {
      occ[kv.first].erase(e);
    }
    active[e] = 0;
    n_active -= 1;
  }
  return true;
}

// clang-format off
#endif  // SRC_SPARSE_SOLVER_SPARSEEXACTSOLVER_H_
// clang-format on
