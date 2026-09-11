// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEROOTSYSTEM_H_
#define SRC_LATT_LATTICEROOTSYSTEM_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_ROOT_SYSTEM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_ROOT_SYSTEM
#endif

/*
  Recognition of the root system of an even positive definite lattice and
  the analytic order of its Weyl group. This is the foundational piece of
  the Vinberg decomposition Aut(L) = W(R) rtimes Aut(L, rho): the factor
  W(R) is a product of the Weyl groups of the irreducible components of
  the root system, whose orders are tabulated by Dynkin type, so it needs
  no group computation at all.

  Roots are the vectors of norm 2. A fundamental system (simple roots) is
  extracted by choosing a generic linear functional, taking the positive
  roots, and keeping those that are not the sum of two positive roots.
  The Dynkin diagram is the graph on the simple roots with an edge
  between s, t whenever <s, t> = -1 (the lattice being even, the only
  off-diagonal products of distinct simple roots are 0 and -1, so the
  system is simply laced: type A, D or E). The type of each connected
  component is read off its shape, and |W| multiplied over components.
 */

// Dynkin type of one irreducible component: the letter and the rank.
struct DynkinComponent {
  char letter;   // 'A', 'D', or 'E'
  int rank;
};

template <typename Tint> struct RootSystemData {
  int n_roots;                            // roots up to sign
  int rank;                               // rank of the root sublattice
  MyMatrix<Tint> SimpleRoots;             // rows: a fundamental system
  std::vector<DynkinComponent> components;
};

// The order of the Weyl group of one irreducible component.
inline void WeylOrderComponent(DynkinComponent const &c,
                               mpz_class &acc) {
  int k = c.rank;
  auto fact = [](int m) -> mpz_class {
    mpz_class r(1);
    for (int i = 2; i <= m; i++) {
      r *= i;
    }
    return r;
  };
  if (c.letter == 'A') {
    // |W(A_k)| = (k+1)!
    acc *= fact(k + 1);
  } else if (c.letter == 'D') {
    // |W(D_k)| = 2^{k-1} k!
    mpz_class w = fact(k);
    for (int i = 0; i < k - 1; i++) {
      w *= 2;
    }
    acc *= w;
  } else if (c.letter == 'E') {
    if (k == 6) {
      acc *= mpz_class(51840);
    } else if (k == 7) {
      acc *= mpz_class(2903040);
    } else if (k == 8) {
      acc *= mpz_class(696729600);
    } else {
      std::cerr << "ROOTSYS: E_" << k << " is not a root system\n";
      throw TerminalException{1};
    }
  } else {
    std::cerr << "ROOTSYS: unknown Dynkin letter " << c.letter << "\n";
    throw TerminalException{1};
  }
}

template <typename Tint>
mpz_class WeylGroupOrder(RootSystemData<Tint> const &rs) {
  mpz_class acc(1);
  for (auto &c : rs.components) {
    WeylOrderComponent(c, acc);
  }
  return acc;
}

// The scalar product of two integer vectors for the form GramMat.
template <typename T, typename Tint>
Tint RootScal(MyMatrix<T> const &GramMat, MyVector<Tint> const &u,
              MyVector<Tint> const &v) {
  int n = u.size();
  T sum(0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      sum += GramMat(i, j) * UniversalScalarConversion<T, Tint>(u(i)) *
             UniversalScalarConversion<T, Tint>(v(j));
    }
  }
  return UniversalScalarConversion<Tint, T>(sum);
}

/*
  Classify one connected component of the Dynkin graph given as an
  adjacency list on its vertices. The lattice is even so the diagram is
  simply laced; the shape determines the type:
    all degrees <= 2                       -> A_k (path; A_1 is one node)
    one degree-3 node, arms (1,1,c)        -> D_{c+3}
    one degree-3 node, arms (1,2,2/3/4)    -> E_6 / E_7 / E_8
 */
inline DynkinComponent
ClassifyComponent(std::vector<int> const &nodes,
                  std::vector<std::vector<int>> const &adj) {
  int k = nodes.size();
  if (k == 1) {
    return {'A', 1};
  }
  // degrees within the component
  int deg3 = -1;
  int n_deg3 = 0, n_deg1 = 0;
  for (int v : nodes) {
    int d = adj[v].size();
    if (d == 1) {
      n_deg1++;
    }
    if (d == 3) {
      n_deg3++;
      deg3 = v;
    }
    if (d > 3) {
      std::cerr << "ROOTSYS: a node of degree " << d
                << " is not simply laced ADE\n";
      throw TerminalException{1};
    }
  }
  if (n_deg3 == 0) {
    // a path: A_k
    return {'A', k};
  }
  if (n_deg3 == 1) {
    // measure the three arm lengths from the degree-3 node
    std::vector<int> arms;
    for (int start : adj[deg3]) {
      int len = 1, prev = deg3, cur = start;
      while (true) {
        int nxt = -1;
        for (int w : adj[cur]) {
          if (w != prev) {
            nxt = w;
            break;
          }
        }
        if (nxt == -1) {
          break;
        }
        prev = cur;
        cur = nxt;
        len++;
      }
      arms.push_back(len);
    }
    std::sort(arms.begin(), arms.end());
    if (arms[0] == 1 && arms[1] == 1) {
      return {'D', arms[2] + 3};
    }
    if (arms[0] == 1 && arms[1] == 2 && arms[2] == 2) {
      return {'E', 6};
    }
    if (arms[0] == 1 && arms[1] == 2 && arms[2] == 3) {
      return {'E', 7};
    }
    if (arms[0] == 1 && arms[1] == 2 && arms[2] == 4) {
      return {'E', 8};
    }
    std::cerr << "ROOTSYS: unrecognised arms (" << arms[0] << "," << arms[1]
              << "," << arms[2] << ")\n";
    throw TerminalException{1};
  }
  std::cerr << "ROOTSYS: " << n_deg3 << " branch nodes in one component\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
RootSystemData<Tint> ComputeRootSystem(MyMatrix<T> const &GramMat,
                                       std::ostream &os) {
  int n = GramMat.rows();
  RootSystemData<Tint> rs;
  // Roots up to sign: the norm-2 vectors.
  MyMatrix<Tint> Roots = T_ShortVector_fixed<T, Tint>(GramMat, T(2), os);
  rs.n_roots = Roots.rows();
  if (rs.n_roots == 0) {
    rs.rank = 0;
    rs.SimpleRoots = MyMatrix<Tint>(0, n);
    return rs;
  }

  // Positive roots via a generic functional f: f(r) != 0 for every root.
  // f = (1, c, c^2, ...) with c large enough that no root cancels it.
  std::vector<MyVector<Tint>> pos;
  {
    // choose the base so that |f(r)| never vanishes: bound by max entry.
    Tint bound(1);
    for (int i = 0; i < rs.n_roots; i++) {
      for (int j = 0; j < n; j++) {
        Tint a = T_abs(Roots(i, j));
        if (a > bound) {
          bound = a;
        }
      }
    }
    Tint base = 2 * bound * Tint(n) + 1;
    MyVector<Tint> f(n);
    Tint p(1);
    for (int j = 0; j < n; j++) {
      f(j) = p;
      p *= base;
    }
    for (int i = 0; i < rs.n_roots; i++) {
      MyVector<Tint> r = GetMatrixRow(Roots, i);
      Tint val(0);
      for (int j = 0; j < n; j++) {
        val += f(j) * r(j);
      }
#ifdef SANITY_CHECK_LATTICE_ROOT_SYSTEM
      if (val == 0) {
        std::cerr << "ROOTSYS: the functional vanished on a root\n";
        throw TerminalException{1};
      }
#endif
      pos.push_back(val > 0 ? r : MyVector<Tint>(-r));
    }
  }

  // Simple roots: positive roots not equal to a sum of two positive roots.
  std::unordered_map<MyVector<Tint>, int> posIndex;
  for (size_t i = 0; i < pos.size(); i++) {
    posIndex[pos[i]] = static_cast<int>(i);
  }
  std::vector<MyVector<Tint>> simple;
  for (size_t i = 0; i < pos.size(); i++) {
    bool is_sum = false;
    for (size_t j = 0; j < pos.size() && !is_sum; j++) {
      if (j == i) {
        continue;
      }
      MyVector<Tint> diff = pos[i] - pos[j];
      if (posIndex.count(diff) == 1) {
        is_sum = true;
      }
    }
    if (!is_sum) {
      simple.push_back(pos[i]);
    }
  }
  rs.SimpleRoots = MatrixFromVectorFamilyDim(n, simple);
  rs.rank = simple.size();
#ifdef DEBUG_LATTICE_ROOT_SYSTEM
  os << "ROOTSYS: n_roots=" << rs.n_roots << " rank=" << rs.rank << "\n";
#endif

  // The Dynkin graph: an edge between simple roots with product -1.
  int k = rs.rank;
  std::vector<std::vector<int>> adj(k);
  for (int a = 0; a < k; a++) {
    for (int b = a + 1; b < k; b++) {
      Tint sc = RootScal(GramMat, simple[a], simple[b]);
      if (sc == Tint(-1)) {
        adj[a].push_back(b);
        adj[b].push_back(a);
      }
#ifdef SANITY_CHECK_LATTICE_ROOT_SYSTEM
      else if (sc != Tint(0)) {
        std::cerr << "ROOTSYS: simple roots with product " << sc
                  << ", not simply laced\n";
        throw TerminalException{1};
      }
#endif
    }
  }

  // Connected components -> Dynkin types.
  std::vector<int> comp(k, -1);
  int n_comp = 0;
  for (int start = 0; start < k; start++) {
    if (comp[start] != -1) {
      continue;
    }
    std::vector<int> nodes;
    std::vector<int> stack{start};
    comp[start] = n_comp;
    while (!stack.empty()) {
      int v = stack.back();
      stack.pop_back();
      nodes.push_back(v);
      for (int w : adj[v]) {
        if (comp[w] == -1) {
          comp[w] = n_comp;
          stack.push_back(w);
        }
      }
    }
    rs.components.push_back(ClassifyComponent(nodes, adj));
    n_comp++;
  }
  return rs;
}

// clang-format off
#endif  // SRC_LATT_LATTICEROOTSYSTEM_H_
// clang-format on
