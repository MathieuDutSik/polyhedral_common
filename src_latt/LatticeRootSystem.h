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
  extracted by choosing a positive system -- the roots whose first nonzero
  coordinate is positive, the order of a generic linear functional -- and
  keeping the positive roots that are not the sum of two positive roots.
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
  MyMatrix<Tint> PositiveRoots;           // rows: the positive roots
  std::vector<DynkinComponent> components;
};

// The order of the Weyl group of one irreducible component.
template <typename Tgroup>
void WeylOrderComponent(DynkinComponent const &c,
                        typename Tgroup::Tint &acc) {
  using TintGroup = typename Tgroup::Tint;
  int k = c.rank;
  auto fact = [](int m) -> TintGroup {
    TintGroup r(1);
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
    TintGroup w = fact(k);
    w *= MyPow<TintGroup>(2, k - 1);
    acc *= w;
  } else if (c.letter == 'E') {
    if (k == 6) {
      acc *= TintGroup(51840);
    } else if (k == 7) {
      acc *= TintGroup(2903040);
    } else if (k == 8) {
      acc *= TintGroup(696729600);
    } else {
      std::cerr << "ROOTSYS: E_" << k << " is not a root system\n";
      throw TerminalException{1};
    }
  } else {
    std::cerr << "ROOTSYS: unknown Dynkin letter " << c.letter << "\n";
    throw TerminalException{1};
  }
}

template <typename Tgroup, typename Tint>
typename Tgroup::Tint WeylGroupOrder(RootSystemData<Tint> const &rs) {
  typename Tgroup::Tint acc(1);
  for (auto &c : rs.components) {
    WeylOrderComponent<Tgroup>(c, acc);
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
  int n_deg3 = 0;
#ifdef SANITY_CHECK_LATTICE_ROOT_SYSTEM
  int n_deg1 = 0;
#endif
  for (int v : nodes) {
    int d = adj[v].size();
#ifdef SANITY_CHECK_LATTICE_ROOT_SYSTEM
    if (d == 1) {
      n_deg1++;
    }
#endif
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
    // a path: A_k. All degrees 2 would be a cycle, the affine diagram of
    // A_{k-1}, which the root system of a definite lattice cannot have.
#ifdef SANITY_CHECK_LATTICE_ROOT_SYSTEM
    if (n_deg1 != 2) {
      std::cerr << "ROOTSYS: a component with no node of degree 3 has "
                << n_deg1 << " ends, so it is not a path\n";
      throw TerminalException{1};
    }
#endif
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
    std::cerr << "ROOTSYS: unrecognized arms (" << arms[0] << "," << arms[1]
              << "," << arms[2] << ")\n";
    throw TerminalException{1};
  }
  std::cerr << "ROOTSYS: " << n_deg3 << " branch nodes in one component\n";
  throw TerminalException{1};
}

// The roots of the lattice, the vectors of norm 2, one per antipodal pair.
template <typename T, typename Tint>
MyMatrix<Tint> GetLatticeRoots(MyMatrix<T> const &GramMat,
                               GramReduction reduction, std::ostream &os) {
  CVPSolver<T, Tint> solver = GetCVPSolver<T, Tint>(GramMat, reduction, os);
  return MatrixFromVectorFamilyDim(GramMat.rows(),
                                   solver.fixed_norm_vectors(T(2)));
}

/*
  The root system from the roots of GetLatticeRoots, for a caller that also
  needs them elsewhere (ComputeRootDecompositionFromRoots) and enumerates
  them once.
 */
template <typename T, typename Tint>
RootSystemData<Tint>
ComputeRootSystemFromRoots(MyMatrix<T> const &GramMat,
                           MyMatrix<Tint> const &Roots,
                           [[maybe_unused]] std::ostream &os) {
  int n = GramMat.rows();
  RootSystemData<Tint> rs;
  rs.n_roots = Roots.rows();
  if (rs.n_roots == 0) {
    rs.rank = 0;
    rs.SimpleRoots = MyMatrix<Tint>(0, n);
    rs.PositiveRoots = MyMatrix<Tint>(0, n);
    return rs;
  }

  /*
    Positive roots: those whose first nonzero coordinate is positive, as
    SignCanonicalizeVector makes them. That is the order of the functional
    f = (c^{n-1}, ..., c, 1) for any c larger than twice the largest
    coordinate of a root, the term of the first nonzero coordinate then
    dominating the others; it is the limit of generic functionals, so the
    positive roots form a positive system. Reading the sign of a coordinate
    gives them without forming f, whose entries would exceed the range of a
    machine integer at the ranks of interest.
   */
  std::vector<MyVector<Tint>> pos;
  for (int i = 0; i < rs.n_roots; i++) {
    pos.push_back(SignCanonicalizeVector(GetMatrixRow(Roots, i)));
  }

  /*
    Simple roots: the positive roots that are not the sum of two positive
    roots. A positive root alpha is such a sum when alpha - beta is a positive
    root for some positive root beta. For two roots, alpha - beta has norm
    4 - 2 <alpha, beta>, so it is a root exactly when <alpha, beta> = 1, and it
    is positive exactly when beta < alpha in the lexicographic order, which
    defines the positive roots and is translation invariant. So only the
    products of the positive roots are needed, one matrix product.
   */
  rs.PositiveRoots = MatrixFromVectorFamilyDim(n, pos);
  MyMatrix<Tint> G = UniversalMatrixConversion<Tint, T>(GramMat);
  MyMatrix<Tint> PosGram = rs.PositiveRoots * G * rs.PositiveRoots.transpose();
  std::vector<MyVector<Tint>> simple;
  size_t n_pos = pos.size();
  for (size_t i = 0; i < n_pos; i++) {
    bool is_sum = false;
    for (size_t j = 0; j < n_pos && !is_sum; j++) {
      if (PosGram(i, j) == 1 && pos[j] < pos[i]) {
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

  // The Dynkin graph: an edge between simple roots with product -1. The
  // products are those of the integral form, one matrix product for all the
  // pairs.
  int k = rs.rank;
  MyMatrix<Tint> GE = G * rs.SimpleRoots.transpose();
  std::vector<std::vector<int>> adj(k);
  for (int a = 0; a < k; a++) {
    for (int b = a + 1; b < k; b++) {
      Tint sc = simple[a].dot(GE.col(b));
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

template <typename T, typename Tint>
RootSystemData<Tint> ComputeRootSystem(MyMatrix<T> const &GramMat,
                                       GramReduction reduction,
                                       std::ostream &os) {
  MyMatrix<Tint> Roots = GetLatticeRoots<T, Tint>(GramMat, reduction, os);
  return ComputeRootSystemFromRoots<T, Tint>(GramMat, Roots, os);
}

// clang-format off
#endif  // SRC_LATT_LATTICEROOTSYSTEM_H_
// clang-format on
