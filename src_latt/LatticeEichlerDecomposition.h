// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEEICHLERDECOMPOSITION_H_
#define SRC_LATT_LATTICEEICHLERDECOMPOSITION_H_

// clang-format off
#include "Shvec_exact.h"
#include "ClassicLLL.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <map>
#include <numeric>
#include <optional>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_EICHLER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_EICHLER
#endif

/*
  The decomposition of a positive definite lattice into indecomposable
  orthogonal summands. By Eichler's theorem it exists and is unique: the
  summands themselves, not only their isometry classes, are determined by
  the lattice. This fails for indefinite lattices (U + <-1> = <1> + <-1> +
  <-1>).

  A vector v of L is decomposable if v = x + y with x, y non-zero in L and
  b(x, y) = 0. Then N(x), N(y) < N(v), so by induction on the norm every
  vector is a sum of indecomposable ones, and an indecomposable vector lies
  in a single summand of any orthogonal decomposition. Hence, for a set S of
  indecomposable vectors generating L, the connected components of the
  graph on S with an edge when b(v, w) != 0 span the summands: their spans
  are orthogonal, they generate L, and each lies in one summand of the
  finest decomposition and fills it, the summand being indecomposable.

  Such an S is built from an LLL-reduced basis by splitting: b(x, v - x) = 0
  iff N(x - v/2) = N(v)/4, so v is decomposable iff the sphere of diameter
  [0, v] holds a lattice point other than 0 and v, a fixed-distance search;
  if so v is replaced by x and v - x, which keeps S generating. A vector of
  norm below 2 min(L) is indecomposable (N(v) = N(x) + N(y) >= 2 min), so
  only the others are searched. The orthogonal splits keep the sum of the
  norms, so S stays of size at most sum_i N(b_i) / min(L). S need not be
  invariant under the isometries: the decomposition does not depend on it.
 */

template <typename T, typename Tint> struct EichlerComponent {
  // the rows: a basis of the summand, in the coordinates of the input
  MyMatrix<Tint> Basis;
  MyMatrix<T> Gram;
};

// Whether v is decomposable; if so, the x of a split v = x + (v - x). The
// enumeration stops at the first such x.
template <typename T, typename Tint>
std::optional<MyVector<Tint>>
EichlerSplitVector(CVPSolver<T, Tint> const &solver, MyMatrix<T> const &G,
                   MyVector<Tint> const &v) {
  int n = G.rows();
  MyVector<T> vT = UniversalVectorConversion<T, Tint>(v);
  T norm = vT.dot(G * vT);
  MyVector<T> center = vT / T(2);
  std::optional<MyVector<Tint>> result;
  auto f_ins = [&](MyVector<Tint> const &x) -> bool {
    bool is_zero = true, is_v = true;
    for (int i = 0; i < n; i++) {
      if (x(i) != 0) {
        is_zero = false;
      }
      if (x(i) != v(i)) {
        is_v = false;
      }
    }
    if (is_zero || is_v) {
      return true;
    }
    result = x;
    return false;
  };
  solver.fixed_dist_vectors_f(center, f_ins, norm / T(4));
  return result;
}

template <typename T, typename Tint>
std::vector<EichlerComponent<T, Tint>>
EichlerDecomposition(MyMatrix<T> const &G, std::ostream &os) {
  int n = G.rows();
  // the starting generating set: an LLL-reduced basis, whose reduction the
  // solver reuses rather than compute its own
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(G, os);
  CVPSolver<T, Tint> solver(G, red, os);
  T min = solver.shortest_vectors().min;
  std::vector<MyVector<Tint>> work;
  for (int i = 0; i < n; i++) {
    work.push_back(GetMatrixRow(red.Pmat, i));
  }
  auto norm_of = [&](MyVector<Tint> const &v) -> T {
    MyVector<T> vT = UniversalVectorConversion<T, Tint>(v);
    return vT.dot(G * vT);
  };
  // the splits, down to indecomposable vectors
  std::vector<MyVector<Tint>> S;
  [[maybe_unused]] size_t n_split = 0;
  while (!work.empty()) {
    MyVector<Tint> v = work.back();
    work.pop_back();
    if (norm_of(v) < T(2) * min) {
      S.push_back(v);
      continue;
    }
    std::optional<MyVector<Tint>> opt = EichlerSplitVector(solver, G, v);
    if (opt) {
      MyVector<Tint> w = v - *opt;
      work.push_back(*opt);
      work.push_back(w);
      n_split++;
    } else {
      S.push_back(v);
    }
  }
#ifdef DEBUG_LATTICE_EICHLER
  os << "EICHLER: n=" << n << " min=" << min << " |S|=" << S.size()
     << " splits=" << n_split << "\n";
#endif
  // the connected components of the non-orthogonality graph on S
  size_t m = S.size();
  std::vector<MyVector<T>> SG;
  for (auto const &v : S) {
    SG.push_back(G * UniversalVectorConversion<T, Tint>(v));
  }
  std::vector<size_t> parent(m);
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](size_t a) -> size_t {
    while (parent[a] != a) {
      parent[a] = parent[parent[a]];
      a = parent[a];
    }
    return a;
  };
  for (size_t i = 0; i < m; i++) {
    for (size_t j = i + 1; j < m; j++) {
      if (find(i) == find(j)) {
        continue;
      }
      T scal(0);
      for (int k = 0; k < n; k++) {
        scal += SG[i](k) * UniversalScalarConversion<T, Tint>(S[j](k));
      }
      if (scal != 0) {
        parent[find(i)] = find(j);
      }
    }
  }
  std::map<size_t, std::vector<size_t>> classes;
  for (size_t i = 0; i < m; i++) {
    classes[find(i)].push_back(i);
  }
  std::vector<EichlerComponent<T, Tint>> result;
  for (auto const &[root, members] : classes) {
    MyMatrix<Tint> Gens(members.size(), n);
    for (size_t r = 0; r < members.size(); r++) {
      AssignMatrixRow(Gens, r, S[members[r]]);
    }
    EichlerComponent<T, Tint> comp;
    comp.Basis = GetZbasis(Gens);
    MyMatrix<T> BT = UniversalMatrixConversion<T, Tint>(comp.Basis);
    comp.Gram = BT * G * BT.transpose();
    result.push_back(std::move(comp));
  }
  // a deterministic order: by rank, then by determinant
  std::sort(result.begin(), result.end(),
            [](EichlerComponent<T, Tint> const &a,
               EichlerComponent<T, Tint> const &b) -> bool {
              if (a.Gram.rows() != b.Gram.rows()) {
                return a.Gram.rows() < b.Gram.rows();
              }
              return DeterminantMat(a.Gram) < DeterminantMat(b.Gram);
            });
#ifdef SANITY_CHECK_LATTICE_EICHLER
  {
    // the summands are orthogonal and their bases form a basis of L
    int rank = 0;
    T prod_det(1);
    for (auto const &comp : result) {
      rank += comp.Gram.rows();
      prod_det *= DeterminantMat(comp.Gram);
    }
    if (rank != n || prod_det != DeterminantMat(G)) {
      std::cerr << "EICHLER: the summands do not form a basis of L, rank="
                << rank << " n=" << n << " prod_det=" << prod_det
                << " det=" << DeterminantMat(G) << "\n";
      throw TerminalException{1};
    }
    for (size_t a = 0; a < result.size(); a++) {
      for (size_t b = a + 1; b < result.size(); b++) {
        MyMatrix<T> Ba = UniversalMatrixConversion<T, Tint>(result[a].Basis);
        MyMatrix<T> Bb = UniversalMatrixConversion<T, Tint>(result[b].Basis);
        if (!IsZeroMatrix(MyMatrix<T>(Ba * G * Bb.transpose()))) {
          std::cerr << "EICHLER: two summands are not orthogonal\n";
          throw TerminalException{1};
        }
      }
    }
  }
#endif
  return result;
}

// clang-format off
#endif  // SRC_LATT_LATTICEEICHLERDECOMPOSITION_H_
// clang-format on
