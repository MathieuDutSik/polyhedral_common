// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEFAMILYPERMUTATION_H_
#define SRC_LATT_LATTICEFAMILYPERMUTATION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

/*
  The automorphisms of a lattice as permutations of an invariant family of
  vectors of full rank: g (acting on the row vectors, g A g^T = A) permutes
  the vectors of the family V (both v and -v), v_i -> v_i g = v_{pi(i)},
  and the map g -> pi is an injective homomorphism (the family spans Q^n),
  with the products in the order of permutalib: pi_{g h} = pi_g pi_h. Back,
  g = B^{-1} C for B the rows of n vectors of V forming a basis of Q^n and C
  the rows of their images.

  This turns the automorphism group of a tuple (A, A P_1^T, ..., A P_k^T),
  P_i automorphisms of A, into the centralizer of the pi_{P_i} in the
  permutation group of O(A), and an isometry of two tuples into a
  conjugation, both by the partition backtrack of permutalib.
 */

template <typename Tint, typename Telt> struct LatticeFamilyPermutation {
  using Tidx = typename Telt::Tidx;
  int n;
  // the vectors, v and -v both
  std::vector<MyVector<Tint>> vectors;
  std::unordered_map<MyVector<Tint>, Tidx> index;
  // n vectors forming a basis, and the inverse of their matrix
  std::vector<Tidx> basis_idx;
  MyMatrix<Tint> basis;
  MyMatrix<Tint> basis_adj;
  Tint basis_det;

  // SHVhalf: one of v, -v per row, of full rank n.
  explicit LatticeFamilyPermutation(MyMatrix<Tint> const &SHVhalf)
      : n(SHVhalf.cols()) {
    int m = SHVhalf.rows();
    for (int i = 0; i < m; i++) {
      MyVector<Tint> v = GetMatrixRow(SHVhalf, i);
      for (int sign = 0; sign < 2; sign++) {
        MyVector<Tint> w = (sign == 0) ? v : MyVector<Tint>(-v);
        index[w] = static_cast<Tidx>(vectors.size());
        vectors.push_back(w);
      }
    }
    MyMatrix<Tint> Mfull(static_cast<int>(vectors.size()), n);
    for (size_t i = 0; i < vectors.size(); i++) {
      AssignMatrixRow(Mfull, static_cast<int>(i), vectors[i]);
    }
    // rows independent modulo p are independent over Q; otherwise over Q
    std::vector<int> rows =
        ColumnRankProfileModPrime(MyMatrix<Tint>(Mfull.transpose()));
    if (static_cast<int>(rows.size()) != n) {
      using Tfield = typename overlying_field<Tint>::field_type;
      rows = TMat_SelectRowCol(UniversalMatrixConversion<Tfield, Tint>(Mfull))
                 .ListRowSelect;
    }
    if (static_cast<int>(rows.size()) != n) {
      std::cerr << "LATT_PERM: the vector family does not have full rank\n";
      throw TerminalException{1};
    }
    basis = MyMatrix<Tint>(n, n);
    for (int i = 0; i < n; i++) {
      basis_idx.push_back(static_cast<Tidx>(rows[i]));
      AssignMatrixRow(basis, i, vectors[rows[i]]);
    }
    std::pair<MyMatrix<Tint>, Tint> pr = AdjugateDeterminant(basis);
    basis_adj = std::move(pr.first);
    basis_det = pr.second;
  }

  size_t degree() const { return vectors.size(); }

  // The permutation of g, or nothing when g does not preserve the family.
  std::optional<Telt> permutation(MyMatrix<Tint> const &g) const {
    std::vector<Tidx> img(vectors.size());
    for (size_t i = 0; i < vectors.size(); i++) {
      MyVector<Tint> w = g.transpose() * vectors[i];
      auto iter = index.find(w);
      if (iter == index.end()) {
        return {};
      }
      img[i] = iter->second;
    }
    return Telt(std::move(img));
  }

  // The matrix of a permutation coming from an automorphism.
  MyMatrix<Tint> matrix(Telt const &p) const {
    MyMatrix<Tint> C(n, n);
    for (int i = 0; i < n; i++) {
      AssignMatrixRow(C, i, vectors[p.at(basis_idx[i])]);
    }
    // B g = C, g = adj(B) C / det(B)
    MyMatrix<Tint> g = basis_adj * C;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        g(i, j) /= basis_det;
      }
    }
#ifdef SANITY_CHECK_LATTICE_FAMILY_PERMUTATION
    if (basis * g != C) {
      std::cerr << "LATT_PERM: the permutation is not of a matrix\n";
      throw TerminalException{1};
    }
#endif
    return g;
  }
};

// clang-format off
#endif  // SRC_LATT_LATTICEFAMILYPERMUTATION_H_
// clang-format on
