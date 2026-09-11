// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEROOTDECOMPOSITION_H_
#define SRC_LATT_LATTICEROOTDECOMPOSITION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include "ClassicLLL.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_ROOT_DECOMPOSITION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
#endif

/*
  Decomposition of a positive definite lattice along its root system, the
  computational skeleton of the canonical form in three steps:

    (1) the root sublattice R, spanned by the vectors of norm 2;
    (2) its orthogonal complement R^perp in L, root free;
    (3) the glue that reconstructs L from R (+) R^perp.

  This mirrors the decomposition of Aut(L) as W(R) rtimes Aut(L, rho)
  (Vinberg; Conway-Sloane): the two orthogonal pieces are canonicalized
  by the graph based ComputeCanonicalForm of src_latt/LatticeStabEquiCan.h,
  and the glue is canonicalized by
  LinPolytopeIntegral_Canonicalization_Subspaces of src_group/MatrixGroup.h.
  No Dynkin classification is required: the canonical form of the root
  sublattice is the canonical form of the small configuration of its
  norm 2 vectors.

  This header provides the decomposition itself, validated, and the
  canonical forms of the two components. The glue assembly into a single
  canonical Gram of L is the remaining step; the function to use is named
  in the comment at the end.
 */

template <typename T, typename Tint> struct RootDecomposition {
  // Basis (rows) of the root sublattice R and of its orthogonal
  // complement R^perp in L, in the coordinates of L.
  MyMatrix<Tint> RootBasis;
  MyMatrix<Tint> PerpBasis;
  // Gram matrices of the two components in those bases.
  MyMatrix<T> RootGram;
  MyMatrix<T> PerpGram;
  // The number of roots (norm 2 vectors up to sign), the rank of the
  // root sublattice, and the glue index [L : R (+) R^perp].
  int n_roots;
  int root_rank;
  Tint glue_index;
};

/*
  The integral orthogonal complement in L of the sublattice spanned by
  the rows of RootBasis, for the form GramMat: the vectors x of L with
  <x, r> = 0 for every root r, i.e. x (G R^T) = 0 over Z.
 */
template <typename T, typename Tint>
MyMatrix<Tint> IntegralOrthogonalComplement(MyMatrix<T> const &GramMat,
                                            MyMatrix<Tint> const &RootBasis) {
  MyMatrix<T> RootBasis_T = UniversalMatrixConversion<T, Tint>(RootBasis);
  // G R^T is n x root_rank; the left integral nullspace are the x with
  // x (G R^T) = 0.
  MyMatrix<T> Prod_T = GramMat * RootBasis_T.transpose();
  MyMatrix<Tint> Prod = UniversalMatrixConversion<Tint, T>(Prod_T);
  return NullspaceIntMat(Prod);
}

template <typename T, typename Tint, typename Tgroup>
RootDecomposition<T, Tint>
ComputeRootDecomposition(MyMatrix<T> const &GramMat, std::ostream &os) {
  int n = GramMat.rows();
  RootDecomposition<T, Tint> dec;

  // (1) The roots: vectors of norm 2, up to sign. For an even lattice of
  // minimum 2 these are the shortest vectors.
  MyMatrix<Tint> Roots = T_ShortVector_fixed<T, Tint>(GramMat, T(2), os);
  dec.n_roots = Roots.rows();
#ifdef DEBUG_LATTICE_ROOT_DECOMPOSITION
  os << "ROOTDEC: n_roots=" << dec.n_roots << "\n";
#endif

  if (dec.n_roots == 0) {
    // Root free: R is trivial, R^perp = L.
    dec.RootBasis = MyMatrix<Tint>(0, n);
    dec.PerpBasis = IdentityMat<Tint>(n);
    dec.RootGram = MyMatrix<T>(0, 0);
    dec.PerpGram = GramMat;
    dec.root_rank = 0;
    dec.glue_index = Tint(1);
    return dec;
  }

  // The root sublattice basis (a Z-basis of the span of the roots).
  dec.RootBasis = GetZbasis(Roots);
  dec.root_rank = dec.RootBasis.rows();

  // (2) The orthogonal complement.
  dec.PerpBasis = IntegralOrthogonalComplement<T, Tint>(GramMat, dec.RootBasis);
#ifdef SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
  if (dec.RootBasis.rows() + dec.PerpBasis.rows() != n) {
    std::cerr << "ROOTDEC: rank(R)=" << dec.RootBasis.rows()
              << " + rank(R^perp)=" << dec.PerpBasis.rows()
              << " != n=" << n << "\n";
    throw TerminalException{1};
  }
#endif

  // Component Gram matrices.
  MyMatrix<T> RootBasis_T = UniversalMatrixConversion<T, Tint>(dec.RootBasis);
  MyMatrix<T> PerpBasis_T = UniversalMatrixConversion<T, Tint>(dec.PerpBasis);
  dec.RootGram = RootBasis_T * GramMat * RootBasis_T.transpose();
  dec.PerpGram = PerpBasis_T * GramMat * PerpBasis_T.transpose();

  // (3) The glue index [L : R (+) R^perp]. The concatenated basis spans
  // R (+) R^perp; its index in L is |det| of the concatenation expressed
  // in the basis of L, which is L = Z^n here (identity basis matrix).
  MyMatrix<Tint> Concat(n, n);
  for (int i = 0; i < dec.RootBasis.rows(); i++) {
    Concat.row(i) = dec.RootBasis.row(i);
  }
  for (int i = 0; i < dec.PerpBasis.rows(); i++) {
    Concat.row(dec.RootBasis.rows() + i) = dec.PerpBasis.row(i);
  }
  dec.glue_index = T_abs(DeterminantMat(Concat));
#ifdef DEBUG_LATTICE_ROOT_DECOMPOSITION
  os << "ROOTDEC: root_rank=" << dec.root_rank
     << " perp_rank=" << dec.PerpBasis.rows()
     << " glue_index=" << dec.glue_index << "\n";
#endif
  return dec;
}

/*
  The canonical forms of the two components (step (1) and step (2) of the
  three step method), through the graph based canonical form. Returns the
  canonical Gram matrices of R and of R^perp; two isometric lattices have
  isometric components, so equal canonical component Grams are a necessary
  condition for isometry, and the pair (canonical RootGram, canonical
  PerpGram) is already a strong isometry invariant.

  Step (3), assembling these together with the glue into a single
  canonical Gram of L, is done by
  LinPolytopeIntegral_Canonicalization_Subspaces
  (src_group/MatrixGroup.h): it canonicalizes the position of L relative
  to the spanning family R (+) R^perp under the product of the component
  automorphism groups. That assembly, and the validation that the whole
  is an isometry invariant agreeing with ComputeCanonicalForm, is the
  remaining work.
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<MyMatrix<T>, MyMatrix<T>>
CanonicalComponentGrams(RootDecomposition<T, Tint> const &dec,
                        std::ostream &os) {
  auto canon = [&](MyMatrix<T> const &Gram) -> MyMatrix<T> {
    if (Gram.rows() == 0) {
      return Gram;
    }
    MyMatrix<Tint> B = ComputeCanonicalForm<T, Tint, Tgroup>(Gram, os);
    MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
    return B_T * Gram * B_T.transpose();
  };
  return {canon(dec.RootGram), canon(dec.PerpGram)};
}

// clang-format off
#endif  // SRC_LATT_LATTICEROOTDECOMPOSITION_H_
// clang-format on
