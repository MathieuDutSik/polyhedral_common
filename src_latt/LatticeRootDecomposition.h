// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEROOTDECOMPOSITION_H_
#define SRC_LATT_LATTICEROOTDECOMPOSITION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include "LatticeRootSystem.h"
#include "InvariantVectorFamily.h"
#include "ClassicLLL.h"
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_ROOT_DECOMPOSITION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_ROOT_DECOMPOSITION
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

/*
  The decomposed characteristic vector family of L: the roots (a
  characteristic family for the root sublattice R) together with a
  characteristic full-rank family of the complement R^perp, lifted from
  the coordinates of R^perp into those of L. The family is invariant
  under Aut(L) --- the roots are, and R^perp with its characteristic
  vectors is canonically determined --- and of full rank, but spans only
  R (+) R^perp, i.e. it has index glue_index in L.

  The point is conditioning and size: for a lattice whose successive
  minima are spread out, the monolithic characteristic family of L climbs
  to high norms, whereas the complement R^perp, being root free and of
  smaller rank, has a small well-conditioned characteristic family of its
  own. On TestData/SlowCanonic this is the difference between a graph the
  canonical labelling can handle and one it cannot.
 */
template <typename T, typename Tint>
MyMatrix<Tint>
DecomposedCharacteristicFamily(MyMatrix<T> const &GramMat,
                               RootDecomposition<T, Tint> const &dec,
                               std::ostream &os) {
  int n = GramMat.rows();
  // The roots, in the coordinates of L, with both signs: fixed_norm_vectors
  // returns one vector per antipodal pair, and which one depends on the
  // basis, so only the full root system is an isometry invariant.
  MyMatrix<Tint> Roots = T_ShortVector_fixed<T, Tint>(GramMat, T(2), os);
  // The characteristic family of R^perp, lifted to L: a vector u in the
  // basis PerpBasis becomes u * PerpBasis in the coordinates of L.
  std::vector<MyVector<Tint>> family;
  std::unordered_set<MyVector<Tint>> seen;
  for (int i = 0; i < Roots.rows(); i++) {
    MyVector<Tint> v = GetMatrixRow(Roots, i);
    MyVector<Tint> v_neg = -v;
    for (MyVector<Tint> const &w : {v, v_neg}) {
      if (seen.insert(w).second) {
        family.push_back(w);
      }
    }
  }
  if (dec.PerpBasis.rows() > 0) {
    CanonicVectorFamily<Tint> perpFam =
        GetCanonicVectorFamily<T, Tint>(dec.PerpGram, os);
    MyMatrix<Tint> perpFull = perpFam.get_full();
    for (int i = 0; i < perpFull.rows(); i++) {
      MyVector<Tint> u = GetMatrixRow(perpFull, i);
      MyVector<Tint> uL = dec.PerpBasis.transpose() * u;
      family.push_back(uL);
    }
  }
  return MatrixFromVectorFamilyDim(n, family);
}

/*
  The canonical form of L assembled along the root decomposition: build
  the decomposed characteristic family and hand it to the full-rank
  canonical form of LatticeStabEquiCan.h, which canonically orders the
  family (steps 1 and 2, the two components) and canonicalizes the
  position of L relative to the sublattice it spans (step 3, the glue,
  through LinPolytopeIntegral_Canonicalization_Subspaces of
  src_group/MatrixGroup.h).

  The returned matrix B lies in GL_n(Z) and B * GramMat * B^T is the
  canonical Gram. As for every characteristic-set canonical form, the
  value depends on the family, so this canonical form may only be
  compared with another produced by this same function, never with
  ComputeCanonicalForm.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
ComputeCanonicalFormRootDecomposed(MyMatrix<T> const &GramMat,
                                   std::ostream &os) {
  RootDecomposition<T, Tint> dec =
      ComputeRootDecomposition<T, Tint, Tgroup>(GramMat, os);
  if (dec.n_roots == 0) {
    // Root free: no decomposition to exploit, fall back to the monolithic
    // characteristic-set canonical form.
    return ComputeCanonicalForm<T, Tint, Tgroup>(GramMat, os);
  }
  MyMatrix<Tint> family =
      DecomposedCharacteristicFamily<T, Tint>(GramMat, dec, os);
#ifdef SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
  if (RankMat(family) != GramMat.rows()) {
    std::cerr << "ROOTDEC: the decomposed family is not of full rank\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> ListMat{GramMat};
  return ComputeCanonicalFormFullRank_family<T, Tint, Tgroup>(ListMat, family,
                                                              os);
}

/*
  The Weyl-reduced family of L (Brandhorst-Hofmann, Automorphisms and
  canonical forms of integer lattices, Prop. 7.6): the simple roots E_rho of
  a positive system, with both signs, together with the characteristic
  family of R^perp lifted to L. It replaces the whole root system of
  DecomposedCharacteristicFamily by 2 * rank(R) vectors, so the graph and
  the automorphism group of the configuration lose the Weyl group W(R).

  The family depends on the chosen Weyl vector rho, but only up to
  isometry: an isometry f from L to L' maps rho to a Weyl vector of L',
  which is w rho' for some w in W(R'), and w fixes R'^perp pointwise, so
  w^{-1} f maps the family of (L, rho) onto that of (L', rho'). The
  canonical form of a lattice together with a family depends only on the
  isometry class of the pair, hence the choice of rho drops out.

  Both signs are needed because CanonicallyReorder_SHV_kernel requires a
  centrally symmetric family. The family has full rank and spans only
  R (+) R^perp; the glue is canonicalized by
  ComputeCanonicalFormFullRank_family as for the decomposed family.
 */
template <typename T, typename Tint>
MyMatrix<Tint>
WeylReducedCharacteristicFamily(RootSystemData<Tint> const &rs,
                                RootDecomposition<T, Tint> const &dec,
                                std::ostream &os) {
  int n = rs.SimpleRoots.cols();
  std::vector<MyVector<Tint>> family;
  for (int i = 0; i < rs.SimpleRoots.rows(); i++) {
    MyVector<Tint> v = GetMatrixRow(rs.SimpleRoots, i);
    family.push_back(v);
    family.push_back(-v);
  }
  if (dec.PerpBasis.rows() > 0) {
    CanonicVectorFamily<Tint> perpFam =
        GetCanonicVectorFamily<T, Tint>(dec.PerpGram, os);
    MyMatrix<Tint> perpFull = perpFam.get_full();
    for (int i = 0; i < perpFull.rows(); i++) {
      MyVector<Tint> u = GetMatrixRow(perpFull, i);
      family.push_back(dec.PerpBasis.transpose() * u);
    }
  }
  return MatrixFromVectorFamilyDim(n, family);
}

/*
  The canonical form of L from the Weyl-reduced family. Its values are not
  comparable with those of ComputeCanonicalFormRootDecomposed or of
  ComputeCanonicalForm, only with other values of this function.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
ComputeCanonicalFormWeylReduced(MyMatrix<T> const &GramMat, std::ostream &os) {
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  MicrosecondTime time;
#endif
  RootDecomposition<T, Tint> dec =
      ComputeRootDecomposition<T, Tint, Tgroup>(GramMat, os);
  if (dec.n_roots == 0) {
    return ComputeCanonicalForm<T, Tint, Tgroup>(GramMat, os);
  }
  RootSystemData<Tint> rs = ComputeRootSystem<T, Tint>(GramMat, os);
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylReduced, decomposition and root system|=" << time
     << "\n";
#endif
  MyMatrix<Tint> family =
      WeylReducedCharacteristicFamily<T, Tint>(rs, dec, os);
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylReduced, family of " << family.rows()
     << " vectors|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
  if (RankMat(family) != GramMat.rows()) {
    std::cerr << "ROOTDEC: the Weyl-reduced family is not of full rank\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> ListMat{GramMat};
  MyMatrix<Tint> B = ComputeCanonicalFormFullRank_family<T, Tint, Tgroup>(
      ListMat, family, os);
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylReduced, canonical form|=" << time << "\n";
#endif
  return B;
}

/*
  The Weyl-chamber family of L (Brandhorst-Hofmann, Prop. 7.6 and Lemma 7.5):
  for a Z-spanning characteristic family V of L, the vectors of V lying in
  the closed fundamental chamber P_rho that are not roots, together with the
  simple roots E_rho, all with both signs.

  V is invariant under Aut(L), hence under the Weyl group, and P_rho is a
  fundamental domain for W(R), so V and (V cap P_rho) span the same lattice
  modulo the root lattice <E_rho>, and the family spans L itself (Lemma
  7.5). The canonical basis is then the Hermite normal form of the
  canonically ordered family, as for ComputeCanonicalForm: neither the
  automorphism group of the configuration nor the glue step is needed. The
  choice of rho drops out for the reason given at
  WeylReducedCharacteristicFamily: an isometry carries chamber to chamber,
  and the Weyl group moves chambers onto each other while preserving V.

  With simple_both_signs the family is centrally symmetric, as the canonical
  form needs; its symmetries can then include -1 and elements of W(R)
  fixing the non-root part. Without it the only roots of the family are
  E_rho, so a symmetry preserves E_rho, hence rho, and the symmetry group
  is exactly Aut(L, rho), which is what the automorphism group and the
  isometry test need. Both versions span L.

  V is the Z-spanning root-iterated-shortest family of
  InvariantVectorFamily.h. V cap P_rho is computed directly (Remark 7.8),
  without forming V: at the top level of the recursion, the closest vector
  problems in the root lattice R and in its saturation are replaced by
  Lemma 4.1. The coset y + R has a unique vector of minimal norm in P_rho,
  the minuscule one, and it is reached from y by reflecting into the chamber
  and subtracting a positive root delta whenever <y, delta> >= 2, which
  shortens y; both moves stay in y + R.
 */
template <typename T, typename Tint>
MyMatrix<Tint> WeylChamberCharacteristicFamily(MyMatrix<T> const &GramMat,
                                               RootSystemData<Tint> const &rs,
                                               bool simple_both_signs,
                                               std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int n = GramMat.rows();
  int k = rs.SimpleRoots.rows();
  int n_pos = rs.PositiveRoots.rows();
  MyMatrix<Tint> G = UniversalMatrixConversion<Tint, T>(GramMat);
  // Columns: the linear forms x -> <x, e> for the simple and positive roots.
  MyMatrix<Tint> GE = G * rs.SimpleRoots.transpose();
  MyMatrix<Tint> GP = G * rs.PositiveRoots.transpose();
  auto scal = [&](MyVector<Tint> const &z, MyMatrix<Tint> const &M,
                  int j) -> Tint {
    Tint sum(0);
    for (int u = 0; u < n; u++) {
      sum += z(u) * M(u, j);
    }
    return sum;
  };
  // Babai rounding modulo R: z - round(c) E, with c the coordinates of the
  // projection of z on R in the basis E of simple roots. It brings the
  // component of z in R into a bounded box, so that the descent below takes
  // few steps whatever the basis; without it the descent subtracts one root
  // at a time from a component that can be huge in a non-reduced basis.
  MyMatrix<Tfield> E_F =
      UniversalMatrixConversion<Tfield, Tint>(rs.SimpleRoots);
  MyMatrix<Tfield> EG_F = E_F * UniversalMatrixConversion<Tfield, Tint>(G);
  MyMatrix<Tfield> EGEt = EG_F * E_F.transpose();
  MyMatrix<Tfield> EGEt_inv = Inverse(EGEt);
  auto babai = [&](MyVector<Tint> const &z) -> MyVector<Tint> {
    MyVector<Tfield> z_F = UniversalVectorConversion<Tfield, Tint>(z);
    MyVector<Tfield> c = EGEt_inv * (EG_F * z_F);
    MyVector<Tint> zr = z;
    for (int i = 0; i < k; i++) {
      Tint ci = UniversalNearestScalarInteger<Tint, Tfield>(c(i));
      if (ci != 0) {
        for (int u = 0; u < n; u++) {
          zr(u) -= ci * rs.SimpleRoots(i, u);
        }
      }
    }
    return zr;
  };
  auto reduce = [&](MyVector<Tint> z) -> MyVector<Tint> {
    z = babai(z);
    while (true) {
      bool moved = true;
      while (moved) {
        moved = false;
        for (int i = 0; i < k; i++) {
          Tint sc = scal(z, GE, i);
          if (sc < 0) {
            // The reflection in e_i: z - <z, e_i> e_i.
            for (int u = 0; u < n; u++) {
              z(u) -= sc * rs.SimpleRoots(i, u);
            }
            moved = true;
          }
        }
      }
      bool shortened = false;
      for (int j = 0; j < n_pos && !shortened; j++) {
        if (scal(z, GP, j) >= 2) {
          for (int u = 0; u < n; u++) {
            z(u) -= rs.PositiveRoots(j, u);
          }
          shortened = true;
        }
      }
      if (!shortened) {
        return z;
      }
    }
  };
  InnerSpanChamber<Tint> chamber =
      [&](MyVector<Tint> const &y, std::vector<MyVector<Tint>> const &reps)
      -> std::vector<MyVector<Tint>> {
    std::vector<MyVector<Tint>> best;
    Tint best_norm(0);
    for (auto &x : reps) {
      MyVector<Tint> z = reduce(y + x);
      Tint norm = z.dot(G * z);
      if (best.empty() || norm < best_norm) {
        best.clear();
        best_norm = norm;
      }
      if (norm == best_norm) {
        best.push_back(z);
      }
    }
    return best;
  };
  MyMatrix<Tfield> GramMat_F = UniversalMatrixConversion<Tfield, T>(GramMat);
  std::vector<MyVector<Tint>> cvL =
      inner_span_iterated_shortest_rec<Tfield, Tint>(GramMat_F, true, os,
                                                     &chamber);
  std::vector<MyVector<Tint>> family;
  std::unordered_set<MyVector<Tint>> seen;
  auto insert_pm = [&](MyVector<Tint> const &v) -> void {
    MyVector<Tint> v_neg = -v;
    for (MyVector<Tint> const &w : {v, v_neg}) {
      if (seen.insert(w).second) {
        family.push_back(w);
      }
    }
  };
  for (int i = 0; i < k; i++) {
    MyVector<Tint> e = GetMatrixRow(rs.SimpleRoots, i);
    if (simple_both_signs) {
      insert_pm(e);
    } else {
      seen.insert(e);
      family.push_back(e);
    }
  }
  for (auto &v : cvL) {
    insert_pm(v);
  }
#ifdef SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION_EXTENSIVE
  std::cerr << "ROOTDEC: SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION_EXTENSIVE, "
               "comparing with the chamber filter of the whole family\n";
  {
    MyMatrix<Tint> V = inner_span_root_iterated_shortest<T, Tint>(GramMat, os);
    std::unordered_set<MyVector<Tint>> filtered;
    for (int i = 0; i < k; i++) {
      MyVector<Tint> e = GetMatrixRow(rs.SimpleRoots, i);
      filtered.insert(e);
      if (simple_both_signs) {
        filtered.insert(-e);
      }
    }
    for (int i = 0; i < V.rows(); i++) {
      MyVector<Tint> v = GetMatrixRow(V, i);
      if (v.dot(G * v) == Tint(2)) {
        continue;
      }
      bool in_chamber = true;
      for (int j = 0; j < k; j++) {
        if (scal(v, GE, j) < 0) {
          in_chamber = false;
        }
      }
      if (in_chamber) {
        filtered.insert(v);
        filtered.insert(-v);
      }
    }
    if (filtered != seen) {
      std::cerr << "ROOTDEC: the direct chamber family has " << seen.size()
                << " vectors, the filtered one " << filtered.size() << "\n";
      throw TerminalException{1};
    }
  }
  std::cerr << "ROOTDEC: SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION_EXTENSIVE, "
               "done\n";
#endif
#ifdef DEBUG_LATTICE_ROOT_DECOMPOSITION
  os << "ROOTDEC: WeylChamber |cvL|=" << cvL.size()
     << " |family|=" << family.size() << "\n";
#endif
  return MatrixFromVectorFamilyDim(n, family);
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> ComputeCanonicalFormWeylChamber(MyMatrix<T> const &GramMat,
                                               std::ostream &os) {
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  MicrosecondTime time;
#endif
  // The family is built in an LLL-reduced basis: the canonical form does not
  // depend on the basis, but the cost of the closest vector problems and of
  // the Hermite normal form does.
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<T> const &GramRed = red.GramMatRed;
  RootSystemData<Tint> rs = ComputeRootSystem<T, Tint>(GramRed, os);
  if (rs.n_roots == 0) {
    return ComputeCanonicalForm<T, Tint, Tgroup>(GramMat, os);
  }
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylChamber, LLL and root system|=" << time << "\n";
#endif
  MyMatrix<Tint> family =
      WeylChamberCharacteristicFamily<T, Tint>(GramRed, rs, true, os);
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylChamber, family of " << family.rows()
     << " vectors|=" << time << "\n";
#endif
#ifdef SANITY_CHECK_LATTICE_ROOT_DECOMPOSITION
  MyMatrix<Tint> Zbasis = GetZbasis(family);
  if (Zbasis.rows() != GramMat.rows() ||
      T_abs(DeterminantMatBareiss(Zbasis)) != Tint(1)) {
    std::cerr << "ROOTDEC: the Weyl-chamber family does not span L\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> ListMat{GramRed};
  MyMatrix<Tint> B_red =
      ComputeCanonicalForm_inner<T, Tint>(ListMat, family, os);
#ifdef TIMINGS_LATTICE_ROOT_DECOMPOSITION
  os << "|ROOTDEC: WeylChamber, canonical form|=" << time << "\n";
#endif
  // B_red GramRed B_red^T is canonical and GramRed = Pmat GramMat Pmat^T.
  return B_red * red.Pmat;
}

// clang-format off
#endif  // SRC_LATT_LATTICEROOTDECOMPOSITION_H_
// clang-format on
