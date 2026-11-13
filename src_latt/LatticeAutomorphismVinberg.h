// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEAUTOMORPHISMVINBERG_H_
#define SRC_LATT_LATTICEAUTOMORPHISMVINBERG_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "LatticeRootSystem.h"
#include "LatticeRootDecomposition.h"
#include "LatticeStabEquiCan.h"
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_AUTOMORPHISM_VINBERG
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_AUTOMORPHISM_VINBERG
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_AUTOMORPHISM_VINBERG
#endif

/*
  The automorphism group of an even positive definite lattice through the
  Vinberg decomposition Aut(L) = W(R) rtimes Aut(L, rho):

    * W(R), the Weyl group of the root system, is not computed as a
      group; its order is read off the Dynkin type (LatticeRootSystem.h)
      and its generators are the reflections in the simple roots.
    * The residual Aut(L, rho), the stabilizer of a Weyl vector, is the
      automorphism group of the small family consisting of the simple
      roots and a characteristic family of the orthogonal complement
      R^perp: an integral form-preserving map that permutes the simple
      roots stabilizes the fundamental chamber, hence fixes rho, and
      conversely. That family is far smaller than the full characteristic
      family of L, which is where the speedup comes from.

  The order is |W(R)| * |Aut(L, rho)| and the generators are the
  simple-root reflections together with the generators of Aut(L, rho);
  the huge factor W(R) is never enumerated.
 */

template <typename Tint, typename Tgroup> struct VinbergAutom {
  typename Tgroup::Tint order;           // |Aut(L)| = |W(R)| * |Aut(L,rho)|
  typename Tgroup::Tint weyl_order;      // |W(R)|
  typename Tgroup::Tint residual_order;  // |Aut(L, rho)|
  std::vector<MyMatrix<Tint>> ListGen;  // generators of Aut(L)
};

// The reflection in a root r, as an integral matrix acting on rows:
// s_r(x) = x - <x, r> r, with <r, r> = 2. As a matrix, g = I - (G r^T) r,
// integral since r and G are integral and <r,r> = 2 cancels the usual 2.
// G_int is the Gram matrix over the ring, which the caller converts once for
// all the roots it reflects in.
template <typename Tint>
MyMatrix<Tint> RootReflection(MyMatrix<Tint> const &G_int,
                              MyVector<Tint> const &r) {
  int n = G_int.rows();
  // w = G r^T (column)
  MyVector<Tint> w = G_int * r;
  MyMatrix<Tint> g = IdentityMat<Tint>(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g(i, j) -= w(i) * r(j);
    }
  }
  return g;
}

/*
  The residual family: the simple roots together with a characteristic
  family of R^perp lifted to the coordinates of L. Invariant under
  Aut(L, rho) and of full rank, but spanning only R (+) R^perp.
 */
template <typename T, typename Tint>
MyMatrix<Tint>
ResidualFamily(RootSystemData<Tint> const &rs,
               RootDecomposition<T, Tint> const &dec, std::ostream &os) {
  int n = rs.SimpleRoots.cols();
  std::vector<MyVector<Tint>> family;
  for (int i = 0; i < rs.SimpleRoots.rows(); i++) {
    family.push_back(GetMatrixRow(rs.SimpleRoots, i));
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
  The automorphism group of a root free lattice: no Weyl factor, Aut(L) is
  the residual itself, computed on the whole characteristic family.
 */
template <typename T, typename Tint, typename Tgroup>
VinbergAutom<Tint, Tgroup>
ComputeAutomorphismRootFree(MyMatrix<T> const &GramMat, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  IntAutomorphismGroup<Tint, Tgroup> aut =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
  VinbergAutom<Tint, Tgroup> res;
  res.weyl_order = TintGroup(1);
  res.residual_order = aut.order;
  res.order = res.residual_order;
  res.ListGen = std::move(aut.ListGen);
  return res;
}

// GramMat is LLL-reduced by the caller (LatticeAutomorphismGroup).
template <typename T, typename Tint, typename Tgroup>
VinbergAutom<Tint, Tgroup>
ComputeAutomorphismVinberg(MyMatrix<T> const &GramMat, std::ostream &os) {
  VinbergAutom<Tint, Tgroup> res;
  // The roots, enumerated once for the root system and the decomposition.
  MyMatrix<Tint> Roots =
      GetLatticeRoots<T, Tint>(GramMat, GramReduction::AlreadyReduced, os);
  RootSystemData<Tint> rs =
      ComputeRootSystemFromRoots<T, Tint>(GramMat, Roots, os);
  if (rs.n_roots == 0) {
    return ComputeAutomorphismRootFree<T, Tint, Tgroup>(GramMat, os);
  }
  res.weyl_order = WeylGroupOrder<Tgroup>(rs);
  RootDecomposition<T, Tint> dec =
      ComputeRootDecompositionFromRoots<T, Tint>(GramMat, Roots, os);
  MyMatrix<Tint> SHV = ResidualFamily<T, Tint>(rs, dec, os);
  std::vector<T> Vdiag(SHV.rows(), T(0));
  std::vector<MyMatrix<T>> ListMat{GramMat};
  // The automorphisms of the residual family = Aut(L, rho).
  IntAutomorphismGroup<Tint, Tgroup> resAut =
      GetIntAutomorphism_ListMat_Vdiag<T, Tint, Tgroup>(SHV, ListMat, Vdiag,
                                                        os);
  res.residual_order = resAut.order;
  res.order = res.weyl_order * res.residual_order;
  // Generators of Aut(L): the simple-root reflections (generating W(R))
  // and the residual generators.
  MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(GramMat);
  for (int i = 0; i < rs.SimpleRoots.rows(); i++) {
    res.ListGen.push_back(
        RootReflection<Tint>(G_int, GetMatrixRow(rs.SimpleRoots, i)));
  }
  for (auto &g : resAut.ListGen) {
    res.ListGen.push_back(g);
  }
#ifdef SANITY_CHECK_LATTICE_AUTOMORPHISM_VINBERG
  for (auto &g : res.ListGen) {
    if (g * G_int * g.transpose() != G_int) {
      std::cerr << "VINBERG: a generator does not preserve the form\n";
      throw TerminalException{1};
    }
  }
#endif
  return res;
}

/*
  The isometry test through the same decomposition. An isometry phi from L1
  to L2 maps the roots to the roots and the fundamental chamber of L1 to a
  chamber of L2, which an element w of W(R2) brings back to the fundamental
  chamber of L2. Then w phi maps the simple roots of L1 to those of L2, and
  the characteristic family of R1^perp to that of R2^perp, the family being
  intrinsic. So L1 and L2 are isometric exactly when the residual families
  are equivalent, which is a test on two small families. The returned P
  satisfies P * GramMat1 * P^T = GramMat2.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
ComputeIsometryVinberg(MyMatrix<T> const &GramMat1,
                       MyMatrix<T> const &GramMat2, std::ostream &os) {
  // The roots of each side, enumerated once for the root system and the
  // decomposition. Both matrices are LLL-reduced by the caller
  // (LatticeIsometry).
  MyMatrix<Tint> Roots1 =
      GetLatticeRoots<T, Tint>(GramMat1, GramReduction::AlreadyReduced, os);
  MyMatrix<Tint> Roots2 =
      GetLatticeRoots<T, Tint>(GramMat2, GramReduction::AlreadyReduced, os);
  RootSystemData<Tint> rs1 =
      ComputeRootSystemFromRoots<T, Tint>(GramMat1, Roots1, os);
  RootSystemData<Tint> rs2 =
      ComputeRootSystemFromRoots<T, Tint>(GramMat2, Roots2, os);
  if (rs1.n_roots != rs2.n_roots) {
    return {};
  }
  std::vector<MyMatrix<T>> ListMat1{GramMat1};
  std::vector<MyMatrix<T>> ListMat2{GramMat2};
  if (rs1.n_roots == 0) {
    // Root free: no Weyl factor, the ordinary test on the full family.
    return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                          os);
  }
  RootDecomposition<T, Tint> dec1 =
      ComputeRootDecompositionFromRoots<T, Tint>(GramMat1, Roots1, os);
  RootDecomposition<T, Tint> dec2 =
      ComputeRootDecompositionFromRoots<T, Tint>(GramMat2, Roots2, os);
  MyMatrix<Tint> SHV1 = ResidualFamily<T, Tint>(rs1, dec1, os);
  MyMatrix<Tint> SHV2 = ResidualFamily<T, Tint>(rs2, dec2, os);
  if (SHV1.rows() != SHV2.rows()) {
    return {};
  }
  std::vector<T> Vdiag(SHV1.rows(), T(0));
  std::optional<MyMatrix<Tint>> opt =
      TestIntEquivalence_ListMat_Vdiag<T, Tint, Tgroup>(
          SHV1, ListMat1, Vdiag, SHV2, ListMat2, Vdiag, os);
#ifdef SANITY_CHECK_LATTICE_AUTOMORPHISM_VINBERG
  if (opt) {
    MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(*opt);
    if (P_T * GramMat1 * P_T.transpose() != GramMat2) {
      std::cerr << "VINBERG: the isometry does not map the forms\n";
      throw TerminalException{1};
    }
  }
#endif
  return opt;
}

/*
  The same decomposition Aut(L) = W(R) rtimes Aut(L, rho), with Aut(L, rho)
  computed on the Weyl-chamber family (Brandhorst-Hofmann, Prop. 7.6) in its
  automorphism version: the simple roots E_rho, one sign, and the
  characteristic vectors of the chamber that are not roots, both signs (see
  WeylChamberCharacteristicFamily). Its symmetry group is exactly
  Aut(L, rho), and unlike ResidualFamily it spans L, so every permutation of
  it preserving the scalar products is realized by an integral matrix: the
  integral lift is direct, with no search over the glue. GramMat is
  LLL-reduced by the caller (LatticeAutomorphismGroup reduces its input).
 */
template <typename T, typename Tint, typename Tgroup>
VinbergAutom<Tint, Tgroup>
ComputeAutomorphismWeylChamber(MyMatrix<T> const &GramMat, std::ostream &os) {
#ifdef TIMINGS_LATTICE_AUTOMORPHISM_VINBERG
  MicrosecondTime time;
#endif
  MyMatrix<T> const &GramRed = GramMat;
  RootSystemData<Tint> rs =
      ComputeRootSystem<T, Tint>(GramRed, GramReduction::AlreadyReduced, os);
  if (rs.n_roots == 0) {
    return ComputeAutomorphismRootFree<T, Tint, Tgroup>(GramMat, os);
  }
  MyMatrix<Tint> SHV =
      WeylChamberCharacteristicFamily<T, Tint>(GramRed, rs, false, os);
#ifdef TIMINGS_LATTICE_AUTOMORPHISM_VINBERG
  os << "|VINBERG: WeylChamber, family of " << SHV.rows() << " vectors|="
     << time << "\n";
#endif
  std::vector<T> Vdiag(SHV.rows(), T(0));
  std::vector<MyMatrix<T>> ListMat{GramRed};
  IntAutomorphismGroup<Tint, Tgroup> resAut =
      GetIntAutomorphism_ListMat_Vdiag<T, Tint, Tgroup>(SHV, ListMat, Vdiag,
                                                        os);
#ifdef TIMINGS_LATTICE_AUTOMORPHISM_VINBERG
  os << "|VINBERG: WeylChamber, Aut(L, rho)|=" << time << "\n";
#endif
  VinbergAutom<Tint, Tgroup> res;
  res.weyl_order = WeylGroupOrder<Tgroup>(rs);
  res.residual_order = resAut.order;
  res.order = res.weyl_order * res.residual_order;
  MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(GramMat);
  for (int i = 0; i < rs.SimpleRoots.rows(); i++) {
    MyVector<Tint> r = GetMatrixRow(rs.SimpleRoots, i);
    res.ListGen.push_back(RootReflection<Tint>(G_int, r));
  }
  for (auto &g : resAut.ListGen) {
    res.ListGen.push_back(g);
  }
#ifdef SANITY_CHECK_LATTICE_AUTOMORPHISM_VINBERG
  for (auto &g : res.ListGen) {
    if (g * G_int * g.transpose() != G_int) {
      std::cerr << "VINBERG: a Weyl-chamber generator does not preserve the "
                   "form\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef TIMINGS_LATTICE_AUTOMORPHISM_VINBERG
  os << "|VINBERG: WeylChamber, order and generators|=" << time << "\n";
#endif
  return res;
}

/*
  The isometry test on the Weyl-chamber families. An isometry phi from L1 to
  L2, composed with the element of W(R2) bringing the image of the chamber
  of rho1 back to that of rho2, maps the family of L1 onto that of L2; and
  conversely a scalar product preserving bijection between the two families
  is an isometry, the families spanning the lattices. The returned P
  satisfies P * GramMat1 * P^T = GramMat2. Both matrices are LLL-reduced by
  the caller (LatticeIsometry reduces its input).
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
ComputeIsometryWeylChamber(MyMatrix<T> const &GramMat1,
                           MyMatrix<T> const &GramMat2, std::ostream &os) {
  if (GramMat1.rows() != GramMat2.rows()) {
    return {};
  }
  RootSystemData<Tint> rs1 =
      ComputeRootSystem<T, Tint>(GramMat1, GramReduction::AlreadyReduced, os);
  RootSystemData<Tint> rs2 =
      ComputeRootSystem<T, Tint>(GramMat2, GramReduction::AlreadyReduced, os);
  if (rs1.n_roots != rs2.n_roots || rs1.rank != rs2.rank) {
    return {};
  }
  if (rs1.n_roots == 0) {
    // Root free on both sides (the counts agree): the test on the whole
    // characteristic families, as ComputeIsometryVinberg does.
    std::vector<MyMatrix<T>> ListMat1{GramMat1};
    std::vector<MyMatrix<T>> ListMat2{GramMat2};
    return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                          os);
  }
  MyMatrix<Tint> SHV1 =
      WeylChamberCharacteristicFamily<T, Tint>(GramMat1, rs1, false, os);
  MyMatrix<Tint> SHV2 =
      WeylChamberCharacteristicFamily<T, Tint>(GramMat2, rs2, false, os);
  if (SHV1.rows() != SHV2.rows()) {
    return {};
  }
  std::vector<T> Vdiag(SHV1.rows(), T(0));
  std::vector<MyMatrix<T>> ListMat1{GramMat1};
  std::vector<MyMatrix<T>> ListMat2{GramMat2};
  std::optional<MyMatrix<Tint>> opt =
      TestIntEquivalence_ListMat_Vdiag<T, Tint, Tgroup>(
          SHV1, ListMat1, Vdiag, SHV2, ListMat2, Vdiag, os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> P = *opt;
#ifdef SANITY_CHECK_LATTICE_AUTOMORPHISM_VINBERG
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  if (P_T * GramMat1 * P_T.transpose() != GramMat2) {
    std::cerr << "VINBERG: the Weyl-chamber isometry does not map the forms\n";
    throw TerminalException{1};
  }
#endif
  return P;
}

// clang-format off
#endif  // SRC_LATT_LATTICEAUTOMORPHISMVINBERG_H_
// clang-format on
