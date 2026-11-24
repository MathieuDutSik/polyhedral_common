// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSCOMPONENTS_H_
#define SRC_GENERA_MASSCOMPONENTS_H_

// clang-format off
#include "ResidueArithmeticLegacy.h"
#include "FiniteMatrixGroupFund.h"
#include "MAT_MatrixPolynomial.h"
#include "ResidueArithmetic.h"
#include "WeightMatrix.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include "ScalarAbelianField.h"
#include "CyclotomicArithmetic.h"
// clang-format on

#ifdef DEBUG
# define DEBUG_ISOMETRY_MASS
#endif

#ifdef SANITY_CHECK
# define SANITY_CHECK_ISOMETRY_MASS
#endif

/*
  The components of the mass of a pair (L, G) (MassIsometry.h): the
  sublattices L_d = L cap ker Phi_d(sigma) of an isometry sigma, or L_chi
  of an abelian group, with their Gram matrices, the actions of the
  generators and the coordinates of their discriminant groups
  (IsometryComponent), and the decompositions of an isometry
  (IsometryComponents) and of an abelian group (AbelianGroupComponents).
  The components of a non-abelian group are in MassGroupComponents.h.
 */

/*
  One component L_d: a basis (rows, coordinates in L), its Gram matrix, the
  action of sigma (sigma(basis_i) = sum_j sigma(i, j) basis_j), and the
  Smith form of the Gram matrix, which coordinatises D(L_d): with ROW G COL
  = diag(s), x in L_d^# (coordinates in the basis) is z ROW with z_i in
  (1/s_i) Z, and its residues are z_i s_i mod s_i. The s_i are unbounded
  integers; only the orders of the primary parts of D(L_d) become machine
  integers, in PrimaryQuadModule, where their size is checked.
 */
/*
  The kind of the centralizer of G on a component, which says how its mass
  is computed (ComponentMass, MassIsometry.h):
  - Orthogonal: d = 1, 2, G acts by +-1, the mass of the genus of L_d;
  - CyclotomicHermitian: d >= 3, the Hermitian lattice L_d over Z[zeta_d];
  - Multiplicity: End_G = Q, the genus of the multiplicity lattice X;
  - AbelianHermitian, AbelianQuadratic: End_G a field K with O_K X = X,
    X Hermitian over K (K CM) or quadratic (K totally real);
  - Precomputed: the mass found with the decomposition, End_G quaternionic
    or an order of K acting on X that is not O_K;
  - Unsupported: no method.
 */
enum class ComponentKind {
  Orthogonal,
  CyclotomicHermitian,
  Multiplicity,
  AbelianHermitian,
  AbelianQuadratic,
  Precomputed,
  Unsupported
};

template <typename T> struct IsometryComponent {
  // d = 0: a component of a non-abelian group on which it does not act
  // through a cyclic quotient
  int d;
  ComponentKind kind = ComponentKind::Unsupported;
  MyMatrix<T> basis;
  MyMatrix<T> gram;
  MyMatrix<T> sigma;
  // the actions of the further generators of the group (non-abelian case)
  std::vector<MyMatrix<T>> extra_sigma;
  // End_G = Q: the lattice X = L cap (l (x) U) of the multiplicity space,
  // whose genus gives the mass up to the local factors at p | |G|
  std::optional<MyMatrix<T>> multiplicity_gram;
  // the product over p | |G| of [O(X_p) : Stab(L_p)], when computed
  std::optional<T> multiplicity_factor;
  // End_G with a center K other than Q (AbelianHermitian,
  // AbelianQuadratic): the multiplicity lattice X with the actions of a
  // Z-basis of O_K
  std::optional<AbelianField> center_field;
  MyMatrix<T> center_gram;
  std::vector<MyMatrix<T>> center_omega;
  // Precomputed: the mass of the component when found (End_G = M_r(D), D
  // a totally definite quaternion algebra, MassQuaternionic.h; or an order
  // of K on X that is not O_K, by comparison with O_K X), to be multiplied
  // by multiplicity_factor
  std::optional<T> precomputed_mass;
  MyMatrix<T> row;
  MyMatrix<T> row_inv;
  std::vector<typename underlying_z_ring<T>::ring_type> elem;
};

// The order of P, finite (IntegralMatrixFiniteOrder: the crystallographic
// bound of GL_n(Z)).
template <typename T> int IsometryOrder(MyMatrix<T> const &P) {
  std::optional<size_t> order = IntegralMatrixFiniteOrder(P);
#ifdef SANITY_CHECK_ISOMETRY_MASS
  if (!order) {
    std::cerr << "ISOM_MASS: the isometry is not of finite order\n";
    throw TerminalException{1};
  }
#endif
  return static_cast<int>(*order);
}

/*
  The component of the saturated sublattice with basis E (rows), on which
  P acts with the d-th cyclotomic polynomial as minimal polynomial.
 */
template <typename T, typename Tint>
IsometryComponent<T> BuildIsometryComponent(MyMatrix<T> const &A,
                                            MyMatrix<T> const &E,
                                            MyMatrix<T> const &P, int d) {
  IsometryComponent<T> comp;
  comp.d = d;
  comp.kind = (d <= 2) ? ComponentKind::Orthogonal
                       : ComponentKind::CyclotomicHermitian;
  comp.basis = E;
  comp.gram = E * A * E.transpose();
  comp.sigma = E * P * A * E.transpose() * Inverse(comp.gram);
#ifdef SANITY_CHECK_ISOMETRY_MASS
  if (!IsIntegralMatrix(comp.sigma) || comp.sigma * E != E * P) {
    std::cerr << "ISOM_MASS: the component of d=" << d << " is not stable\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> G_i = UniversalMatrixConversion<Tint, T>(comp.gram);
  ResultSmithNormalForm<Tint> snf = SmithNormalForm(G_i);
  MyMatrix<Tint> D_i = snf.ROW * G_i * snf.COL;
  comp.row = UniversalMatrixConversion<T, Tint>(snf.ROW);
  comp.row_inv = Inverse(comp.row);
  for (int i = 0; i < D_i.rows(); i++) {
    Tint s = T_abs(D_i(i, i));
    comp.elem.push_back(UniversalScalarConversion<
                        typename underlying_z_ring<T>::ring_type, Tint>(s));
  }
  return comp;
}

/*
  The rows x of the saturated lattice spanned by the rows of E with
  x Phi_d(P) = 0, in the coordinates of L.
 */
template <typename T, typename Tint>
MyMatrix<T> CyclotomicKernel(MyMatrix<T> const &E, MyMatrix<T> const &P,
                             int d) {
  MyMatrix<T> Phi = EvaluatePolynomialMatrix(CyclotomicPolynomial(d), P);
  MyMatrix<T> EPhi = E * Phi;
  MyMatrix<Tint> EPhi_i = UniversalMatrixConversion<Tint, T>(EPhi);
  // NullspaceIntTrMat gives the x with M x^T = 0
  MyMatrix<Tint> K_i = NullspaceIntTrMat(MyMatrix<Tint>(EPhi_i.transpose()));
  if (K_i.rows() == 0) {
    return MyMatrix<T>(0, E.cols());
  }
  MyMatrix<T> K = UniversalMatrixConversion<T, Tint>(K_i) * E;
#ifdef SANITY_CHECK_ISOMETRY_MASS
  if (!IsZeroMatrix(MyMatrix<T>(K * Phi))) {
    std::cerr << "ISOM_MASS: the kernel of Phi_" << d << " is wrong\n";
    throw TerminalException{1};
  }
#endif
  return K;
}

template <typename T, typename Tint>
std::vector<IsometryComponent<T>> IsometryComponents(MyMatrix<T> const &A,
                                                     MyMatrix<T> const &P,
                                                     int m) {
  std::vector<IsometryComponent<T>> res;
  int n = A.rows();
  MyMatrix<T> Id = IdentityMat<T>(n);
  int tot = 0;
  for (int d = 1; d <= m; d++) {
    if (m % d != 0) {
      continue;
    }
    MyMatrix<T> E = CyclotomicKernel<T, Tint>(Id, P, d);
    if (E.rows() == 0) {
      continue;
    }
    tot += E.rows();
    res.push_back(BuildIsometryComponent<T, Tint>(A, E, P, d));
  }
#ifdef SANITY_CHECK_ISOMETRY_MASS
  if (tot != n) {
    std::cerr << "ISOM_MASS: the components have total rank " << tot
              << " instead of " << n << "\n";
    throw TerminalException{1};
  }
#endif
  return res;
}

/*
  The elements of the group generated by the ListP (acting on the rows),
  by closure; nothing beyond max_group elements.
 */
template <typename T>
std::vector<MyMatrix<T>> GroupElementsClosure(std::vector<MyMatrix<T>> const &ListP,
                                              size_t max_group) {
  int n = ListP[0].rows();
  auto flatten = [&](MyMatrix<T> const &M) -> std::vector<T> {
    std::vector<T> v;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v.push_back(M(i, j));
      }
    }
    return v;
  };
  std::vector<MyMatrix<T>> elts{IdentityMat<T>(n)};
  std::set<std::vector<T>> seen{flatten(elts[0])};
  for (size_t pos = 0; pos < elts.size(); pos++) {
    for (auto const &h : ListP) {
      MyMatrix<T> x = elts[pos] * h;
      if (seen.insert(flatten(x)).second) {
        elts.push_back(x);
        if (elts.size() > max_group) {
          std::cerr << "ISOM_MASS: the group has more than " << max_group
                    << " elements\n";
          throw TerminalException{1};
        }
      }
    }
  }
  return elts;
}

/*
  The components of a finite abelian group G. A rational irreducible
  representation of G is the Galois orbit of a character chi, which factors
  through the cyclic quotient G / ker chi of order e, faithfully; the
  orbit is determined by the orders of the chi(h), and V_chi is the
  intersection over h in G of the kernels of Phi_{ord chi(h)}(h). On V_chi
  the group acts through the powers of one element g_chi with chi(g_chi) of
  order e, so the centralizer of G in O(V_chi) is that of g_chi, and the
  component is that of the cyclic case for (L_chi, g_chi) and d = e.
 */
template <typename T, typename Tint>
std::vector<IsometryComponent<T>>
AbelianGroupComponents(MyMatrix<T> const &A,
                       std::vector<MyMatrix<T>> const &elts) {
  int n = A.rows();
  struct Block {
    MyMatrix<T> E;
    std::vector<int> sig;
  };
  std::vector<Block> blocks{{IdentityMat<T>(n), {}}};
  for (auto const &h : elts) {
    int o = IsometryOrder(h);
    std::vector<Block> next;
    for (auto const &blk : blocks) {
      for (int d = 1; d <= o; d++) {
        if (o % d != 0) {
          continue;
        }
        MyMatrix<T> K = CyclotomicKernel<T, Tint>(blk.E, h, d);
        if (K.rows() == 0) {
          continue;
        }
        std::vector<int> sig = blk.sig;
        sig.push_back(d);
        next.push_back({K, sig});
      }
    }
    blocks = std::move(next);
  }
  std::vector<IsometryComponent<T>> res;
  int tot = 0;
  for (auto const &blk : blocks) {
    size_t i_max = 0;
    for (size_t i = 0; i < blk.sig.size(); i++) {
      if (blk.sig[i] > blk.sig[i_max]) {
        i_max = i;
      }
    }
    IsometryComponent<T> comp = BuildIsometryComponent<T, Tint>(
        A, blk.E, elts[i_max], blk.sig[i_max]);
#ifdef SANITY_CHECK_ISOMETRY_MASS
    // every element acts on the component as a power of sigma
    for (auto const &h : elts) {
      MyMatrix<T> h_c = blk.E * h * A * blk.E.transpose() * Inverse(comp.gram);
      MyMatrix<T> pw = IdentityMat<T>(blk.E.rows());
      bool found = false;
      for (int k = 0; k < blk.sig[i_max] && !found; k++) {
        found = (pw == h_c);
        pw = pw * comp.sigma;
      }
      if (!found) {
        std::cerr << "ISOM_MASS: the group is not cyclic on a component\n";
        throw TerminalException{1};
      }
    }
#endif
    tot += blk.E.rows();
    res.push_back(comp);
  }
#ifdef SANITY_CHECK_ISOMETRY_MASS
  if (tot != n) {
    std::cerr << "ISOM_MASS: the components have total rank " << tot
              << " instead of " << n << "\n";
    throw TerminalException{1};
  }
#endif
  return res;
}

/*
  The mass of the genus of (L, sigma). herm_mass gives the masses of the
  Hermitian components by d; the total mass is returned only when every
  component has one.
 */
template <typename T>
void CheckIsometryInput(MyMatrix<T> const &A, MyMatrix<T> const &g) {
  int n = A.rows();
  if (!IsSymmetricMatrix(A) || !IsIntegralMatrix(A) || !IsIntegralMatrix(g)) {
    std::cerr << "ISOM_MASS: A must be symmetric integral and g integral\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++) {
    if (!IsInteger(A(i, i) / T(2))) {
      std::cerr << "ISOM_MASS: the lattice must be even\n";
      throw TerminalException{1};
    }
  }
  if (g.transpose() * A * g != A) {
    std::cerr << "ISOM_MASS: g^T A g is not A\n";
    throw TerminalException{1};
  }
  if (T_abs(DeterminantMat(g)) != T(1)) {
    std::cerr << "ISOM_MASS: g is not in GL_n(Z)\n";
    throw TerminalException{1};
  }
}

// clang-format off
#endif  // SRC_GENERA_MASSCOMPONENTS_H_
// clang-format on
