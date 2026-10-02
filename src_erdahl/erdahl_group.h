// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_GROUP_H_
#define SRC_ERDAHL_ERDAHL_GROUP_H_

// clang-format off
#include "erdahl_polytope.h"
#include <optional>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_GROUP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ERDAHL_GROUP
#endif

/*
  Symmetries of Delaunay polyhedra with a nontrivial isotropy lattice.

  For a Delaunay polyhedron X = P + L(X), by Theorem "AutomAut_D" of the
  paper Aut(X) = Aff(X) x| Aut(P(X)) with
  * Aff(X) the transformations that are the identity modulo L(X): in the
    adapted coordinates (1, y, z) of erdahl_polyhedron.h they are
    (1, y, z) -> (1, y, z A + y B + b) with A in GL_d(Z).
  * Aut(P(X)) finite: the integral automorphisms of the reduced
    configuration of the rows (1, y) of the representatives.
  An element of Aut(P(X)) is lifted to Aut(X) by extending it by the
  identity on the z coordinates. If X is contained in Y then
  Aff(X) subset Aut(Y), so Stab(X, Y) = Aff(X) x| G_1 with G_1 the finite
  group of the elements of Aut(P(X)) whose lift preserves Y. The elements
  of G_1 preserve the scalar products of X and of Y on the reduced
  configuration of X (Theorem "ScalProdPreservation"), which is how they
  are computed, as in the GAP functions InfDel_PairStabilizer and
  InfDel_TripleEquivalence.

  The space W. For the full space every transformation is allowed. For the
  space of functions of center c, the allowed transformations are those
  fixing c, and all the polyhedra occurring are symmetric around c. If g
  maps X1 to X2, both symmetric around c, then g(c) and c are two centers
  of symmetry of X2, which differ by an element of L(X2) tensor Q, and an
  element a of Aff(X2) maps g(c) to c (see erdahl_correct_for_space). Since
  a preserves X2 and every polyhedron containing X2, equivalence under the
  full group implies equivalence under the stabilizer of c, and it is
  enough to correct the transformations obtained from the full group.
  Other spaces are not supported for polyhedra with a nontrivial isotropy
  lattice.
 */

template <typename T>
void erdahl_check_supported_space(ErdahlFunctionSpace<T> const &W) {
  if (W.center) {
    return;
  }
  if (static_cast<int>(W.basis.size()) == (W.n + 1) * (W.n + 2) / 2) {
    return;
  }
  std::cerr << "ERDAHL: the symmetry groups of Delaunay polyhedra with a "
               "nontrivial isotropy lattice are only available for the full "
               "space and the centered spaces, not for space="
            << W.name << "\n";
  throw TerminalException{1};
}

// Whether e -> e g maps Y onto Y, ab being the adapted basis of Y.
template <typename T, typename Tint>
bool erdahl_transformation_preserves(DelaunayPolyhedron<T, Tint> const &Y,
                                     ErdahlAdaptedBasis<Tint> const &ab,
                                     MyMatrix<Tint> const &g) {
  int n = erdahl_dimension(Y);
  for (int i = 0; i < Y.L.rows(); i++) {
    MyVector<Tint> dir = erdahl_direction_row<Tint>(GetMatrixRow(Y.L, i));
    MyVector<Tint> img = g.transpose() * dir;
    MyVector<Tint> img_red(n);
    for (int j = 0; j < n; j++) {
      img_red(j) = img(j + 1);
    }
    if (!SolutionIntMat(Y.L, img_red)) {
      return false;
    }
  }
  for (int i = 0; i < Y.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(Y.EXT, i);
    MyVector<Tint> img = g.transpose() * e;
    if (!erdahl_contains_point(Y, ab, img)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename Tint>
bool erdahl_transformation_preserves(DelaunayPolyhedron<T, Tint> const &Y,
                                     MyMatrix<Tint> const &g) {
  return erdahl_transformation_preserves(Y, erdahl_adapted_basis(Y), g);
}

// The image g(S) = { e g : e in S }, whose function is F o g^{-1}.
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_apply_transformation(DelaunayPolyhedron<T, Tint> const &S,
                            MyMatrix<Tint> const &g) {
  int n = erdahl_dimension(S);
  MyMatrix<Tint> EXT = S.EXT * g;
  MyMatrix<Tint> L(S.L.rows(), n);
  for (int i = 0; i < S.L.rows(); i++) {
    MyVector<Tint> dir = erdahl_direction_row<Tint>(GetMatrixRow(S.L, i));
    MyVector<Tint> img = g.transpose() * dir;
    for (int j = 0; j < n; j++) {
      L(i, j) = img(j + 1);
    }
  }
  ErdahlLatticeSet<Tint> can = erdahl_canonical_lattice_set<Tint>({EXT, L});
  MyMatrix<Tint> ginv = Inverse(g);
  MyMatrix<T> F = RemoveFractionMatrix(erdahl_transform_function(S.F, ginv));
  DelaunayPolyhedron<T, Tint> Simg{can.EXT, can.L, F};
#ifdef SANITY_CHECK_ERDAHL_GROUP
  for (int i = 0; i < Simg.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(Simg.EXT, i);
    if (EvaluationQuadForm<T, Tint>(F, e) != 0) {
      std::cerr << "ERDAHL: the transformed function does not vanish on the "
                   "transformed polyhedron\n";
      throw TerminalException{1};
    }
  }
#endif
  return Simg;
}

/*
  The reduced configuration of X relative to a chain of polyhedra
  containing it: the rows (1, y) of the representatives of X in the
  adapted coordinates, the scalar products of X and the scalar products of
  each of the super polyhedra, expressed on those coordinates.
 */
template <typename T, typename Tint> struct ErdahlChainConfig {
  ErdahlAdaptedBasis<Tint> ab;
  MyMatrix<Tint> EXTred;
  std::vector<MyMatrix<T>> ListMat;
};

/*
  The data of a super polyhedron Y used by the chain configurations, the
  same for all the polyhedra X it contains: it is computed once by the
  callers testing many X against the same supers.
 */
template <typename T, typename Tint> struct ErdahlSuperData {
  DelaunayPolyhedron<T, Tint> Y;
  ErdahlAdaptedBasis<Tint> abY;
  // The discriminant matrix of Y in its adapted coordinates.
  MyMatrix<T> MY;
  // The first p_Y + 1 columns of the inverse of the adapted basis of Y.
  MyMatrix<T> InvY_red;
};

template <typename T, typename Tint>
ErdahlSuperData<T, Tint> erdahl_super_data(DelaunayPolyhedron<T, Tint> const &Y) {
  ErdahlAdaptedBasis<Tint> abY = erdahl_adapted_basis(Y);
  int np1 = abY.n + 1;
  MyMatrix<Tint> EXTredY = erdahl_reduced_ext(Y, abY);
  MyMatrix<T> MY = erdahl_discriminant_matrix<T, Tint>(EXTredY);
  MyMatrix<T> InvY_T = UniversalMatrixConversion<T, Tint>(abY.AffBasisInv);
  MyMatrix<T> InvY_red(np1, abY.p + 1);
  for (int i = 0; i < np1; i++) {
    for (int j = 0; j <= abY.p; j++) {
      InvY_red(i, j) = InvY_T(i, j);
    }
  }
  return {Y, abY, MY, InvY_red};
}

template <typename T, typename Tint>
std::vector<ErdahlSuperData<T, Tint>>
erdahl_list_super_data(std::vector<DelaunayPolyhedron<T, Tint>> const &supers) {
  std::vector<ErdahlSuperData<T, Tint>> l_sd;
  for (auto &Y : supers) {
    l_sd.push_back(erdahl_super_data(Y));
  }
  return l_sd;
}

template <typename T, typename Tint>
bool erdahl_preserves_all(std::vector<ErdahlSuperData<T, Tint>> const &l_sd,
                          MyMatrix<Tint> const &g) {
  for (auto &sd : l_sd) {
    if (!erdahl_transformation_preserves(sd.Y, sd.abY, g)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename Tint>
ErdahlChainConfig<T, Tint>
erdahl_chain_config(DelaunayPolyhedron<T, Tint> const &X,
                    std::vector<ErdahlSuperData<T, Tint>> const &l_sd) {
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(X);
  MyMatrix<Tint> EXTred = erdahl_reduced_ext(X, ab);
  std::vector<MyMatrix<T>> ListMat{erdahl_discriminant_matrix<T, Tint>(EXTred)};
  int np1 = ab.n + 1;
  MyMatrix<T> AffBasis_T = UniversalMatrixConversion<T, Tint>(ab.AffBasis);
  MyMatrix<T> topX(ab.p + 1, np1);
  for (int i = 0; i <= ab.p; i++) {
    for (int j = 0; j < np1; j++) {
      topX(i, j) = AffBasis_T(i, j);
    }
  }
  for (auto &sd : l_sd) {
    MyMatrix<T> S = topX * sd.InvY_red;
    ListMat.push_back(S * sd.MY * S.transpose());
  }
  return {ab, EXTred, ListMat};
}

template <typename T, typename Tint>
ErdahlChainConfig<T, Tint>
erdahl_chain_config(DelaunayPolyhedron<T, Tint> const &X,
                    std::vector<DelaunayPolyhedron<T, Tint>> const &supers) {
  return erdahl_chain_config(X, erdahl_list_super_data(supers));
}

/*
  An invariant of X under the transformations preserving the super
  polyhedra: the hash of the matrices of scalar products of the chain
  configuration, with the rows put in the canonical order of
  Canonicalization_ListMat_Vdiag. Equivalent polyhedra (for
  erdahl_equivalence with the same supers) have the same invariant, and
  equal invariants mean equivalence up to a rational transformation
  preserving the scalar products, so the integral equivalence test is
  almost only run on actually equivalent polyhedra.
 */
template <typename T, typename Tint>
size_t erdahl_invariant_hash(ErdahlChainConfig<T, Tint> const &cfg,
                             std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(cfg.EXTred);
  int n_ext = EXT_T.rows();
  std::vector<T> Vdiag(n_ext, T(0));
  std::vector<uint32_t> ord = Canonicalization_ListMat_Vdiag<T, Tfield, uint32_t>(
      EXT_T, cfg.ListMat, Vdiag, THRESHOLD_USE_SUBSET_SCHEME_CANONIC, os);
  std::vector<T> l_val{T(cfg.ab.d), T(n_ext)};
  for (auto &M : cfg.ListMat) {
    for (int i = 0; i < n_ext; i++) {
      MyVector<T> Vi = GetMatrixRow(EXT_T, ord[i]);
      MyVector<T> MVi = M * Vi;
      for (int j = i; j < n_ext; j++) {
        MyVector<T> Vj = GetMatrixRow(EXT_T, ord[j]);
        l_val.push_back(MVi.dot(Vj));
      }
    }
  }
  return std::hash<std::vector<T>>()(l_val);
}

template <typename T, typename Tint>
size_t erdahl_invariant_hash(DelaunayPolyhedron<T, Tint> const &X,
                             std::vector<DelaunayPolyhedron<T, Tint>> const &supers,
                             std::ostream &os) {
  return erdahl_invariant_hash(erdahl_chain_config(X, supers), os);
}

// The lift of the reduced transformation h from X1 to X2 (same d).
template <typename Tint>
MyMatrix<Tint> erdahl_lift(MyMatrix<Tint> const &h,
                           ErdahlAdaptedBasis<Tint> const &ab1,
                           ErdahlAdaptedBasis<Tint> const &ab2) {
  int np1 = ab1.n + 1;
  MyMatrix<Tint> big = IdentityMat<Tint>(np1);
  for (int i = 0; i <= ab1.p; i++) {
    for (int j = 0; j <= ab1.p; j++) {
      big(i, j) = h(i, j);
    }
  }
  return ab1.AffBasisInv * big * ab2.AffBasis;
}

template <typename T, typename Tint>
bool erdahl_preserves_all(std::vector<DelaunayPolyhedron<T, Tint>> const &supers,
                          MyMatrix<Tint> const &g) {
  for (auto &Y : supers) {
    if (!erdahl_transformation_preserves(Y, g)) {
      return false;
    }
  }
  return true;
}

// A unimodular matrix whose first row is the primitive vector V.
template <typename Tint>
MyMatrix<Tint> erdahl_unimodular_with_first_row(MyVector<Tint> const &V) {
  int d = V.size();
  MyMatrix<Tint> M1(1, d);
  AssignMatrixRow(M1, 0, V);
  MyMatrix<Tint> Compl = SubspaceCompletionInt(M1, d);
  return Concatenate(M1, Compl);
}

/*
  A transformation g mapping X1 onto X2 (and preserving the super
  polyhedra) is changed into a g a with a in Aff(X2) that preserves W.
  For the centered spaces, a maps g(c) to c: in the adapted coordinates of
  X2, g(c) and c have the same y coordinates and a is
  (1, y, z) -> (1, y, z A + y B + b).
  * If c_y is not integral, A = Id and B, b are chosen such that
    c_y B + b = c_z - q_z, which lies in (1/2) Z^d.
  * If c_y is integral, then c_z and q_z are not (c is not a lattice
    point), and A in GL_d(Z) maps 2 q_z to 2 c_z modulo 2, which is
    possible since GL_d(Z) is transitive on the nonzero vectors of
    (Z/2Z)^d.
 */
template <typename T, typename Tint>
MyMatrix<Tint> erdahl_correct_for_space(ErdahlFunctionSpace<T> const &W,
                                        MyMatrix<Tint> const &g,
                                        DelaunayPolyhedron<T, Tint> const &X2) {
  if (!W.center) {
    erdahl_check_supported_space(W);
    return g;
  }
  int n = W.n;
  MyVector<T> e_c(n + 1);
  e_c(0) = 1;
  for (int i = 0; i < n; i++) {
    e_c(i + 1) = (*W.center)(i);
  }
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
  MyVector<T> e_q = g_T.transpose() * e_c;
  if (e_q == e_c) {
    return g;
  }
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(X2);
  int p = ab.p;
  int d = ab.d;
  MyMatrix<T> Inv_T = UniversalMatrixConversion<T, Tint>(ab.AffBasisInv);
  MyVector<T> u_q = Inv_T.transpose() * e_q;
  MyVector<T> u_c = Inv_T.transpose() * e_c;
  for (int i = 0; i <= p; i++) {
    if (u_q(i) != u_c(i)) {
      std::cerr << "ERDAHL: the image of the center differs from the center "
                   "outside of the isotropy lattice\n";
      throw TerminalException{1};
    }
  }
  MyVector<T> qz(d), cz(d), cy(p);
  for (int i = 0; i < p; i++) {
    cy(i) = u_c(1 + i);
  }
  for (int k = 0; k < d; k++) {
    qz(k) = u_q(1 + p + k);
    cz(k) = u_c(1 + p + k);
  }
  MyMatrix<Tint> A = IdentityMat<Tint>(d);
  MyMatrix<Tint> B = ZeroMatrix<Tint>(p, d);
  MyVector<T> b_T(d);
  int j_half = -1;
  for (int i = 0; i < p; i++) {
    if (!IsInteger(cy(i))) {
      j_half = i;
    }
  }
  if (j_half >= 0) {
    MyVector<T> delta = cz - qz;
    for (int k = 0; k < d; k++) {
      B(j_half, k) = UniversalScalarConversion<Tint, T>(2 * delta(k));
    }
    for (int k = 0; k < d; k++) {
      b_T(k) = delta(k) - cy(j_half) * 2 * delta(k);
    }
  } else {
    auto mod2 = [&](MyVector<T> const &V) -> MyVector<Tint> {
      MyVector<Tint> W2(d);
      for (int k = 0; k < d; k++) {
        Tint val = UniversalScalarConversion<Tint, T>(2 * V(k));
        Tint res = ResInt(val, Tint(2));
        W2(k) = res;
      }
      return W2;
    };
    MyVector<Tint> u = mod2(qz);
    MyVector<Tint> v = mod2(cz);
    if (IsZeroVector(u) || IsZeroVector(v)) {
      std::cerr << "ERDAHL: the center should not be a lattice point\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> Mu = erdahl_unimodular_with_first_row(u);
    MyMatrix<Tint> Mv = erdahl_unimodular_with_first_row(v);
    A = Inverse(Mu) * Mv;
    MyMatrix<T> A_T = UniversalMatrixConversion<T, Tint>(A);
    MyVector<T> qzA = A_T.transpose() * qz;
    b_T = cz - qzA;
  }
  MyMatrix<Tint> big = IdentityMat<Tint>(n + 1);
  for (int k = 0; k < d; k++) {
    big(0, 1 + p + k) = UniversalScalarConversion<Tint, T>(b_T(k));
    for (int i = 0; i < p; i++) {
      big(1 + i, 1 + p + k) = B(i, k);
    }
    for (int l = 0; l < d; l++) {
      big(1 + p + l, 1 + p + k) = A(l, k);
    }
  }
  MyMatrix<Tint> a = ab.AffBasisInv * big * ab.AffBasis;
  MyMatrix<Tint> gc = g * a;
#ifdef SANITY_CHECK_ERDAHL_GROUP
  MyMatrix<T> gc_T = UniversalMatrixConversion<T, Tint>(gc);
  MyVector<T> e_img = gc_T.transpose() * e_c;
  if (e_img != e_c) {
    std::cerr << "ERDAHL: the corrected transformation does not fix the "
                 "center\n";
    throw TerminalException{1};
  }
  if (!erdahl_transformation_preserves(X2, a)) {
    std::cerr << "ERDAHL: the correction does not preserve X2\n";
    throw TerminalException{1};
  }
#endif
  return gc;
}

/*
  The finite group of the reduced automorphisms of X whose lifts preserve
  all the super polyhedra, as a permutation group on the rows of EXTred,
  together with the lifts of its generators.
 */
template <typename T, typename Tint, typename Tgroup> struct ErdahlFiniteGroup {
  ErdahlChainConfig<T, Tint> cfg;
  Tgroup grp;
};

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> erdahl_lift_permutation(ErdahlChainConfig<T, Tint> const &cfg,
                                       typename Tgroup::Telt const &elt) {
  MyMatrix<T> EXTred_T = UniversalMatrixConversion<T, Tint>(cfg.EXTred);
  MyMatrix<T> h_T = FindTransformation(EXTred_T, EXTred_T, elt);
  MyMatrix<Tint> h = UniversalMatrixConversion<Tint, T>(h_T);
  return erdahl_lift(h, cfg.ab, cfg.ab);
}

template <typename T, typename Tint, typename Tgroup>
ErdahlFiniteGroup<T, Tint, Tgroup>
erdahl_finite_group(DelaunayPolyhedron<T, Tint> const &X,
                    std::vector<DelaunayPolyhedron<T, Tint>> const &supers,
                    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  ErdahlChainConfig<T, Tint> cfg = erdahl_chain_config(X, supers);
  int n_ext = cfg.EXTred.rows();
  std::vector<T> Vdiag(n_ext, T(0));
  std::vector<MyMatrix<Tint>> l_red =
      GetIntAutomorphism_ListMat_Vdiag<T, Tint, Tgroup>(cfg.EXTred, cfg.ListMat,
                                                        Vdiag, os);
  std::vector<Telt> l_perm;
  bool all_ok = true;
  for (auto &h : l_red) {
    std::optional<Telt> opt =
        erdahl_permutation_of_rows<Tint, Telt>(cfg.EXTred, cfg.EXTred, h);
    if (!opt) {
      std::cerr << "ERDAHL: a reduced automorphism does not preserve the "
                   "reduced configuration\n";
      throw TerminalException{1};
    }
    l_perm.push_back(*opt);
    MyMatrix<Tint> g = erdahl_lift(h, cfg.ab, cfg.ab);
    if (!erdahl_preserves_all(supers, g)) {
      all_ok = false;
    }
  }
  Telt id(static_cast<Tidx>(n_ext));
  Tgroup grp_all(l_perm, id);
  if (all_ok) {
    return {cfg, grp_all};
  }
  // The scalar products only give a supergroup: filter the elements.
#ifdef DEBUG_ERDAHL_GROUP
  os << "ERDAHL: filtering a group of order " << grp_all.size()
     << " for the preservation of the super polyhedra\n";
#endif
  std::vector<Telt> l_sel;
  for (auto &elt : grp_all) {
    MyMatrix<Tint> g = erdahl_lift_permutation<T, Tint, Tgroup>(cfg, elt);
    if (erdahl_preserves_all(supers, g)) {
      l_sel.push_back(elt);
    }
  }
  Tgroup grp(l_sel, id);
  return {cfg, grp};
}

/*
  A transformation g preserving W and the super polyhedra, mapping X1
  onto X2, if one exists.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
erdahl_equivalence(ErdahlFunctionSpace<T> const &W,
                   DelaunayPolyhedron<T, Tint> const &X1,
                   ErdahlChainConfig<T, Tint> const &cfg1,
                   DelaunayPolyhedron<T, Tint> const &X2,
                   ErdahlChainConfig<T, Tint> const &cfg2,
                   std::vector<ErdahlSuperData<T, Tint>> const &l_sd,
                   std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  if (X1.L.rows() != X2.L.rows() || X1.EXT.rows() != X2.EXT.rows()) {
    return {};
  }
  int n_ext = cfg1.EXTred.rows();
  std::vector<T> Vdiag(n_ext, T(0));
  std::optional<MyMatrix<Tint>> opt =
      TestIntEquivalence_ListMat_Vdiag<T, Tint, Tgroup>(
          cfg1.EXTred, cfg1.ListMat, Vdiag, cfg2.EXTred, cfg2.ListMat, Vdiag,
          os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> h = Inverse(*opt);
  if (!erdahl_permutation_of_rows<Tint, Telt>(cfg1.EXTred, cfg2.EXTred, h)) {
    std::cerr << "ERDAHL: the reduced equivalence does not map the "
                 "configurations\n";
    throw TerminalException{1};
  }
  auto finalize = [&](MyMatrix<Tint> const &g) -> MyMatrix<Tint> {
    MyMatrix<Tint> gc = erdahl_correct_for_space(W, g, X2);
#ifdef SANITY_CHECK_ERDAHL_GROUP
    DelaunayPolyhedron<T, Tint> X1img = erdahl_apply_transformation(X1, gc);
    if (X1img.EXT != X2.EXT || X1img.L != X2.L) {
      std::cerr << "ERDAHL: the equivalence does not map X1 to X2\n";
      throw TerminalException{1};
    }
    if (!erdahl_preserves_all(l_sd, gc)) {
      std::cerr << "ERDAHL: the equivalence does not preserve the supers\n";
      throw TerminalException{1};
    }
#endif
    return gc;
  };
  MyMatrix<Tint> g = erdahl_lift(h, cfg1.ab, cfg2.ab);
  if (erdahl_preserves_all(l_sd, g)) {
    return finalize(g);
  }
  // The scalar products are preserved but not the super polyhedra: search
  // the other equivalences h k with k in the automorphisms of X2.
  ErdahlFiniteGroup<T, Tint, Tgroup> fg2 =
      erdahl_finite_group<T, Tint, Tgroup>(X2, {}, os);
  MyMatrix<T> EXTred2_T = UniversalMatrixConversion<T, Tint>(cfg2.EXTred);
  for (auto &elt : fg2.grp) {
    MyMatrix<T> k_T = FindTransformation(EXTred2_T, EXTred2_T, elt);
    MyMatrix<Tint> k = UniversalMatrixConversion<Tint, T>(k_T);
    MyMatrix<Tint> gk = erdahl_lift(MyMatrix<Tint>(h * k), cfg1.ab, cfg2.ab);
    if (erdahl_preserves_all(l_sd, gk)) {
      return finalize(gk);
    }
  }
  return {};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
erdahl_equivalence(ErdahlFunctionSpace<T> const &W,
                   DelaunayPolyhedron<T, Tint> const &X1,
                   DelaunayPolyhedron<T, Tint> const &X2,
                   std::vector<DelaunayPolyhedron<T, Tint>> const &supers,
                   std::ostream &os) {
  if (X1.L.rows() != X2.L.rows() || X1.EXT.rows() != X2.EXT.rows()) {
    return {};
  }
  std::vector<ErdahlSuperData<T, Tint>> l_sd = erdahl_list_super_data(supers);
  ErdahlChainConfig<T, Tint> cfg1 = erdahl_chain_config(X1, l_sd);
  ErdahlChainConfig<T, Tint> cfg2 = erdahl_chain_config(X2, l_sd);
  return erdahl_equivalence<T, Tint, Tgroup>(W, X1, cfg1, X2, cfg2, l_sd, os);
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_GROUP_H_
// clang-format on
