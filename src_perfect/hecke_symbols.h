// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_SYMBOLS_H_
#define SRC_PERFECT_HECKE_SYMBOLS_H_

// clang-format off
#include "hecke_common.h"
#include "ImagQuadRing.h"
// clang-format on

/*
  The symbols of the sharbly method: the frames of n vectors of Z^n, or of
  m lines over O_K for a Hermitian T-space (SymbolGeometry), their level
  in the complex and their orientation.
 */

// One vector per pair {v, -v} of the rows, each with its first non-zero
// entry positive. The EXT of a face holds both v and -v.
template <typename Tint>
MyMatrix<Tint> half_vector_set(MyMatrix<Tint> const &EXT) {
  std::set<MyVector<Tint>> set;
  for (int i_row = 0; i_row < EXT.rows(); i_row++) {
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    set.insert(CanonicalizeVectorToInvertible(V));
  }
  MyMatrix<Tint> M(set.size(), EXT.cols());
  int pos = 0;
  for (auto &V : set) {
    AssignMatrixRow(M, pos, V);
    pos += 1;
  }
  return M;
}

/*
  The symbols of a Hermitian T-space over an imaginary quadratic field
  K = Q(omega), omega^2 = s omega - p. The vectors of Z^{2m} are the
  elements of O_K^m with the layout (x_1, ..., x_m, y_1, ..., y_m) for
  z_j = x_j + y_j omega, and the multiplication by omega is the matrix J
  of ListComm acting on the right (GetCommQuadratic). A symbol is a family
  of m lines over K: its cell has the dimension m, so the symbol level is
  dim T - m, and the EXT of the cell holds all the unit multiples of the
  generators of the lines.

  The Ash-Rudolph reduction works the same over a norm-Euclidean ring
  O_K: for a symbol [v_1, ..., v_m] with |N(det)| > 1 and a lattice point
  w = sum_i q_i v_i with all N(q_i) < 1,
    [v_1, ..., v_m] = sum_{q_i <> 0} [v_1, .., v_{i-1}, w, v_{i+1}, .., v_m]
  and each term has the norm of determinant N(q_i) |N(det)|. The lattice
  points modulo the O_K-span of the v_i are the classes of Z^{2m} modulo
  the Z-span of the rows of the realification R(W) = (W ; W J), and the
  coordinate q_i of the class f is (f_i, f_{m+i}) modulo O_K.
 */
template <typename T> struct ImagQuadStructure {
  int m;
  T eSum;
  T eProd;
  MyMatrix<T> J;
};

template <typename T>
std::optional<ImagQuadStructure<T>> get_imag_quad_structure(LinSpaceMatrix<T> const &LinSpa) {
  if (LinSpa.ListComm.size() != 1) {
    return {};
  }
  int n = LinSpa.n;
  if (n % 2 != 0) {
    return {};
  }
  int m = n / 2;
  MyMatrix<T> const &J = LinSpa.ListComm[0];
  T eSum = J(m, m);
  T eProd = -J(m, 0);
  T disc = eSum * eSum - 4 * eProd;
  if (disc >= 0) {
    return {};
  }
  if (J != GetCommQuadratic(m, eSum, eProd)) {
    return {};
  }
  return ImagQuadStructure<T>{m, eSum, eProd, J};
}

/*
  The geometry of the symbols: the classic T-space (the symbols are n
  vectors of Z^n) or a Hermitian T-space (m lines over O_K). All the
  matrices of vectors are of rational integers, a line of O_K^m being
  given by one of its generators.
 */
template <typename T, typename Tint> struct SymbolGeometry {
  std::optional<ImagQuadStructure<T>> field;
  // The arithmetic of O_K (Euclidean division, gcd, units).
  std::optional<ImagQuadRing<T>> ring;
  // The number of rows of a symbol: n, or m.
  int dim;
  int n;
  // The units of O_K, as pairs (a, b) for a + b omega.
  std::vector<std::pair<T, T>> l_unit;
  using Elt = std::pair<T, T>;
  explicit SymbolGeometry(LinSpaceMatrix<T> const &LinSpa)
      : field(get_imag_quad_structure(LinSpa)), n(LinSpa.n) {
    if (field) {
      dim = field->m;
      ring = ImagQuadRing<T>(field->eSum, field->eProd);
      l_unit = ring->units();
    } else {
      dim = n;
    }
  }
  bool is_hermitian() const { return field.has_value(); }
  Elt get_entry(MyVector<T> const &V, int const &j) const { return {V(j), V(field->m + j)}; }
  void set_entry(MyVector<T> &V, int const &j, Elt const &x) const {
    V(j) = x.first;
    V(field->m + j) = x.second;
  }
  MyVector<T> scale(Elt const &u, MyVector<T> const &V) const {
    MyVector<T> W(n);
    for (int j = 0; j < field->m; j++) {
      set_entry(W, j, ring->mul(u, get_entry(V, j)));
    }
    return W;
  }
  // The generator of the line of V in O_K^m (or Z^n) that is the
  // largest in the lexicographic order among its unit multiples.
  MyVector<Tint> primitive_vector(MyVector<T> const &V) const {
    MyVector<T> Vred = RemoveFractionVector(V);
    if (!field) {
      return UniversalVectorConversion<Tint, T>(Vred);
    }
    std::vector<Elt> l_entry;
    for (int j = 0; j < field->m; j++) {
      l_entry.push_back(get_entry(Vred, j));
    }
    std::optional<Elt> opt_g = ring->gcd(l_entry);
    if (!opt_g) {
      throw HeckeFailure{"sharbly: a Euclidean division failed, O_K is not norm-Euclidean"};
    }
    MyVector<T> Vprim(n);
    for (int j = 0; j < field->m; j++) {
      set_entry(Vprim, j, ring->div(l_entry[j], *opt_g));
    }
    std::optional<MyVector<T>> best;
    for (auto &u : l_unit) {
      MyVector<T> Vu = scale(u, Vprim);
      if (!best || *best < Vu) {
        best = Vu;
      }
    }
    return UniversalVectorConversion<Tint, T>(*best);
  }
  MyMatrix<Tint> primitive_rows(MyMatrix<T> const &W_T) const {
    MyMatrix<Tint> W(W_T.rows(), W_T.cols());
    for (int i_row = 0; i_row < W_T.rows(); i_row++) {
      AssignMatrixRow(W, i_row, primitive_vector(GetMatrixRow(W_T, i_row)));
    }
    return W;
  }
  // The vectors of a cell, one per line (one per pair {v, -v} for the
  // classic T-space), in a canonical order.
  MyMatrix<Tint> symbol_rows(MyMatrix<Tint> const &EXT) const {
    if (!field) {
      return half_vector_set(EXT);
    }
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    std::set<MyVector<Tint>> set;
    for (int i_row = 0; i_row < EXT.rows(); i_row++) {
      set.insert(primitive_vector(GetMatrixRow(EXT_T, i_row)));
    }
    MyMatrix<Tint> M(set.size(), EXT.cols());
    int pos = 0;
    for (auto &V : set) {
      AssignMatrixRow(M, pos, V);
      pos += 1;
    }
    return M;
  }
  // The rows of W and of W J.
  MyMatrix<T> realify(MyMatrix<T> const &W) const {
    if (!field) {
      return W;
    }
    MyMatrix<T> WJ = W * field->J;
    MyMatrix<T> R(2 * W.rows(), W.cols());
    for (int i = 0; i < W.rows(); i++) {
      for (int j = 0; j < W.cols(); j++) {
        R(i, j) = W(i, j);
        R(W.rows() + i, j) = WJ(i, j);
      }
    }
    return R;
  }
  MyMatrix<Tint> realify_int(MyMatrix<Tint> const &W) const {
    return UniversalMatrixConversion<Tint, T>(realify(UniversalMatrixConversion<T, Tint>(W)));
  }
  // |det W| for the classic T-space, |N(det W)| for the Hermitian one.
  Tint det_norm(MyMatrix<Tint> const &W) const {
    return T_abs(DeterminantMat(realify_int(W)));
  }
  // Whether the rows of D span the whole space (over K).
  bool is_spanning(MyMatrix<Tint> const &D) const {
    return RankMat(realify_int(D)) == n;
  }
};

// The symbol level: the level of the cones of dimension n (m lines for a
// Hermitian T-space), which is dim T - n since level k holds the cones of
// codimension k.
template <typename T, typename Tint, typename Tgroup>
int get_symbol_level(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                     SymbolGeometry<T, Tint> const &geom) {
  int dimT = fce.pctdi.LinSpa.ListMat.size();
  return dimT - geom.dim;
}

/*
  The orientation of a simplicial cell of the symbol level, given with its
  rays in the order of the rows of W0, relative to the reference
  orientation of the cell (or_info): the cell is symbol_orientation times
  the oriented simplex [W0]. The rank one forms of v and u v (u a unit)
  are the same, so any generator of the lines can be used.
 */
template <typename T, typename Tint, typename Tgroup>
int symbol_orientation(FacePerfectComplex<T, Tint, Tgroup> const &face, MyMatrix<Tint> const &W0,
                       std::vector<MyMatrix<T>> const &ListMat) {
  OrientationInfo const &or_info = face.or_info;
  int dim = or_info.ListColSelect.size();
  if (W0.rows() != dim) {
    std::cerr << "HECKE: symbol_orientation, the symbol has " << W0.rows()
              << " rows but the cell is of dimension " << dim << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> M(dim, dim);
  for (int i = 0; i < dim; i++) {
    MyVector<Tint> V = GetMatrixRow(W0, i);
    for (int j = 0; j < dim; j++) {
      M(i, j) = EvaluationQuadForm<T, Tint>(ListMat[or_info.ListColSelect[j]], V);
    }
  }
  T det = DeterminantMat(M);
  if (det == 0) {
    std::cerr << "HECKE: symbol_orientation, the rays should be independent\n";
    throw TerminalException{1};
  }
  return (det > 0) ? or_info.sign : -or_info.sign;
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_SYMBOLS_H_
// clang-format on
