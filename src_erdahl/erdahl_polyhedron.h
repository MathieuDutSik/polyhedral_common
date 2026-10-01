// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_POLYHEDRON_H_
#define SRC_ERDAHL_ERDAHL_POLYHEDRON_H_

// clang-format off
#include "erdahl_result.h"
#include "erdahl_space.h"
#include "Positivity.h"
#include "Shvec_exact.h"
#include "POLY_LinearProgramming.h"
#include <cstdint>
#include <optional>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_POLYHEDRON
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ERDAHL_POLYHEDRON
#endif

#ifdef TIMINGS
#define TIMINGS_ERDAHL_POLYHEDRON
#endif

/*
  Zero sets of nonnegative degree 2 functions and Delaunay polyhedra.

  By Erdahl (1992), if f >= 0 on Z^n then its zero set Z(f) is either empty
  or of the form P + L with L a saturated sublattice (the isotropy lattice)
  and P a finite set of representatives modulo L. When P + L affinely
  spans R^n, it is a Delaunay polyhedron. It is stored as
  * EXT: the representatives modulo L, as rows (1, x),
  * L: a Z-basis of L, as rows,
  * F: a function f in Erdahl(n) (and in the space W of the computation)
    with Z(f) = P + L.
  See erdahl_space.h for the conventions on functions and transformations.
 */

// Points (rows (1, x)) where a function is strictly negative. When the
// function is unbounded below, also the directions (rows (0, v)) along which
// it is: v with Quad(f)[v] < 0, or v in the kernel of Quad(f) on which the
// linear part does not vanish.
template <typename Tint> struct ErdahlNegativeWitness {
  std::string reason;
  std::vector<MyVector<Tint>> points;
  std::vector<MyVector<Tint>> directions;
};

// The set EXT + L: EXT the representatives as rows (1, x), L a basis.
template <typename Tint> struct ErdahlLatticeSet {
  MyMatrix<Tint> EXT;
  MyMatrix<Tint> L;
};

template <typename T, typename Tint> struct ErdahlLatticeMinimum {
  T min_value;
  ErdahlLatticeSet<Tint> argmin;
};

template <typename Tint> MyVector<Tint> erdahl_direction_row(MyVector<Tint> const &z) {
  int n = z.size();
  MyVector<Tint> e(n + 1);
  e(0) = 0;
  for (int i = 0; i < n; i++) {
    e(i + 1) = z(i);
  }
  return e;
}

/*
  Along the line t -> e0 + t (0, v) the function is a + 2 b t + q t^2. When
  q < 0, or q = 0 and b != 0, it takes negative values at integers t, and
  one such lattice point is returned.
 */
template <typename T, typename Tint>
MyVector<Tint> erdahl_negative_point_on_line(MyMatrix<T> const &F,
                                             MyVector<Tint> const &e0,
                                             MyVector<Tint> const &v) {
  MyVector<Tint> dir = erdahl_direction_row(v);
  T a = erdahl_bilinear(F, e0, e0);
  T b = erdahl_bilinear(F, e0, dir);
  T q = erdahl_bilinear(F, dir, dir);
  if (q > 0 || (q == 0 && b == 0)) {
    std::cerr << "ERDAHL: the function is bounded below on the line\n";
    throw TerminalException{1};
  }
  Tint t(1);
  while (true) {
    for (int sign = -1; sign <= 1; sign += 2) {
      T t_T = UniversalScalarConversion<T, Tint>(t) * sign;
      T val = a + 2 * b * t_T + q * t_T * t_T;
      if (val < 0) {
        Tint t_s = t * sign;
        MyVector<Tint> e = e0 + t_s * dir;
        return e;
      }
    }
    t *= 2;
  }
}

/*
  The minimum of the degree 2 function F (of size m+1) over Z^m.
  * If F is unbounded below, an error with points of negative value.
  * Otherwise the minimum value and the set of minimizers, which is of the
    form EXT + L with L the integral kernel of Quad(F).
 */
template <typename T, typename Tint>
ErdahlResult<ErdahlLatticeMinimum<T, Tint>, ErdahlNegativeWitness<Tint>>
erdahl_lattice_minimum(MyMatrix<T> const &F, std::ostream &os) {
  using Tres =
      ErdahlResult<ErdahlLatticeMinimum<T, Tint>, ErdahlNegativeWitness<Tint>>;
  int m = F.rows() - 1;
  T cst = F(0, 0);
  MyVector<Tint> origin = ZeroVector<Tint>(m + 1);
  origin(0) = 1;
  if (m == 0) {
    MyMatrix<Tint> EXT(1, 1);
    EXT(0, 0) = 1;
    MyMatrix<Tint> L(0, 0);
    ErdahlLatticeSet<Tint> argmin{EXT, L};
    return Tres::ok({cst, argmin});
  }
  MyMatrix<T> Q = erdahl_get_quad(F);
  MyVector<T> lin = erdahl_get_lin(F);
  if (!IsPositiveSemiDefinite(Q, os)) {
    MyVector<Tint> v = GetShortVector<T, Tint>(Q, T(0), os);
    MyVector<Tint> e = erdahl_negative_point_on_line(F, origin, v);
    return Tres::err({"the quadratic form is not positive semidefinite",
                      {e},
                      {erdahl_direction_row(v)}});
  }
  MyMatrix<T> Q_red = RemoveFractionMatrix(Q);
  MyMatrix<Tint> Q_int = UniversalMatrixConversion<Tint, T>(Q_red);
  MyMatrix<Tint> NSP = NullspaceIntMat(Q_int);
  int dim_nsp = NSP.rows();
  for (int i = 0; i < dim_nsp; i++) {
    MyVector<Tint> z = GetMatrixRow(NSP, i);
    MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
    if (lin.dot(z_T) != 0) {
      MyVector<Tint> e = erdahl_negative_point_on_line(F, origin, z);
      return Tres::err({"the linear part does not vanish on the kernel of the "
                        "quadratic form",
                        {e},
                        {erdahl_direction_row(z)}});
    }
  }
  int p = m - dim_nsp;
  if (p == 0) {
    // Q = 0 and lin = 0: the function is the constant cst.
    MyMatrix<Tint> EXT(1, m + 1);
    AssignMatrixRow(EXT, 0, origin);
    ErdahlLatticeSet<Tint> argmin{EXT, NSP};
    return Tres::ok({cst, argmin});
  }
  MyMatrix<Tint> C = SubspaceCompletionInt(NSP, m);
  MyMatrix<T> C_T = UniversalMatrixConversion<T, Tint>(C);
  MyMatrix<T> Q_c = C_T * Q * C_T.transpose();
  MyVector<T> lin_c = C_T * lin;
  // f(y C) = cst + 2 lin_c.y + Q_c[y] = cst - Q_c[y0] + Q_c[y - y0]
  MyMatrix<T> Q_c_inv = Inverse(Q_c);
  MyVector<T> y0 = -(Q_c_inv * lin_c);
  T min_cont = cst + lin_c.dot(y0);
  resultCVP<T, Tint> cvp = NearestVectors<T, Tint>(Q_c, y0, os);
  T min_value = min_cont + cvp.TheNorm;
  MyMatrix<Tint> Xred = cvp.ListVect * C;
  int n_ext = Xred.rows();
  MyMatrix<Tint> EXT(n_ext, m + 1);
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    EXT(i_ext, 0) = 1;
    for (int i = 0; i < m; i++) {
      EXT(i_ext, i + 1) = Xred(i_ext, i);
    }
  }
#ifdef SANITY_CHECK_ERDAHL_POLYHEDRON
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    MyVector<Tint> e = GetMatrixRow(EXT, i_ext);
    T val = EvaluationQuadForm<T, Tint>(F, e);
    if (val != min_value) {
      std::cerr << "ERDAHL: the minimizer does not attain the minimum\n";
      throw TerminalException{1};
    }
  }
#endif
  ErdahlLatticeSet<Tint> argmin{EXT, NSP};
  return Tres::ok({min_value, argmin});
}

/*
  The zero set of F, when F is nonnegative on Z^n. Otherwise an error with
  lattice points where F is negative. An empty zero set is a success with
  an EXT of zero rows.
 */
template <typename T, typename Tint>
ErdahlResult<ErdahlLatticeSet<Tint>, ErdahlNegativeWitness<Tint>>
erdahl_zero_set(MyMatrix<T> const &F, std::ostream &os) {
  using Tres = ErdahlResult<ErdahlLatticeSet<Tint>, ErdahlNegativeWitness<Tint>>;
  auto res = erdahl_lattice_minimum<T, Tint>(F, os);
  if (res.is_err()) {
    return Tres::err(res.get_err());
  }
  ErdahlLatticeMinimum<T, Tint> const &min = res.get_ok();
  if (min.min_value < 0) {
    std::vector<MyVector<Tint>> points;
    for (int i = 0; i < min.argmin.EXT.rows(); i++) {
      points.push_back(GetMatrixRow(min.argmin.EXT, i));
    }
    return Tres::err({"the minimum over the lattice is negative", points, {}});
  }
  if (min.min_value > 0) {
    int np1 = F.rows();
    MyMatrix<Tint> EXT(0, np1);
    MyMatrix<Tint> L(0, np1 - 1);
    return Tres::ok({EXT, L});
  }
  return Tres::ok(min.argmin);
}

/*
  The coordinates adapted to a saturated lattice L of rank d in Z^n:
  a complement C (p = n - d rows) such that the rows of C and L form a basis
  of Z^n, and the affine basis
      AffBasis = [ 1 0 ]
                 [ 0 C ]
                 [ 0 L ]
  A row e is written as u AffBasis with u = (1, y, z), y of length p.
 */
template <typename Tint> struct ErdahlAdaptedBasis {
  int n;
  int p;
  int d;
  MyMatrix<Tint> C;
  MyMatrix<Tint> AffBasis;
  MyMatrix<Tint> AffBasisInv;
};

template <typename Tint>
ErdahlAdaptedBasis<Tint> erdahl_adapted_basis(MyMatrix<Tint> const &L, int n) {
  int d = L.rows();
  int p = n - d;
  MyMatrix<Tint> C = SubspaceCompletionInt(L, n);
  MyMatrix<Tint> AffBasis = ZeroMatrix<Tint>(n + 1, n + 1);
  AffBasis(0, 0) = 1;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < n; j++) {
      AffBasis(1 + i, 1 + j) = C(i, j);
    }
  }
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < n; j++) {
      AffBasis(1 + p + i, 1 + j) = L(i, j);
    }
  }
  MyMatrix<Tint> AffBasisInv = Inverse(AffBasis);
#ifdef SANITY_CHECK_ERDAHL_POLYHEDRON
  if (!IsIntegralMatrix(AffBasisInv)) {
    std::cerr << "ERDAHL: the adapted basis should be unimodular\n";
    throw TerminalException{1};
  }
#endif
  return {n, p, d, C, AffBasis, AffBasisInv};
}

// The representative of e + L with vanishing L-coordinates.
template <typename Tint>
MyVector<Tint> erdahl_canonical_representative(ErdahlAdaptedBasis<Tint> const &ab,
                                               MyVector<Tint> const &e) {
  MyVector<Tint> u = ab.AffBasisInv.transpose() * e;
  for (int i = 0; i < ab.d; i++) {
    u(1 + ab.p + i) = 0;
  }
  return ab.AffBasis.transpose() * u;
}

// The first p+1 coordinates of e in the adapted basis: (1, y).
template <typename Tint>
MyVector<Tint> erdahl_reduced_coordinates(ErdahlAdaptedBasis<Tint> const &ab,
                                          MyVector<Tint> const &e) {
  MyVector<Tint> u = ab.AffBasisInv.transpose() * e;
  MyVector<Tint> y(ab.p + 1);
  for (int i = 0; i <= ab.p; i++) {
    y(i) = u(i);
  }
  return y;
}

// The lattice spanned by the rows of L in Hermite normal form.
template <typename Tint> MyMatrix<Tint> erdahl_canonical_lattice(MyMatrix<Tint> const &L) {
  int n = L.cols();
  if (L.rows() == 0) {
    return MyMatrix<Tint>(0, n);
  }
  MyMatrix<Tint> H = ComputeRowHermiteNormalForm(L).second;
  std::vector<MyVector<Tint>> l_row;
  for (int i = 0; i < H.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(H, i);
    if (!IsZeroVector(V)) {
      l_row.push_back(V);
    }
  }
  return MatrixFromVectorFamilyDim(n, l_row);
}

template <typename Tint>
MyMatrix<Tint> erdahl_sorted_rows(std::vector<MyVector<Tint>> l_row, int len) {
  std::sort(l_row.begin(), l_row.end());
  l_row.erase(std::unique(l_row.begin(), l_row.end()), l_row.end());
  return MatrixFromVectorFamilyDim(len, l_row);
}

// The canonical form of EXT + L: HNF lattice, canonical sorted representatives.
template <typename Tint>
ErdahlLatticeSet<Tint> erdahl_canonical_lattice_set(ErdahlLatticeSet<Tint> const &ls) {
  int n = ls.EXT.cols() - 1;
  MyMatrix<Tint> L = erdahl_canonical_lattice(ls.L);
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(L, n);
  std::vector<MyVector<Tint>> l_row;
  for (int i = 0; i < ls.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i);
    l_row.push_back(erdahl_canonical_representative(ab, e));
  }
  MyMatrix<Tint> EXT = erdahl_sorted_rows(l_row, n + 1);
  return {EXT, L};
}

// Whether EXT + L affinely spans R^n.
template <typename Tint>
bool erdahl_is_full_dimensional(ErdahlLatticeSet<Tint> const &ls) {
  int n = ls.EXT.cols() - 1;
  if (ls.EXT.rows() == 0) {
    return false;
  }
  MyMatrix<Tint> M(ls.EXT.rows() + ls.L.rows(), n + 1);
  for (int i = 0; i < ls.EXT.rows(); i++) {
    AssignMatrixRow(M, i, GetMatrixRow(ls.EXT, i));
  }
  for (int i = 0; i < ls.L.rows(); i++) {
    AssignMatrixRow(M, ls.EXT.rows() + i,
                    erdahl_direction_row(GetMatrixRow(ls.L, i)));
  }
  return RankMat(M) == n + 1;
}

template <typename T, typename Tint> struct DelaunayPolyhedron {
  // Canonical representatives modulo L, rows (1, x), sorted.
  MyMatrix<Tint> EXT;
  // The isotropy lattice L(D) in Hermite normal form.
  MyMatrix<Tint> L;
  // A function of Erdahl(n) with Z(F) = EXT + L.
  MyMatrix<T> F;
};

template <typename T, typename Tint>
int erdahl_dimension(DelaunayPolyhedron<T, Tint> const &D) {
  return D.EXT.cols() - 1;
}

template <typename T, typename Tint>
ErdahlAdaptedBasis<Tint> erdahl_adapted_basis(DelaunayPolyhedron<T, Tint> const &D) {
  return erdahl_adapted_basis(D.L, erdahl_dimension(D));
}

// The representatives in the adapted coordinates, rows (1, y).
template <typename T, typename Tint>
MyMatrix<Tint> erdahl_reduced_ext(DelaunayPolyhedron<T, Tint> const &D,
                                  ErdahlAdaptedBasis<Tint> const &ab) {
  MyMatrix<Tint> EXTred(D.EXT.rows(), ab.p + 1);
  for (int i = 0; i < D.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(D.EXT, i);
    AssignMatrixRow(EXTred, i, erdahl_reduced_coordinates(ab, e));
  }
  return EXTred;
}

/*
  The Delaunay polyhedron of a function F of Erdahl(n). The zero set of F
  has to be nonempty and full dimensional; the computations only produce
  such functions, so anything else is a programming or input error.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_polyhedron_from_function(MyMatrix<T> const &F, std::ostream &os) {
  auto res = erdahl_zero_set<T, Tint>(F, os);
  if (res.is_err()) {
    std::cerr << "ERDAHL: the function is negative on some lattice points, "
                 "reason="
              << res.get_err().reason << "\n";
    throw TerminalException{1};
  }
  ErdahlLatticeSet<Tint> const &zs = res.get_ok();
  if (!erdahl_is_full_dimensional(zs)) {
    std::cerr << "ERDAHL: the zero set is empty or not full dimensional, "
                 "|EXT|="
              << zs.EXT.rows() << " dim(L)=" << zs.L.rows() << "\n";
    throw TerminalException{1};
  }
  ErdahlLatticeSet<Tint> can = erdahl_canonical_lattice_set(zs);
  MyMatrix<T> Fred = RemoveFractionMatrix(F);
  return {can.EXT, can.L, Fred};
}

template <typename T, typename Tint>
bool erdahl_contains_point(DelaunayPolyhedron<T, Tint> const &D,
                           ErdahlAdaptedBasis<Tint> const &ab,
                           MyVector<Tint> const &e) {
  MyVector<Tint> e_can = erdahl_canonical_representative(ab, e);
  for (int i = 0; i < D.EXT.rows(); i++) {
    if (GetMatrixRow(D.EXT, i) == e_can) {
      return true;
    }
  }
  return false;
}

template <typename T, typename Tint>
bool erdahl_contains_point(DelaunayPolyhedron<T, Tint> const &D,
                           MyVector<Tint> const &e) {
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(D);
  return erdahl_contains_point(D, ab, e);
}

// Whether the lattice spanned by L1 is contained in the one spanned by L2.
template <typename Tint>
bool erdahl_is_sublattice(MyMatrix<Tint> const &L1, MyMatrix<Tint> const &L2) {
  for (int i = 0; i < L1.rows(); i++) {
    MyVector<Tint> V = GetMatrixRow(L1, i);
    if (L2.rows() == 0) {
      return false;
    }
    if (!SolutionIntMat(L2, V)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename Tint>
bool erdahl_is_subset(DelaunayPolyhedron<T, Tint> const &D1,
                      DelaunayPolyhedron<T, Tint> const &D2) {
  if (!erdahl_is_sublattice(D1.L, D2.L)) {
    return false;
  }
  ErdahlAdaptedBasis<Tint> ab2 = erdahl_adapted_basis(D2);
  for (int i = 0; i < D1.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(D1.EXT, i);
    if (!erdahl_contains_point(D2, ab2, e)) {
      return false;
    }
  }
  return true;
}

// The polyhedra are in canonical form, so equality is equality of data.
template <typename T, typename Tint>
bool erdahl_is_equal(DelaunayPolyhedron<T, Tint> const &D1,
                     DelaunayPolyhedron<T, Tint> const &D2) {
  bool test = D1.L == D2.L && D1.EXT == D2.EXT;
#ifdef SANITY_CHECK_ERDAHL_POLYHEDRON
  bool test2 = erdahl_is_subset(D1, D2) && erdahl_is_subset(D2, D1);
  if (test != test2) {
    std::cerr << "ERDAHL: the canonical forms are inconsistent\n";
    throw TerminalException{1};
  }
#endif
  return test;
}

/*
  The linear conditions for a function to vanish on EXT + L, expressed on
  the coefficients in the basis of W: one row per condition.
  * f(x) = 0 for x in EXT,
  * e F (0, z)^T = 0 for e in EXT and z in the basis of L,
  * (0, z1) F (0, z2)^T = 0 for z1, z2 in the basis of L.
  Together they say that t -> f(x + t z) vanishes identically.
 */
template <typename T, typename Tint>
MyMatrix<T> erdahl_vanishing_equations(ErdahlFunctionSpace<T> const &W,
                                       ErdahlLatticeSet<Tint> const &ls) {
  std::vector<MyVector<T>> l_equa;
  int n_ext = ls.EXT.rows();
  int d = ls.L.rows();
  std::vector<MyVector<Tint>> l_dir;
  for (int i = 0; i < d; i++) {
    l_dir.push_back(erdahl_direction_row(GetMatrixRow(ls.L, i)));
  }
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i_ext);
    l_equa.push_back(erdahl_bilinear_vector(W, e, e));
    for (auto &dir : l_dir) {
      l_equa.push_back(erdahl_bilinear_vector(W, e, dir));
    }
  }
  for (int i = 0; i < d; i++) {
    for (int j = 0; j <= i; j++) {
      l_equa.push_back(erdahl_bilinear_vector(W, l_dir[i], l_dir[j]));
    }
  }
  return MatrixFromVectorFamilyDim(W.basis.size(), l_equa);
}

/*
  Space_W(D): the functions of W vanishing on D. Returned as coefficient
  vectors in the basis of W (rows).
 */
template <typename T, typename Tint>
MyMatrix<T> erdahl_vanishing_space_coeff(ErdahlFunctionSpace<T> const &W,
                                         ErdahlLatticeSet<Tint> const &ls) {
  MyMatrix<T> Equa = erdahl_vanishing_equations(W, ls);
  if (Equa.rows() == 0) {
    return IdentityMat<T>(W.basis.size());
  }
  return NullspaceTrMat(Equa);
}

template <typename T, typename Tint>
std::vector<MyMatrix<T>>
erdahl_vanishing_space(ErdahlFunctionSpace<T> const &W,
                       ErdahlLatticeSet<Tint> const &ls) {
  MyMatrix<T> NSP = erdahl_vanishing_space_coeff(W, ls);
  std::vector<MyMatrix<T>> l_func;
  for (int i = 0; i < NSP.rows(); i++) {
    MyVector<T> coeff = GetMatrixRow(NSP, i);
    l_func.push_back(
        RemoveFractionMatrix(erdahl_function_from_coefficients(W, coeff)));
  }
  return l_func;
}

template <typename T, typename Tint>
ErdahlLatticeSet<Tint> erdahl_lattice_set(DelaunayPolyhedron<T, Tint> const &D) {
  return {D.EXT, D.L};
}

// The perfection rank of D relative to W: dim Space_W(D).
template <typename T, typename Tint>
int erdahl_perfection_rank(ErdahlFunctionSpace<T> const &W,
                           DelaunayPolyhedron<T, Tint> const &D) {
  return erdahl_vanishing_space_coeff(W, erdahl_lattice_set(D)).rows();
}

/*
  The restriction of G to the coset e + L, as a function on Z^d:
  z -> G(e + z L).
 */
template <typename T, typename Tint>
MyMatrix<T> erdahl_restrict_to_coset(MyMatrix<T> const &G,
                                     MyVector<Tint> const &e,
                                     MyMatrix<Tint> const &L) {
  int d = L.rows();
  std::vector<MyVector<Tint>> l_row{e};
  for (int i = 0; i < d; i++) {
    l_row.push_back(erdahl_direction_row(GetMatrixRow(L, i)));
  }
  MyMatrix<T> Fd(d + 1, d + 1);
  for (int i = 0; i <= d; i++) {
    for (int j = 0; j <= d; j++) {
      Fd(i, j) = erdahl_bilinear(G, l_row[i], l_row[j]);
    }
  }
  return Fd;
}

// The point e + z L from the row (1, z).
template <typename Tint>
MyVector<Tint> erdahl_point_in_coset(MyVector<Tint> const &e,
                                     MyMatrix<Tint> const &L,
                                     MyVector<Tint> const &ez) {
  MyVector<Tint> ret = e;
  int n = e.size() - 1;
  for (int i = 0; i < L.rows(); i++) {
    if (ez(1 + i) != 0) {
      for (int j = 0; j < n; j++) {
        ret(1 + j) += ez(1 + i) * L(i, j);
      }
    }
  }
  return ret;
}

/*
  Points of the lattice set EXT + L where G is negative. The returned list
  is empty if and only if G >= 0 on EXT + L. For each representative, the
  restriction of G to the coset is minimized by a closest vector problem.
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>>
erdahl_negative_points_on_set(MyMatrix<T> const &G,
                              ErdahlLatticeSet<Tint> const &ls,
                              std::ostream &os) {
  std::vector<MyVector<Tint>> l_neg;
  for (int i_ext = 0; i_ext < ls.EXT.rows(); i_ext++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i_ext);
    if (ls.L.rows() == 0) {
      if (EvaluationQuadForm<T, Tint>(G, e) < 0) {
        l_neg.push_back(e);
      }
      continue;
    }
    MyMatrix<T> Fd = erdahl_restrict_to_coset(G, e, ls.L);
    auto res = erdahl_lattice_minimum<T, Tint>(Fd, os);
    if (res.is_err()) {
      for (auto &ez : res.get_err().points) {
        l_neg.push_back(erdahl_point_in_coset(e, ls.L, ez));
      }
      continue;
    }
    ErdahlLatticeMinimum<T, Tint> const &min = res.get_ok();
    if (min.min_value < 0) {
      for (int i = 0; i < min.argmin.EXT.rows(); i++) {
        MyVector<Tint> ez = GetMatrixRow(min.argmin.EXT, i);
        l_neg.push_back(erdahl_point_in_coset(e, ls.L, ez));
      }
    }
  }
#ifdef SANITY_CHECK_ERDAHL_POLYHEDRON
  for (auto &e : l_neg) {
    if (EvaluationQuadForm<T, Tint>(G, e) >= 0) {
      std::cerr << "ERDAHL: a returned point should have negative value\n";
      throw TerminalException{1};
    }
  }
#endif
  return l_neg;
}

/*
  The zero set of G on the lattice set EXT + L, when G >= 0 on it: on each
  coset e + L the restriction of G is minimized (a closest vector problem
  of dimension dim L, involving G only).
 */
template <typename T, typename Tint>
ErdahlLatticeSet<Tint> erdahl_zero_set_on_set(MyMatrix<T> const &G,
                                              ErdahlLatticeSet<Tint> const &ls,
                                              std::ostream &os) {
  int np1 = ls.EXT.cols();
  int n = np1 - 1;
  std::vector<MyVector<Tint>> l_pt;
  std::optional<MyMatrix<Tint>> opt_L;
  for (int i_ext = 0; i_ext < ls.EXT.rows(); i_ext++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i_ext);
    if (ls.L.rows() == 0) {
      T val = EvaluationQuadForm<T, Tint>(G, e);
      if (val < 0) {
        std::cerr << "ERDAHL: the function should be nonnegative on the set\n";
        throw TerminalException{1};
      }
      if (val == 0) {
        l_pt.push_back(e);
      }
      continue;
    }
    MyMatrix<T> Fd = erdahl_restrict_to_coset(G, e, ls.L);
    auto res = erdahl_lattice_minimum<T, Tint>(Fd, os);
    if (res.is_err() || res.get_ok().min_value < 0) {
      std::cerr << "ERDAHL: the function should be nonnegative on the set\n";
      throw TerminalException{1};
    }
    ErdahlLatticeMinimum<T, Tint> const &min = res.get_ok();
    if (min.min_value > 0) {
      continue;
    }
    for (int i = 0; i < min.argmin.EXT.rows(); i++) {
      MyVector<Tint> ez = GetMatrixRow(min.argmin.EXT, i);
      l_pt.push_back(erdahl_point_in_coset(e, ls.L, ez));
    }
    // The isotropy lattice: the kernel of the quadratic part of G on L,
    // the same for all the cosets.
    MyMatrix<Tint> Lzero = min.argmin.L * ls.L;
    if (!opt_L) {
      opt_L = Lzero;
    }
#ifdef SANITY_CHECK_ERDAHL_POLYHEDRON
    if (erdahl_canonical_lattice(*opt_L) != erdahl_canonical_lattice(Lzero)) {
      std::cerr << "ERDAHL: the isotropy lattices of the cosets differ\n";
      throw TerminalException{1};
    }
#endif
  }
  MyMatrix<Tint> EXT = MatrixFromVectorFamilyDim(np1, l_pt);
  MyMatrix<Tint> L = opt_L ? *opt_L : MyMatrix<Tint>(0, n);
  return {EXT, L};
}

/*
  The points x_i + x_j - x_k outside of EXT + L, for x_i, x_j, x_k among the
  representatives and their translates by the basis of L, as in the GAP
  function ExtendByTriangleIneq: the vertices of the adjacent Delaunay
  polytopes are typically among them. The points are given on machine
  integers; the empty optional means that the coordinates are too large.
 */
template <typename Tint>
std::optional<std::vector<std::vector<int64_t>>>
erdahl_triangle_pool(ErdahlLatticeSet<Tint> const &ls) {
  using Tpt = std::vector<int64_t>;
  int np1 = ls.EXT.cols();
  int n = np1 - 1;
  int64_t const coord_max = int64_t(1) << 20;
  Tint b = UniversalScalarConversion<Tint, int64_t>(coord_max);
  auto to_pt = [&](MyVector<Tint> const &e) -> std::optional<Tpt> {
    Tpt pt(np1);
    for (int a = 0; a < np1; a++) {
      if (!(e(a) < b && e(a) > -b)) {
        return {};
      }
      pt[a] = UniversalScalarConversion<int64_t, Tint>(e(a));
    }
    return pt;
  };
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(ls.L, n);
  std::vector<std::vector<int64_t>> Ainv(np1, std::vector<int64_t>(np1));
  std::vector<std::vector<int64_t>> Aff(np1, std::vector<int64_t>(np1));
  for (int a = 0; a < np1; a++) {
    for (int c = 0; c < np1; c++) {
      if (!(ab.AffBasisInv(a, c) < b && ab.AffBasisInv(a, c) > -b) ||
          !(ab.AffBasis(a, c) < b && ab.AffBasis(a, c) > -b)) {
        return {};
      }
      Ainv[a][c] = UniversalScalarConversion<int64_t, Tint>(ab.AffBasisInv(a, c));
      Aff[a][c] = UniversalScalarConversion<int64_t, Tint>(ab.AffBasis(a, c));
    }
  }
  std::unordered_set<Tpt> set_ext;
  std::vector<Tpt> l_vert;
  for (int i_ext = 0; i_ext < ls.EXT.rows(); i_ext++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i_ext);
    std::optional<Tpt> opt = to_pt(e);
    if (!opt) {
      return {};
    }
    set_ext.insert(*opt);
    l_vert.push_back(*opt);
    for (int j = 0; j < ls.L.rows(); j++) {
      Tpt f = *opt;
      for (int k = 0; k < n; k++) {
        f[1 + k] += UniversalScalarConversion<int64_t, Tint>(ls.L(j, k));
      }
      l_vert.push_back(f);
    }
  }
  int p_ab = ab.p;
  bool is_polytope = ls.L.rows() == 0;
  auto is_in_set = [&](Tpt const &pt) -> bool {
    if (is_polytope) {
      return set_ext.count(pt) > 0;
    }
    Tpt u(np1, 0);
    for (int c = 0; c <= p_ab; c++) {
      int64_t sum = 0;
      for (int a = 0; a < np1; a++) {
        sum += pt[a] * Ainv[a][c];
      }
      u[c] = sum;
    }
    Tpt img(np1, 0);
    for (int c = 0; c < np1; c++) {
      int64_t sum = 0;
      for (int a = 0; a <= p_ab; a++) {
        sum += u[a] * Aff[a][c];
      }
      img[c] = sum;
    }
    return set_ext.count(img) > 0;
  };
  // The differences first, deduplicated, then the translates, then the
  // membership test, once per candidate.
  std::unordered_set<Tpt> set_diff;
  for (auto &vj : l_vert) {
    for (auto &vk : l_vert) {
      if (vj != vk) {
        Tpt delta(np1);
        for (int a = 0; a < np1; a++) {
          delta[a] = vj[a] - vk[a];
        }
        set_diff.insert(delta);
      }
    }
  }
  std::unordered_set<Tpt> set_cand;
  for (auto &vi : l_vert) {
    for (auto &delta : set_diff) {
      Tpt f(np1);
      for (int a = 0; a < np1; a++) {
        f[a] = vi[a] + delta[a];
      }
      set_cand.insert(f);
    }
  }
  std::vector<Tpt> l_pool;
  for (auto &f : set_cand) {
    if (!is_in_set(f)) {
      l_pool.push_back(f);
    }
  }
  // A deterministic order.
  std::sort(l_pool.begin(), l_pool.end());
  return l_pool;
}

/*
  What the canonical function of a sub-polyhedron F of a Delaunay polyhedron
  D can reuse from D: the pool of D (which contains the pool of F, since
  the points of F are points of D, and lies outside of D hence of F) and
  the points of D, which are relevant constraints for F.
 */
template <typename Tint> struct ErdahlCanonicalHint {
  std::vector<MyVector<Tint>> extra_points;
  std::optional<std::vector<std::vector<int64_t>>> pool;
};

template <typename Tint>
ErdahlCanonicalHint<Tint> erdahl_canonical_hint(ErdahlLatticeSet<Tint> const &ls) {
  std::vector<MyVector<Tint>> extra_points;
  for (int i = 0; i < ls.EXT.rows(); i++) {
    extra_points.push_back(GetMatrixRow(ls.EXT, i));
  }
  return {extra_points, erdahl_triangle_pool(ls)};
}

/*
  A function of Erdahl(n) cap W whose zero set is exactly the lattice set
  EXT + L (which has to be a Delaunay polyhedron relative to W), with small
  coefficients. This is the analogue of the GAP function
  InfDel_GetCanonicalPosQuadFuncExpression and of Algorithm
  "TestRealizability" of the paper. The function is sought in Space_W(D),
  the functions of W vanishing on D, by the linear program
      minimize   sum_{v in S} f(v)
      subject to f(v) >= 1 for v in S
  for a finite set S of lattice points outside of D, and of directions
  (0, v) for which f((0, v)) = Quad(f)[v]. The objective is
  bounded below by |S|. The solution is checked by computing its zero set:
  the points where it is negative, or zero outside of D, are added to S and
  the process repeated. It terminates since the functions of W with zero
  set D form, up to scaling, the relative interior of a polyhedral face.

  The function depends only on D and the order of the points in S, not on
  the history of the computation: taking instead the function produced by
  a flip, f1 - lambda f2 + k f3, makes the coefficients grow at every
  generation of flips, and the closest vector problems with them.
 */
template <typename T, typename Tint>
MyMatrix<T> erdahl_canonical_function(ErdahlFunctionSpace<T> const &W,
                                      ErdahlLatticeSet<Tint> const &ls_in,
                                      ErdahlCanonicalHint<Tint> const *hint,
                                      std::ostream &os) {
#ifdef TIMINGS_ERDAHL_POLYHEDRON
  MicrosecondTime time;
  size_t n_iter = 0;
  int64_t time_lp = 0, time_zero = 0, time_pool = 0;
  MicrosecondTime time_step;
#endif
  ErdahlLatticeSet<Tint> ls = erdahl_canonical_lattice_set(ls_in);
  int np1 = ls.EXT.cols();
  int n = np1 - 1;
  std::vector<MyMatrix<T>> basis = erdahl_vanishing_space(W, ls);
  int r = basis.size();
  if (r == 0) {
    std::cerr << "ERDAHL: no function of W vanishes on the set\n";
    throw TerminalException{1};
  }
  /*
    The bookkeeping (evaluations of the basis at the points, membership in
    D, generation of the pool) is done on machine integers: the points are
    near D and the basis of Space_W(D) is integral with small coefficients.
    It is checked that no overflow can occur, and the exact arithmetic is
    used otherwise. The linear programs and the zero sets are exact.
   */
  using Tpt = std::vector<int64_t>;
  ErdahlAdaptedBasis<Tint> ab = erdahl_adapted_basis(ls.L, n);
  int64_t const coord_max = int64_t(1) << 20;
  int64_t const coeff_max = int64_t(1) << 18;
  auto fits = [&](Tint const &x, int64_t bound) -> bool {
    Tint b = UniversalScalarConversion<Tint, int64_t>(bound);
    return x < b && x > -b;
  };
  // The basis as int64 matrices if possible.
  bool basis_fast = true;
  std::vector<std::vector<int64_t>> basis_i64(r, std::vector<int64_t>(np1 * np1));
  for (int i = 0; i < r && basis_fast; i++) {
    for (int a = 0; a < np1 && basis_fast; a++) {
      for (int b = 0; b < np1; b++) {
        T const &val = basis[i](a, b);
        if (!IsInteger(val)) {
          basis_fast = false;
          break;
        }
        Tint val_i = UniversalScalarConversion<Tint, T>(val);
        if (!fits(val_i, coeff_max)) {
          basis_fast = false;
          break;
        }
        basis_i64[i][a * np1 + b] = UniversalScalarConversion<int64_t, Tint>(val_i);
      }
    }
  }
  auto to_pt = [&](MyVector<Tint> const &e) -> std::optional<Tpt> {
    Tpt pt(np1);
    for (int a = 0; a < np1; a++) {
      if (!fits(e(a), coord_max)) {
        return {};
      }
      pt[a] = UniversalScalarConversion<int64_t, Tint>(e(a));
    }
    return pt;
  };
  auto to_vec = [&](Tpt const &pt) -> MyVector<Tint> {
    MyVector<Tint> e(np1);
    for (int a = 0; a < np1; a++) {
      e(a) = UniversalScalarConversion<Tint, int64_t>(pt[a]);
    }
    return e;
  };
  // The evaluation row (f_i(e))_i, with |f_i(e)| < coeff_max (n+1)^2
  // coord_max^2 < 2^63 when both bounds hold.
  auto eval_row = [&](MyVector<Tint> const &e) -> MyVector<T> {
    MyVector<T> V(r);
    std::optional<Tpt> opt = to_pt(e);
    if (basis_fast && opt && np1 <= 16) {
      Tpt const &pt = *opt;
      for (int i = 0; i < r; i++) {
        std::vector<int64_t> const &B = basis_i64[i];
        int64_t sum = 0;
        for (int a = 0; a < np1; a++) {
          if (pt[a] == 0) {
            continue;
          }
          int64_t part = 0;
          for (int b = 0; b < np1; b++) {
            part += B[a * np1 + b] * pt[b];
          }
          sum += part * pt[a];
        }
        V(i) = UniversalScalarConversion<T, int64_t>(sum);
      }
      return V;
    }
    for (int i = 0; i < r; i++) {
      V(i) = EvaluationQuadForm<T, Tint>(basis[i], e);
    }
    return V;
  };
  // Membership in D, on machine integers.
  MyMatrix<Tint> const &Ainv = ab.AffBasisInv;
  MyMatrix<Tint> const &Aff = ab.AffBasis;
  bool basis_ab_fast = true;
  for (int a = 0; a < np1; a++) {
    for (int b = 0; b < np1; b++) {
      if (!fits(Ainv(a, b), coeff_max) || !fits(Aff(a, b), coeff_max)) {
        basis_ab_fast = false;
      }
    }
  }
  std::unordered_set<Tpt> set_ext;
  std::unordered_set<MyVector<Tint>> set_ext_vec;
  for (int i = 0; i < ls.EXT.rows(); i++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i);
    set_ext_vec.insert(e);
    std::optional<Tpt> opt = to_pt(e);
    if (opt) {
      set_ext.insert(*opt);
    }
  }
  int p_ab = ab.p;
  bool is_polytope = ls.L.rows() == 0;
  auto is_in_set_pt = [&](Tpt const &pt) -> bool {
    if (is_polytope) {
      return set_ext.count(pt) > 0;
    }
    // u = pt Ainv with the L coordinates zeroed, then u Aff.
    Tpt u(np1, 0);
    for (int b = 0; b <= p_ab; b++) {
      int64_t sum = 0;
      for (int a = 0; a < np1; a++) {
        sum += pt[a] * UniversalScalarConversion<int64_t, Tint>(Ainv(a, b));
      }
      u[b] = sum;
    }
    Tpt img(np1, 0);
    for (int b = 0; b < np1; b++) {
      int64_t sum = 0;
      for (int a = 0; a <= p_ab; a++) {
        sum += u[a] * UniversalScalarConversion<int64_t, Tint>(Aff(a, b));
      }
      img[b] = sum;
    }
    return set_ext.count(img) > 0;
  };
  auto is_in_set = [&](MyVector<Tint> const &e) -> bool {
    std::optional<Tpt> opt = to_pt(e);
    if (opt && basis_ab_fast) {
      return is_in_set_pt(*opt);
    }
    if (is_polytope) {
      return set_ext_vec.count(e) > 0;
    }
    return set_ext_vec.count(erdahl_canonical_representative(ab, e)) > 0;
  };
  std::unordered_set<MyVector<Tint>> set_S;
  std::vector<MyVector<Tint>> l_S;
  // The evaluation rows of the points of S, computed once.
  std::vector<MyVector<T>> l_S_eval;
  // The pool, with the evaluation rows on machine integers (only used when
  // the basis is on machine integers).
  std::vector<Tpt> l_pool;
  std::vector<std::vector<int64_t>> l_pool_eval;
  auto eval_row_i64 = [&](Tpt const &pt) -> std::vector<int64_t> {
    std::vector<int64_t> V(r);
    for (int i = 0; i < r; i++) {
      std::vector<int64_t> const &B = basis_i64[i];
      int64_t sum = 0;
      for (int a = 0; a < np1; a++) {
        if (pt[a] == 0) {
          continue;
        }
        int64_t part = 0;
        for (int b = 0; b < np1; b++) {
          part += B[a * np1 + b] * pt[b];
        }
        sum += part * pt[a];
      }
      V[i] = sum;
    }
    return V;
  };
  auto insert = [&](MyVector<Tint> const &e) -> void {
    // Directions (0, v) are never in the set.
    if (e(0) != 0 && is_in_set(e)) {
      return;
    }
    if (set_S.count(e) > 0) {
      return;
    }
    set_S.insert(e);
    l_S.push_back(e);
    l_S_eval.push_back(eval_row(e));
  };
  // The initial points: the neighbors e +- e_i of the representatives. The
  // points x_i + x_j - x_k for representatives x_i, x_j, x_k (translated by
  // the isotropy basis), as in the GAP function ExtendByTriangleIneq, are
  // where the vertices of the adjacent Delaunay polytopes typically are: they
  // form a pool of lazy constraints, added when violated. The linear programs
  // stay small and the zero set, a closest vector problem, is computed only
  // when the pool is satisfied.
  for (int i_ext = 0; i_ext < ls.EXT.rows(); i_ext++) {
    MyVector<Tint> e = GetMatrixRow(ls.EXT, i_ext);
    for (int i = 0; i < n; i++) {
      for (int sign = -1; sign <= 1; sign += 2) {
        MyVector<Tint> f = e;
        f(1 + i) += sign;
        insert(f);
      }
    }
  }
  if (hint) {
    for (auto &e : hint->extra_points) {
      insert(e);
    }
  }
  {
    // The pool: the one of the super polyhedron if given, else its own.
    // Without machine integers it is skipped: the cuts from the zero sets
    // alone also converge, only more slowly.
    std::optional<std::vector<Tpt>> opt_pool =
        hint ? hint->pool : erdahl_triangle_pool(ls);
    if (opt_pool && basis_fast && np1 <= 16) {
      for (auto &f : *opt_pool) {
        bool ok = true;
        for (int a = 0; a < np1; a++) {
          if (f[a] >= coord_max || f[a] <= -coord_max) {
            ok = false;
          }
        }
        if (!ok) {
          continue;
        }
        l_pool_eval.push_back(eval_row_i64(f));
        l_pool.push_back(f);
      }
    }
  }
  // The objective is fixed by the initial points: the points added later
  // are constraints only, otherwise the far away points pull the optimum
  // away from D.
  size_t n_objective = l_S.size();
#ifdef TIMINGS_ERDAHL_POLYHEDRON
  time_pool = time_step.eval_int64();
#endif
  auto get_function = [&](MyVector<T> const &c) -> MyMatrix<T> {
    MyMatrix<T> F = ZeroMatrix<T>(np1, np1);
    for (int i = 0; i < r; i++) {
      if (c(i) != 0) {
        F += c(i) * basis[i];
      }
    }
    return F;
  };
  while (true) {
#ifdef TIMINGS_ERDAHL_POLYHEDRON
    n_iter++;
#endif
    int n_S = l_S.size();
    MyMatrix<T> ListIneq(n_S, r + 1);
    MyVector<T> ToBeMinimized = ZeroVector<T>(r + 1);
    for (int i_S = 0; i_S < n_S; i_S++) {
      ListIneq(i_S, 0) = -1;
      for (int i = 0; i < r; i++) {
        T const &val = l_S_eval[i_S](i);
        ListIneq(i_S, 1 + i) = val;
        if (static_cast<size_t>(i_S) < n_objective) {
          ToBeMinimized(1 + i) += val;
        }
      }
    }
#ifdef TIMINGS_ERDAHL_POLYHEDRON
    (void)time_step.eval_int64();
#endif
    LpSolution<T> eSol = SIMPLEX_LinearProgramming(ListIneq, ToBeMinimized, os);
#ifdef TIMINGS_ERDAHL_POLYHEDRON
    time_lp += time_step.eval_int64();
#endif
    if (!eSol.DirectSolution || !eSol.DualSolution) {
      std::cerr << "ERDAHL: the linear program for the function should have "
                   "an optimal solution\n";
      throw TerminalException{1};
    }
    MyMatrix<T> F = RemoveFractionMatrix(get_function(*eSol.DirectSolution));
    {
      // The violated constraints of the pool, the most violated first.
      // The values are filtered in floating point: only those within the
      // rounding tolerance of the threshold 1 are computed exactly.
      MyVector<T> const &c = *eSol.DirectSolution;
      std::vector<double> c_d(r);
      for (int i = 0; i < r; i++) {
        c_d[i] = UniversalScalarConversion<double, T>(c(i));
      }
      std::vector<std::pair<T, size_t>> l_viol;
      for (size_t i_pool = 0; i_pool < l_pool.size(); i_pool++) {
        std::vector<int64_t> const &ev = l_pool_eval[i_pool];
        double val_d = 0;
        double mag = 0;
        for (int i = 0; i < r; i++) {
          double term = c_d[i] * static_cast<double>(ev[i]);
          val_d += term;
          mag += std::abs(term);
        }
        double tol = 1e-9 * (1 + mag);
        if (val_d > 1 + tol) {
          continue;
        }
        T val(0);
        for (int i = 0; i < r; i++) {
          val += c(i) * UniversalScalarConversion<T, int64_t>(ev[i]);
        }
        if (val < 1) {
          l_viol.push_back({val, i_pool});
        }
      }
      if (!l_viol.empty()) {
        std::sort(l_viol.begin(), l_viol.end());
        size_t n_add = std::min(l_viol.size(), static_cast<size_t>(4 * r));
        for (size_t u = 0; u < n_add; u++) {
          insert(to_vec(l_pool[l_viol[u].second]));
        }
        continue;
      }
    }
#ifdef TIMINGS_ERDAHL_POLYHEDRON
    (void)time_step.eval_int64();
#endif
    auto res = erdahl_zero_set<T, Tint>(F, os);
#ifdef TIMINGS_ERDAHL_POLYHEDRON
    time_zero += time_step.eval_int64();
#endif
#ifdef DEBUG_ERDAHL_POLYHEDRON
    os << "ERDAHL: canonical_function |S|=" << l_S.size() << " r=" << r
       << " d=" << ls.L.rows() << " |EXT|=" << ls.EXT.rows();
    if (res.is_err()) {
      os << " err=" << res.get_err().reason
         << " |points|=" << res.get_err().points.size()
         << " |directions|=" << res.get_err().directions.size();
      for (auto &e : res.get_err().points) {
        os << " pt=" << StringVector(e) << " val=" << EvaluationQuadForm<T, Tint>(F, e);
      }
    } else {
      os << " zero set |EXT|=" << res.get_ok().EXT.rows()
         << " d=" << res.get_ok().L.rows();
    }
    os << "\n";
#endif
    if (res.is_err()) {
      // As in Algorithm "TestRealizability", a direction v along which F is
      // unbounded below gives the constraint Quad(f)[v] >= 1: it is valid
      // since the kernel of the quadratic part of a function with zero set
      // D is L(D) tensor R. The row (0, v) expresses it with the same
      // evaluation as a point. A far away negative point alone would only
      // cut a sliver off the feasible set.
      for (auto &e : res.get_err().directions) {
        insert(e);
      }
      for (auto &e : res.get_err().points) {
        insert(e);
      }
      continue;
    }
    ErdahlLatticeSet<Tint> zs = erdahl_canonical_lattice_set(res.get_ok());
    if (zs.EXT == ls.EXT && zs.L == ls.L) {
#ifdef TIMINGS_ERDAHL_POLYHEDRON
      os << "|ERDAHL: canonical_function n_iter=" << n_iter
         << " |S|=" << l_S.size() << " r=" << r << " |pool|=" << l_pool.size()
         << " |EXT|=" << ls.EXT.rows() << " d=" << ls.L.rows()
         << " time_pool=" << time_pool << " time_lp=" << time_lp
         << " time_zero=" << time_zero << "|=" << time << "\n";
#endif
      return F;
    }
    // Zeros outside of the set. An isotropy direction of the zero set that
    // is not in L(D) gives the constraint Quad(f)[v] >= 1, then the
    // representatives and their translates are added.
    size_t n_before = l_S.size();
    for (int j = 0; j < zs.L.rows(); j++) {
      MyVector<Tint> v = GetMatrixRow(zs.L, j);
      if (ls.L.rows() == 0 || !SolutionIntMat(ls.L, v)) {
        insert(erdahl_direction_row(v));
      }
    }
    for (int i_ext = 0; i_ext < zs.EXT.rows(); i_ext++) {
      MyVector<Tint> e = GetMatrixRow(zs.EXT, i_ext);
      insert(e);
      for (int j = 0; j < zs.L.rows(); j++) {
        for (int sign = -1; sign <= 1; sign += 2) {
          MyVector<Tint> f = e;
          for (int k = 0; k < n; k++) {
            f(1 + k) += sign * zs.L(j, k);
          }
          insert(f);
        }
      }
    }
    if (l_S.size() == n_before) {
      std::cerr << "ERDAHL: the zero set differs from the set but no new "
                   "point was found\n";
      throw TerminalException{1};
    }
  }
}

/*
  The Delaunay polyhedron Z(G) cap D3, for G >= 0 on D3, with its canonical
  function. The zero set is computed on the cosets of D3. With
  compute_function = false the function is left empty (0 x 0), for a
  caller that may discard the polyhedron: see erdahl_ensure_function.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_polyhedron_extension(ErdahlFunctionSpace<T> const &W,
                            MyMatrix<T> const &G,
                            DelaunayPolyhedron<T, Tint> const &D3,
                            bool compute_function,
                            ErdahlCanonicalHint<Tint> const *hint,
                            std::ostream &os) {
  ErdahlLatticeSet<Tint> zs =
      erdahl_zero_set_on_set(G, erdahl_lattice_set(D3), os);
  if (!erdahl_is_full_dimensional(zs)) {
    std::cerr << "ERDAHL: the zero set is empty or not full dimensional\n";
    throw TerminalException{1};
  }
  ErdahlLatticeSet<Tint> can = erdahl_canonical_lattice_set(zs);
  if (!compute_function) {
    return {can.EXT, can.L, MyMatrix<T>(0, 0)};
  }
  MyMatrix<T> F = erdahl_canonical_function(W, can, hint, os);
  return {can.EXT, can.L, F};
}

// Computes the canonical function of D if it was left empty.
template <typename T, typename Tint>
void erdahl_ensure_function(ErdahlFunctionSpace<T> const &W,
                            DelaunayPolyhedron<T, Tint> &D, std::ostream &os) {
  if (D.F.rows() == 0) {
    D.F = erdahl_canonical_function<T, Tint>(W, erdahl_lattice_set(D),
                                             nullptr, os);
  }
}

template <typename T, typename Tint>
void erdahl_check_polyhedron(ErdahlFunctionSpace<T> const &W,
                             DelaunayPolyhedron<T, Tint> const &D,
                             std::ostream &os) {
  if (!erdahl_is_in_space(W, D.F)) {
    std::cerr << "ERDAHL: the function is not in the space W\n";
    throw TerminalException{1};
  }
  DelaunayPolyhedron<T, Tint> D2 =
      erdahl_polyhedron_from_function<T, Tint>(D.F, os);
  if (D2.L != D.L || D2.EXT != D.EXT) {
    std::cerr << "ERDAHL: the zero set of the function is not the polyhedron\n";
    throw TerminalException{1};
  }
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_POLYHEDRON_H_
// clang-format on
