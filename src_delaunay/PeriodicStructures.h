// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_PERIODICSTRUCTURES_H_
#define SRC_DELAUNAY_PERIODICSTRUCTURES_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "TransformTraits.h"
#include <utility>
#include <vector>
// clang-format on

/*
  Structures for periodic point sets Z^n + {c_1, ..., c_m} and their affine
  transformations. The design decisions (see the iso-Delaunay work):
  --- The cosets c_k have a common denominator N that is a feature of the
  point set: every transformation preserving the point set permutes the
  cosets modulo Z^n, so its translation part lies in (1/N) Z^n as well.
  --- A transformation is x -> x A + w / d with A in GL_n(Z) and w integral.
  Composition of two transformations with the same denominator keeps that
  denominator, and inversion keeps it as well (the linear part inverts over
  the integers). Therefore, as long as no reduction of w / d is attempted,
  the denominators are stable under all the group operations: the whole
  transformation algebra is GCD-free.
 */

#ifdef DEBUG
#define DEBUG_PERIODIC_STRUCTURES
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERIODIC_STRUCTURES
#endif

template <typename Tring> struct PeriodicAffineTransform {
  // The map is x -> x A + w / d for a row vector x.
  MyMatrix<Tring> A;
  MyVector<Tring> w;
  Tring d;
};

template <typename Tring>
bool operator==(PeriodicAffineTransform<Tring> const &x,
                PeriodicAffineTransform<Tring> const &y) {
  // Compared without reduction: equality of maps requires cross products
  // since the denominators are never reduced.
  if (!TestEqualityMatrix(x.A, y.A)) {
    return false;
  }
  int n = x.w.size();
  for (int i = 0; i < n; i++) {
    if (x.w(i) * y.d != y.w(i) * x.d) {
      return false;
    }
  }
  return true;
}

namespace boost::serialization {
template <class Archive, typename Tring>
inline void serialize(Archive &ar, PeriodicAffineTransform<Tring> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("A", eRec.A);
  ar &make_nvp("w", eRec.w);
  ar &make_nvp("d", eRec.d);
}
} // namespace boost::serialization

template <typename Tring>
PeriodicAffineTransform<Tring>
IdentityPeriodicAffineTransform(int const &n, Tring const &d) {
  MyMatrix<Tring> A = IdentityMat<Tring>(n);
  MyVector<Tring> w = ZeroVector<Tring>(n);
  return {std::move(A), std::move(w), d};
}

/*
  Composition, in the same order as the matrix product of the corresponding
  (n+1) x (n+1) matrices: (x * M1) * M2, that is x -> x A1 A2 + (w1 A2)/d1
  + w2/d2. When the denominators agree (the invariant of the periodic
  setting) the common denominator is kept; otherwise the product of the
  denominators is used, without reduction.
 */
template <typename Tring>
PeriodicAffineTransform<Tring>
operator*(PeriodicAffineTransform<Tring> const &x,
          PeriodicAffineTransform<Tring> const &y) {
  MyMatrix<Tring> A = x.A * y.A;
  MyVector<Tring> w1_img = y.A.transpose() * x.w;
  if (x.d == y.d) {
    MyVector<Tring> w = w1_img + y.w;
    return {std::move(A), std::move(w), x.d};
  }
  MyVector<Tring> w = y.d * w1_img + x.d * y.w;
  Tring d = x.d * y.d;
  return {std::move(A), std::move(w), std::move(d)};
}

template <typename Tring>
PeriodicAffineTransform<Tring>
Inverse(PeriodicAffineTransform<Tring> const &x) {
  // A is in GL_n(Z), so the ring inversion is exact; the translation
  // becomes -(w Ainv)/d with the denominator unchanged.
  MyMatrix<Tring> Ainv = Inverse(x.A);
  MyVector<Tring> w = -(Ainv.transpose() * x.w);
  return {std::move(Ainv), std::move(w), x.d};
}

/*
  The transform_traits bridge to the flipping machinery: the field matrix of
  x -> x A + w/d is the (n+1) x (n+1) matrix M with M(0,0)=1, M(i,0)=0,
  M(0,j) = w_j/d and M(i,j) = A(i,j) for 1 <= i, j <= n, acting on rows
  (1, v).
 */
template <typename Tring>
struct transform_traits<PeriodicAffineTransform<Tring>> {
  // The denominator of the point set, which a matrix of the acting frame
  // does not carry: its translation is integral there.
  struct Tcontext {
    Tring N;
  };
  template <typename T>
  static PeriodicAffineTransform<Tring>
  from_field_acting(MyMatrix<T> const &M, Tcontext const &ctx) {
    int n = M.rows() - 1;
    MyMatrix<Tring> A(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A(i, j) = UniversalScalarConversion<Tring, T>(M(i + 1, j + 1));
      }
    }
    MyVector<Tring> w(n);
    for (int j = 0; j < n; j++) {
      w(j) = UniversalScalarConversion<Tring, T>(M(0, j + 1));
    }
    return {std::move(A), std::move(w), ctx.N};
  }
  template <typename T>
  static PeriodicAffineTransform<Tring> from_field(MyMatrix<T> const &M) {
    int n = M.rows() - 1;
#ifdef SANITY_CHECK_PERIODIC_STRUCTURES
    if (M.cols() != n + 1 || M(0, 0) != T(1)) {
      std::cerr << "PERIODIC: from_field: M is not affine normalized\n";
      throw TerminalException{1};
    }
    for (int i = 1; i <= n; i++) {
      if (M(i, 0) != T(0)) {
        std::cerr << "PERIODIC: from_field: nonzero first column\n";
        throw TerminalException{1};
      }
    }
#endif
    MyMatrix<Tring> A(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A(i, j) = UniversalScalarConversion<Tring, T>(M(i + 1, j + 1));
      }
    }
    // The translation is scaled together with a trailing 1: the scaled
    // vector is integral, so its first n entries are the numerators and its
    // last entry is the common denominator, both over the ring. Scaling the
    // translation alone would not do, since the scaling coefficient of
    // RemoveFractionVectorPlusCoeff is itself a fraction in general.
    MyVector<T> trans_ext(n + 1);
    for (int j = 0; j < n; j++) {
      trans_ext(j) = M(0, j + 1);
    }
    trans_ext(n) = T(1);
    FractionVector<T> fr = RemoveFractionVectorPlusCoeff(trans_ext);
    MyVector<Tring> w(n);
    for (int j = 0; j < n; j++) {
      w(j) = UniversalScalarConversion<Tring, T>(fr.TheVect(j));
    }
    Tring d = UniversalScalarConversion<Tring, T>(fr.TheVect(n));
    PeriodicAffineTransform<Tring> ret{std::move(A), std::move(w),
                                       std::move(d)};
#ifdef SANITY_CHECK_PERIODIC_STRUCTURES
    MyMatrix<T> M_back = to_field<T>(ret);
    if (!TestEqualityMatrix(M_back, M)) {
      std::cerr << "PERIODIC: from_field: round trip failed\n";
      throw TerminalException{1};
    }
#endif
    return ret;
  }
  /*
    The matrix acting on the coordinates scaled by the denominator, in
    which the vertices of a periodic tessellation are stored: with u = N x
    the map x -> x A + w / d of denominator d = N reads u -> u A + w, of
    integral translation. That is the form the flipping needs, its
    products of vertex matrices being in the scaled coordinates; the
    inequalities it derives are then scaled by N^2, a positive multiple
    which leaves the iso-Delaunay cone unchanged.
   */
  template <typename T>
  static MyMatrix<T> to_field_acting(PeriodicAffineTransform<Tring> const &x) {
    int n = x.A.rows();
    MyMatrix<T> M = ZeroMatrix<T>(n + 1, n + 1);
    M(0, 0) = T(1);
    for (int j = 0; j < n; j++) {
      M(0, j + 1) = UniversalScalarConversion<T, Tring>(x.w(j));
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        M(i + 1, j + 1) = UniversalScalarConversion<T, Tring>(x.A(i, j));
      }
    }
    return M;
  }
  template <typename T>
  static MyMatrix<T> to_field(PeriodicAffineTransform<Tring> const &x) {
    int n = x.A.rows();
    MyMatrix<T> M = ZeroMatrix<T>(n + 1, n + 1);
    M(0, 0) = T(1);
    T d_T = UniversalScalarConversion<T, Tring>(x.d);
    for (int j = 0; j < n; j++) {
      T w_T = UniversalScalarConversion<T, Tring>(x.w(j));
      M(0, j + 1) = w_T / d_T;
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        M(i + 1, j + 1) = UniversalScalarConversion<T, Tring>(x.A(i, j));
      }
    }
    return M;
  }
};

/*
  A periodic point set Z^n + {c_1, ..., c_m}: the cosets are stored as
  integral numerators with the common denominator N, reduced to [0, N)^n.
 */
template <typename Tring> struct PeriodicPointSet {
  MyMatrix<Tring> cosets_num;
  Tring N;
};

// Entrywise representative in [0, N)^n of a scaled point.
template <typename Tring>
MyVector<Tring> PeriodicVectorMod(MyVector<Tring> const &u, Tring const &N) {
  int n = u.size();
  MyVector<Tring> ret(n);
  for (int i = 0; i < n; i++) {
    ret(i) = ResInt(u(i), N);
  }
  return ret;
}

/*
  Whether the cosets form a group modulo Z^n, in which case the point set
  is a lattice containing Z^n rather than a genuinely periodic set.

  The coset formalism is not the right one for that case -- the set should
  be described as the lattice it is -- so it is rejected at construction.
  Z^2 + {(0,0), (1/2,1/2)} is such a set, the cosets being closed under
  addition; Z^2 + {(0,0), (1/3,0)} is not, since (1/3,0) + (1/3,0) is
  (2/3,0), which is not a coset.

  The set contains the origin after normalization and is finite, so being
  closed under addition is all that has to be checked.
 */
template <typename Tring>
bool IsPeriodicPointSetLattice(PeriodicPointSet<Tring> const &pps) {
  int n_coset = pps.cosets_num.rows();
  int n = pps.cosets_num.cols();
  for (int k1 = 0; k1 < n_coset; k1++) {
    for (int k2 = k1; k2 < n_coset; k2++) {
      MyVector<Tring> eSum(n);
      for (int j = 0; j < n; j++) {
        eSum(j) = pps.cosets_num(k1, j) + pps.cosets_num(k2, j);
      }
      if (!GetCosetIndex(pps, eSum)) {
        return false;
      }
    }
  }
  return true;
}

/*
  Build the periodic point set from the rational coset matrix: N is the
  common denominator of all the entries and the numerators are reduced
  modulo N.

  The set is translated so that its first coset is the origin. A periodic
  point set and its translates have the same Delaunay decomposition up to
  that translation, and the Delaunay geometry needs a point of the set at
  the origin: its defining inequalities read "at least as close to 0 as to
  v", so the origin has to play the role of a point of the set. Rather
  than carrying a base point through the geometry, one coset is made to be
  the origin here, once.
 */
template <typename Tring, typename T>
PeriodicPointSet<Tring>
PeriodicPointSetFromRational_Kernel(MyMatrix<T> const &Cosets) {
  int n_coset = Cosets.rows();
  int n = Cosets.cols();
  // Same trailing-1 trick as in transform_traits::from_field: all the coset
  // entries are scaled together with a 1, so that the scaled entries are
  // the numerators and the scaled 1 is the common denominator N.
  MyVector<T> V(n_coset * n + 1);
  for (int k = 0; k < n_coset; k++) {
    for (int j = 0; j < n; j++) {
      V(k * n + j) = Cosets(k, j);
    }
  }
  V(n_coset * n) = T(1);
  FractionVector<T> fr = RemoveFractionVectorPlusCoeff(V);
  Tring N = UniversalScalarConversion<Tring, T>(fr.TheVect(n_coset * n));
#ifdef SANITY_CHECK_PERIODIC_STRUCTURES
  if (N <= 0) {
    std::cerr << "PERIODIC: PeriodicPointSetFromRational: N should be "
                 "positive\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tring> cosets_num(n_coset, n);
  for (int k = 0; k < n_coset; k++) {
    for (int j = 0; j < n; j++) {
      Tring val =
          UniversalScalarConversion<Tring, T>(fr.TheVect(k * n + j));
      cosets_num(k, j) = ResInt(val, N);
    }
  }
  // The translation putting the first coset at the origin.
  MyVector<Tring> shift = GetMatrixRow(cosets_num, 0);
  for (int k = 0; k < n_coset; k++) {
    for (int j = 0; j < n; j++) {
      Tring val = cosets_num(k, j) - shift(j);
      cosets_num(k, j) = ResInt(val, N);
    }
  }
  // The translation can leave a denominator larger than necessary, and
  // the cost of the coset machinery is in N^n, so it is brought back down
  // to the denominator the translated cosets actually need.
  Tring g = N;
  for (int k = 0; k < n_coset; k++) {
    for (int j = 0; j < n; j++) {
      g = GcdPair(g, cosets_num(k, j));
    }
  }
  if (g != 1) {
    for (int k = 0; k < n_coset; k++) {
      for (int j = 0; j < n; j++) {
        cosets_num(k, j) = QuoInt(cosets_num(k, j), g);
      }
    }
    N = QuoInt(N, g);
  }
  return {std::move(cosets_num), std::move(N)};
}

// The point set, or nothing when the cosets form a group and the set is
// therefore a lattice. For a caller that has to decide rather than to
// require, a random generator of point sets among others.
template <typename Tring, typename T>
std::optional<PeriodicPointSet<Tring>>
PeriodicPointSetFromRational_Opt(MyMatrix<T> const &Cosets) {
  PeriodicPointSet<Tring> ret =
      PeriodicPointSetFromRational_Kernel<Tring, T>(Cosets);
  if (IsPeriodicPointSetLattice(ret)) {
    return {};
  }
  return ret;
}

template <typename Tring, typename T>
PeriodicPointSet<Tring> PeriodicPointSetFromRational(MyMatrix<T> const &Cosets) {
  PeriodicPointSet<Tring> ret =
      PeriodicPointSetFromRational_Kernel<Tring, T>(Cosets);
#ifdef SANITY_CHECK_PERIODIC_STRUCTURES
  if (IsPeriodicPointSetLattice(ret)) {
    std::cerr << "PERIODIC: PeriodicPointSetFromRational: the cosets form a "
                 "group modulo Z^n, so the point set is a lattice and has to "
                 "be described as one rather than by cosets\n";
    throw TerminalException{1};
  }
#endif
  return ret;
}

// Index of the coset matching the scaled point u modulo N, if any.
template <typename Tring>
std::optional<size_t> GetCosetIndex(PeriodicPointSet<Tring> const &pps,
                                    MyVector<Tring> const &u) {
  MyVector<Tring> u_red = PeriodicVectorMod(u, pps.N);
  int n_coset = pps.cosets_num.rows();
  int n = pps.cosets_num.cols();
  for (int k = 0; k < n_coset; k++) {
    bool is_match = true;
    for (int j = 0; j < n; j++) {
      if (pps.cosets_num(k, j) != u_red(j)) {
        is_match = false;
        break;
      }
    }
    if (is_match) {
      return k;
    }
  }
  return {};
}

// clang-format off
#endif  // SRC_DELAUNAY_PERIODICSTRUCTURES_H_
// clang-format on
