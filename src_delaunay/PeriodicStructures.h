// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_PERIODICSTRUCTURES_H_
#define SRC_DELAUNAY_PERIODICSTRUCTURES_H_

// clang-format off
#include "MAT_MatrixInt.h"
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
  // The coset entries are scaled together with a 1, so that the scaled
  // entries are the numerators and the scaled 1 is the denominator N.
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
