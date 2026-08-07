// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_TRANSFORMTRAITS_H_
#define SRC_DELAUNAY_TRANSFORMTRAITS_H_

// clang-format off
#include "MAT_Matrix.h"
#include "TypeConversion.h"
// clang-format on

/*
  The transformations of a Delaunay tessellation are kept in a dedicated
  storage type Ttransform, while the flipping machinery computes its
  transformation chains over the field T. The bridge between the two is
  transform_traits, whose from_field must reject (by termination) a field
  matrix that is not representable in the storage type: that is what makes
  the storage invariant an enforced one.

  The instantiations are:
  --- MyMatrix<Tring> for the lattice case, the transformations being
  integral affine matrices (specialization below),
  --- PeriodicAffineTransform<Tring> for the periodic case, the translation
  part being fractional with a denominator that is a feature of the
  periodic point set (specialization in PeriodicStructures.h).
 */
template <typename Ttransform> struct transform_traits;

template <typename Tring> struct transform_traits<MyMatrix<Tring>> {
  template <typename T>
  static MyMatrix<Tring> from_field(MyMatrix<T> const &M) {
    return UniversalMatrixConversion<Tring, T>(M);
  }
  template <typename T>
  static MyMatrix<T> to_field(MyMatrix<Tring> const &x) {
    return UniversalMatrixConversion<T, Tring>(x);
  }
};

template <typename Ttransform, typename T>
Ttransform TransformFromFieldMatrix(MyMatrix<T> const &M) {
  return transform_traits<Ttransform>::from_field(M);
}

template <typename T, typename Ttransform>
MyMatrix<T> TransformToFieldMatrix(Ttransform const &x) {
  return transform_traits<Ttransform>::template to_field<T>(x);
}

// clang-format off
#endif  // SRC_DELAUNAY_TRANSFORMTRAITS_H_
// clang-format on
