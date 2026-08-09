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
  // What reading a transformation off a matrix of the acting frame needs
  // besides the matrix. Nothing for a lattice.
  struct Tcontext {};
  template <typename T>
  static MyMatrix<Tring> from_field(MyMatrix<T> const &M) {
    return UniversalMatrixConversion<Tring, T>(M);
  }
  template <typename T>
  static MyMatrix<Tring> from_field_acting(MyMatrix<T> const &M,
                                           [[maybe_unused]] Tcontext const &ctx) {
    return UniversalMatrixConversion<Tring, T>(M);
  }
  template <typename T>
  static MyMatrix<T> to_field(MyMatrix<Tring> const &x) {
    return UniversalMatrixConversion<T, Tring>(x);
  }
  // The matrix acting on the coordinates in which the vertices are
  // stored. For a lattice those are the coordinates themselves, so this
  // is the matrix; for a periodic point set the vertices are stored
  // scaled and the acting matrix differs from the one above.
  template <typename T>
  static MyMatrix<T> to_field_acting(MyMatrix<Tring> const &x) {
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

// The matrix acting on the stored vertex coordinates, which is what the
// flipping needs since it forms products of the vertex matrices with it.
template <typename T, typename Ttransform>
MyMatrix<T> TransformToActingMatrix(Ttransform const &x) {
  return transform_traits<Ttransform>::template to_field_acting<T>(x);
}

/*
  The transformation of a matrix of the acting frame. It cannot be read
  off the matrix alone in the periodic case: on the coordinates scaled by
  the denominator the translation is integral, so nothing in the matrix
  says what the denominator is. The context carries it.

  It is also why the denominator is a feature of the point set and not of
  each transformation: were it to vary, the product and the inverse would
  have to reconcile the denominators of their operands, which is exactly
  the reduction the algebra is built to avoid.
 */
template <typename Ttransform, typename T>
Ttransform TransformFromActingMatrix(
    MyMatrix<T> const &M,
    typename transform_traits<Ttransform>::Tcontext const &ctx) {
  return transform_traits<Ttransform>::from_field_acting(M, ctx);
}

// clang-format off
#endif  // SRC_DELAUNAY_TRANSFORMTRAITS_H_
// clang-format on
