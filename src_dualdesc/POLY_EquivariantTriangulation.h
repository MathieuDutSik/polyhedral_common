// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
#define SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_

// clang-format off
#include "POLY_RecursiveDualDesc.h"
// clang-format on

// The computation of several things require using triangulations.
//
// We have essentially two methods for computing them:
// * The LRS can produce a triangulation.
// * For a symmetric polytope, compute the barycenter.
//   Then from the list of orbit of facets, compute the
//   decomposition from the isobarycenter.
//
// We cannot compute an equivariant triangulation in a
// simple way but fortunately, we do not really neeed that.
// We do not store the triangulations, that can be done
// by the functions but it is not a requirement.
//
// Possible usage for this:
// * Computing the vector inside a cone of fixed norm in a
//   Polyhedral cone using Copositive programming.
// * Computing an integral over a polytope.
//
// The templatized functions must look like:
// * The f_lrs that takes a triangle in lrs (A mapping of the function
//     f used for Kernel_Simplices_cond.
// * A type T_merge that is updated.
// * A function f_compilation that compiles everything that has
//   been computed.
//   --For the integral of a polynomial, it is just an average over
//     a group.
//   --For the fixed norm vectors, it is just the compilation of the
//     database, so nothing particular.
// * The function is then

template <typename T> struct EquivariantTriangulation {
  MyMatrix<T> RelVectors;
  vectface vf_trig;
};

template <typename T> struct ExtendibleEquivariantTriangulation {
  std::vector<MyVector<T>> ListVert;
  std::unordered_map<MyVector<T>, size_t> MapVert;
  std::vector<std::vector<size_t>> ListTrig;
};

// clang-format off
#endif  // SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
// clang-format on
