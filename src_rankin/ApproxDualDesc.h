// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_APPROXDUALDESC_H_
#define SRC_RANKIN_APPROXDUALDESC_H_

/*
  We want to deal with the dual descriptions computed over the
  double number and we want to reduce them.
  ---
  Several phenomena can occur:
  * An artificial facet shows up that does not actually correspond to a real facet.
  * A facet that actually merged into a bigger one.
  The first phenomenon is detected by the tolDet and the second by the tolFacet.

  */
template<typename T>
vectface ReduceVectfaceApproximate(MyMatrix<T> const& EXT, vectface const& vf, T const& tolDet, T const& tolFacet) {

  return ProjMat;
}







// clang-format off
#endif  // SRC_RANKIN_APPROXDUALDESC_H_
// clang-format on
