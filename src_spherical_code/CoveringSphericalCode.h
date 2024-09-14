// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_COVERINGSPHERICALCODE_H_
#define SRC_DUALDESC_COVERINGSPHERICALCODE_H_

/*
  We want to enumerate the covering radius of
  a spherical code.
  ---
  A direct approach is simply to enumerate the
  centers of the empty circumradius.
  ---
  Let us use some notations and write the
  vectors of the spherical code as
  {x_1, ...., x_M}. Then we want to find sphercial
  points e such that for some subset S we
  have C = || x_s - e || for s \in S and
  C <= || x_i - e || for all 1 <= i <= M.
  Passing to the squares the expression || x - y||
  expands as x^2 - 2x.y + y^2 which is equal to
  2 - 2x.y since both x and y belong to the sphere.
  Thus the inequality gets expressed as
  C = x_s.e for s \in S and
  C >= x_i.e for all 1 <= i <= M.
  (we are pretty indifferent about the constant C
  in question)
  ---
  The sets as defining uniquely gets us vertices
  of the polyhedral cone defined by the rays
  (1, x_i).
  So, the method is simply to enumerate all the
  facets of this cone. Each facet gives the
  incidence and from the incidence we can by
  solving a linear system obtain the center e
  up to some positive scalar. That scalar is then
  determined by the norm equal to 1 for the spherical
  code.
 */

// clang-format off
#endif  // SRC_DUALDESC_COVERINGSPHERICALCODE_H_
// clang-format on
