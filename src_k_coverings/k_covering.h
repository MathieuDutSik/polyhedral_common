// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_K_COVERING_K_COVERINGS_H_
#define SRC_K_COVERING_K_COVERINGS_H_

// clang-format off
#include "Shvec_exact.h"
#include "LatticeDelaunay.h"
// clang-format on


/*
  Given a lattice L, and an integer k, the k-covering radius r
  is such that for each x in R^n, there are k points of L in
  the sphere of center x and radius r.
  ----
  Now, how to make that work with an analysis similar to the
  Delaunay / Voronoi polytopes.
  What we can do:
  * We can clearly define for a fixed set of points v_1, ...., v_k
    the set of constraints
    || v - x || >= || v_i - x || for 1 <= i <= k and v != v_i.
    That is analogous to the Voronoi polytope. And it is a polytope.
    The set of vertices is named S={v_1, ...., v_k} and the
    Voronoi is named Vor(S).
  * It is not just an intersection of Voronoi polytopes, since we
    have to exclude inequalities with v = v_j. But it looks
    computable:
    + Take some vectors from the Voronoi relevant ones, consider them
      as candidate.
    + Then compute the vertices of those polytopes.
    + Check that the vertices satisfy those inequalities. If not
      add more vertices.
  * What geometry we can extract from this:
    + The facets of the polytope indicate how to switch from
      {v_1, ...., v_k} to something else.
      Note that the switch may be for more than of the vertices
      of the set S.
    + Altogether we should be able to get all the possible
      sets S up to equivalence.
    + The vertices of Vor(S) and their adjacencies should indicate
      us how to get the defining inequalities on the matrix
      coefficients in order to get the (L,k)-types domains
  * 

 */





// clang-format off
#endif  // SRC_K_COVERING_K_COVERINGS_H_
// clang-format on
