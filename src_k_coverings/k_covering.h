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
    + Together, they make a nice tiling of the space by polytopes.
      By itself, this is quite nice.
    + The vertices of Vor(S) and their adjacencies should indicate
      us how to get the defining inequalities on the matrix
      coefficients in order to get the (L,k)-types domains
  * Now, the vertices of those polytopes would be the analogue
    of the Delaunay polytopes. But can we get a nice structure
    about them?
    + We can get an initial vertex by the same technique that we
      use for finding an initial Delaunay polytope.
    + The mystery of the story is how to define the vertices of that
      Delaunay polytope? The natural way would be to select the
      vectors v showing up in the definition. However, if we are
      doing so, we are certainly not going to get a tiling. But
      that might be the play here.
    + Pushing in the direction of the (L,k)-polytopes makes
      sense if we want to define the (L,k)-types.
    + The theory by Aurenhammer 1990 is what we need. It should
      allow to find adjacent (L,k)-polytopes.
  * The theory of the polytopes Q_k. For each of the Voronoi
    polytopes, we form the points
      (sum(S), p(S)) in R^(d+1) with
      sum(S) = sum_{x\in S} x
      and p(S) = sum_{x\in S} ||x||^2
    + It is quite obvious just from stating it that this is
      the right construction for us.
    + For the case of the A2, it gives the kagome tiling
      by hexagon and triangles. Very nice. This is what we want
    + But it requires some work since we need to find the
      adjacent (L,k)-polytopes.
  * Theory of (L,k)-types:
    + We need that theory because otherwise, we cannot do much
      with respect to publishing.
    + The theory of Aurenhammer 1990 of (L,k)-polytopes seems
      to be a direct generalization of what we have for k=1.
      So, likely, all the constructions of T-spaces of Delaunay
      polytopes should apply here.

    References:
    * Lee 1982a
    * Chazelle and Edelsbrunner 1987:
    * Clarkson 1987
    * Herbert Edelsbrunner, Alexey Garber, Teresa Heiss, Morteza
      Saghafian, On Angles in Higher Order Brillouin Tessellations and
      Related Tilings in the Plane.
      A recent paper on Order-k Voronoi diagrams. Mentions the Delaunay
      version.
    * Franz Aurenhammer, Rolf Klein, Der-Tsai Lee, Voronoi Diagrams
      and Delaunay Triangulations.
      Section 6.5 has some nice info, but no Delaunay diagrams.
    * Franz Aurenhammer, A new duality result concerning Voronoi
      diagrams, 1990.
      Introduce a polytope in Q_{k+1} in R^(d+1) for computing it.
      The Delaunay_k polytopes and the Voronoi_k polytopes. This is
      the thing to pursue.
 */





// clang-format off
#endif  // SRC_K_COVERING_K_COVERINGS_H_
// clang-format on
