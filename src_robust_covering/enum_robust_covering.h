// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "generalized_polytopes.h"
#include "parall_search.h"
#include "LatticeDelaunay.h"
// clang-format on

/*
  The motivation is explained in
  "New upper bound for lattice covering by spheres"
  https://arxiv.org/pdf/2508.06446
*/

/*
  Given a lattice L, and a point x in L \otimes R,
  we want to find the smallest radius r such that
  the set of lattice points contains a parallelepiped
  of volume 1.
  The enumeration is done directly and that part is
  fine.
  That allow us to get the distance at a point x.
  ----
  What we want to do next is to compute the local
  maximum. This is for sure going to be something
  polyhedral.
  If we have a random point in the lattice, then
  there is a unique parallelohedron of volume 1
  closest to the point x. And among all the vertices
  of the polyhedron, there is just 1 that attains
  the minimum.
  ----
  Now looking at the modelization:
  + The basic combinatorial data is a parallelepiped P
    and a vertex v of P that realize the maximum.
    -
    It is a very precise combinatorial data. This
    is bad for the enumeration as it gets us a
    combinatorial explosion. But there is no
    alternative.
  + A parallelepiped P and a vertex v defines a
    set Ineq(P, v) by the inequalities
    || x - v || <= || x - w || for w vertex of P.
    -
    The Ineq(P, v) is unbounded. The directions can
    be obtained by the vertices adjacent to v.
  + If the point is not random, then there are
    several paralelepipeds that are contained.
    -
    And this creates some equalities.
  + The definition of vertices has to take into
    account both sources of inequalities.
    -
    When doing the switch, we of course have to
    take into account the nature of the inequality.
    Hard combinatorics ahead.
    -
    However, this should not impact the definition
    of the vertices. All the inequalities should have
    equal weight there. Otherwise, we will reach a
    contradiction at some point.
  + A mystery is why we need that, while for Delaunay
    polytopes, it is nowhere near as complex. But,
    I guess this is what it is. We encountered similar
    situation when computing the Voronoi polytope for
    polyhedral norms.
  + The algorithm is to simply take a random point
    and to determine in this point the corresponding cell.
    x The cell data gives a set of defining inequalities.
    x It also gives the objective function since we know
      which point is realizing the maximum. But it is not
      linear.
  + The main object is thus the cell in question. It is
    a nice combinatorial object to compute with.
  + So, what are the needed functions?
    x Finding one initial cell.
    x From one cell, finding the adjacent cells.
    x Finding stabilizer of cell / equivalence.
      We can use the geometrically unique center for computing
      the center and testing equivalence.
 */

/*
  The data structure for the P-polytope:
  (A) The parallelepiped realizing the minimum.
  (B) The list of other parallelepipeds being
  used to define the structure.
  ---- FUNDAMENTAL PROBLEM
  The problems is that if we add some more parallelepipeds
  to the structure and they define additional inequalities
  Therefore, the construction of the P-polytope is not
  canonical. In particular, this makes the decomposition
  not face-to-face. It is not a purely linear theory like
  L-type.
  ----
  If a polytope is correct, that is the function is well defined
  over it and the vertices are fine, then we know that
  we can use it to compute the robust covering density.
  And that works whether it is good or whether it is very refined
  for no particular reason. It just works.
  ---- CANONICITY ----
  Canonicality question: Is there a notion of canonical maximal
  decomposition?
  * In the A2 case, we have something that we clearly
  expect to be the canonical P-polytope: The inner triangle
  formed by two vertices and a center of a Delaunay polytope.
  (Figure 1 of https://arxiv.org/pdf/2508.06446).
  * One facet comes from setting the farthest vertex of the
  parallelepiped. Two facets come from having two parallelepiped
  giving a better parallelepiped. Both of those inequalities
  are perfectly normal and we expect them to define facets.
  * Could it be that we have just those inequalities? What may
  be true in A2 may not be in general.
  * So, only the inequalities of parallelepipeds which are not
  minimal specifying which vertex are the farthest actually
  creates our troubles.
  * This implies at a minimum that we should keep track of
  the origin of the facets.
  * At a minimum we can define a canonical domain to be one
  for which (P, v) is the optimal solution. This is canonic.
  It just might not be convex or connected and can be difficult
  to compute.
  ---- CONSIDERATION
  If we do not have canonicality, can we look for boundaries
  that are not real ones? That is the ones that after removal
  will still guarantee that we are ok, that is that the
  vertices have their norm determinedby the minimal
  parallelepiped.
  Sure, that is possible. But, in a canonical way? Not
  guaranteed at all.
  If we drop a facet, we should actually drop the whole
  parallelepiped. That is the only way that makes sense.
  So, the only facet that are candidates for removal are
  the ones that are inner facet of cells for which the outer
  inequalities do not define a facet.
  But of course, removal is not guaranteed to work:
  (A) The outer inequality might become a facet defining
      inequality.
  (B) The vertices might not be at the right distance.
  So inequality removal is absolutely not guarenteed to be
  possible. It is also possible to have inequality removal
  that are path dependent. No reason to think it cannot
  occur.
  ----
  If we do not have canonicity, can we ensure that if we have
  a facet, then what we have on the other side is still just
  one P-domain?
  If we accept that there can be facet from inner inequalities
  which are not defining inequalities. If we do not have
  canonicity, then a lot of very bad things can occur:
  * We could have face-to-face and that when an over covering
  because we go over the limits and never see that we match
  with another object.
  * Face-to-face is absolutely not guaranteed. Because basically
  everything can happen.
  * We can define some
  ---- THE GENERALIZED POLYTOPE ----
  [ Generalized polytopes are unions of polytopes see the
    other file ]
  Define for a parallelepiped P, the function
  f(x,y) = || x - y ||^2
  Then the generalized polytope Can(P,v) in question will be:
  phi_P(x) = max_{y vertex of P} f(x,y)
  So phi_P(x) >= f(x,y)
  So if the optimal configuration is (P,v) then the domain
  that we want is
  f(x,v) < phi_Q(x) for all parallelepipeds Q != P.
  If we select for each of the point in Q a specific
  point Q_v then the set of inequalities
  f(x,v) < f(x, Q_v).
  What about this object:
  * It will be smaller because f(x, Q_v) <= phi_Q(x).
  * It will be bigger because we select a finite set
    of point.
  But it is a polytope. But Can(P, v) is not necessarily
  connected. It will be an union of polytopes.

  But some features need to be addressed:
  * Can we find from the start enough vertices so that
    the generated P-polytope will actually be smaller.
  * Test that a P-polytope is indeed correct, that is
    all the points inside of it have (P, v) as minimum.
    That does not say that it could not be extended.
  * Testing for each point on the boundary whether
    it can be extended.

  Starting with enough points:
  * Can we get an upper bound on the relevant distances?
    + If we select a parallelepiped P, compute the maximum
      pairwise distance between the vertices. Then that give
      us an upper bound.
    + For the A2 the distance is about (2/3)^2. So, that seems
      overall reasonable.
  * But we need to determine a diameter of the P-polytope.
    + The covering radius could help.
    + Does that really matter?
      - We could compute an initial set of points so that it
        is bounded.
      - From the list of vertices, we get some upper bound.
      - Then use that upper bound for getting a global
        upper bound.

  Testing correctness of the P-polytope:
  * Checking that the vertices are correct. But it is
    only indicative.
  * It seems, that we cannot avoid having the upper bound
    as a fundamental piece of the proof that it is correct.

  What algorithms can we write:
  * The P-polytope as defined by a best (P, v) but also
    some other parallelepipeds (P_i, v_i) with the best
    are kind of useful. But the big problem is that we
    cannot use all of the P_i in the choice since that
    break the P-polytope further and further. That is
    the original source of the problem.
  * But the P-polytopes have the advantage of being
    polytopes and so computable. So, somehow we have to
    use them.
  * So, if we take a generic point x, we find the
    minimum (P_min, v_min). Then we list the neighboring
    parallelepipeds P_i and vertices v_i.
    + If the P_i are far away, then all those
      inequalities are redundant.
  * Now what other vertices w_i of the polytope of P_i
    could so:
    + For a far away P_i the inequality between w_i
      and v_i being the best can split the domain, even
      if very far.
    + So, over the polytope P we define the additional
      f(x, v) <= f(x, w_i) and f(x, v_i) < f(x, w_i),
      and f(x, v_i) < f(x,v) [Because otherwise this
                            is covered by the preceding]
    + If the additional polytope is non-empty, then
      it has to be added to the structure. But it is
      strictly an union. No convexity here. It is a
      generalized polytope.
    + So, we can potentially generate the object.

  How to build the generelized polytope:
  * Compute the upper bound, generate a correct
    initial P-polytope.
  * If a parallelepiped has vertices v_i, w_i and
    f(x,v) >= f(x,v_i) with f(x,v) <= f(w_i) is meaningful.
    Then this implies that the inequality
    f(x,v) <= f(x,v_i) is a facet defining inequality.
  * The inequality from the (P,v) are hard one that cannot
    be changed. They are not something we can pass over.
    Those are the hard facets.
  * Possible algorithm:
    + Generate the boundary of the polytope.
    + Determine the hard facets. Have it ready to be
      removed from insertion, at the start or later in the
      process.
    + For each open facet, consider all the possibilities w_i
      of replacing equalities by
      f(x,v_i) <= f(x,v) <= f(x,w_i)
    + The boundary keeps track of what remains. To avoid
      abstraction leakage, we could have a function that takes
      the boundary and returns an interior point. Then the
      interior point can give you the i_polytope and the i_face.
    + All the inequalities being generated are of the form
      f(x,v) <=/>= f(x,w)
      This is somewhat important as it means that if an
      inequality in two different forms, then they correspond
      to the same vector w.
  * What could we work out to make it work:
    + Keeping track of hard and soft facets is mandatory.
    + If there is no extension from a soft facet, then we
      need to convert it to hard.
    + 


    need to switch to

  ---- THE L-TYPE THEORY SIDE ----
  Here is what we have:
  * Each vertex is defined like for Delaunay.
  * The adjacent vertices should give the defining equalities.
  * What is the structure of the cells?
    x It contains a bunch of inequalities.
    x They can be reduced by linear programming (Clarkson).
    x The flipping is done via computing an adjacent point
      checking if it shares a facet, or not.
    x Add a facetness check for sanity check.
    x There is iteration. But we have iteration for
      Delaunay as well. So, maybe this is what to expect.
  ----
 */

#ifdef DEBUG
#define DEBUG_ENUM_P_POLYTOPES
#endif


#ifdef SANITY_CHECK
#define SANITY_CHECK_ENUM_P_POLYTOPES
#endif

#ifdef DISABLE_DEBUG_ENUM_P_POLYTOPES
#undef DEBUG_ENUM_P_POLYTOPES
#endif





// A family of vectors with the index which is the farthest.
template <typename Tint> struct GenericRobustM {
  int index;
  MyMatrix<Tint> M;
  MyVector<Tint> v_long() const { return GetMatrixRow(M, index); }
  std::vector<MyVector<Tint>> get_short_vertors() const {
    std::vector<MyVector<Tint>> l_v;
    for (int i_row=0; i_row<M.rows(); i_row++) {
      if (i_row != index) {
        MyVector<Tint> V = GetMatrixRow(M, i_row);
        l_v.emplace_back(std::move(V));
      }
    }
    return l_v;
  }
};

template <typename T, typename Tint> struct ExtendedGenericRobustM {
  T max;
  bool is_correct;
  GenericRobustM<Tint> robust_m;
};


template<typename T, typename Tint>
struct FullIneq {
  MyVector<T> eIneq;
  // If none, then the inequality is a hard one.
  // If some, then it is a soft one and that indicates what
  // can change.
  std::optional<GenericRobustM<Tint>> opt_soft;
};

// The convex block of the construction
template <typename T, typename Tint>
struct ConvexBlock {
  std::vector<GenericRobustM<Tint>> list_robust_m;
  SinglePolytope<T> sp;
};

template<typename T>
struct HardConvexBoundary {
  int index_cb; // The corresponding face;
  ConvexBoundary<T> sp;
};

template<typename T, typename Tint>
struct SoftConvexBoundary {
  int index_cb; // The corresponding face;
  ConvexBoundary<T> cb;
  std::vector<MyVector<Tint>> l_excluded_max; // The excluded vectors. Cannot use the ConvexBlock since they can vary.
  GenericRobustM<Tint> robust_m;
};


/*
  The robust_m_min is defining the P-polytope.
  This is what we are after in the end.
  ----
  It is a full enumeration result.
 */
template <typename T, typename Tint> struct PVoronoi {
  GenericRobustM<Tint> robust_m_min;
  std::vector<ConvexBlock<T,Tint>> l_cb; // The list of convex blocks.
  std::vector<HardConvexBoundary<T>> l_hcb;
  GeneralizedPolytope<T> gp;
  MyMatrix<T> EXT; // The list of vertices as defined from the generalized polytope.
};

/*
  The robust_m_min is defining the P-polytope.
  This is what we are after in the end.
  ----
  It is a partial enumeration result.
 */
template <typename T, typename Tint> struct PVoronoiPart {
  GenericRobustM<Tint> robust_m_min;
  std::vector<ConvexBlock<T,Tint>> l_cb; // The list of convex blocks.
  std::vector<HardConvexBoundary<T>> l_hcb;
  std::vector<SoftConvexBoundary<T,Tint>> l_scb;
};








template <typename T, typename Tint>
ExtendedGenericRobustM<T, Tint>
get_generic_robust_m(MyMatrix<Tint> const &M, MyMatrix<T> const &G,
                     MyVector<T> const &eV, [[maybe_unused]] std::ostream &os) {
  T max(0);
  int best_index = 0;
  int n_ineq = M.rows();
  size_t n_att = 0;
  for (int index = 0; index < n_ineq; index++) {
    MyVector<Tint> fV = GetMatrixRow(M, index);
    MyVector<T> diff = UniversalVectorConversion<T, Tint>(fV) - eV;
    T norm = EvaluationQuadForm(G, diff);
#ifdef DEBUG_ENUM_P_POLYTOPES
    double norm_d = UniversalScalarConversion<double, T>(norm);
    os << "ROBUST:   get_generic_robust_m, index=" << index
       << " fV=" << StringVector(fV) << " norm=" << norm << " norm_d=" << norm_d
       << "\n";
#endif
    if (index == 0) {
      max = norm;
      best_index = index;
      n_att = 1;
    } else {
      if (norm == max) {
        n_att += 1;
      } else {
        if (norm > max) {
          max = norm;
          best_index = index;
          n_att = 1;
        }
      }
    }
  }
  bool is_correct = true;
  if (n_att > 1) {
    is_correct = false;
  }
#ifdef DEBUG_ENUM_P_POLYTOPES
  double max_d = UniversalScalarConversion<double, T>(max);
  os << "ROBUST:   get_generic_robust_m, best_index=" << best_index
     << " max=" << max << " max_d=" << max_d << "\n";
#endif
  GenericRobustM<Tint> robust_m{best_index, M};
  return {max, is_correct, robust_m};
};

/*
  We have G[x - v_short] <= G[x - v_long]
  So,
  -2 x G v_short + G[v_short] <= -2 x G v_long + G[v_long]
  which gets
  0 <= G[v_long] - G[v_short] + x ( 2 G (v_short - v_long))
 */
template <typename T, typename Tint>
MyVector<T> get_ineq(MyMatrix<T> const &G, MyVector<Tint> const &v_short,
                     MyVector<Tint> const &v_long) {
  int dim = G.rows();
  T norm_short = EvaluationQuadForm(G, v_short);
  T norm_long = EvaluationQuadForm(G, v_long);
  MyVector<Tint> diff = v_short - v_long;
  MyVector<T> diff_T = UniversalVectorConversion<T, Tint>(diff);
  MyVector<T> G_v = G * diff_T;
  MyVector<T> ineq(1 + dim);
  ineq(0) = norm_long - norm_short;
  for (int i = 0; i < dim; i++) {
    ineq(1 + i) = 2 * G_v(i);
  }
  return ScalarCanonicalizationVector(ineq);
}

// In a robust structure robust_m, the longest vector
template <typename T, typename Tint>
void insert_inner_ineqs_parallelepiped(GenericRobustM<Tint> const &robust_m,
                                       MyMatrix<T> const &G,
                                       std::vector<FullIneq<T,Tint>> & list_full_ineq,
                                       [[maybe_unused]] std::ostream &os) {
  int n_row = robust_m.M.rows();
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST:   insert_inner_ineqs_parallelepiped, M=\n";
  WriteMatrix(os, robust_m.M);
#endif
  MyVector<Tint> v_long = robust_m.v_long();
  for (int i = 0; i < n_row; i++) {
    if (i != robust_m.index) {
      MyVector<Tint> v_short = GetMatrixRow(robust_m.M, i);
      MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   insert_inner_ineqs_parallelepiped, eIneq="
         << StringVector(eIneq) << " i=" << i << " index=" << robust_m.index
         << "\n";
#endif
      FullIneq<T,Tint> full_ineq{eIneq, {}};
      list_full_ineq.push_back(full_ineq);
    }
  }
}

// The farthest vector of the robust structure has to be farther from the
// v_short (which is the farthest of the best one)
template <typename T, typename Tint>
void insert_outer_ineqs_parallelepiped(GenericRobustM<Tint> const &robust_m,
                                       MyMatrix<T> const &G,
                                       MyVector<Tint> const &v_short,
                                       std::vector<FullIneq<T, Tint>> &list_full_ineq,
                                       [[maybe_unused]] std::ostream &os) {
  MyVector<Tint> v_long = robust_m.v_long();
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_short="
     << StringVector(v_short) << "\n";
  os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_long="
     << StringVector(v_long) << "\n";
#endif
  if (v_short != v_long) {
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   insert_outer_ineqs_parallelepiped, eIneq="
       << StringVector(eIneq) << "\n";
#endif
    FullIneq<T,Tint> full_ineq{eIneq, robust_m};
    list_full_ineq.push_back(full_ineq);
  } else {
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_short="
       << StringVector(v_short) << "\n";
#endif
  }
}

template <typename T, typename Tint>
void insert_excluded_max(GenericRobustM<Tint> const &robust_m,
                         MyMatrix<T> const &G,
                         std::vector<MyVector<Tint>> const& list_excluded_max,
                         std::vector<FullIneq<T,Tint>> & list_full_ineq,
                         [[maybe_unused]] std::ostream &os) {
  MyVector<Tint> v_long = robust_m.v_long();
  for (auto & v_short: list_excluded_max) {
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
    FullIneq<T,Tint> full_ineq{eIneq, {}};
    list_full_ineq.push_back(full_ineq);
  }
}

template <typename T, typename Tint>
SinglePolytope<T> get_single_p_polytope([[maybe_unused]] CVPSolver<T, Tint> const &solver,
                                        MyMatrix<T> const &FAC,
                                        [[maybe_unused]] MyVector<Tint> const &v_short,
                                        std::ostream &os) {
  std::string ansProg = "lrs";
  MyMatrix<T> EXTbig = DirectFacetComputationInequalities(FAC, ansProg, os);
  int n_ext = EXTbig.rows();
  int dim = EXTbig.cols() - 1;
  MyMatrix<T> EXT(n_ext, dim + 1);
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    if (EXTbig(i_ext, 0) <= 0) {
      std::cerr
          << "We should have v(0) > 0, because we expect to have a polytope\n";
      throw TerminalException{1};
    }
    MyVector<T> eEXT(dim);
    for (int i = 0; i < dim; i++) {
      eEXT(i) = EXTbig(i_ext, i + 1) / EXTbig(i_ext, 0);
    }
    MyMatrix<T> const &GramMat = solver.GramMat;
    MyVector<T> v_short_T = UniversalVectorConversion<T, Tint>(v_short);
    MyVector<T> diff = eEXT - v_short_T;
    T norm = EvaluationQuadForm(GramMat, diff);
    std::vector<MyVector<Tint>> ListVect =
        solver.at_most_dist_vectors(eEXT, norm);
    resultCVP<T, Tint> res_cvp{norm, MatrixFromVectorFamily(ListVect)};
    std::optional<ResultDirectEnumeration<T, Tint>> opt_rde =
        compute_and_enumerate_structures(solver.GramMat, res_cvp, eEXT, os);
    if (!opt_rde) {
      std::cerr << "ROBUST: get_p_polytope_vertices_and_test_them failed. We "
                   "should have one\n";
      throw TerminalException{1};
    }
#endif
    EXT(i_ext, 0) = 1;
    for (int i = 0; i < dim; i++) {
      EXT(i_ext, i+1) = EXTbig(i_ext, i + 1) / EXTbig(i_ext, 0);
    }
  }
  return get_single_polytope(FAC, EXT);
}

template <typename T, typename Tint>
MyMatrix<T> get_list_ineq(std::vector<FullIneq<T, Tint>> &list_full_ineq) {
  int n_ineq = list_full_ineq.size();
  if (n_ineq == 0) {
    return ZeroMatrix<T>(0,0);
  }
  int dim = list_full_ineq[0].eIneq.size();
  MyMatrix<T> M(n_ineq, dim);
  for (int i_ineq=0; i_ineq<n_ineq; i_ineq++) {
    for (int i=0; i<dim; i++) {
      M(i_ineq, i) = list_full_ineq[i].eIneq(i);
    }
  }
  return M;
}


// Compute
template <typename T, typename Tint>
T get_upper_bound_covering(CVPSolver<T, Tint> const &solver, std::ostream &os) {
  int denom = 2;
  int dim = solver.GramMat.rows();
  T upper_bound(0);
  int iter = 0;
  int max_iter = 10;
  while (true) {
    MyVector<T> eV = get_random_vector<T>(denom, dim);
    if (!IsIntegralVector(eV)) {
      ResultRobustClosest<T, Tint> rrc =
          compute_robust_closest<T, Tint>(solver, eV, os);
      T value = compute_upper_bound_rrc(solver.GramMat, rrc);
      if (upper_bound == 0) {
        // First step
        upper_bound = value;
      } else {
        if (value < upper_bound) {
          upper_bound = value;
          iter += 1;
        }
      }
    }
    if (iter == max_iter) {
      break;
    }
  }
  return upper_bound;
}

template <typename T>
T get_upper_bound_ext(MyMatrix<T> const &GramMat, MyMatrix<T> const &EXT) {
  int n = GramMat.rows();
  int n_ext = EXT.rows();
  MyVector<T> diff(n);
  T upper_bound(0);
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    for (int j_ext = i_ext + 1; j_ext < n_ext; j_ext++) {
      for (int i = 0; i < n; i++) {
        diff(i) = EXT(i_ext, i + 1) - EXT(j_ext, i + 1);
        T norm = EvaluationQuadForm(GramMat, diff);
        if (norm > upper_bound) {
          upper_bound = norm;
        }
      }
    }
  }
  return upper_bound;
}



// Find the defining inequalities of a polytope.
// It should fail and return None if the point eV is not generic enough.
// Which should lead to an increase in randomness.
template <typename T, typename Tint>
std::optional<PVoronoiPart<T, Tint>>
kernel_initial_p_polytope_part(CVPSolver<T, Tint> const &solver,
                               std::vector<MyVector<Tint>> const& list_excluded_max,
                               MyVector<T> const &eV, std::ostream &os) {
  MyMatrix<T> const &G = solver.GramMat;
  if (IsIntegralVector(eV)) {
    // Nothing can be done here
    return {};
  }
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: initial_vertex_data_test_ev eV=" << StringVectorGAP(eV)
     << "\n";
#endif
  // Working variables
  bool is_correct = true;
  PVoronoiPart<T, Tint> ppoly;
  // The lambda function.
  auto f_insert = [&](ResultDirectEnumeration<T, Tint> const &rde) -> bool {
    //
    // Building the set of inequalities from the definition.
    //
    T const &min = rde.min;
    std::vector<MyMatrix<Tint>> const &list_min_parallelepipeds =
        rde.list_min_parallelepipeds;
    std::vector<MyMatrix<Tint>> const &tot_list_parallelepipeds =
        rde.tot_list_parallelepipeds;
    std::vector<FullIneq<T,Tint>> list_full_ineq;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "--------------------------------------------------------------------"
          "---------------------------\n";
    os << "ROBUST:   initial_vertex_data_test_ev, |list_min_parallelepipeds|="
       << list_min_parallelepipeds.size()
       << " |tot_list_parallelepipeds|=" << tot_list_parallelepipeds.size()
       << " min=" << min << "\n";
#endif
    if (list_min_parallelepipeds.size() > 1) {
      // Terminate the enumeration
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by "
            "|list_min_parallelepipeds| > 1\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, pass 1\n";
#endif
    std::vector<MyVector<T>> ListIneq;
    MyMatrix<Tint> const &min_m = list_min_parallelepipeds[0];
    ExtendedGenericRobustM<T, Tint> ext_robust_m_min =
        get_generic_robust_m(min_m, G, eV, os);
    if (!ext_robust_m_min.is_correct) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by "
            "!ext_robust_m_min.is_correct\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, pass 2\n";
#endif
    if (min == 0) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by "
            "min=0\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3\n";
#endif
    GenericRobustM<Tint> const &robust_m_min = ext_robust_m_min.robust_m;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, robust_m_min, index="
       << robust_m_min.index << " M=\n";
    WriteMatrix(os, robust_m_min.M);
#endif
    insert_inner_ineqs_parallelepiped(robust_m_min, G, list_full_ineq, os);
    MyVector<Tint> v_short =
        robust_m_min.v_long(); // It is the shortest for the other structures!
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, v_short="
       << StringVectorGAP(v_short) << "\n";
#endif
    std::vector<GenericRobustM<Tint>> list_robust_m;
    insert_excluded_max(robust_m_min,
                        G,
                        list_excluded_max,
                        list_full_ineq,
                        os);

#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3, step 1\n";
    size_t i_m = 0;
#endif
    for (auto &eM : tot_list_parallelepipeds) {
      if (eM != min_m) {
#ifdef DEBUG_ENUM_P_POLYTOPES
        os << "ROBUST:   --------- " << i_m << " ------------\n";
        i_m += 1;
#endif
        ExtendedGenericRobustM<T, Tint> ext_robust_m =
            get_generic_robust_m(eM, G, eV, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
        os << "ROBUST:   initial_vertex_data_test_ev, ext_robust_m.robust_m, "
              "index="
           << ext_robust_m.robust_m.index << " M=\n";
        WriteMatrix(os, ext_robust_m.robust_m.M);
#endif
        if (!ext_robust_m.is_correct) {
#ifdef DEBUG_ENUM_P_POLYTOPES
          os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by "
                "!ext_robust_m.is_correct\n";
#endif
          is_correct = false;
          return true;
        }
        GenericRobustM<Tint> const &robust_m = ext_robust_m.robust_m;
        if (ext_robust_m.max <= min) {
          std::cerr << "ROBUST: ext_robust_m.min=" << ext_robust_m.max
                    << " min=" << min << "\n";
          std::cerr << "ROBUST: The parallelepiped has an even lower minimum\n";
          throw TerminalException{1};
        }
        insert_outer_ineqs_parallelepiped(robust_m, G, v_short, list_full_ineq, os);
        list_robust_m.push_back(robust_m);
      }
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3, step 2\n";
#endif
    //
    // Testing definition of the polytopes.
    //
    MyMatrix<T> list_ineq = get_list_ineq(list_full_ineq);
    bool test = is_full_dimensional_bounded_polytope(list_ineq, os);
    if (!test) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: initial_vertex_data_test_ev, failing by "
            "is_full_dimensional_bounded_polytope\n";
#endif
      return false;
    }
    //
    // Now doing the processing.
    //
    std::vector<int> list_irred = cdd::RedundancyReductionClarkson(list_ineq, os);
    MyMatrix<T> FAC = SelectRow(list_ineq, list_irred);
    SinglePolytope<T> sp = get_single_p_polytope(solver, FAC, v_short, os);
    std::vector<HardConvexBoundary<T>> l_hcb;
    std::vector<SoftConvexBoundary<T,Tint>> l_scb;
    for (size_t i_irred=0; i_irred<list_irred.size(); i_irred++) {
      int j_irred = list_irred[i_irred];
      FullIneq<T,Tint> const& full_ineq = list_full_ineq[j_irred];
      ConvexBoundary<T> c_bnd = get_convex_boundary(sp, i_irred);
      if (full_ineq.opt_soft) {
        SoftConvexBoundary<T,Tint> scb{0, c_bnd, list_excluded_max, *full_ineq.opt_soft};
        l_scb.push_back(scb);
      } else {
        HardConvexBoundary<T> hcb{0, c_bnd};
        l_hcb.push_back(hcb);
      }
    }
    ConvexBlock<T,Tint> c_bl{list_robust_m, sp};
    ppoly = PVoronoiPart<T,Tint>{robust_m_min, {c_bl}, l_hcb, l_scb};
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
  if (is_correct) {
    return ppoly;
  } else {
    return {};
  }
}

/*
  Get a possible vector to consider
 */
template<typename T, typename Tint>
std::optional<MyVector<Tint>> get_next_side_vector(MyMatrix<T> const& G,
                                                   MyVector<Tint> const& v_long,
                                                   std::vector<MyVector<Tint>> const& l_cand,
                                                   std::vector<MyVector<T>> const& l_vertices) {
  T delta(0);
  MyVector<T> v_long_T = UniversalVectorConversion<T,Tint>(v_long);
  auto f_max=[&](MyVector<Tint> const& cand) -> T {
    MyVector<T> cand_T = UniversalVectorConversion<T,Tint>(cand);
    T delta(0);
    for (auto & vert: l_vertices) {
      MyVector<T> diff1 = vert - v_long_T;
      MyVector<T> diff2 = vert - cand_T;
      T norm1 = EvaluationQuadForm(G, diff1);
      T norm2 = EvaluationQuadForm(G, diff2);
      if (norm2 > norm1) {
        delta = norm2 - norm1;
      }
    }
    return delta;
  };
  std::optional<MyVector<Tint>> opt;
  for (auto & cand: l_cand) {
    T del = f_max(cand);
    if (del > delta) {
      delta = del;
      opt = cand;
    }
  }
  return opt;
}



template <typename T, typename Tint>
std::optional<PVoronoi<T, Tint>>
find_p_voronoi(CVPSolver<T, Tint> const &solver, MyVector<T> const &eV, std::ostream &os) {
  MyMatrix<T> const& G = solver.GramMat;
  MyMatrix<T> G_inv = Inverse(G);
  std::optional<PVoronoiPart<T,Tint>> opt = kernel_initial_p_polytope_part<T,Tint>(solver, {}, eV, os);
  if (!opt) {
    return {};
  }
  PVoronoiPart<T, Tint> p_voronoi = *opt;
  MyVector<Tint> v_crit = p_voronoi.robust_m_min.v_long();

  auto f_process_entry=[&]() -> bool {
    size_t len = p_voronoi.l_scb.size();
    if (len == 0) {
      return true;
    }
    SoftConvexBoundary<T,Tint> scb = p_voronoi.l_scb[len-1];
    p_voronoi.l_scb.pop_back();
    std::vector<MyVector<Tint>> l_excluded_max = scb.l_excluded_max;
    MyVector<Tint> v_long = scb.robust_m.v_long();
    l_excluded_max.push_back(v_long);
    MyVector<T> eIneq = get_ineq(G, v_crit, v_long);
    MyVector<T> eIneq_op = - eIneq;
    std::vector<MyVector<T>> l_vertices = scb.cb.get_list_vertices();
    std::vector<MyVector<Tint>> l_cand = DifferenceVect(scb.robust_m.get_short_vertors(), l_excluded_max);
    std::optional<MyVector<Tint>> opt1 = get_next_side_vector(G, v_long, l_cand, l_vertices);
    if (!opt1) {
      HardConvexBoundary<T> hcb{scb.index_cb, scb.cb};
      p_voronoi.l_hcb.push_back(hcb);
      return false;
    }
    MyVector<Tint> const& v_long_new = *opt1;
    MyVector<T> eIneq_new = get_ineq(G, v_crit, v_long_new);
    std::optional<ConvexBoundary<T>> opt2 = convexboundary_halfspace_int(scb.cb, eIneq_new, os);
    if (!opt2) {
      std::cerr << "ROBUST: The opt2 should be nontrivial\n";
      throw TerminalException{1};
    }
    ConvexBoundary<T> const& cb_new = *opt2;
    MyVector<T> eIso = cb_new.get_isobarycenter();
    MyVector<T> dir = G_inv * cb_new.V;
    T shift(1);
    while(true) {
      MyVector<T> fV = eIso + shift * dir;
      std::optional<PVoronoiPart<T,Tint>> opt3 = kernel_initial_p_polytope_part<T,Tint>(solver, l_excluded_max, fV, os);
      if (!opt3) {
        continue;
      }
      PVoronoiPart<T,Tint> const& p_poly_vor_part = *opt3;
      int pos = get_position_vec_in_mat(p_poly_vor_part.l_cb[0].sp.FAC, eIneq_op);
      if (pos == -1) {
        continue;
      }
      SinglePolytope<T> const& sp = p_poly_vor_part.l_cb[0].sp;
      std::vector<ConvexBoundary<T>> l_cb = convec_boundary_minus_sp(scb.cb, sp, os);
      for (auto& cb2: l_cb) {
        SoftConvexBoundary<T,Tint> scb_new{scb.index_cb, cb2, l_excluded_max, scb.robust_m};
        p_voronoi.l_scb.push_back(scb_new);
      }
      for (auto& hcb: p_poly_vor_part.l_hcb) {
        p_voronoi.l_hcb.push_back(hcb);
      }
      p_voronoi.l_cb.push_back(p_poly_vor_part.l_cb[0]);
      int index = p_voronoi.l_cb.size() - 1;
      // That part needs to be improved, since the inserted faces could match 
      for (auto& scb: p_poly_vor_part.l_scb) {
        if (scb.cb.V != eIneq_op) {
          SoftConvexBoundary<T,Tint> scb_new = scb;
          scb_new.index_cb = index;
          p_voronoi.l_scb.push_back(scb_new);
        }
      }
    }
    return false;
  };


  while(true) {
    bool test = f_process_entry();
    if (test) {
      break;
    }
  }
  std::cerr << "ROBUST: find_p_voronoi should never reach that stage\n";
  throw TerminalException{1};
}





template <typename T, typename Tint>
PVoronoi<T, Tint>
initial_p_polytope(CVPSolver<T, Tint> const &solver, std::ostream &os) {
  int dim = solver.GramMat.rows();
  int denom = 2;
  while (true) {
    MyVector<T> eV = get_random_vector<T>(denom, dim);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: initial_vertex_data, before initial_vertex_data_test_ev, eV="
       << StringVectorGAP(eV) << " denom=" << denom << "\n";
#endif
    std::optional<PVoronoi<T, Tint>> opt = find_p_voronoi(solver, eV, os);
#ifdef DEBUG_ENUM_P_POLYTOPES_DISABLE
    std::cerr << "ROBUST: Stopping the execution\n";
    throw TerminalException{1};
#endif
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: initial_vertex_data, after initial_vertex_data_test_ev\n";
#endif
    if (opt) {
      return *opt;
    }
    denom += 1;
  }
}





template <typename T, typename Tint, typename Tgroup>
std::vector<PVoronoi<T, Tint>>
find_list_adjacent_p_voronoi(DataLattice<T, Tint, Tgroup> &eData,
                             PVoronoi<T, Tint> const &pvd) {
  std::ostream &os = eData.rddo.os;
  CVPSolver<T, Tint> const &solver = eData.solver;
  int dim = eData.solver.GramMat.rows();
  auto get_adj_p_polytope = [&](MyVector<T> const &TestFAC,
                                MyVector<T> const &x)
      -> std::optional<PVoronoi<T, Tint>> {
    MyVector<T> x_red(dim);
    for (int i = 0; i < dim; i++) {
      x_red(i) = x(i + 1);
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: get_adj_p_polytope x_red=" << StringVector(x_red) << "\n";
#endif
    std::optional<PVoronoi<T, Tint>> opt = find_p_voronoi(solver, x_red, os);
    if (!opt) {
      return {};
    }
    PVoronoi<T, Tint> const &ppoly = *opt;
    size_t m_facet = ppoly.sp.FAC.rows();
    int m_ext = ppoly.EXT.rows();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: get_adj_p_polytope m_facet=" << m_facet << " m_ext=" << m_ext
       << "\n";
#endif
    for (size_t j_facet = 0; j_facet < m_facet; j_facet++) {
      if (GetMatrixRow(ppoly.sp.FAC, j_facet) == TestFAC) {
        return ppoly;
      }
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: get_adj_p_polytope, Failed to find a correct adjacent "
          "cone\n";
#endif
    return {};
  };
  auto get_adj = [&](Face const &f, MyVector<T> const &eFAC,
                     MyVector<T> const &eIso) -> PVoronoi<T, Tint> {
    std::unordered_set<MyVector<T>> set;
    int n_ext = pvd.EXT.rows();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      if (f[i_ext] == 1) {
        MyVector<T> eEXT = GetMatrixRow(pvd.EXT, i_ext);
        set.insert(eEXT);
      }
    }
    MyVector<T> TestFAC = -eFAC;
    MyVector<T> delta_x = eIso - pvd.eIso;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: get_adj, eIso=" << StringVector(eIso)
       << " pvd.eIso=" << StringVector(pvd.eIso) << " delta_x=" << delta_x
       << "\n";
#endif
    T factor(1);
#ifdef DEBUG_ENUM_P_POLYTOPES
    size_t n_iter = 0;
#endif
    while (true) {
      MyVector<T> x = eIso + factor * delta_x;
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: get_adj, n_iter=" << n_iter << " factor=" << factor
         << "\n";
      n_iter += 1;
#endif
      std::optional<PVoronoi<T, Tint>> opt =
          get_adj_p_polytope(TestFAC, x);
      if (opt) {
        return *opt;
      }
      factor = factor / 2;
    }
  };
  std::vector<PVoronoi<T, Tint>> l_adj;
  size_t n_facet = pvd.sp.FAC.rows();
  for (size_t i_facet = 0; i_facet < n_facet; i_facet++) {
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: find_list_adjacent_p_voronoi i_facet=" << i_facet << " / "
       << n_facet << "\n";
#endif
    Face facet = pvd.sp.facets[i_facet];
    MyVector<T> eFAC = GetMatrixRow(pvd.sp.FAC, i_facet);
    MyVector<T> eIso = get_interior_facet_pt(pvd.sp, i_facet);
    PVoronoi<T, Tint> eAdj = get_adj(facet, eFAC, eIso);
    l_adj.push_back(eAdj);
  }
  return l_adj;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<PVoronoi<T, Tint>>
compute_all_p_polytopes(DataLattice<T, Tint, Tgroup> &eData) {
  std::ostream &os = eData.rddo.os;
  CVPSolver<T, Tint> const &solver = eData.solver;
  std::vector<PVoronoi<T, Tint>> l_ppoly;
  PVoronoi<T, Tint> ppoly = initial_p_polytope(solver, os);
  l_ppoly.push_back(ppoly);
  auto f_insert = [&](PVoronoi<T, Tint> const &f_ppoly) -> void {
    for (auto &e_ppoly : l_ppoly) {
      std::optional<MyMatrix<T>> opt_equiv =
          Polytope_TestEquivalence<T, Tint, Tgroup>(eData, e_ppoly.EXT,
                                                    f_ppoly.EXT);
      if (opt_equiv.has_value()) {
        return;
      }
    }
    l_ppoly.push_back(f_ppoly);
  };
  size_t start = 0;
  while (true) {
    size_t len = l_ppoly.size();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: compute_all_p_polytopes, start=" << start << " len=" << len
       << "\n";
#endif
    for (size_t i_ppoly = start; i_ppoly < len; i_ppoly++) {
      std::vector<PVoronoi<T, Tint>> l_adj =
          find_list_adjacent_p_voronoi(eData, l_ppoly[i_ppoly]);
      for (auto &eAdj : l_adj) {
        f_insert(eAdj);
      }
    }
    if (len == l_ppoly.size()) {
      break;
    }
    start = len;
  }
  return l_ppoly;
}

template <typename T, typename Tint, typename Tgroup>
T compute_square_robust_covering_radius(DataLattice<T, Tint, Tgroup> &eData) {
  std::vector<PVoronoi<T, Tint>> l_ppoly =
      compute_all_p_polytopes(eData);
  T max_sqr_radius(0);
  MyMatrix<T> const &GramMat = eData.solver.GramMat;
  int dim = GramMat.rows();
  MyVector<T> diff(dim);
  for (auto &ppoly : l_ppoly) {
    MyVector<Tint> v_long = ppoly.robust_m_min.v_long();
    MyVector<T> v_long_T = UniversalVectorConversion<T, Tint>(v_long);
    int n_ext = ppoly.EXT.rows();
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      for (int i = 0; i < dim; i++) {
        diff(i) = v_long_T(i) - ppoly.EXT(i_ext, i + 1);
      }
      T norm = EvaluationQuadForm(GramMat, diff);
      if (norm > max_sqr_radius) {
        max_sqr_radius = norm;
      }
    }
  }
  return max_sqr_radius;
}

// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
