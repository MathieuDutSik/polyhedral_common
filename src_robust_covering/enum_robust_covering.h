// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "generalized_polytopes.h"
#include "parall_search.h"
#include "LatticeDelaunay.h"
#include "subgroup_algorithms.h"
#include "POLY_RedundancyElimination.h"
#include "COMB_Combinatorics.h"
#include "boost_serialization.h"
// clang-format on

/*
  The motivation is explained in
  "New upper bound for lattice covering by spheres"
  https://arxiv.org/pdf/2508.06446
  ---
  The quantity of interest in that case is
  log_2(Theta_n) / n
  ---
  In dimension 2, we have for A2
  Theta_2 = 8 pi / (3 Sqrt(3)).
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
#define DEBUG_ROBUST_VERTEX_ENUM
#define DEBUG_GET_INEQ_P_POLYTOPES
#define DEBUG_P_VORONOI_STABILIZER
#endif

#ifdef DEBUG
#define PRINT_ENUM_P_POLYTOPES
#endif



#ifdef SANITY_CHECK
#define SANITY_CHECK_ENUM_P_POLYTOPES
#define SANITY_CHECK_ROBUST_VERTEX_ENUM
#endif

#ifdef DISABLE_DEBUG_ENUM_P_POLYTOPES
#undef DEBUG_ENUM_P_POLYTOPES
#endif





// A family of vectors with the index which is the farthest.
template <typename Tint> struct GenericRobustM {
  int index;
  MyMatrix<Tint> M;
  MyVector<Tint> v_long() const { return GetMatrixRow(M, index); }
  std::vector<MyVector<Tint>> get_short_vectors() const {
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

template <typename Tint>
std::ostream &operator<<(std::ostream &os, GenericRobustM<Tint> const &grm) {
  os << "GenericRobustM(index=" << grm.index << " M=\n";
  WriteMatrix(os, grm.M);
  os << ")";
  return os;
}

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out, GenericRobustM<Tint> const &grm) {
  os_out << "rec(index:=" << grm.index << ", M:=";
  WriteMatrixGAP(os_out, grm.M);
  os_out << ")";
}

template <typename T, typename Tint> struct ExtendedGenericRobustM {
  T max;
  bool is_correct;
  GenericRobustM<Tint> robust_m;
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, ExtendedGenericRobustM<T, Tint> const &egrm) {
  os << "ExtendedGenericRobustM(max=" << egrm.max << " is_correct=" << egrm.is_correct << " robust_m=" << egrm.robust_m << ")";
  return os;
}


template<typename T, typename Tint>
struct MapFullIneq {
  // Map from the inequalities to the parallelepipeds
  std::map<MyVector<T>, std::vector<GenericRobustM<Tint>>> m;
  // If std::vector<_> is empty, then the inequality is a hard one.
  // If std::vector<_>, then it is a soft one and that indicates what
  void insert_hard_ineq(MyVector<T> const& eIneq) {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    if (m.contains(eIneq)) {
      std::cerr << "ROBUST: A hard inequality should be present only once. eIneq=" << StringVectorGAP(eIneq) << "\n";
      throw TerminalException{1};
    }
#endif
    m[eIneq] = {};
  }
  void insert_soft_ineq(MyVector<T> const& eIneq, GenericRobustM<Tint> const& grm) {
    m[eIneq].push_back(grm);
  }
  MyMatrix<T> get_list_ineq([[maybe_unused]] std::ostream &os) {
    int n_ineq = m.size();
    if (n_ineq == 0) {
      return ZeroMatrix<T>(0,0);
    }
    int dim = m.begin()->first.size();
    MyMatrix<T> M(n_ineq, dim);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: get_list_ineq n_ineq=" << n_ineq << " dim=" << dim << "\n";
#endif
    int i_ineq = 0;
    for (auto & kv: m) {
      for (int i=0; i<dim; i++) {
        M(i_ineq, i) = kv.first(i);
      }
      i_ineq += 1;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES_DISABLE
    os << "ROBUST: get_list_ineq M=\n";
    WriteMatrix(os, M);
#endif
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES_DISABLE
    bool test = no_duplicated_scalar_multiple(M);
    if (!test) {
      std::cerr << "ROBUST: The matrix M has duplication\n";
      throw TerminalException{1};
    }
#endif
    return M;
  }
  std::vector<GenericRobustM<Tint>> get_list_paralls(MyVector<T> const& eIneq) const {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    if (!m.contains(eIneq)) {
      std::cerr << "ROBUST: get_generic_robust_paralls, The inequality is missing\n";
      throw TerminalException{1};
    }
#endif
    return m.at(eIneq);
  }
};

// The convex block of the construction
template <typename T, typename Tint>
struct ConvexBlock {
  std::vector<GenericRobustM<Tint>> list_robust_m;
  SinglePolytope<T> sp;
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, ConvexBlock<T, Tint> const &cb) {
  os << "ConvexBlock(|list_robust_m|=" << cb.list_robust_m.size() << " FAC=\n";
  WriteMatrix(os, cb.sp.FAC);
  os << "EXT=\n";
  WriteMatrix(os, cb.sp.EXT);
  os << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, ConvexBlock<T, Tint> const &cb) {
  os_out << "rec(list_robust_m:=[";
  bool IsFirst = true;
  for (auto &grm : cb.list_robust_m) {
    if (!IsFirst) {
      os_out << ",";
    }
    IsFirst = false;
    WriteEntryGAP(os_out, grm);
  }
  os_out << "], sp:=";
  WriteEntryGAP(os_out, cb.sp);
  os_out << ")";
}

template<typename T>
struct HardConvexBoundary {
  int index_cb; // The corresponding face;
  ConvexBoundary<T> sp;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, HardConvexBoundary<T> const &hcb) {
  os << "HardConvexBoundary(index_cb=" << hcb.index_cb << " V=" << StringVectorGAP(hcb.sp.V) << ")";
  return os;
}

template <typename T>
void WriteEntryGAP(std::ostream &os_out, HardConvexBoundary<T> const &hcb) {
  os_out << "rec(index_cb:=" << hcb.index_cb << ", sp:=";
  WriteEntryGAP(os_out, hcb.sp);
  os_out << ")";
}

template<typename T, typename Tint>
struct SoftConvexBoundary {
  int index_cb; // The corresponding face;
  ConvexBoundary<T> cb;
  std::vector<MyVector<Tint>> l_excluded_max; // The excluded vectors. Cannot use the ConvexBlock since they can vary.
  std::vector<GenericRobustM<Tint>> l_robust_m;
  MyVector<Tint> v_long() const {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    if (l_robust_m.empty()) {
      std::cerr << "ROBUST: l_robust_m should be non-empty for getting the v_long\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> v_long = l_robust_m[0].v_long();
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
    size_t len = l_robust_m.size();
    for (size_t i=1; i<len; i++) {
      MyVector<Tint> v = l_robust_m[i].v_long();
      if (v != v_long) {
        std::cerr << "ROBUST: Incoherent v_long in the structure\n";
        throw TerminalException{1};
      }
    }
#endif
    return v_long;
  }
};

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SoftConvexBoundary<T, Tint> const &scb) {
  os << "SoftConvexBoundary(index_cb=" << scb.index_cb
     << " V=" << StringVectorGAP(scb.cb.V)
     << " |l_excluded_max|=" << scb.l_excluded_max.size()
     << " |l_robust_m|=" << scb.l_robust_m.size() << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, SoftConvexBoundary<T, Tint> const &scb) {
  os_out << "rec(index_cb:=" << scb.index_cb << ", cb:=";
  WriteEntryGAP(os_out, scb.cb);
  os_out << ", l_excluded_max:=[";
  for (size_t i = 0; i < scb.l_excluded_max.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    os_out << StringVectorGAP(scb.l_excluded_max[i]);
  }
  os_out << "], l_robust_m:=[";
  for (size_t i = 0; i < scb.l_robust_m.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, scb.l_robust_m[i]);
  }
  os_out << "])";
}

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

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, PVoronoi<T, Tint> const &pv) {
  os << "PVoronoi(\n  robust_m_min=" << pv.robust_m_min << "\n";
  os << "  |l_cb|=" << pv.l_cb.size() << "\n";
  for (size_t i = 0; i < pv.l_cb.size(); i++) {
    os << "  l_cb[" << i << "]=" << pv.l_cb[i] << "\n";
  }
  os << "  |l_hcb|=" << pv.l_hcb.size() << "\n";
  for (size_t i = 0; i < pv.l_hcb.size(); i++) {
    os << "  l_hcb[" << i << "]=" << pv.l_hcb[i] << "\n";
  }
  os << "  EXT=\n";
  WriteMatrix(os, pv.EXT);
  os << ")";
  return os;
}

template <typename T, typename Tint>
void WriteEntryGAP(std::ostream &os_out, PVoronoi<T, Tint> const &pv) {
  os_out << "rec(robust_m_min:=";
  WriteEntryGAP(os_out, pv.robust_m_min);
  os_out << ", l_cb:=[";
  for (size_t i = 0; i < pv.l_cb.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, pv.l_cb[i]);
  }
  os_out << "], l_hcb:=[";
  for (size_t i = 0; i < pv.l_hcb.size(); i++) {
    if (i > 0) {
      os_out << ",";
    }
    WriteEntryGAP(os_out, pv.l_hcb[i]);
  }
  os_out << "], gp:=";
  WriteEntryGAP(os_out, pv.gp);
  os_out << ", EXT:=";
  WriteMatrixGAP(os_out, pv.EXT);
  os_out << ")";
}

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
std::ostream &operator<<(std::ostream &os, PVoronoiPart<T, Tint> const &pvp) {
  os << "PVoronoiPart(\n  robust_m_min=" << pvp.robust_m_min << "\n";
  os << "  |l_cb|=" << pvp.l_cb.size() << "\n";
  for (size_t i = 0; i < pvp.l_cb.size(); i++) {
    os << "  l_cb[" << i << "]=" << pvp.l_cb[i] << "\n";
  }
  os << "  |l_hcb|=" << pvp.l_hcb.size() << "\n";
  for (size_t i = 0; i < pvp.l_hcb.size(); i++) {
    os << "  l_hcb[" << i << "]=" << pvp.l_hcb[i] << "\n";
  }
  os << "  |l_scb|=" << pvp.l_scb.size() << "\n";
  for (size_t i = 0; i < pvp.l_scb.size(); i++) {
    os << "  l_scb[" << i << "]=" << pvp.l_scb[i] << "\n";
  }
  os << ")";
  return os;
}

template<typename T>
T min_pairwise_norm(MyMatrix<T> const& EXT, MyMatrix<T> const& G) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  int dim = n_col - 1;
  T min_norm(0);
  bool is_first=true;
  MyVector<T> diff(dim);
  for (int i_row=0; i_row<n_row; i_row++) {
    for (int j_row=i_row+1; j_row<n_row; j_row++) {
      for (int i=0; i<dim; i++) {
        diff(i) = EXT(i_row, i+1) - EXT(j_row, i+1);
      }
      T norm = EvaluationQuadForm(G, diff);
      if (is_first) {
        min_norm = norm;
        is_first = false;
      } else {
        if (norm < min_norm) {
          min_norm = norm;
        }
      }
    }
  }
  return min_norm;
}

template <typename T, typename Tint>
PVoronoi<T,Tint> convert_p_voronoi_part(PVoronoiPart<T,Tint> const& pvp, std::ostream& os) {
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
  if (!pvp.l_scb.empty()) {
    std::cerr << "ROBUST: We still have soft boundaries\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 1\n";
#endif
  std::vector<SinglePolytope<T>> polytopes;
  for (auto & cb: pvp.l_cb) {
    polytopes.push_back(cb.sp);
  }
  GeneralizedPolytope<T> gp{polytopes};
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 2\n";
#endif
  BoundaryGeneralizedPolytope<T> bnd = find_generalized_polytope_boundary(gp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 3\n";
#endif
  std::vector<MyVector<T>> l_vert = get_vertices(gp, bnd, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 4\n";
#endif
  MyMatrix<T> EXT = MatrixFromVectorFamily(l_vert);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: convert_p_voronoi_part, step 5\n";
#endif
  return {pvp.robust_m_min,
          pvp.l_cb,
          pvp.l_hcb,
          gp,
          EXT};
}






template <typename T, typename Tint>
ExtendedGenericRobustM<T, Tint>
get_generic_robust_m(MyMatrix<Tint> const &M, MyMatrix<T> const &G,
                     MyVector<T> const &eV, [[maybe_unused]] std::ostream &os) {
  T max(0);
  int best_index = 0;
  int n_ineq = M.rows();
  size_t n_att = 0;
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
  if (G.rows() != eV.size()) {
    std::cerr << "ROBUST: G and eV have inconsistent size\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
  os << "ROBUST:   ggrm, eV=" << StringVectorGAP(eV) << "\n";
#endif
  for (int index = 0; index < n_ineq; index++) {
    MyVector<Tint> fV = GetMatrixRow(M, index);
    MyVector<T> diff = UniversalVectorConversion<T, Tint>(fV) - eV;
    T norm = EvaluationQuadForm(G, diff);
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
    double norm_d = UniversalScalarConversion<double, T>(norm);
    os << "ROBUST:   ggrm, index=" << index
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
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
  double max_d = UniversalScalarConversion<double, T>(max);
  os << "ROBUST:   ggrm, best_index=" << best_index
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
                                       MapFullIneq<T,Tint> & m_full_ineq,
                                       [[maybe_unused]] std::ostream &os) {
  int n_row = robust_m.M.rows();
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
  os << "ROBUST:   iiip, M=\n";
  WriteMatrix(os, robust_m.M);
#endif
  MyVector<Tint> v_long = robust_m.v_long();
  for (int i = 0; i < n_row; i++) {
    if (i != robust_m.index) {
      MyVector<Tint> v_short = GetMatrixRow(robust_m.M, i);
      MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
      os << "ROBUST:   iiip, eIneq="
         << StringVector(eIneq) << " i=" << i << " index=" << robust_m.index
         << "\n";
#endif
      m_full_ineq.insert_hard_ineq(eIneq);
    }
  }
}

// The farthest vector of the robust structure has to be farther from the
// v_short (which is the farthest of the best one)
template <typename T, typename Tint>
void insert_outer_ineqs_parallelepiped(GenericRobustM<Tint> const &robust_m,
                                       MyMatrix<T> const &G,
                                       MyVector<Tint> const &v_short,
                                       MapFullIneq<T, Tint> &m_full_ineq,
                                       [[maybe_unused]] std::ostream &os) {
  MyVector<Tint> v_long = robust_m.v_long();
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
  os << "ROBUST:   ioip, v_short=" << StringVector(v_short) << "\n";
  os << "ROBUST:   ioip, v_long=" << StringVector(v_long) << "\n";
#endif
  if (v_short != v_long) {
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
    os << "ROBUST:   ioip, eIneq=" << StringVector(eIneq) << "\n";
#endif
    m_full_ineq.insert_soft_ineq(eIneq, robust_m);
  } else {
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
    os << "ROBUST:   ioip, v_short=" << StringVector(v_short) << "\n";
#endif
  }
}

template <typename T, typename Tint>
void insert_excluded_max(GenericRobustM<Tint> const &robust_m,
                         MyMatrix<T> const &G,
                         std::vector<MyVector<Tint>> const& l_excluded_max,
                         MapFullIneq<T,Tint> & m_full_ineq,
                         [[maybe_unused]] std::ostream &os) {
  MyVector<Tint> v_long = robust_m.v_long();
  for (auto & v_short: l_excluded_max) {
#ifdef DEBUG_GET_INEQ_P_POLYTOPES
    os << "ROBUST:   iem, v_short=" << StringVector(v_short) << "\n";
#endif
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
    m_full_ineq.insert_hard_ineq(eIneq);
  }
}

template <typename T, typename Tint>
SinglePolytope<T> get_single_p_polytope([[maybe_unused]] CVPSolver<T, Tint> const &solver,
                                        MyMatrix<T> const &FAC,
                                        [[maybe_unused]] MyVector<Tint> const &v_short,
                                        std::ostream &os) {
  MyMatrix<T> EXTbig = DirectDualDescription_mat(FAC, os);
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

// When we have a parallelotope, we can compute an upper bound
// on the robust covering density.
// Here we iterate 10 times to find better and better upper bounds.
template <typename T, typename Tint>
T get_upper_bound_covering(CVPSolver<T, Tint> const &solver, std::ostream &os) {
  int denom = 10000;
  int dim = solver.GramMat.rows();
  T upper_bound(0);
  int max_iter = 10;
  for (int iter=0; iter<max_iter; iter++) {
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
        }
      }
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
      }
      T norm = EvaluationQuadForm(GramMat, diff);
      if (norm > upper_bound) {
        upper_bound = norm;
      }
    }
  }
  return upper_bound;
}


template <typename T, typename Tint>
bool are_vertices_correct(CVPSolver<T, Tint> const &solver,
                          MyVector<Tint> const& v_short,
                          SinglePolytope<T> const& sp,
                          std::ostream &os) {
  MyMatrix<T> const& G = solver.GramMat;
  int dim = G.rows();
  MyVector<T> v_short_T = UniversalVectorConversion<T,Tint>(v_short);
  int n_ext=sp.EXT.rows();
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: are_vertices_correct, n_ext=" << n_ext << "\n";
#endif
  MyVector<T> diff(dim);
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<T> eV = GetMatrixRow(sp.EXT, i_ext);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: are_vertices_correct, i_ext=" << i_ext << " eV=" << StringVectorGAP(eV) << "\n";
#endif
    ResultRobustClosest<T, Tint> rrc =
      compute_robust_closest<T, Tint>(solver, eV, os);
    for (int i=0; i<dim; i++) {
      diff(i) = sp.EXT(i_ext, i+1) - v_short_T(i);
    }
    T norm = EvaluationQuadForm(G, diff);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: are_vertices_correct, i_ext=" << i_ext << " norm=" << norm << " robust_minimum=" << rrc.robust_minimum << "\n";
#endif
    if (norm > rrc.robust_minimum) {
      return false;
    }
  }
  return true;
}



// Find the defining inequalities of a polytope.
// It should fail and return None if the point eV is not generic enough.
// Which should lead to an increase in randomness.
template <typename T, typename Tint>
std::optional<PVoronoiPart<T, Tint>>
kernel_l2_p_polytope_part(CVPSolver<T, Tint> const &solver,
                               std::vector<MyVector<Tint>> const& l_excluded_max,
                               MyVector<T> const &eV, std::ostream &os) {
  MyMatrix<T> const &G = solver.GramMat;
  if (IsIntegralVector(eV)) {
    // Nothing can be done here
    return {};
  }
  int dim = G.rows();
  MyVector<T> eV_red(dim);
  for (int i=0; i<dim; i++) {
    eV_red(i) = eV(i + 1);
  }
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: kippp eV=" << StringVectorGAP(eV) << "\n";
#endif
  // The is_correct variable indicate whether the initial point eV is correct.
  // In the course of the enumeration, we encounter some errors. Some merely
  // indicate that the enumeration is incomplete. That is we keep is_correct=true.
  // But others indicate that no matter the level of the enumeration, the
  // point in incorrect.
  bool is_correct = true;
  PVoronoiPart<T, Tint> ppoly;
  // The lambda function.
  auto f_insert = [&](ResultDirectEnumeration<T, Tint> const &rde, [[maybe_unused]] T const& TheNorm) -> bool {
    //
    // Building the set of inequalities from the definition.
    //
    T const &min = rde.min;
    std::vector<MyMatrix<Tint>> const &list_min_parallelepipeds =
        rde.list_min_parallelepipeds;
    std::vector<MyMatrix<Tint>> const &tot_list_parallelepipeds =
        rde.tot_list_parallelepipeds;
    MapFullIneq<T,Tint> m_full_ineq;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "-----------------------------------------------------------\n";
    os << "ROBUST:   kippp, |list_min_parallelepipeds|="
       << list_min_parallelepipeds.size()
       << " |tot_list_parallelepipeds|=" << tot_list_parallelepipeds.size()
       << " min=" << min << "\n";
#endif
    if (list_min_parallelepipeds.size() > 1) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   kippp, is_correct=false by n_list_min_parallelepipeds > 1\n";
#endif
      // We mark is_correct=false since there is no chance that this degeneracy will be
      // resolved by further enumeration. The degeneracy is there and will remain.
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 1\n";
#endif
    MyMatrix<Tint> const &min_m = list_min_parallelepipeds[0];
    T upper_bound = compute_upper_bound_mat(G, min_m);
    ExtendedGenericRobustM<T, Tint> ext_robust_m_min =
        get_generic_robust_m(min_m, G, eV_red, os);
    if (!ext_robust_m_min.is_correct) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   kippp, is_correct=false by "
            "!ext_robust_m_min.is_correct\n";
#endif
      // We mark is_correct=false since the minimal has at least two extremal
      // vectors. This problem will remain no matter how many other parallelepiped
      // we have
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 2\n";
#endif
    if (min == 0) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST:   kippp, is_correct=false by min=0\n";
#endif
      // Wrong value of "min". This will remain no matter what. So we mark
      // with is_correct=false
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 3\n";
#endif
    GenericRobustM<Tint> const &robust_m_min = ext_robust_m_min.robust_m;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, robust_m_min, index="
       << robust_m_min.index << " M=\n";
    WriteMatrix(os, robust_m_min.M);
#endif
    insert_inner_ineqs_parallelepiped(robust_m_min, G, m_full_ineq, os);
    MyVector<Tint> v_short =
        robust_m_min.v_long(); // It is the shortest for the other structures!
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, v_short="
       << StringVectorGAP(v_short) << "\n";
    os << "ROBUST:   kippp, |l_excluded_max|=" << l_excluded_max.size() << "\n";
#endif
    std::vector<GenericRobustM<Tint>> list_robust_m;
    insert_excluded_max(robust_m_min,
                        G,
                        l_excluded_max,
                        m_full_ineq,
                        os);

#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 3, step 1\n";
#endif
#ifdef DEBUG_ENUM_P_POLYTOPES_DISABLE
    size_t i_m = 0;
#endif
    for (auto &eM : tot_list_parallelepipeds) {
      if (eM != min_m) {
#ifdef DEBUG_ENUM_P_POLYTOPES_DISABLE
        os << "ROBUST:   --------- "
           << i_m
           << "/"
           << tot_list_parallelepipeds.size()
           << " ------------\n";
        i_m += 1;
#endif
        T val = compute_upper_bound_mat(G, eM);
        if (val < upper_bound) {
          upper_bound = val;
        }
        ExtendedGenericRobustM<T, Tint> ext_robust_m =
            get_generic_robust_m(eM, G, eV_red, os);
#ifdef DEBUG_ENUM_P_POLYTOPES_DISABLE
        os << "ROBUST:   kippp, ext_robust_m.robust_m, index="
           << ext_robust_m.robust_m.index << " M=\n";
        WriteMatrix(os, ext_robust_m.robust_m.M);
#endif
        if (!ext_robust_m.is_correct) {
#ifdef DEBUG_ENUM_P_POLYTOPES
          os << "ROBUST:   kippp, is_correct=false by "
                "!ext_robust_m.is_correct\n";
#endif
          // Wrong parallelepipeds occur with ambiguous input. This will
          // remain if we have more parallelepipeds. So is_correct=false
          is_correct = false;
          return true;
        }
        GenericRobustM<Tint> const &robust_m = ext_robust_m.robust_m;
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
        if (ext_robust_m.max <= min) {
          std::cerr << "ROBUST: ext_robust_m.min=" << ext_robust_m.max
                    << " min=" << min << "\n";
          std::cerr << "ROBUST: The parallelepiped has an even lower minimum\n";
          throw TerminalException{1};
        }
#endif
        insert_outer_ineqs_parallelepiped(robust_m, G, v_short, m_full_ineq, os);
        list_robust_m.push_back(robust_m);
      }
    }
    if (upper_bound > TheNorm) {
      // It is unfortunate that we have to compute that upper bound.
      // Which is frankly not that good. We need a better upper bound.
      //
      // We are below the upper bound, which indicates that we could have some
      // additional parallelepipeds missing. So, we do not set is_correct=false.
      return false;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 3, step 2\n";
#endif
    //
    // Testing definition of the polytopes.
    //
    MyMatrix<T> list_ineq = m_full_ineq.get_list_ineq(os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST:   kippp, pass 3, step 2(B)\n";
    {
      bool test = IsFullDimensional(list_ineq, os);
      os << "ROBUST:   kippp, is_full_dim, test=" << test << "\n";
      if (test) {
        std::vector<int> list_irred = get_non_redundant_indices(list_ineq, os);
        os << "ROBUST:   kippp, |list_irred|=" << list_irred.size() << "\n";
      }
    }
#endif
    bool test = is_full_dimensional_bounded_polytope(list_ineq, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, is_full_dim_bounded_p, test=" << test << "\n";
#endif
    if (!test) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: kippp, failing by "
            "is_full_dimensional_bounded_polytope\n";
#endif
      return false;
    }
    //
    // Now doing the processing.
    //
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 1\n";
#endif
    std::vector<int> list_irred = get_non_redundant_indices(list_ineq, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 2\n";
#endif
    MyMatrix<T> FAC = SelectRow(list_ineq, list_irred);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 3\n";
#endif
    SinglePolytope<T> sp = get_single_p_polytope(solver, FAC, v_short, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 4\n";
#endif
    bool test_vert = are_vertices_correct(solver, v_short, sp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 5, test_vert=" << test_vert << "\n";
#endif
    if (!test_vert) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: kippp, final, is_correct=false because of test_vert\n";
#endif
      // We do NOT mark is_correct=false since the problem could be addressed
      // by further enumeration of the parallelepipeds.
      return false;
    }
    std::vector<HardConvexBoundary<T>> l_hcb;
    std::vector<SoftConvexBoundary<T,Tint>> l_scb;
    for (size_t i_irred=0; i_irred<list_irred.size(); i_irred++) {
      ConvexBoundary<T> c_bnd = get_convex_boundary(sp, i_irred, os);
      MyVector<T> const& V = c_bnd.V;
      std::vector<GenericRobustM<Tint>> l_grm = m_full_ineq.get_list_paralls(V);
      if (!l_grm.empty()) {
        SoftConvexBoundary<T,Tint> scb{0, c_bnd, l_excluded_max, l_grm};
        l_scb.push_back(scb);
      } else {
        HardConvexBoundary<T> hcb{0, c_bnd};
        l_hcb.push_back(hcb);
      }
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 5\n";
#endif
    ConvexBlock<T,Tint> c_bl{list_robust_m, sp};
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 6\n";
#endif
    ppoly = PVoronoiPart<T,Tint>{robust_m_min, {c_bl}, l_hcb, l_scb};
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: kippp, final, step 7 TheNorm=" << TheNorm << "\n";
    os << "ROBUST: kippp, final, step 7 upper_bound=" << upper_bound << "\n";
#endif
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: kippp, is_correct=" << is_correct << "\n";
#endif
  if (is_correct) {
    return ppoly;
  } else {
    return {};
  }
}

template <typename T, typename Tint>
std::optional<PVoronoiPart<T, Tint>>
kernel_l1_p_polytope_part(CVPSolver<T, Tint> const &solver,
                               std::vector<MyVector<Tint>> const& l_excluded_max,
                               MyVector<T> const &eV, std::ostream &os) {
  std::optional<PVoronoiPart<T, Tint>> opt =
      kernel_l2_p_polytope_part(solver, l_excluded_max, eV, os);
#ifdef SANITY_CHECK_EXTENSIVE_ENUM_P_POLYTOPES
  if (opt) {
    PVoronoiPart<T, Tint> const& pvp1 = *opt;
    MyMatrix<T> const& GramMat = solver.GramMat;
    int dim = GramMat.rows();
    MyVector<T> eV_red(dim);
    for (int i=0; i<dim; i++) {
      eV_red(i) = eV(i + 1);
    }
    // Find the vertex farthest from v_long
    MyVector<Tint> v_long = pvp1.robust_m_min.v_long();
    MyMatrix<Tint> M1 = reorder_matrix(pvp1.robust_m_min.M);
    MyVector<T> v_long_T = UniversalVectorConversion<T, Tint>(v_long);
    std::cerr << "ROBUST: kernel_l1_p_polytope_part v_long=" << StringVectorGAP(v_long) << "\n";
    MyMatrix<T> EXT1 = reorder_matrix(pvp1.l_cb[0].sp.EXT);
    int n_row = EXT1.rows();
    T max_norm(0);
    MyVector<T> diff(dim);
    for (int i_row = 0; i_row < n_row; i_row++) {
      for (int i = 0; i < dim; i++) {
        diff(i) = v_long_T(i) - EXT1(i_row, i + 1);
      }
      T norm = EvaluationQuadForm(GramMat, diff);
      if (norm > max_norm) {
        max_norm = norm;
      }
    }
    std::unordered_set<MyMatrix<Tint>> set;
    size_t n_robust = pvp1.l_cb[0].list_robust_m.size();
    for (size_t i_robust=0; i_robust<n_robust; i_robust++) {
      MyMatrix<Tint> M = reorder_matrix(pvp1.l_cb[0].list_robust_m[i_robust].M);
      set.insert(M);
    }
    // Generate random points inside the polytope and check consistency
    int n_test = 20;
    int N = 10;
    for (int i_test = 0; i_test < n_test; i_test++) {
      std::cerr << "ROBUST: kernel_l1_p_polytope_part i_test=" << i_test << "/" << n_test << "\n";
      MyVector<T> fV = random_interior_pt(EXT1, N, os);
      MyVector<T> fV_red(dim);
      for (int i=0; i<dim; i++) {
        fV_red(i) = fV(i + 1);
      }
      // Check with compute_robust_closest
      ResultRobustClosest<T, Tint> rrc =
          compute_robust_closest<T, Tint>(solver, fV, os);
      if (rrc.robust_minimum > max_norm) {
        std::cerr << "ROBUST: rrc failed to work fV=" << StringVectorGAP(fV) << "\n";
        std::cerr << "ROBUST: max_norm=" << max_norm << " robust_minimum=" << rrc.robust_minimum << "\n";
        throw TerminalException{1};
      }
      // The polytope kernel_l2_p_polytope_part(solver, l_excluded_max, fV, os);
      // will not necessarily be the same depending on the chosen point.
      // So, we cannot make that check. And of course it could return None.
      if (rrc.list_parallelepipeds.size() > 1) {
        size_t n_parall = rrc.list_parallelepipeds.size();
        std::cerr << "ROBUST: |rrc.list_parallelepipeds|=" << n_parall << " fV=" << StringVectorGAP(fV) << " M1=" << StringMatrixGAP_line(M1) << "\n";
        for (size_t i_parall=0; i_parall<n_parall; i_parall++) {
          MyMatrix<Tint> M_paral = reorder_matrix(rrc.list_parallelepipeds[i_parall]);
          bool in_list = set.contains(M_paral);
          ExtendedGenericRobustM<T, Tint> egr = get_generic_robust_m(M_paral, GramMat, fV_red, os);
          MyVector<Tint> v_l = egr.robust_m.v_long();
          MyVector<T> diff = eV_red - UniversalVectorConversion<T,Tint>(v_l);
          T norm = EvaluationQuadForm(GramMat, diff);
          std::cerr << "ROBUST i_parall=" << i_parall << "/" << n_parall << " norm=" << egr.max << " v_long=" << StringVectorGAP(v_l) << " in_list=" << in_list << " norm=" << norm << " M=" << StringMatrixGAP_line(M_paral) << "\n";
        }
        std::cerr << "ROBUST: That is not what is expected for an inner point\n";
        throw TerminalException{1};
      }

      MyMatrix<Tint> M2 = reorder_matrix(rrc.list_parallelepipeds[0]);
      ExtendedGenericRobustM<T, Tint> egr = get_generic_robust_m(M2, GramMat, fV_red, os);
      if (v_long != egr.robust_m.v_long()) {
        std::cerr << "ROBUST: v_long is inconsistent\n";
        throw TerminalException{1};
      }
      if (M1 != M2) {
        std::cerr << "ROBUST: M1=\n";
        WriteMatrix(std::cerr, M1);
        std::cerr << "ROBUST: M2=\n";
        WriteMatrix(std::cerr, M2);
        std::cerr << "ROBUST: M1 should be equal to M2\n";
        std::cerr << "ROBUST: That is not what is expected for an inner point\n";
        throw TerminalException{1};
      }
      N += 1;
    }
  }
#endif
  return opt;
}

/*
  Get a possible vector to consider
 */
template<typename T, typename Tint>
std::optional<ConvexBoundary<T>> get_next_side_vector(SoftConvexBoundary<T,Tint> const& scb,
                                                      MyMatrix<T> const& G,
                                                      MyVector<Tint> const& v_crit,
                                                      MyVector<Tint> const& v_long,
                                                      std::vector<MyVector<Tint>> const& l_excluded_max,
                                                      std::ostream& os) {

  std::vector<MyVector<T>> l_vertices = scb.cb.get_list_vertices();
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: gnsv v_long=" << StringVectorGAP(v_long) << "\n";
  os << "ROBUST: gnsv l_vertices=\n";
  WriteMatrix(os, MatrixFromVectorFamily(l_vertices));
#endif
  MyVector<T> v_long_T = UniversalVectorConversion<T,Tint>(v_long);
  int dim = G.rows();
  MyVector<T> diff1(dim);
  MyVector<T> diff2(dim);
  auto is_relevant=[&](MyVector<Tint> const& cand) -> bool {
    MyVector<T> cand_T = UniversalVectorConversion<T,Tint>(cand);
    for (auto & vert: l_vertices) {
      for (int i=0; i<dim; i++) {
        diff1(i) = vert(i + 1) - v_long_T(i);
        diff2(i) = vert(i + 1) - cand_T(i);
      }
      T norm1 = EvaluationQuadForm(G, diff1);
      T norm2 = EvaluationQuadForm(G, diff2);
      if (norm2 > norm1) {
        return true;
      }
    }
    return false;
  };

  std::vector<std::vector<MyVector<Tint>>> ll_cand_rel;
  std::vector<int> VectSiz;
#ifdef DEBUG_ENUM_P_POLYTOPES
  int iter = 0;
#endif
  for (auto &robust_m: scb.l_robust_m) {
    std::vector<MyVector<Tint>> l_cand = DifferenceVect(robust_m.get_short_vectors(), l_excluded_max);
    std::vector<MyVector<Tint>> l_cand_rel;
    for (auto & cand: l_cand) {
      if (is_relevant(cand)) {
        l_cand_rel.push_back(cand);
      }
    }
    int n_choice = l_cand_rel.size();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: gnsv iter=" << iter << " n_choice=" << n_choice << "\n";
    iter += 1;
#endif
    if (n_choice == 0) {
      return {};
    }
    VectSiz.push_back(n_choice);
    ll_cand_rel.push_back(l_cand_rel);
  }
  BlockIterationMultiple bim(VectSiz);
  for (auto & choices: bim) {
    std::vector<MyVector<T>> l_ineq;
    for (size_t u=0; u<ll_cand_rel.size(); u++) {
      int ch = choices[u];
      MyVector<Tint> const& v_long_new = ll_cand_rel[u][ch];
      MyVector<T> eIneq_new = get_ineq(G, v_crit, v_long_new);
      l_ineq.push_back(eIneq_new);
    }
    MyMatrix<T> FAC_inp = MatrixFromVectorFamily(l_ineq);
    std::optional<ConvexBoundary<T>> opt1 = convexboundary_halfspaces_int(scb.cb, FAC_inp, os);
    if (opt1) {
      ConvexBoundary<T> const& cb_new = *opt1;
      return cb_new;
    }
  }
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: gnsv nothing was found\n";
#endif
  return {};
}



template <typename T, typename Tint>
std::optional<PVoronoi<T, Tint>>
find_p_voronoi(CVPSolver<T, Tint> const &solver, MyVector<T> const &eV, std::ostream &os) {
  MyMatrix<T> const& G = solver.GramMat;
  std::optional<PVoronoiPart<T,Tint>> opt = kernel_l1_p_polytope_part<T,Tint>(solver, {}, eV, os);
  if (!opt) {
    return {};
  }
  PVoronoiPart<T, Tint> pvp = *opt;
  T min_norm = min_pairwise_norm(pvp.l_cb[0].sp.EXT, G);
  MyVector<Tint> v_crit = pvp.robust_m_min.v_long();
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: find_p_voronoi, step 1\n";
#endif

  auto f_process_entry=[&]() -> bool {
    size_t len = pvp.l_scb.size();
    if (len == 0) {
      return true;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 1\n";
#endif
    SoftConvexBoundary<T,Tint> scb = pvp.l_scb[len-1];
    pvp.l_scb.pop_back();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 2\n";
#endif
    std::vector<MyVector<Tint>> l_excluded_max = scb.l_excluded_max;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 3\n";
#endif
    MyVector<Tint> v_long = scb.v_long();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 4\n";
#endif
    l_excluded_max.push_back(v_long);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 5\n";
#endif
    MyVector<T> eIneq = get_ineq(G, v_crit, v_long);
    MyVector<T> eIneq_op = - eIneq;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 6\n";
#endif
    std::optional<ConvexBoundary<T>> opt1 = get_next_side_vector(scb, G, v_crit, v_long, l_excluded_max, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 7\n";
#endif
    if (!opt1) {
      HardConvexBoundary<T> hcb{scb.index_cb, scb.cb};
      pvp.l_hcb.push_back(hcb);
      return false;
    }
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 8\n";
#endif
    ConvexBoundary<T> const& cb_new = *opt1;
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: fpe, step 9 min_norm=" << min_norm << "\n";
#endif
    T shift = min_norm / T(10);
    int N = 1;
    while(true) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_1, N=" << N << " shift=" << shift << "\n";
#endif
      MyVector<T> eIso = cb_new.random_interior_point(N, os);
      InteriorPtDir<T> ipd_new{eIso, cb_new.V};
      MyVector<T> fV = ipd_new.get_point(shift);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_2, shift=" << shift << "\n";
#endif
      std::optional<PVoronoiPart<T,Tint>> opt3 = kernel_l1_p_polytope_part<T,Tint>(solver, l_excluded_max, fV, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_3\n";
#endif
      if (!opt3) {
#ifdef DEBUG_ENUM_P_POLYTOPES
        os << "ROBUST: fpe, op3 exit\n";
#endif
        shift = shift / 2;
        N += 1;
        continue;
      }
      PVoronoiPart<T,Tint> const& p_poly_vor_part = *opt3;
      int pos = get_position_vec_in_mat(p_poly_vor_part.l_cb[0].sp.FAC, eIneq_op);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_4\n";
#endif
      if (pos == -1) {
#ifdef DEBUG_ENUM_P_POLYTOPES
        os << "ROBUST: fpe, p_poly_vor_part.l_cb[0].sp.FAC=\n";
        WriteMatrix(os, p_poly_vor_part.l_cb[0].sp.FAC);
        os << "ROBUST: fpe, eIneq_op=" << StringVectorGAP(eIneq_op) << "\n";
        os << "ROBUST: fpe, pos = -1 exit\n";
#endif
        shift = shift / 2;
        N += 1;
        continue;
      }
      SinglePolytope<T> const& sp = p_poly_vor_part.l_cb[0].sp;
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_5\n";
#endif
      std::vector<ConvexBoundary<T>> l_cb = convec_boundary_minus_sp(scb.cb, sp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_6\n";
#endif
      for (auto& cb2: l_cb) {
        SoftConvexBoundary<T,Tint> scb_new{scb.index_cb, cb2, l_excluded_max, scb.l_robust_m};
        pvp.l_scb.push_back(scb_new);
      }
      for (auto& hcb: p_poly_vor_part.l_hcb) {
        pvp.l_hcb.push_back(hcb);
      }
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_7\n";
#endif
      pvp.l_cb.push_back(p_poly_vor_part.l_cb[0]);
      int index = pvp.l_cb.size() - 1;
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: fpe, step 10_8\n";
#endif
      // That part needs to be improved, since the inserted faces could match
      for (auto& scb: p_poly_vor_part.l_scb) {
        if (scb.cb.V != eIneq_op) {
          SoftConvexBoundary<T,Tint> scb_new = scb;
          scb_new.index_cb = index;
          pvp.l_scb.push_back(scb_new);
        }
      }
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: f_process_entry, step 10_9\n";
#endif
      return false;
    }
  };


#ifdef DEBUG_ENUM_P_POLYTOPES
  size_t iter = 0;
#endif
  while(true) {
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: find_p_voronoi, |l_hcb|=" << pvp.l_hcb.size()
       << " |l_scb|=" << pvp.l_scb.size() << " iter=" << iter << "\n";
    iter += 1;
#endif
    bool test = f_process_entry();
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: find_p_voronoi, test=" << test << "\n";
#endif
    if (test) {
      break;
    }
  }
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: find_p_voronoi, Before convert_p_voronoi_part\n";
#endif
  PVoronoi<T,Tint> p_voronoi = convert_p_voronoi_part(pvp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: find_p_voronoi, EXT=\n";
  WriteMatrix(os, p_voronoi.EXT);
  MyMatrix<T> FAC = DirectDualDescription_mat(p_voronoi.EXT, os);
  os << "ROBUST: find_p_voronoi, FAC=\n";
  WriteMatrix(os, FAC);
#endif
  return p_voronoi;
}





template <typename T, typename Tint>
PVoronoi<T, Tint>
initial_p_polytope(CVPSolver<T, Tint> const &solver, std::ostream &os) {
  int dim = solver.GramMat.rows();
  int denom = 20;
  while (true) {
    MyVector<T> eV = get_random_vector<T>(denom, dim);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: ipp, before find_p_voronoi, eV="
       << StringVectorGAP(eV) << " denom=" << denom << "\n";
#endif
    std::optional<PVoronoi<T, Tint>> opt = find_p_voronoi(solver, eV, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: ipp, after find_p_voronoi\n";
#endif
    if (opt) {
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: ipp, return a PVoronoi\n";
      os << "ROBUST: ipp, Pvoronoi=" << *opt << "\n";
#endif
      return *opt;
    }
    denom += 1;
  }
}

template<typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> get_p_voronoi_stabilizer(DataLattice<T, Tint, Tgroup> &eData,
                                                     PVoronoi<T, Tint> const &pv) {
  std::ostream &os = eData.rddo.os;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using TgroupOut = std::pair<std::vector<Telt>, Tgroup>;
  auto f_group_out=[&](MyMatrix<T> const& EXTin) -> TgroupOut {
    MyMatrix<T> const& GramMat = eData.solver.GramMat;
    MyMatrix<T> const& SHV = eData.SHV;
#ifdef DEBUG_P_VORONOI_STABILIZER
    os << "ROBUST: f_group_out, step 1\n";
#endif
    Tgroup GRPbig = Polytope_StabilizerKernel<T,Tint,Tgroup>(GramMat, SHV, EXTin, os);
#ifdef DEBUG_P_VORONOI_STABILIZER
    os << "ROBUST: f_group_out, step 2\n";
#endif
    Tidx n_act = EXTin.rows();
    std::vector<Telt> LGenBig = GRPbig.SmallGeneratingSet();
    std::vector<Telt> LGenSma;
    auto f_correct=[&](Telt const& x) -> bool {
#ifdef DEBUG_P_VORONOI_STABILIZER
      os << "ROBUST: f_group_out / f_correct, step 1\n";
#endif
      MyMatrix<T> P = RepresentVertexPermutation(EXTin, EXTin, x);
#ifdef DEBUG_P_VORONOI_STABILIZER
      os << "ROBUST: f_group_out / f_correct, step 2\n";
#endif
      GeneralizedPolytope<T> gp_img = mat_product(pv.gp, P);
#ifdef DEBUG_P_VORONOI_STABILIZER
      os << "ROBUST: f_group_out / f_correct, step 3\n";
#endif
      bool test = is_equal(pv.gp, gp_img, os);
#ifdef DEBUG_P_VORONOI_STABILIZER
      os << "ROBUST: f_group_out / f_correct, step 4\n";
#endif
      return test;
    };
#ifdef DEBUG_P_VORONOI_STABILIZER
    os << "ROBUST: f_group_out, step 3\n";
#endif
    return get_intermediate_group<Tgroup,decltype(f_correct)>(n_act, LGenSma, LGenBig, f_correct, os);
  };
  MyMatrix<T> const& EXT_T = pv.EXT;
  std::pair<std::vector<Telt>, Tgroup> pair1 = f_group_out(EXT_T);
  std::vector<MyMatrix<Tint>> l_gens;
  for (auto & eGen: pair1.first) {
    MyMatrix<T> P_T = RepresentVertexPermutation(EXT_T, EXT_T, eGen);
    MyMatrix<Tint> P = UniversalMatrixConversion<Tint,T>(P_T);
    l_gens.emplace_back(std::move(P));
  }
#ifdef DEBUG_P_VORONOI_STABILIZER
  int dim = eData.solver.GramMat.rows();
  os << "ROBUST: get_p_voronoi_stabilizer |G|=" << pair1.second.size() << "\n";
  std::vector<MyVector<Tint>> EXT1 = pv.robust_m_min.get_short_vectors();
  int n_vect = EXT1.size();
  MyMatrix<T> EXT2(n_vect, dim + 1);
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    EXT2(i_vect, 0) = 1;
    MyVector<Tint> const& V = EXT1[i_vect];
    for (int i=0; i<dim; i++) {
      EXT2(i_vect, i+1) = UniversalScalarConversion<T,Tint>(V(i));
    }
  }
  os << "ROBUST: We have EXT2=\n";
  WriteMatrix(os, EXT2);
  std::pair<std::vector<Telt>, Tgroup> pair2 = f_group_out(EXT2);
  os << "ROBUST: get_p_voronoi_stabilizer |G(pv)|=" << pair1.second.size()
     << " |G(parall)|=" << pair2.second.size() << "\n";
#endif
  return l_gens;
}






template <typename T, typename Tint, typename Tgroup>
std::vector<PVoronoi<T, Tint>>
find_list_adjacent_p_voronoi(DataLattice<T, Tint, Tgroup> &eData,
                             PVoronoi<T, Tint> const &pv) {
  std::ostream &os = eData.rddo.os;
  CVPSolver<T, Tint> const &solver = eData.solver;
  MyMatrix<T> const& G = solver.GramMat;
  T min_norm = min_pairwise_norm(pv.EXT, G);
  BoundaryGeneralizedPolytope<T> bnd = find_generalized_polytope_boundary(pv.gp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: flapv, We have bnd\n";
  //  print_raw_boundary(bnd, os);
#endif
  auto get_adj_p_polytope = [&](InteriorPtDir<T> const& ipd_test,
                                MyVector<T> const &x)
      -> std::optional<PVoronoi<T, Tint>> {
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: flapv, gapp x=" << StringVector(x) << "\n";
#endif
    std::optional<PVoronoi<T, Tint>> opt = find_p_voronoi(solver, x, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: flapv, opt.has_value()=" << opt.has_value() << "\n";
#endif
    if (!opt) {
      return {};
    }
    PVoronoi<T, Tint> const &ppoly_adj = *opt;
    BoundaryGeneralizedPolytope<T> bnd_adj = find_generalized_polytope_boundary(ppoly_adj.gp, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: flapv, before is_boundary_point\n";
#endif
    bool test = is_boundary_point(ipd_test, bnd_adj, os);
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: flapv, ipd_test=" << ipd_test.to_string() << " test=" << test << "\n";
#endif
    if (!test) {
      return {};
    }
    return opt;
  };
  auto get_adj = [&]() -> PVoronoi<T, Tint> {
    T factor = min_norm;
    int N = 1;
    while (true) {
      std::optional<InteriorPtDir<T>> opt1 = get_random_interior_point_bnd(bnd, N, os);
#ifdef SANITY_CHECK_ENUM_P_POLYTOPES
      if (!opt1) {
        std::cerr << "ROBUST: Failed to find the random interior point\n";
        throw TerminalException{1};
      }
#endif
      InteriorPtDir<T> const& ipd = *opt1;
      InteriorPtDir<T> ipd_opp = ipd_opposite(ipd);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: flapv, ga, N=" << N << "     ipd=" << ipd.to_string() << "\n";
      os << "ROBUST: flapv, ga, N=" << N << " ipd_opp=" << ipd_opp.to_string() << "\n";
#endif
      MyVector<T>  x = ipd_opp.get_point(factor);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: flapv, ga, factor=" << factor << "\n";
#endif
      std::optional<PVoronoi<T, Tint>> opt2 =
          get_adj_p_polytope(ipd_opp, x);
      if (opt2) {
        return *opt2;
      }
      factor = factor / 2;
      N += 1;
    }
  };
  std::vector<PVoronoi<T, Tint>> l_adj;
  while(true) {
#ifdef DEBUG_ENUM_P_POLYTOPES
    os << "ROBUST: flapv |l_adj|=" << l_adj.size() << "\n";
#endif
    bool test = bnd.is_empty();
    if (test) {
      break;
     }
    PVoronoi<T, Tint> eadj = get_adj();
    reduce_boundary_generalized_polytope(bnd, eadj.gp, os);
    l_adj.push_back(eadj);
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
#ifdef DEBUG_ENUM_P_POLYTOPES
  os << "ROBUST: capp, After initial_p_polytope\n";
#endif
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
    os << "ROBUST: capp, start=" << start << " len=" << len << "\n";
#endif
    for (size_t i_p = start; i_p < len; i_p++) {
      std::vector<PVoronoi<T, Tint>> l_adj =
          find_list_adjacent_p_voronoi(eData, l_ppoly[i_p]);
#ifdef DEBUG_ENUM_P_POLYTOPES
      os << "ROBUST: capp, i_p=" << i_p << " |l_adj|=" << l_adj.size() << "\n";
#endif
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

template <typename T, typename Tint>
T random_vertex_estimation_robust_covering(MyMatrix<T> const &GramMat, size_t n_iter,
                                           std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  int dim = GramMat.rows();
  T max_cov(0);
  std::vector<MyVector<Tint>> l_excluded_max;
  MyVector<T> diff(dim);
  for (size_t iter = 0; iter < n_iter; iter++) {
    int denom = random() % 1000000 + 1;
    MyVector<T> eV = get_random_vector<T>(denom, dim);
#ifdef PRINT_ENUM_P_POLYTOPES
    os << "ROBUST: robust vertex iter=" << iter << "/" << n_iter << "\n";
#endif
#ifdef DEBUG_ROBUST_VERTEX_ENUM
    os << "ROBUST: iter=" << iter << " kernel_l1_p_polytope_part eV=" << StringVectorGAP(eV)
       << " denom=" << denom << "\n";
#endif
    std::optional<PVoronoiPart<T, Tint>> opt =
      kernel_l1_p_polytope_part(solver, l_excluded_max, eV, os);
#ifdef DEBUG_ROBUST_VERTEX_ENUM
    os << "ROBUST: After kernel_l1_p_polytope_part\n";
#endif
    if (opt) {
      PVoronoiPart<T, Tint> const& pvp = *opt;
      MyVector<Tint> v_long = pvp.robust_m_min.v_long();
      MyVector<T> v_long_T = UniversalVectorConversion<T,Tint>(v_long);
      int n_row = pvp.l_cb[0].sp.EXT.rows();
      T max_local(0);
#ifdef SANITY_CHECK_ROBUST_VERTEX_ENUM
      MyVector<T> V_test(dim+1);
#endif
      for (int i_row=0; i_row<n_row; i_row++) {
        for (int i=0; i<dim; i++) {
          diff(i) = v_long_T(i) - pvp.l_cb[0].sp.EXT(i_row, i+1);
        }
        T norm = EvaluationQuadForm(GramMat, diff);
        if (norm > max_local) {
          max_local = norm;
#ifdef SANITY_CHECK_ROBUST_VERTEX_ENUM
          for (int i=0; i<=dim; i++) {
            V_test(i) = pvp.l_cb[0].sp.EXT(i_row, i);
          }
#endif
        }
      }
      if (max_local > max_cov) {
        max_cov = max_local;
#ifdef SANITY_CHECK_ROBUST_VERTEX_ENUM
        ResultRobustClosest<T, Tint> rrc =
          compute_robust_closest<T, Tint>(solver, V_test, os);
        if (rrc.robust_minimum != max_local) {
          std::cerr << "ROBUST: max_local=" << max_local << " robust_minimum=" << rrc.robust_minimum << "\n";
          std::cerr << "ROBUST: inconsistent norms\n";
          throw TerminalException{1};
        }
#endif
      }
    }
  }
  return max_cov;
}





template <typename Tint>
void WriteEntryCPP(std::ostream &os, GenericRobustM<Tint> const &grm) {
  os << grm.index << "\n";
  WriteMatrix(os, grm.M);
}

template <typename Tint>
GenericRobustM<Tint> ReadEntryCPP_GenericRobustM(std::istream &is) {
  int index;
  is >> index;
  MyMatrix<Tint> M = ReadMatrix<Tint>(is);
  return {index, M};
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, ConvexBlock<T, Tint> const &cb) {
  size_t n = cb.list_robust_m.size();
  os << n << "\n";
  for (size_t i = 0; i < n; i++) {
    WriteEntryCPP(os, cb.list_robust_m[i]);
  }
  WriteEntryCPP(os, cb.sp);
}

template <typename T, typename Tint>
ConvexBlock<T, Tint> ReadEntryCPP_ConvexBlock(std::istream &is) {
  size_t n;
  is >> n;
  std::vector<GenericRobustM<Tint>> list_robust_m;
  for (size_t i = 0; i < n; i++) {
    list_robust_m.push_back(ReadEntryCPP_GenericRobustM<Tint>(is));
  }
  SinglePolytope<T> sp = ReadEntryCPP_SinglePolytope<T>(is);
  return {list_robust_m, sp};
}

template <typename T>
void WriteEntryCPP(std::ostream &os, HardConvexBoundary<T> const &hcb) {
  os << hcb.index_cb << "\n";
  WriteEntryCPP(os, hcb.sp);
}

template <typename T>
HardConvexBoundary<T> ReadEntryCPP_HardConvexBoundary(std::istream &is) {
  int index_cb;
  is >> index_cb;
  ConvexBoundary<T> sp = ReadEntryCPP_ConvexBoundary<T>(is);
  return {index_cb, sp};
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, SoftConvexBoundary<T, Tint> const &scb) {
  os << scb.index_cb << "\n";
  WriteEntryCPP(os, scb.cb);
  size_t n_excl = scb.l_excluded_max.size();
  os << n_excl << "\n";
  for (size_t i = 0; i < n_excl; i++) {
    WriteVector(os, scb.l_excluded_max[i]);
  }
  size_t n_robust = scb.l_robust_m.size();
  os << n_robust << "\n";
  for (size_t i = 0; i < n_robust; i++) {
    WriteEntryCPP(os, scb.l_robust_m[i]);
  }
}

template <typename T, typename Tint>
SoftConvexBoundary<T, Tint> ReadEntryCPP_SoftConvexBoundary(std::istream &is) {
  int index_cb;
  is >> index_cb;
  ConvexBoundary<T> cb = ReadEntryCPP_ConvexBoundary<T>(is);
  size_t n_excl;
  is >> n_excl;
  std::vector<MyVector<Tint>> l_excluded_max;
  for (size_t i = 0; i < n_excl; i++) {
    l_excluded_max.push_back(ReadVector<Tint>(is));
  }
  size_t n_robust;
  is >> n_robust;
  std::vector<GenericRobustM<Tint>> l_robust_m;
  for (size_t i = 0; i < n_robust; i++) {
    l_robust_m.push_back(ReadEntryCPP_GenericRobustM<Tint>(is));
  }
  return {index_cb, cb, l_excluded_max, l_robust_m};
}

template <typename T, typename Tint>
void WriteEntryCPP(std::ostream &os, PVoronoi<T, Tint> const &pv) {
  WriteEntryCPP(os, pv.robust_m_min);
  size_t n_cb = pv.l_cb.size();
  os << n_cb << "\n";
  for (size_t i = 0; i < n_cb; i++) {
    WriteEntryCPP(os, pv.l_cb[i]);
  }
  size_t n_hcb = pv.l_hcb.size();
  os << n_hcb << "\n";
  for (size_t i = 0; i < n_hcb; i++) {
    WriteEntryCPP(os, pv.l_hcb[i]);
  }
  WriteEntryCPP(os, pv.gp);
  WriteMatrix(os, pv.EXT);
}

template <typename T, typename Tint>
PVoronoi<T, Tint> ReadEntryCPP_PVoronoi(std::istream &is) {
  GenericRobustM<Tint> robust_m_min = ReadEntryCPP_GenericRobustM<Tint>(is);
  size_t n_cb;
  is >> n_cb;
  std::vector<ConvexBlock<T, Tint>> l_cb;
  for (size_t i = 0; i < n_cb; i++) {
    l_cb.push_back(ReadEntryCPP_ConvexBlock<T, Tint>(is));
  }
  size_t n_hcb;
  is >> n_hcb;
  std::vector<HardConvexBoundary<T>> l_hcb;
  for (size_t i = 0; i < n_hcb; i++) {
    l_hcb.push_back(ReadEntryCPP_HardConvexBoundary<T>(is));
  }
  GeneralizedPolytope<T> gp = ReadEntryCPP_GeneralizedPolytope<T>(is);
  MyMatrix<T> EXT = ReadMatrix<T>(is);
  return {robust_m_min, l_cb, l_hcb, gp, EXT};
}

template <typename T, typename Tint>
PVoronoi<T, Tint> ReadEntryCPP_PVoronoi_File(std::string const& file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadMatrixFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadEntryCPP_PVoronoi<T,Tint>(is);
}

namespace boost::serialization {

template <class Archive, typename Tint>
inline void serialize(Archive &ar, GenericRobustM<Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index", val.index);
  ar &make_nvp("M", val.M);
}

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, ConvexBlock<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("list_robust_m", val.list_robust_m);
  ar &make_nvp("sp", val.sp);
}

template <class Archive, typename T>
inline void serialize(Archive &ar, HardConvexBoundary<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index_cb", val.index_cb);
  ar &make_nvp("sp", val.sp);
}

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, SoftConvexBoundary<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("index_cb", val.index_cb);
  ar &make_nvp("cb", val.cb);
  ar &make_nvp("l_excluded_max", val.l_excluded_max);
  ar &make_nvp("l_robust_m", val.l_robust_m);
}

template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, PVoronoi<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("robust_m_min", val.robust_m_min);
  ar &make_nvp("l_cb", val.l_cb);
  ar &make_nvp("l_hcb", val.l_hcb);
  ar &make_nvp("gp", val.gp);
  ar &make_nvp("EXT", val.EXT);
}

} // namespace boost::serialization

// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
