// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "Shvec_exact.h"
#include "LatticeDelaunay.h"
// clang-format on

/*
  The motivation is explained in
  "New upper bound for lattice covering by spheres"
  https://arxiv.org/pdf/2508.06446
  ----
  It is about parallelepipeds of determinant 1.
  But this is a very specific case.
  Other case to consider are Order-k Voronoi
  polytopes.
  ---
  The technology applies as well.
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
  + What is the structure of the cells?
    x It contains a bunch of inequalities.
    x They can be reduced by linear programming (Clarkson).
    x The flipping is done via computing an adjacent point
      checking if it shares a facet, or not.
    x Add a facetness check for sanity check.
    x There is iteration. But we have iteration for
      Delaunay as well. So, maybe this is what to expect.
  + Can we have a L-type theory?
    x The vertices of the cells are defined in the same way
      as the Delaunay center. Expressing that some distances
      are equal lead to some linear inequalities just
      as in this context.

 */

#ifdef DEBUG
#define DEBUG_ENUM_ROBUST_COVERING
#endif

struct PartSolution {
  int vert;
  std::vector<int> l_dir;
  Face full_set;
};

template<typename Tint>
struct DataVect {
  int n_vect;
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
};

template<typename Tint>
DataVect<Tint> get_data_vect(MyMatrix<Tint> const& M) {
  int n_vect = M.rows();
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(M, i_vect);
    map[V] = i_vect;
    ListV.push_back(V);
  }
  return {n_vect, std::move(ListV), std::move(map)};
}

/*
  We want to enumerate the possible parallelepipeds
  of fixed dimension.
  It is a simple tree search
 */
template<typename Tint, typename Finsert>
void kernel_enumerate_parallelepiped(DataVect<Tint> const& dv, int const& p, Finsert f_insert, [[maybe_unused]] std::ostream& os) {
  int n_vect = dv.n_vect;
  int miss_val = std::numeric_limits<int>::max();
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   kernel_enumerate_parallelepiped, n_vect=" << n_vect << "\n";
#endif

  auto span_new_solution=[&](PartSolution const& psol, int const& newdir) -> std::optional<PartSolution> {
    Face new_set = psol.full_set;
    MyVector<Tint> trans = dv.ListV[newdir] - dv.ListV[psol.vert];
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      if (psol.full_set[i_vect] == 1) {
        MyVector<Tint> newV = dv.ListV[i_vect] + trans;
        auto iter = dv.map.find(newV);
        if (iter == dv.map.end()) {
          return {};
        }
        int pos = iter->second;
        new_set[pos] = 1;
      }
    }
    std::vector<int> l_dir = psol.l_dir;
    l_dir.push_back(newdir);
    PartSolution newsol{psol.vert, std::move(l_dir), std::move(new_set)};
    return newsol;
  };

  auto span_part_solution=[&](PartSolution const& psol) -> std::vector<PartSolution> {
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    os << "ROBUST:   span_part_solution |full_set|=" << psol.full_set.size() << " / " << psol.full_set.count() << "\n";
#endif
    std::vector<PartSolution> list_sol;
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      if (psol.full_set[i_vect] == 0) {
        std::optional<PartSolution> opt = span_new_solution(psol, i_vect);
        if (opt) {
          list_sol.push_back(*opt);
        }
      }
    }
    return list_sol;
  };
  auto get_empty=[&]() -> PartSolution {
    return {miss_val, {}, {}};
  };
  auto get_all_starts=[&]() -> std::vector<PartSolution> {
    std::vector<PartSolution> l_sol;
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      Face full_set(n_vect);
      full_set[i_vect] = 1;
      PartSolution esol{i_vect, {}, full_set};
      l_sol.push_back(esol);
    }
    return l_sol;
  };
  struct OneLevel {
    PartSolution prev_sol;
    std::vector<PartSolution> l_sol;
    size_t choice;
  };
  auto get_initial=[&]() -> OneLevel {
    PartSolution prev_sol = get_empty();
    std::vector<PartSolution> l_sol = get_all_starts();
    size_t choice = 0;
    return {prev_sol, l_sol, choice};
  };
  std::vector<OneLevel> l_levels{get_initial()};
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
  os << "ROBUST:   kernel_enumerate_parallelepiped, l_levels\n";
#endif
  int i_level = 0;
  auto GoUpNextInTree=[&]() -> bool {
    while(true) {
      OneLevel & level = l_levels[i_level];
      if (level.choice < level.l_sol.size() - 1) {
        level.choice += 1;
        return true;
      }
      if (i_level == 0) {
        return false;
      }
      i_level -= 1;
    }
  };
  auto NextInTree=[&]() -> bool {
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    os << "ROBUST:   NextInTree, i_level=" << i_level << "\n";
#endif
    int choice = l_levels[i_level].choice;
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    os << "ROBUST:   NextInTree, choice=" << choice << "\n";
#endif
    PartSolution const& psol = l_levels[i_level].l_sol[choice];
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    os << "ROBUST:   NextInTree, we have psol\n";
#endif
    if (i_level == p) {
      f_insert(psol);
      return GoUpNextInTree();
    } else {
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
      os << "ROBUST:   NextInTree, before span_part_solution\n";
#endif
      std::vector<PartSolution> new_sols = span_part_solution(psol);
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
      os << "ROBUST:   NextInTree, after span_part_solution |new_sols|=" << new_sols.size() << "\n";
#endif
      if (new_sols.size() == 0) {
        return GoUpNextInTree();
      }
      size_t new_choice = 0;
      OneLevel new_level{psol, new_sols, new_choice};
      int new_i_level = i_level + 1;
      if (l_levels.size() >= static_cast<size_t>(new_i_level + 1)) {
        l_levels[new_i_level] = new_level;
      } else {
        l_levels.push_back(new_level);
      }
      i_level = new_i_level;
      return true;
    }
  };
  while(true) {
    bool test = NextInTree();
    if (!test) {
      break;
    }
  }
}

int pow_two(int dim) {
  int pow = 1;
  for (int u=0; u<dim; u++) {
    pow *= 2;
  }
  return pow;
}

template<typename Tint>
std::vector<Face> enumerate_parallelepiped(MyMatrix<Tint> const& M, std::ostream& os) {
  DataVect<Tint> dv = get_data_vect(M);
  std::unordered_set<Face> set_face;
  int dim = M.cols();
  int pow = pow_two(dim);
  if (dv.n_vect < pow) {
    // No point trying to enumerate when there are no solutions.
    return {};
  }
  MyMatrix<Tint> Mdet(dim,dim);
  auto f_insert=[&](PartSolution const& psol) -> void {
    int e_vert = psol.vert;
    for (int i=0; i<dim; i++) {
      int f_vert = psol.l_dir[i];
      for (int j=0; j<dim; j++) {
        Mdet(i,j) = M(f_vert, j) - M(e_vert, j);
      }
    }
    Tint det = DeterminantMat(Mdet);
    if (T_abs(det) == 1) {
      set_face.insert(psol.full_set);
    }
  };
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   Before kernel_enumerate_parallelepiped\n";
#endif
  kernel_enumerate_parallelepiped(dv, dim, f_insert, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   After kernel_enumerate_parallelepiped\n";
#endif
  std::vector<Face> l_face;
  for (auto & eFace: set_face) {
    l_face.push_back(eFace);
  }
  return l_face;
}

template<typename T, typename Tint>
struct ResultRobustClosest {
  T robust_minimum;
  std::vector<MyMatrix<Tint>> list_parallelepipeds;
};

template<typename T, typename Tint>
struct ResultDirectEnumeration {
  T min;
  std::vector<MyMatrix<Tint>> list_min_parallelepipeds;
  std::vector<MyMatrix<Tint>> tot_list_parallelepipeds;
};

template<typename T, typename Tint>
std::optional<ResultDirectEnumeration<T,Tint>> compute_and_enumerate_structures(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, std::optional<T> & opt, std::ostream& os) {
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
  os << "ROBUST:   compute_and_enumerate_structures, beginning, eV=" << StringVector(eV) << "\n";
#endif
  int dim = eV.size();
  int pow = pow_two(dim);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
  os << "ROBUST:   compute_and_enumerate_structures, dim=" << dim << " pow=" << pow << "\n";
#endif
  MyMatrix<Tint> M_sol(pow, dim);
  T factor = T(3) / T(2);
  auto get_msol=[&](MyMatrix<Tint> const& Min, Face const& eFace) -> MyMatrix<Tint> {
    int pos = 0;
    for (int& vert: FaceToVector<int>(eFace)) {
      for (int i=0; i<dim; i++) {
        M_sol(pos, i) = Min(vert, i);
      }
      pos += 1;
    }
    return M_sol;
  };
  if (opt) {
    T const& min_search = *opt;
    T new_val = min_search * factor;
    opt = new_val;
    std::vector<MyVector<Tint>> elist = solver.at_most_dist_vectors(eV, min_search);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
    os << "ROBUST:   After solver.at_most_dist_vectors min_search=" << min_search << " |elist|=" << elist.size() << "\n";
#endif
    std::vector<T> l_norm;
    for (auto & fV: elist) {
      MyVector<T> diff = UniversalVectorConversion<T,Tint>(fV) - eV;
      T norm = EvaluationQuadForm(solver.GramMat, diff);
      l_norm.push_back(norm);
    }
    MyMatrix<Tint> M = MatrixFromVectorFamily(elist);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
    os << "ROBUST:   Before enumerate_parallelepiped\n";
#endif
    std::vector<Face> l_face = enumerate_parallelepiped(M, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
    os << "ROBUST:   After enumerate_parallelepiped |l_face|=" << l_face.size() << "\n";
#endif
    if (l_face.size() > 0) {
      T eff_min = min_search + T(1);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
      os << "ROBUST:   enumerating, eff_min=" << eff_min << "\n";
      int i_face = 0;
#endif
      std::vector<MyMatrix<Tint>> list_min_parallelepipeds;
      std::vector<MyMatrix<Tint>> tot_list_parallelepipeds;
      for (auto & eFace: l_face) {
        T local_max_norm(0);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
        os << "ROBUST:   i_face=" << i_face << " eFace=" << eFace << "\n";
#endif
        for (int& vert: FaceToVector<int>(eFace)) {
          if (l_norm[vert] > local_max_norm) {
            local_max_norm = l_norm[vert];
          }
        }
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
        os << "ROBUST:   i_face=" << i_face << " local_max_norm=" << local_max_norm << "\n";
        i_face += 1;
#endif
        MyMatrix<Tint> Mparall = get_msol(M, eFace);
        tot_list_parallelepipeds.push_back(Mparall);
        if (local_max_norm < eff_min) {
          list_min_parallelepipeds.clear();
          eff_min = local_max_norm;
          list_min_parallelepipeds.push_back(Mparall);
        } else {
          if (local_max_norm == eff_min) {
            list_min_parallelepipeds.push_back(Mparall);
          }
        }
      }
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
      os << "ROBUST:   eff_min=" << eff_min << "\n";
#endif
      ResultDirectEnumeration<T,Tint> rde{eff_min, list_min_parallelepipeds, tot_list_parallelepipeds};
      return rde;
    } else {
      return {};
    }
  } else {
    resultCVP<T, Tint> res_cvp = solver.nearest_vectors(eV);
#ifdef DEBUG_ENUM_ROBUST_COVERING_PARALL_ENUM
    os << "ROBUST:   After solver.nearest_vectors\n";
#endif
    T min = res_cvp.TheNorm;
    opt = min;
    std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect, os);
    if (l_face.size() > 0) {
      std::vector<MyMatrix<Tint>> list_parallelepipeds;
      for (auto & eFace: l_face) {
        list_parallelepipeds.push_back(get_msol(res_cvp.ListVect, eFace));
      }
      ResultDirectEnumeration<T,Tint> rde{min, list_parallelepipeds, list_parallelepipeds};
      return rde;
    } else {
      return {};
    }
  }
}









template<typename T, typename Tint, typename Finsert>
void compute_robust_close_f(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, Finsert f_insert, std::ostream& os) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST: compute_robust_close_f, step 1\n";
  int n_iter = 0;
#endif
  std::optional<T> min_search;
  while(true) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: compute_robust_close_f, step 2, n_iter=" << n_iter << "\n";
    n_iter += 1;
#endif
    std::optional<ResultDirectEnumeration<T,Tint>> opt = compute_and_enumerate_structures(solver, eV, min_search, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: compute_robust_close_f, After compute_and_enumerate_structures\n";
#endif
    if (opt) {
      ResultDirectEnumeration<T,Tint> const& rde = *opt;
      bool test = f_insert(rde);
      if (test) {
        return;
      }
    }
  }
}

// Find the robust closest minimum with the lambda expression.
template<typename T, typename Tint>
ResultRobustClosest<T,Tint> compute_robust_closest(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, std::ostream& os) {
  ResultRobustClosest<T,Tint> result;
  auto f_insert=[&](ResultDirectEnumeration<T,Tint> const& rde) {
    result = {rde.min, rde.list_min_parallelepipeds};
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
  return result;
}

template<typename T, typename Tint>
T random_estimation_robust_covering(MyMatrix<T> const& GramMat, size_t n_iter, std::ostream & os) {
  CVPSolver<T,Tint> solver(GramMat, os);
  int dim = GramMat.rows();
  T max_cov(0);
  MyVector<T> eV(dim);
  for (size_t iter=0; iter<n_iter; iter++) {
    int denom = random() % 1000000000000000;
    T denom_T(denom);
    for (int i=0; i<dim; i++) {
      int val = random() % denom;
      T val_T(val);
      eV(i) = val_T / denom_T;
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: Before compute_robust_closest eV=" << StringVectorGAP(eV) << " denom=" << denom << "\n";
#endif
    ResultRobustClosest<T,Tint> rrc = compute_robust_closest<T,Tint>(solver, eV, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: After compute_robust_closest\n";
#endif
    if (rrc.robust_minimum > max_cov) {
      max_cov = rrc.robust_minimum;
    }
  }
  return max_cov;
}





// A family of vectors with the index which is the farthest.
template<typename Tint>
struct GenericRobustM {
  int index;
  MyMatrix<Tint> M;
  MyVector<Tint> v_long() const {
    return GetMatrixRow(M, index);
  }
};





template<typename T, typename Tint>
struct ExtendedGenericRobustM {
  T max;
  bool is_correct;
  GenericRobustM<Tint> robust_m;
};



template<typename T, typename Tint>
ExtendedGenericRobustM<T,Tint> get_generic_robust_m(MyMatrix<Tint> const& M, MyMatrix<T> const&G, MyVector<T> const& eV, [[maybe_unused]] std::ostream& os) {
  T max(0);
  int best_index = 0;
  int n_ineq = M.rows();
  size_t n_att = 0;
  for (int index=0; index<n_ineq; index++) {
    MyVector<Tint> fV = GetMatrixRow(M, index);
    MyVector<T> diff = UniversalVectorConversion<T,Tint>(fV) - eV;
    T norm = EvaluationQuadForm(G, diff);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    double norm_d = UniversalScalarConversion<double,T>(norm);
    os << "ROBUST:   get_generic_robust_m, index=" << index << " fV=" << StringVector(fV) << " norm=" << norm << " norm_d=" << norm_d << "\n";
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
  double max_d = UniversalScalarConversion<double,T>(max);
  os << "ROBUST:   get_generic_robust_m, best_index=" << best_index << " max=" << max << " max_d=" << max_d << "\n";
#endif
  GenericRobustM<Tint> robust_m{best_index, M};
  return {max, is_correct, robust_m};
};


template<typename T>
struct PpolytopeFacetIncidence {
  std::vector<Face> l_face;
  std::vector<MyVector<T>> l_FAC;
  std::vector<MyVector<T>> l_Iso;
};

/*
  The data structure for the P-polytope:
  (A) The parallelepiped realizing the minimum.
  (B) The list of other parallelepipeds being
  used to define the structure.
  (C) Other bureaucratic stuff: EXT, FAC, Iso, incidence, ...
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
  Define for a parallelepiped P, the function
  Define f(x,y) = || x - y ||^2
  Then the generalized polytope Can(P,v) in question will be:
  phi_P(x) = max_{y vertex of P} f(x,y)
  So phi_P(x) >= f(x,y)
  So if the optimal configuration is (P,v) then the domain
  that we want is
  f(x,v) < phi_Q(x) for all parallelepipeds Q != P.
  If we select for each of the point in Q a specific
  point Q_v then the set of inequalities
  f(x,v) < f(x, Q_v).
  This is a smaller object because f(x, Q_v) <= phi_Q(x).
  But it is a polytope. But Can(P, v) is not necessarily
  connected.
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
    strictly an union. No convexity here.
    + So, we have a non-deterministic algorithm
  ---- HOW TO MANAGE THE POLYTOPE UNIONS ----
  The decomposition as an union of polytope is not
  unique. So, it cannot be used for testing equivalence
  What can we do in the database:
  * Store all the vertices.
  * Store the facets that are opened to the outer world.
  * The family of polytopes making the decomposition.
  [ How does that help? unsure ]
  * The adjacencies between vertices.
  ---- THE L-TYPE THEORY SIDE ----
  Here is what we have:
  * Each vertex is defined like for Delaunay.
  * The adjacent vertices should 
  ----
 */
template<typename T, typename Tint>
struct PpolytopeVoronoiData {
  GenericRobustM<Tint> robust_m_min;
  std::vector<GenericRobustM<Tint>> list_robust_m;
  MyMatrix<T> FAC;
  MyMatrix<T> EXT;
  MyVector<T> eIso;
  PpolytopeFacetIncidence<T> ppfi;
};

/*
  We have G[x - v_short] <= G[x - v_long]
  So,
  -2 x G v_short + G[v_short] <= -2 x G v_long + G[v_long]
  which gets
  0 <= G[v_long] - G[v_short] + x ( 2 G (v_short - v_long))
 */
template<typename T, typename Tint>
MyVector<T> get_ineq(MyMatrix<T> const& G, MyVector<Tint> const& v_short, MyVector<Tint> const& v_long) {
  int dim = G.rows();
  T norm_short = EvaluationQuadForm(G, v_short);
  T norm_long = EvaluationQuadForm(G, v_long);
  MyVector<Tint> diff = v_short - v_long;
  MyVector<T> diff_T = UniversalVectorConversion<T,Tint>(diff);
  MyVector<T> G_v = G * diff_T;
  MyVector<T> ineq(1 + dim);
  ineq(0) = norm_long - norm_short;
  for (int i=0; i<dim; i++) {
    ineq(1 + i) = 2 * G_v(i);
  }
  return ScalarCanonicalizationVector(ineq);
}

// In a robust structure robust_m, the longest vector
template<typename T, typename Tint>
void insert_inner_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, std::vector<MyVector<T>> & ListIneq, [[maybe_unused]] std::ostream& os) {
  int n_row = robust_m.M.rows();
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   insert_inner_ineqs_parallelepiped, M=\n";
  WriteMatrix(os, robust_m.M);
#endif
  MyVector<Tint> v_long = robust_m.v_long();
  for (int i=0; i<n_row; i++) {
    if (i != robust_m.index) {
      MyVector<Tint> v_short = GetMatrixRow(robust_m.M, i);
      MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   insert_inner_ineqs_parallelepiped, eIneq=" << StringVector(eIneq) << " i=" << i << " index=" << robust_m.index << "\n";
#endif
      ListIneq.push_back(eIneq);
    }
  }
}

// The farthest vector of the robust structure has to be farther from the v_short (which is the farthest of
// the best one)
template<typename T, typename Tint>
void insert_outer_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, MyVector<Tint> const& v_short, std::vector<MyVector<T>> & ListIneq, [[maybe_unused]] std::ostream& os) {
  MyVector<Tint> v_long = robust_m.v_long();
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_short=" << StringVector(v_short) << "\n";
  os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_long=" << StringVector(v_long) << "\n";
#endif
  if (v_short != v_long) {
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   insert_outer_ineqs_parallelepiped, eIneq=" << StringVector(eIneq) << "\n";
#endif
    ListIneq.push_back(eIneq);
  } else {
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   insert_outer_ineqs_parallelepiped, v_short=" << StringVector(v_short) << "\n";
#endif
  }
}


template<typename T>
std::vector<MyVector<T>> get_p_polytope_vertices(MyMatrix<T> const& FAC, std::ostream& os) {
  std::string ansProg = "lrs";
  MyMatrix<T> EXTbig = DirectFacetComputationInequalities(FAC, ansProg, os);
  int n_ext = EXTbig.rows();
  int dim = EXTbig.cols() - 1;
  std::vector<MyVector<T>> EXT;
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<T> v = GetMatrixRow(EXTbig, i_ext);
    if (v(0) <= 0) {
      std::cerr << "We should have v(0) > 0, because we expect to have a polytope\n";
      throw TerminalException{1};
    }
    MyVector<T> eEXT(dim);
    for (int i=0; i<dim; i++) {
      eEXT(i) = v(i+1) / v(0);
    }
    EXT.push_back(eEXT);
  }
  return EXT;
}

template<typename T, typename Tint>
std::optional<MyMatrix<T>> get_p_polytope_vertices_and_test_them(CVPSolver<T,Tint> const& solver, MyMatrix<T> const& FAC, MyVector<Tint> const& v_short, std::ostream& os) {
  std::vector<MyVector<T>> EXT = get_p_polytope_vertices(FAC, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST: get_p_polytope_vertices_and_test_them |EXT|=" << EXT.size() << "\n";
#endif
  MyMatrix<T> const& GramMat = solver.GramMat;
  MyVector<T> v_short_T = UniversalVectorConversion<T,Tint>(v_short);
  for (auto & eEXT: EXT) {
    MyVector<T> diff = eEXT - v_short_T;
    T norm = EvaluationQuadForm(GramMat, diff);
    std::optional<T> opt = norm;
    std::optional<ResultDirectEnumeration<T,Tint>> opt_rde = compute_and_enumerate_structures(solver, eEXT, opt, os);
    if (!opt_rde) {
      std::cerr << "ROBUST: get_p_polytope_vertices_and_test_them failed. We should have one\n";
      throw TerminalException{1};
    }
    ResultDirectEnumeration<T,Tint> const& rde = *opt_rde;
    if (rde.min < norm) {
      return {};
    }
  }
  int n_ext = EXT.size();
  int dim = GramMat.cols();
  MyMatrix<T> EXTret(n_ext, dim+1);
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    EXTret(i_ext, 0) = 1;
    for (int i=0; i<dim; i++) {
      EXTret(i_ext, i+1) = EXT[i_ext](i);
    }
  }
  return EXTret;
}

template<typename T>
PpolytopeFacetIncidence<T> get_p_polytope_incidence(MyMatrix<T> const& FAC, MyMatrix<T> const& EXT, [[maybe_unused]] std::ostream& os) {
  int n_ext = EXT.rows();
  int dim = EXT.cols();
  int dim_ext = dim - 1;
#ifdef DEBUG_ENUM_ROBUST_COVERING_INCIDENCE
  os << "ROBUST: get_p_polytope_incidence n_ext=" << n_ext << " dim=" << dim << " dim_ext=" << dim_ext << " EXT=\n";
  WriteMatrix(os, EXT);
  os << "ROBUST: get_p_polytope_incidence rank=" << RankMat(EXT) << "\n";
#endif
  auto get_face=[&](MyVector<T> const& eFAC) -> Face {
    Face f(n_ext);
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      T scal(0);
      for (int i=0; i<dim; i++) {
        scal += eFAC(i) * EXT(i_ext, i);
      }
      if (scal == 0) {
        f[i_ext] = 1;
      }
    }
    return f;
  };
  auto get_iso_facet=[&](Face const& f) -> std::optional<MyVector<T>> {
    int cnt = f.count();
    if (cnt < dim_ext) {
      return {};
    }
    MyMatrix<T> EXTincd(cnt, dim);
    int pos = 0;
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      if (f[i_ext] == 1) {
        for (int i=0; i<dim; i++) {
          EXTincd(pos, i) = EXT(i_ext, i);
        }
        pos += 1;
      }
    }
    if (RankMat(EXTincd) != dim_ext) {
      return {};
    }
    MyVector<T> eIso = Isobarycenter(EXTincd);
    return eIso;
  };
  int n_fac = FAC.rows();
#ifdef DEBUG_ENUM_ROBUST_COVERING_INCIDENCE
  os << "ROBUST: get_p_polytope_incidence n_fac=" << n_fac << "\n";
#endif
  using Tpair = std::pair<MyVector<T>, MyVector<T>>;
  std::unordered_map<Face, Tpair> map;
  for (int i_fac=0; i_fac<n_fac; i_fac++) {
    MyVector<T> eFAC = GetMatrixRow(FAC, i_fac);
    MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
    Face f = get_face(fFAC);
#ifdef DEBUG_ENUM_ROBUST_COVERING_INCIDENCE
    os << "ROBUST: get_p_polytope_incidence i_fac=" << i_fac << " eFAC=" << StringVector(eFAC) << " fFAC=" << StringVector(fFAC) << " f=" << f << "\n";
#endif
    std::optional<MyVector<T>> opt = get_iso_facet(f);
    if (opt) {
      MyVector<T> const& eIso = *opt;
#ifdef DEBUG_ENUM_ROBUST_COVERING_INCIDENCE
      os << "ROBUST: get_p_polytope_incidence eIso=" << StringVector(eIso) << "\n";
#endif
      Tpair pair{fFAC, eIso};
      map[f] = pair;
    }
  }
  std::vector<Face> l_face;
  std::vector<MyVector<T>> l_FAC;
  std::vector<MyVector<T>> l_Iso;
  for (auto & kv: map) {
    l_face.push_back(kv.first);
    l_FAC.push_back(kv.second.first);
    l_Iso.push_back(kv.second.second);
  }
  return {l_face, l_FAC, l_Iso};
}

// Find the defining inequalities of a polytope.
// It should fail and return None if the point eV is not generic enough.
// Which should lead to an increase in randomness.
template<typename T, typename Tint>
std::optional<PpolytopeVoronoiData<T,Tint>> initial_vertex_data_test_ev(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, std::ostream& os) {
  MyMatrix<T> const& G = solver.GramMat;
  if (IsIntegralVector(eV)) {
    // Nothing cam be done here
    return {};
  }
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST: initial_vertex_data_test_ev eV=" << StringVectorGAP(eV) << "\n";
#endif
  // Working variables
  bool is_correct = true;
  PpolytopeVoronoiData<T,Tint> ppoly;
  // The lambda function.
  auto f_insert=[&](ResultDirectEnumeration<T,Tint> const& rde) {
    T const& min = rde.min;
    std::vector<MyMatrix<Tint>> const& list_min_parallelepipeds = rde.list_min_parallelepipeds;
    std::vector<MyMatrix<Tint>> const& tot_list_parallelepipeds = rde.tot_list_parallelepipeds;
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "-----------------------------------------------------------------------------------------------\n";
    os << "ROBUST:   initial_vertex_data_test_ev, |list_min_parallelepipeds|=" << list_min_parallelepipeds.size() << " |tot_list_parallelepipeds|=" << tot_list_parallelepipeds.size() << " min=" << min << "\n";
#endif
    if (list_min_parallelepipeds.size() > 1) {
      // Terminate the enumeration
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by |list_min_parallelepipeds| > 1\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, pass 1\n";
#endif
    std::vector<MyVector<T>> ListIneq;
    MyMatrix<Tint> const& min_m = list_min_parallelepipeds[0];
    ExtendedGenericRobustM<T, Tint> ext_robust_m_min = get_generic_robust_m(min_m, G, eV, os);
    if (!ext_robust_m_min.is_correct) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by !ext_robust_m_min.is_correct\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, pass 2\n";
#endif
    if (min == 0) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by min=0\n";
#endif
      is_correct = false;
      return true;
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3\n";
#endif
    GenericRobustM<Tint> const& robust_m_min = ext_robust_m_min.robust_m;
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, robust_m_min, index=" << robust_m_min.index << " M=\n";
    WriteMatrix(os, robust_m_min.M);
#endif
    insert_inner_ineqs_parallelepiped(robust_m_min, G, ListIneq, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, Initial FAC=\n";
    WriteMatrix(os, MatrixFromVectorFamily(ListIneq));
#endif
    MyVector<Tint> v_short = robust_m_min.v_long(); // It is the shortest for the other structures!
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, v_short=" << StringVectorGAP(v_short) << "\n";
#endif
    std::vector<GenericRobustM<Tint>> list_robust_m;
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3, step 1\n";
    size_t i_m = 0;
#endif
    for (auto& eM: tot_list_parallelepipeds) {
      if (eM != min_m) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
        os << "ROBUST:   ------------------------- " << i_m << " ------------------------\n";
        i_m += 1;
#endif
        ExtendedGenericRobustM<T,Tint> ext_robust_m = get_generic_robust_m(eM, G, eV, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
        os << "ROBUST:   initial_vertex_data_test_ev, ext_robust_m.robust_m, index=" << ext_robust_m.robust_m.index << " M=\n";
        WriteMatrix(os, ext_robust_m.robust_m.M);
#endif
        if (!ext_robust_m.is_correct) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
          os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by !ext_robust_m.is_correct\n";
#endif
          is_correct = false;
          return true;
        }
        GenericRobustM<Tint> const& robust_m = ext_robust_m.robust_m;
        if (ext_robust_m.max <= min) {
          std::cerr << "ROBUST: ext_robust_m.min=" << ext_robust_m.max << " min=" << min << "\n";
          std::cerr << "ROBUST: The parallelepiped has an even lower minimum\n";
          throw TerminalException{1};
        }
        insert_inner_ineqs_parallelepiped(robust_m, G, ListIneq, os);
        insert_outer_ineqs_parallelepiped(robust_m, G, v_short, ListIneq, os);
        list_robust_m.push_back(robust_m);
      }
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   initial_vertex_data_test_ev, pass 3, step 2\n";
#endif
    MyMatrix<T> FAC = MatrixFromVectorFamily(ListIneq);
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    os << "ROBUST:   initial_vertex_data_test_ev, |FAC|=" << FAC.rows() << " / " << FAC.cols() << "\n";
    MyMatrix<T> EXTbig = DirectFacetComputationInequalities(FAC, "lrs", os);
    os << "ROBUST:   initial_vertex_data_test_ev, EXTbig=\n";
    WriteMatrix(os, EXTbig);
    int n_fac = FAC.rows();
    int n_ext_big = EXTbig.rows();
    for (int i_fac=0; i_fac<n_fac; i_fac++) {
      for (int i_ext_big=0; i_ext_big<n_ext_big; i_ext_big++) {
        T scal(0);
        for (int i=0; i<FAC.cols(); i++) {
          scal += FAC(i_fac, i) * EXTbig(i_ext_big, i);
        }
        if (scal < 0) {
          std::cerr << "Incorrect vertices\n";
          throw TerminalException{1};
        }
      }
    }
#endif
    bool test = is_full_dimensional_bounded_polytope(FAC, os);
    if (!test) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST: initial_vertex_data_test_ev, failing by is_full_dimensional_bounded_polytope\n";
#endif
      return false;
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: initial_vertex_data_test_ev, before get_p_polytope_vertices_and_test_them\n";
#endif
    std::optional<MyMatrix<T>> opt_ext = get_p_polytope_vertices_and_test_them<T,Tint>(solver, FAC, v_short, os);
    if (!opt_ext) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST: initial_vertex_data_test_ev, failing by get_p_polytope_vertices_and_test_them\n";
#endif
      return false;
    }
    MyMatrix<T> const& EXT = *opt_ext;
    MyVector<T> eIso = Isobarycenter(EXT);
    PpolytopeFacetIncidence<T> ppfi = get_p_polytope_incidence(FAC, EXT, os);
    ppoly = PpolytopeVoronoiData<T,Tint>{robust_m_min, list_robust_m, FAC, EXT, eIso, ppfi};
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: initial_vertex_data_test_ev, successful end\n";
    os << "ROBUST: initial_vertex_data_test_ev, EXT=\n";
    WriteMatrix(os, EXT);
    os << "ROBUST: initial_vertex_data_test_ev, FAC=\n";
    WriteMatrix(os, MatrixFromVectorFamily(ppfi.l_FAC));
#endif
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
  if (is_correct) {
    return ppoly;
  } else {
    return {};
  }
}






template<typename T, typename Tint>
PpolytopeVoronoiData<T,Tint> initial_vertex_data(CVPSolver<T,Tint> const& solver, std::ostream& os) {
  int dim = solver.GramMat.rows();
  int denom = 2;
  MyVector<T> eV(dim);
  while(true) {
    T denom_T(denom);
    for (int i=0; i<dim; i++) {
      int val1 = random();
      int val = val1 % denom;
      T val_T(val);
      T quot = val_T / denom_T;
      eV(i) = quot;
    }
    eV(0) = T(7) / T(12);
    eV(1) = T(7) / T(13);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: initial_vertex_data, before initial_vertex_data_test_ev, eV=" << StringVectorGAP(eV) << " denom=" << denom << "\n";
#endif
    std::optional<PpolytopeVoronoiData<T,Tint>> opt = initial_vertex_data_test_ev(solver, eV, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
    std::cerr << "ROBUST: Stopping the execution\n";
    throw TerminalException{1};
#endif
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: initial_vertex_data, after initial_vertex_data_test_ev\n";
#endif
    if (opt) {
      return *opt;
    }
    denom += 1;
  }
}

template<typename T, typename Tint, typename Tgroup>
std::vector<PpolytopeVoronoiData<T,Tint>> find_adjacent_p_polytopes(DataLattice<T, Tint, Tgroup> &eData, PpolytopeVoronoiData<T,Tint> const& pvd) {
  std::ostream &os = eData.rddo.os;
  CVPSolver<T,Tint> const& solver = eData.solver;
  int dim = eData.GramMat.rows();
  auto get_adj_p_polytope=[&](MyVector<T> const& TestFAC, MyVector<T> const& x, std::unordered_set<MyVector<T>> const& set) -> std::optional<PpolytopeVoronoiData<T,Tint>> {
    MyVector<T> x_red(dim);
    for (int i=0; i<dim; i++) {
      x_red(i) = x(i + 1);
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: get_adj_p_polytope x_red=" << StringVector(x_red) << "\n";
#endif
    std::optional<PpolytopeVoronoiData<T,Tint>> opt = initial_vertex_data_test_ev(solver, x_red, os);
    if (!opt) {
      return {};
    }
    PpolytopeVoronoiData<T,Tint> const& ppoly = *opt;
    size_t m_facet = ppoly.ppfi.l_FAC.size();
    int m_ext = ppoly.EXT.rows();
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: get_adj_p_polytope m_facet=" << m_facet << " m_ext=" << m_ext << "\n";
#endif
    for (size_t j_facet=0; j_facet<m_facet; j_facet++) {
      if (ppoly.ppfi.l_FAC[j_facet] == TestFAC) {
        if (set.size() != ppoly.ppfi.l_face[j_facet].count()) {
          std::cerr << "ROBUST: Different incidence, so this is not what we expect\n";
          std::cerr << "ROBUST: Non-facetness is a big problem. Rethink the code if not a bug\n";
          throw TerminalException{1};
        }
        for (int j_ext=0; j_ext<m_ext; j_ext++) {
          if (ppoly.ppfi.l_face[j_facet][j_ext] == 1) {
            MyVector<T> fEXT = GetMatrixRow(ppoly.EXT, j_ext);
            if (set.count(fEXT) == 0) {
              std::cerr << "ROBUST: Vertex is not contained in the facet.\n";
              std::cerr << "ROBUST: Non-facetness is a big problem. Rethink the code if not a bug\n";
              throw TerminalException{1};
            }
          }
        }
        return ppoly;
      }
    }
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: get_adj_p_polytope, Failed to find a correct adjacent cone\n";
#endif
    return {};
  };
  auto get_adj=[&](Face const& f, MyVector<T> const& eFAC, MyVector<T> const& eIso) -> PpolytopeVoronoiData<T,Tint> {
    std::unordered_set<MyVector<T>> set;
    int n_ext = pvd.EXT.rows();
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      if (f[i_ext] == 1) {
        MyVector<T> eEXT = GetMatrixRow(pvd.EXT, i_ext);
        set.insert(eEXT);
      }
    }
    MyVector<T> TestFAC = - eFAC;
    MyVector<T> delta_x = eIso - pvd.eIso;
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: get_adj, eIso=" << StringVector(eIso) << " pvd.eIso=" << StringVector(pvd.eIso) << " delta_x=" << delta_x << "\n";
#endif
    T factor(1);
#ifdef DEBUG_ENUM_ROBUST_COVERING
    size_t n_iter = 0;
#endif
    while(true) {
      MyVector<T> x = eIso + factor * delta_x;
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST: get_adj, n_iter=" << n_iter << " factor=" << factor << "\n";
      n_iter += 1;
#endif
      std::optional<PpolytopeVoronoiData<T,Tint>> opt = get_adj_p_polytope(TestFAC, x, set);
      if (opt) {
        return *opt;
      }
      factor = factor / 2;
    }
  };
  std::vector<PpolytopeVoronoiData<T,Tint>> l_adj;
  size_t n_facet = pvd.ppfi.l_FAC.size();
  for (size_t i_facet=0; i_facet<n_facet; i_facet++) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: find_adjacent_p_polytopes i_facet=" << i_facet << " / " << n_facet << "\n";
#endif
    PpolytopeVoronoiData<T,Tint> eAdj = get_adj(pvd.ppfi.l_face[i_facet], pvd.ppfi.l_FAC[i_facet], pvd.ppfi.l_Iso[i_facet]);
    l_adj.push_back(eAdj);
  }
  return l_adj;
}







template<typename T, typename Tint, typename Tgroup>
std::vector<PpolytopeVoronoiData<T,Tint>> compute_all_p_polytopes(DataLattice<T, Tint, Tgroup> &eData) {
  std::ostream &os = eData.rddo.os;
  CVPSolver<T,Tint> const& solver = eData.solver;
  std::vector<PpolytopeVoronoiData<T,Tint>> l_ppoly;
  PpolytopeVoronoiData<T,Tint> ppoly = initial_vertex_data(solver, os);
  l_ppoly.push_back(ppoly);
  auto f_insert=[&](PpolytopeVoronoiData<T,Tint> const& f_ppoly) -> void {
    for (auto & e_ppoly: l_ppoly) {
      std::optional<MyMatrix<T>> opt_equiv = Polytope_TestEquivalence<T,Tint,Tgroup>(eData, e_ppoly.EXT, f_ppoly.EXT);
      if (opt_equiv.has_value()) {
        return;
      }
    }
    l_ppoly.push_back(f_ppoly);
  };
  size_t start = 0;
  while(true) {
    size_t len = l_ppoly.size();
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: compute_all_p_polytopes, start=" << start << " len=" << len << "\n";
#endif
    for (size_t i_ppoly=start; i_ppoly<len; i_ppoly++) {
      std::vector<PpolytopeVoronoiData<T,Tint>> l_adj = find_adjacent_p_polytopes(eData, l_ppoly[i_ppoly]);
      for (auto & eAdj: l_adj) {
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

template<typename T, typename Tint, typename Tgroup>
T compute_square_robust_covering_radius(DataLattice<T, Tint, Tgroup> &eData) {
  std::vector<PpolytopeVoronoiData<T,Tint>> l_ppoly = compute_all_p_polytopes(eData);
  T max_sqr_radius(0);
  MyMatrix<T> const& GramMat = eData.GramMat;
  int dim = GramMat.rows();
  MyVector<T> diff(dim);
  for (auto & ppoly: l_ppoly) {
    MyVector<Tint> v_long = ppoly.robust_m_min.v_long();
    MyVector<T> v_long_T = UniversalVectorConversion<T,Tint>(v_long);
    int n_ext = ppoly.EXT.rows();
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      for (int i=0; i<dim; i++) {
        diff(i) = v_long_T(i) - ppoly.EXT(i_ext, i+1);
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
