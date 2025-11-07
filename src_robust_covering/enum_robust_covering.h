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
      are equal lead to some numerical 

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
  int dim = eV.size();
  int pow = pow_two(dim);
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   enumerating, eff_min=" << eff_min << "\n";
      int i_face = 0;
#endif
      std::vector<MyMatrix<Tint>> list_min_parallelepipeds;
      std::vector<MyMatrix<Tint>> tot_list_parallelepipeds;
      for (auto & eFace: l_face) {
        T local_max_norm(0);
#ifdef DEBUG_ENUM_ROBUST_COVERING
        os << "ROBUST:   i_face=" << i_face << " eFace=" << eFace << "\n";
#endif
        for (int& vert: FaceToVector<int>(eFace)) {
          if (l_norm[vert] > local_max_norm) {
            local_max_norm = l_norm[vert];
          }
        }
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
  os << "ROBUST: compute_robust_closest, step 1\n";
  int n_iter = 0;
#endif
  std::optional<T> min_search;
  while(true) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: compute_robust_closest, step 2, n_iter=" << n_iter << "\n";
    n_iter += 1;
#endif
    std::optional<ResultDirectEnumeration<T,Tint>> opt = compute_and_enumerate_structures(solver, eV, min_search, os);
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
  T min;
  bool is_correct;
  GenericRobustM<Tint> robust_m;
};



template<typename T, typename Tint>
ExtendedGenericRobustM<T,Tint> get_generic_robust_m(MyMatrix<Tint> const& M, MyMatrix<T> const&G, MyVector<T> const& eV) {
  T min(0);
  int best_index = 0;
  int n_ineq = M.rows();
  size_t n_att = 0;
  for (int index=0; index<n_ineq; index++) {
    MyVector<Tint> fV = GetMatrixRow(M, index);
    MyVector<T> diff = UniversalVectorConversion<T,Tint>(fV) - eV;
    T norm = EvaluationQuadForm(G, diff);
    if (index == 0) {
      min = norm;
      best_index = index;
      n_att = 1;
    } else {
      if (norm == min) {
        n_att += 1;
      } else {
        if (norm < min) {
          min = norm;
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
  GenericRobustM<Tint> robust_m{best_index, M};
  return {min, is_correct, robust_m};
};


template<typename T>
struct PpolytopeFacetIncidence {
  std::vector<Face> l_face;
  std::vector<MyVector<T>> l_FAC;
  std::vector<MyVector<T>> l_Iso;
};

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
    ineq(1+i) = 2 * G_v(i);
  }
  return ineq;
}

// In a robust structure robust_m, the longest vector
template<typename T, typename Tint>
void insert_inner_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, std::vector<MyVector<T>> & ListIneq) {
  int n_row = robust_m.M.rows();
  MyVector<Tint> v_long = robust_m.v_long();
  for (int i=0; i<n_row; i++) {
    if (i != robust_m.index) {
      MyVector<Tint> v_short = GetMatrixRow(robust_m.M, i);
      MyVector<T> eIneq = get_ineq(G, v_short, v_long);
      ListIneq.push_back(eIneq);
    }
  }
}

// The farthest vector of the robust structure has to be farther from the v_short (which is the farthest of
// the best one)
template<typename T, typename Tint>
void insert_outer_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, MyVector<Tint> const& v_short, std::vector<MyVector<T>> & ListIneq) {
  MyVector<Tint> v_long = robust_m.v_long();
  MyVector<T> eIneq = get_ineq(G, v_short, v_long);
  ListIneq.push_back(eIneq);
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
PpolytopeFacetIncidence<T> get_p_polytope_incidence(MyMatrix<T> const& FAC, MyMatrix<T> const& EXT) {
  int n_ext = EXT.rows();
  int dim = EXT.cols();
  int dim_ext = dim - 1;
  auto get_face=[&](MyVector<T> const& eFAC) -> Face {
    Face f(n_ext);
    for (int i_ext=0; i_ext<n_ext; i_ext++) {
      T scal(0);
      for (int i=0; i<dim; i++) {
        scal += eFAC(i) * EXT(i_ext, i);
      }
      if (scal != 0) {
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
    MyVector<T> eIso = Isobarycenter(EXT);
    return eIso;
  };
  int n_fac = FAC.rows();
  using Tpair = std::pair<MyVector<T>, MyVector<T>>;
  std::unordered_map<Face, Tpair> map;
  for (int i_fac=0; i_fac<n_fac; i_fac++) {
    MyVector<T> eFAC = GetMatrixRow(FAC, i_fac);
    MyVector<T> fFAC = ScalarCanonicalizationVector(eFAC);
    Face f = get_face(fFAC);
    std::optional<MyVector<T>> opt = get_iso_facet(f);
    if (opt) {
      MyVector<T> const& eIso = *opt;
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
    std::vector<MyVector<T>> ListIneq;
    MyMatrix<Tint> const& min_m = list_min_parallelepipeds[0];
    ExtendedGenericRobustM<T, Tint> ext_robust_m_min = get_generic_robust_m(min_m, G, eV);
    if (!ext_robust_m_min.is_correct) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by !ext_robust_m_min.is_correct\n";
#endif
      is_correct = false;
      return true;
    }
    if (min == 0) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by min=0\n";
#endif
      is_correct = false;
      return true;
    }
    GenericRobustM<Tint> const& robust_m_min = ext_robust_m_min.robust_m;
    insert_inner_ineqs_parallelepiped(robust_m_min, G, ListIneq);
    MyVector<Tint> v_short = robust_m_min.v_long(); // It is the shortest for the other structures!
    std::vector<GenericRobustM<Tint>> list_robust_m;
    for (auto& eM: tot_list_parallelepipeds) {
      if (eM != min_m) {
        ExtendedGenericRobustM<T,Tint> ext_robust_m = get_generic_robust_m(eM, G, eV);
        if (!ext_robust_m.is_correct) {
#ifdef DEBUG_ENUM_ROBUST_COVERING_DISABLE
          os << "ROBUST:   initial_vertex_data_test_ev, is_correct=false by !ext_robust_m.is_correct\n";
#endif
          is_correct = false;
          return true;
        }
        GenericRobustM<Tint> const& robust_m = ext_robust_m_min.robust_m;
        if (ext_robust_m.min <= min) {
          std::cerr << "ROBUST: The parallelepiped has an even lower minimum\n";
          throw TerminalException{1};
        }
        insert_inner_ineqs_parallelepiped(robust_m, G, ListIneq);
        insert_outer_ineqs_parallelepiped(robust_m, G, v_short, ListIneq);
        list_robust_m.push_back(robust_m);
      }
    }
    MyMatrix<T> FAC = MatrixFromVectorFamily(ListIneq);
    bool test = is_full_dimensional_bounded_polytope(FAC, os);
    if (!test) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST: initial_vertex_data_test_ev, failing by is_full_dimensional_bounded_polytope\n";
#endif
      return false;
    }
    std::optional<MyMatrix<T>> opt_ext = get_p_polytope_vertices_and_test_them<T,Tint>(solver, FAC, v_short, os);
    if (!opt_ext) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST: initial_vertex_data_test_ev, failing by get_p_polytope_vertices_and_test_them\n";
#endif
      return false;
    }
    MyMatrix<T> const& EXT = *opt_ext;
    MyVector<T> eIso = Isobarycenter(EXT);
    PpolytopeFacetIncidence<T> ppfi = get_p_polytope_incidence(FAC, EXT);
    ppoly = PpolytopeVoronoiData<T,Tint>{robust_m_min, list_robust_m, FAC, EXT, eIso, ppfi};
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: initial_vertex_data, before initial_vertex_data_test_ev, eV=" << StringVectorGAP(eV) << " denom=" << denom << "\n";
#endif
    std::optional<PpolytopeVoronoiData<T,Tint>> opt = initial_vertex_data_test_ev(solver, eV, os);
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
  auto get_adj_p_polytope=[&](MyVector<T> const& TestFAC, MyVector<T> const& x, std::unordered_set<MyVector<T>> const& set) -> std::optional<PpolytopeVoronoiData<T,Tint>> {
    std::optional<PpolytopeVoronoiData<T,Tint>> opt = initial_vertex_data_test_ev(solver, x, os);
    if (!opt) {
      return {};
    }
    PpolytopeVoronoiData<T,Tint> const& ppoly = *opt;
    size_t m_facet = ppoly.ppfi.l_FAC.size();
    int m_ext = ppoly.EXT.rows();
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
    T factor(1);
    while(true) {
      MyVector<T> x = eIso + factor * delta_x;
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
