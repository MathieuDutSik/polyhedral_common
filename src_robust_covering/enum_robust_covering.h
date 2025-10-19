// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "Shvec_exact.h"
#include "FundamentalDelaunay.h"
// clang-format on

/*
  The scheme is explained in
  "New upper bound for lattice covering by spheres"
  https://arxiv.org/pdf/2508.06446
  ----
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
  + So, what are the needed functions?
    x Finding one initial vertex.
    x From one vertex, finding the adjacent vertices.
    x Finding stabilizer of vertices.
    x Testing equivalence of vertices.
  + What is the structure of the vertices?
    x It has to contain all the parallelepiped in
      which the vertex is contained.
    x Together with the vertices that are at this distance.
    x The collection of those vertices define

    make the 

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
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   NextInTree, i_level=" << i_level << "\n";
#endif
    int choice = l_levels[i_level].choice;
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   NextInTree, choice=" << choice << "\n";
#endif
    PartSolution const& psol = l_levels[i_level].l_sol[choice];
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST:   NextInTree, we have psol\n";
#endif
    if (i_level == p) {
      f_insert(psol);
      return GoUpNextInTree();
    } else {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   NextInTree, before span_part_solution\n";
#endif
      std::vector<PartSolution> new_sols = span_part_solution(psol);
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
  std::vector<Face> l_face;
  int dim = M.cols();
  int pow = pow_two(dim);
  if (dv.n_vect < pow) {
    // No point trying to enumerate when there are no solutions.
    return l_face;
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
      l_face.push_back(psol.full_set);
    }
  };
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   Before kernel_enumerate_parallelepiped\n";
#endif
  kernel_enumerate_parallelepiped(dv, dim, f_insert, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST:   After kernel_enumerate_parallelepiped\n";
#endif
  return l_face;
}

template<typename T, typename Tint>
struct ResultRobustClosest {
  T robust_minimum;
  std::vector<MyMatrix<Tint>> list_parallelepipeds;
};



template<typename T, typename Tint, typename Finsert>
void compute_robust_close_f(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, Finsert f_insert, [[maybe_unused]] std::ostream& os) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST: compute_robust_closest, step 1\n";
#endif
  T min_search(0);
  int n_iter = 0;
  int dim = eV.size();
  int pow = pow_two(dim);
  MyMatrix<Tint> M_sol(pow, dim);
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
  os << "ROBUST: compute_robust_closest, step 2\n";
#endif
  while(true) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
    os << "ROBUST: compute_robust_closest, n_iter = " << n_iter << "\n";
#endif
    if (n_iter == 0) {
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   Before solver.nearest_vectors\n";
#endif
      resultCVP<T, Tint> res_cvp = solver.nearest_vectors(eV);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   After solver.nearest_vectors\n";
#endif
      min_search = res_cvp.TheNorm;
      std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect, os);
      if (l_face.size() > 0) {
        std::vector<MyMatrix<Tint>> list_parallelepipeds;
        for (auto & eFace: l_face) {
          list_parallelepipeds.push_back(get_msol(res_cvp.ListVect, eFace));
        }
        bool test = f_insert(min_search, list_parallelepipeds, list_parallelepipeds);
        if (test) {
          return;
        }
      }
    } else {
      min_search = (min_search * T(3)) / T(2);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   Before solver.at_most_dist_vectors\n";
#endif
      std::vector<MyVector<Tint>> elist = solver.at_most_dist_vectors(eV, min_search);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   After solver.at_most_dist_vectors |elist|=" << elist.size() << "\n";
      int i_fv = 0;
#endif
      std::vector<T> l_norm;
      for (auto & fV: elist) {
        MyVector<T> diff = UniversalVectorConversion<T,Tint>(fV) - eV;
        T norm = EvaluationQuadForm(solver.GramMat, diff);
#ifdef DEBUG_ENUM_ROBUST_COVERING
        os << "ROBUST:   i_fv=" << i_fv << " norm=" << norm << "\n";
        i_fv += 1;
#endif
        l_norm.push_back(norm);
      }
      MyMatrix<Tint> M = MatrixFromVectorFamily(elist);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   Before enumerate_parallelepiped\n";
#endif
      std::vector<Face> l_face = enumerate_parallelepiped(M, os);
#ifdef DEBUG_ENUM_ROBUST_COVERING
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
        bool test = f_insert(eff_min, list_min_parallelepipeds, tot_list_parallelepipeds);
        if (test) {
          return;
        }
      }
    }
    n_iter += 1;
  }
}

// Find the robust closest minimum with the lambda expression.
template<typename T, typename Tint>
ResultRobustClosest<T,Tint> compute_robust_closest(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, [[maybe_unused]] std::ostream& os) {
  ResultRobustClosest<T,Tint> result;
  auto f_insert=[&](T const& min, std::vector<MyMatrix<Tint>> const& list_min_parallelepipeds, [[maybe_unused]] std::vector<MyMatrix<Tint>> const& tot_list_parallelepipeds) {
    result = {min, list_min_parallelepipeds};
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
    os << "ROBUST: Before compute_robust_closest eV=" << StringVectorGAP(eV) << "\n";
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




template<typename Tint>
struct GenericRobustM {
  int index;
  MyMatrix<Tint> M;
  MyVector<Tint> v_short() const {
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

template<typename T, typename Tint>
struct InitialVoronoiData {
  GenericRobustM<Tint> robust_m_min;
  std::vector<GenericRobustM<Tint>> list_robust_m;
  MyMatrix<T> FAC;
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

template<typename T, typename Tint>
void insert_inner_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, std::vector<MyVector<T>> & ListIneq) {
  int n_row = robust_m.M.rows();
  MyVector<Tint> v_short = robust_m.v_short();
  for (int i=0; i<n_row; i++) {
    if (i != robust_m.index) {
      MyVector<Tint> v_long = GetMatrixRow(robust_m.M, i);
      MyVector<T> eIneq = get_ineq(G, v_short, v_long);
      ListIneq.push_back(eIneq);
    }
  }
}

template<typename T, typename Tint>
void insert_outer_ineqs_parallelepiped(GenericRobustM<Tint> const& robust_m, MyMatrix<T> const& G, MyVector<Tint> const& v_short, std::vector<MyVector<T>> & ListIneq) {
  int n_row = robust_m.M.rows();
  for (int i=0; i<n_row; i++) {
    MyVector<Tint> v_long = GetMatrixRow(robust_m.M, i);
    MyVector<T> eIneq = get_ineq(G, v_short, v_long);
    ListIneq.push_back(eIneq);
  }
}





// Find the defining inequalities of a polytope.
// It should fail and return None if the point eV is not generic enough.
// Which should lead to an increase in randomness.
template<typename T, typename Tint>
std::optional<InitialVoronoiData<T,Tint>> initial_vertex_data_test_ev(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, std::ostream& os) {
  MyMatrix<T> const& G = solver.GramMat;
  // Working variables
  bool is_correct = true;
  InitialVoronoiData<T,Tint> ivd;
  // The lambda function.
  auto f_insert=[&](T const& min, std::vector<MyMatrix<Tint>> const& list_min_parallelepipeds, std::vector<MyMatrix<Tint>> const& tot_list_parallelepipeds) {
    if (list_min_parallelepipeds.size() > 1) {
      // Terminate the enumeration
      is_correct = false;
      return true;
    }
    std::vector<MyVector<T>> ListIneq;
    MyMatrix<Tint> const& min_m = list_min_parallelepipeds[0];
    ExtendedGenericRobustM<T, Tint> ext_robust_m_min = get_generic_robust_m(min_m, G, eV);
    if (!ext_robust_m_min.is_correct) {
      is_correct = false;
      return true;
    }
    if (min == 0) {
      is_correct = false;
      return true;
    }
    GenericRobustM<Tint> const& robust_m_min = ext_robust_m_min.robust_m;
    insert_inner_ineqs_parallelepiped(robust_m_min, G, ListIneq);
    MyVector<Tint> v_short = robust_m_min.v_short();
    std::vector<GenericRobustM<Tint>> list_robust_m;
    for (auto& eM: tot_list_parallelepipeds) {
      if (eM != min_m) {
        ExtendedGenericRobustM<T,Tint> ext_robust_m = get_generic_robust_m(eM, G, eV);
        if (!ext_robust_m.is_correct) {
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
    if (test) {
      ivd = InitialVoronoiData<T,Tint>{robust_m_min, list_robust_m, FAC};
      return true;
    } else {
      return false;
    }
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
  if (is_correct) {
    return ivd;
  } else {
    return {};
  }
}

template<typename T, typename Tint>
InitialVoronoiData<T,Tint> initial_vertex_data(CVPSolver<T,Tint> const& solver, std::ostream& os) {
  int dim = solver.GramMat.rows();
  int denom = 2;
  MyVector<T> eV(dim);
  while(true) {
    T denom_T(denom);
    for (int i=0; i<dim; i++) {
      int val = random() % denom;
      T val_T(val);
      eV(i) = val_T / denom_T;
    }
    std::optional<InitialVoronoiData<T,Tint>> opt = initial_vertex_data_test_ev(solver, eV, os);
    if (opt) {
      return *opt;
    }
    denom += 1;
  }
}

// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
