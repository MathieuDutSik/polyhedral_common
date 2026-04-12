// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_PARALL_SEARCH_H_
#define SRC_ROBUST_COVERING_PARALL_SEARCH_H_

// clang-format off
#include "Shvec_exact.h"
#include "LatticeDelaunay.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_ENUM_PARALL_SEARCH
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ENUM_PARALL_SEARCH
#endif

#ifdef PRINT
#define PRINT_ENUM_PARALL_SEARCH
#endif

#ifdef DISABLE_DEBUG_ENUM_PARALL_SEARCH
#undef DEBUG_ENUM_PARALL_SEARCH
#endif

struct PartSolution {
  int vert;
  std::vector<int> l_dir;
  Face full_set;
};

template <typename Tint> struct DataVect {
  int n_vect;
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
};

template <typename Tint> DataVect<Tint> get_data_vect(MyMatrix<Tint> const &M) {
  int n_vect = M.rows();
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
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
template <typename Tint, typename Finsert>
void kernel_enumerate_parallelepiped(DataVect<Tint> const &dv, int const &p,
                                     Finsert f_insert,
                                     [[maybe_unused]] std::ostream &os) {
  int n_vect = dv.n_vect;
  int miss_val = std::numeric_limits<int>::max();
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   kernel_enumerate_parallelepiped, n_vect=" << n_vect << "\n";
#endif

  auto span_new_solution =
      [&](PartSolution const &psol,
          int const &newdir) -> std::optional<PartSolution> {
    Face new_set = psol.full_set;
    MyVector<Tint> trans = dv.ListV[newdir] - dv.ListV[psol.vert];
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
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

  auto span_part_solution =
      [&](PartSolution const &psol) -> std::vector<PartSolution> {
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
    os << "PARALL:   span_part_solution |full_set|=" << psol.full_set.size()
       << " / " << psol.full_set.count() << "\n";
#endif
    std::vector<PartSolution> list_sol;
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      if (psol.full_set[i_vect] == 0) {
        std::optional<PartSolution> opt = span_new_solution(psol, i_vect);
        if (opt) {
          list_sol.push_back(*opt);
        }
      }
    }
    return list_sol;
  };
  auto get_empty = [&]() -> PartSolution { return {miss_val, {}, {}}; };
  auto get_all_starts = [&]() -> std::vector<PartSolution> {
    std::vector<PartSolution> l_sol;
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
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
  auto get_initial = [&]() -> OneLevel {
    PartSolution prev_sol = get_empty();
    std::vector<PartSolution> l_sol = get_all_starts();
    size_t choice = 0;
    return {prev_sol, l_sol, choice};
  };
  std::vector<OneLevel> l_levels{get_initial()};
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
  os << "PARALL:   kernel_enumerate_parallelepiped, l_levels\n";
#endif
  int i_level = 0;
  auto GoUpNextInTree = [&]() -> bool {
    while (true) {
      OneLevel &level = l_levels[i_level];
      if (level.choice + 1 < level.l_sol.size()) {
        level.choice += 1;
        return true;
      }
      if (i_level == 0) {
        return false;
      }
      i_level -= 1;
    }
  };
  auto NextInTree = [&]() -> bool {
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
    os << "PARALL:   NextInTree, i_level=" << i_level << "\n";
#endif
    int choice = l_levels[i_level].choice;
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
    os << "PARALL:   NextInTree, choice=" << choice << "\n";
#endif
    PartSolution const &psol = l_levels[i_level].l_sol[choice];
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
    os << "PARALL:   NextInTree, we have psol\n";
#endif
    if (i_level == p) {
      f_insert(psol);
      return GoUpNextInTree();
    } else {
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
      os << "PARALL:   NextInTree, before span_part_solution\n";
#endif
      std::vector<PartSolution> new_sols = span_part_solution(psol);
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
      os << "PARALL:   NextInTree, after span_part_solution |new_sols|="
         << new_sols.size() << "\n";
#endif
      if (new_sols.empty()) {
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
  while (true) {
    bool test = NextInTree();
    if (!test) {
      break;
    }
  }
}

inline int pow_two(int dim) {
  int pow = 1;
  for (int u = 0; u < dim; u++) {
    pow *= 2;
  }
  return pow;
}

template <typename Tint>
std::vector<Face> enumerate_parallelepiped(MyMatrix<Tint> const &M,
                                           std::ostream &os) {
  DataVect<Tint> dv = get_data_vect(M);
  std::unordered_set<Face> set_face;
  int dim = M.cols();
  int pow = pow_two(dim);
  if (dv.n_vect < pow) {
    // No point trying to enumerate when there are no solutions.
    return {};
  }
  MyMatrix<Tint> Mdet(dim, dim);
  auto f_insert = [&](PartSolution const &psol) -> void {
    int e_vert = psol.vert;
    for (int i = 0; i < dim; i++) {
      int f_vert = psol.l_dir[i];
      for (int j = 0; j < dim; j++) {
        Mdet(i, j) = M(f_vert, j) - M(e_vert, j);
      }
    }
    Tint det = DeterminantMat(Mdet);
    if (T_abs(det) == 1) {
      set_face.insert(psol.full_set);
    }
  };
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   Before kernel_enumerate_parallelepiped\n";
#endif
  kernel_enumerate_parallelepiped(dv, dim, f_insert, os);
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   After kernel_enumerate_parallelepiped\n";
#endif
  std::vector<Face> l_face;
  for (auto &eFace : set_face) {
    l_face.push_back(eFace);
  }
  return l_face;
}

template <typename T, typename Tint> struct ResultRobustClosest {
  T robust_minimum;
  std::vector<MyMatrix<Tint>> list_parallelepipeds;
};

template <typename T, typename Tint>
T compute_upper_bound_mat(MyMatrix<T> const &GramMat, MyMatrix<Tint> const &M) {
  int n_ent = M.rows();
  T upper_value(0);
  for (int i_ent = 0; i_ent < n_ent; i_ent++) {
    for (int j_ent = i_ent + 1; j_ent < n_ent; j_ent++) {
      MyVector<Tint> v1 = GetMatrixRow(M, i_ent);
      MyVector<Tint> v2 = GetMatrixRow(M, j_ent);
      MyVector<Tint> diff = v1 - v2;
      T norm = EvaluationQuadForm(GramMat, diff);
      if (norm > upper_value) {
        upper_value = norm;
      }
    }
  }
  return upper_value;
}

template <typename T, typename Tint>
T compute_upper_bound_rrc(MyMatrix<T> const &GramMat,
                          ResultRobustClosest<T, Tint> const &rrc) {
  T upper_value(0);
  for (auto &M : rrc.list_parallelepipeds) {
    T value = compute_upper_bound_mat(GramMat, M);
    if (value == 0) {
      std::cerr << "PARALL: The value should be non-zero\n";
      throw TerminalException{1};
    }
    if (upper_value == 0) {
      // Done at first step
      upper_value = value;
    } else {
      if (value < upper_value) {
        upper_value = value;
      }
    }
  }
  return upper_value;
}

template <typename T, typename Tint> struct ResultDirectEnumeration {
  T min;
  std::vector<MyMatrix<Tint>> list_min_parallelepipeds;
  std::vector<MyMatrix<Tint>> tot_list_parallelepipeds;
};

template <typename T, typename Tint>
std::optional<ResultDirectEnumeration<T, Tint>>
compute_and_enumerate_structures(MyMatrix<T> const &GramMat,
                                 resultCVP<T, Tint> const &res_cvp,
                                 MyVector<T> const &eV, std::ostream &os) {
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   compute_and_enumerate_structures, beginning, eV="
     << StringVector(eV) << "\n";
#endif
  int dim = eV.size();
  int pow = pow_two(dim);
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   compute_and_enumerate_structures, dim=" << dim
     << " pow=" << pow << "\n";
#endif
  MyMatrix<Tint> M_sol(pow, dim);
  auto get_msol = [&](MyMatrix<Tint> const &Min,
                      Face const &eFace) -> MyMatrix<Tint> {
    int pos = 0;
    for (int &vert : FaceToVector<int>(eFace)) {
      for (int i = 0; i < dim; i++) {
        M_sol(pos, i) = Min(vert, i);
      }
      pos += 1;
    }
    return M_sol;
  };
  int n_vect = res_cvp.ListVect.rows();
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   After solver.at_most_dist_vectors res_cvp.TheNorm="
     << res_cvp.TheNorm << " |ListVect|=" << n_vect << "\n";
#endif
  std::vector<T> l_norm;
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<Tint> fV = GetMatrixRow(res_cvp.ListVect, i_vect);
    MyVector<T> diff = UniversalVectorConversion<T, Tint>(fV) - eV;
    T norm = EvaluationQuadForm(GramMat, diff);
    l_norm.push_back(norm);
  }
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   Before enumerate_parallelepiped\n";
#endif
  std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect, os);
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   After enumerate_parallelepiped |l_face|=" << l_face.size()
     << "\n";
#endif
  if (!l_face.empty()) {
    T eff_min = res_cvp.TheNorm + T(1);
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL:   enumerating, eff_min=" << eff_min << "\n";
    int i_face = 0;
#endif
    std::vector<MyMatrix<Tint>> list_min_parallelepipeds;
    std::vector<MyMatrix<Tint>> tot_list_parallelepipeds;
    for (auto &eFace : l_face) {
      T local_max_norm(0);
#ifdef DEBUG_ENUM_PARALL_SEARCH
      os << "PARALL:   i_face=" << i_face << " eFace=" << eFace << "\n";
#endif
      for (int &vert : FaceToVector<int>(eFace)) {
        if (l_norm[vert] > local_max_norm) {
          local_max_norm = l_norm[vert];
        }
      }
#ifdef DEBUG_ENUM_PARALL_SEARCH
      os << "PARALL:   i_face=" << i_face
         << " local_max_norm=" << local_max_norm << "\n";
      i_face += 1;
#endif
      MyMatrix<Tint> Mparall = get_msol(res_cvp.ListVect, eFace);
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
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL:   eff_min=" << eff_min << "\n";
#endif
    ResultDirectEnumeration<T, Tint> rde{eff_min, list_min_parallelepipeds,
                                         tot_list_parallelepipeds};
    return rde;
  } else {
    return {};
  }
}

template <typename T, typename Tint, typename Finsert>
void compute_robust_close_f(CVPSolver<T, Tint> const &solver,
                            MyVector<T> const &eV, Finsert f_insert,
                            std::ostream &os) {
  int dim = solver.GramMat.rows();
#ifdef SANITY_CHECK_ENUM_PARALL_SEARCH
  int dim2 = eV.size() - 1;
  if (dim != dim2) {
    std::cerr << "We have dim=" << dim << " and dim2=" << dim2 << "\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> eV_red(dim);
  for (int i=0; i<dim; i++) {
    eV_red(i) = eV(i + 1);
  }
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL: compute_robust_close_f, step 1\n";
  int n_iter = 0;
#endif
  std::optional<T> opt_norm;
  while (true) {
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL: compute_robust_close_f, step 2, n_iter=" << n_iter << "\n";
    n_iter += 1;
#endif
    resultCVP<T, Tint> res_cvp = solver.increase_distance_vectors(eV_red, opt_norm);
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL: compute_robust_close_f, we have res_cvp\n";
#endif
    std::optional<ResultDirectEnumeration<T, Tint>> opt_rde =
        compute_and_enumerate_structures(solver.GramMat, res_cvp, eV_red, os);
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL: compute_robust_close_f, After "
          "compute_and_enumerate_structures\n";
#endif
    opt_norm = res_cvp.TheNorm;
    if (opt_rde) {
#ifdef DEBUG_ENUM_PARALL_SEARCH
      os << "PARALL: compute_robust_close_f, opt_red.is_some()\n";
#endif
      ResultDirectEnumeration<T, Tint> const &rde = *opt_rde;
      bool test = f_insert(rde);
      if (test) {
        return;
      }
    } else {
#ifdef DEBUG_ENUM_PARALL_SEARCH
      os << "PARALL: compute_robust_close_f, opt_red.is_none()\n";
#endif
    }
  }
}

// Find the robust closest minimum with the lambda expression.
template <typename T, typename Tint>
ResultRobustClosest<T, Tint>
compute_robust_closest(CVPSolver<T, Tint> const &solver, MyVector<T> const &eV,
                       std::ostream &os) {
  ResultRobustClosest<T, Tint> result;
  auto f_insert = [&](ResultDirectEnumeration<T, Tint> const &rde) {
    result = {rde.min, rde.list_min_parallelepipeds};
    return true;
  };
  compute_robust_close_f(solver, eV, f_insert, os);
  return result;
}

template <typename T> MyVector<T> get_random_vector(int denom, int dim) {
  MyVector<T> eV(1 + dim);
  eV(0) = T(1);
  T denom_T(denom);
  for (int i = 0; i < dim; i++) {
    int val1 = random();
    int val = val1 % denom;
    T val_T(val);
    T quot = val_T / denom_T;
    eV(i + 1) = quot;
  }
  return eV;
}

template <typename T, typename Tint>
T random_estimation_robust_covering(MyMatrix<T> const &GramMat, size_t n_iter,
                                    std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  int dim = GramMat.rows();
  T max_cov(0);
  MyVector<T> eV_best;
  auto f_update=[&](MyVector<T> const& fV) -> void {
    ResultRobustClosest<T, Tint> rrc =
        compute_robust_closest<T, Tint>(solver, fV, os);
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL: After compute_robust_closest\n";
#endif
    if (rrc.robust_minimum > max_cov) {
      eV_best = fV;
      max_cov = rrc.robust_minimum;
    }
  };
  for (size_t iter = 0; iter < n_iter; iter++) {
    int denom = random() % 1000000 + 1;
    MyVector<T> eV = get_random_vector<T>(denom, dim);
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL: Before compute_robust_closest eV=" << StringVectorGAP(eV)
       << " denom=" << denom << "\n";
#endif
    f_update(eV);
  }
  MyVector<T> eV = ZeroVector<T>(1 + dim);
  eV(0) = 1;
  f_update(eV);
#ifdef PRINT_ENUM_PARALL_SEARCH
  os << "PARALL: random_estimation_robust_covering eV_best=" << StringVectorGAP(eV_best) << "\n";
#endif
  return max_cov;
}

// clang-format off
#endif  // SRC_ROBUST_COVERING_PARALL_SEARCH_H_
// clang-format on


