// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "ShortestUniversal.h"
// clang-format on



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
}

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
void kernel_enumerate_parallelepiped(DataVect<Tint> const& dv, int const& p, Finsert f_insert) {
  int miss_val = std::numeric_limits<int>::max();

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
        int pos = *iter;
        new_set[pos] = 1;
      }
    }
    std::vector<int> l_dir = psol.l_dir;
    l_dir.push_back(newdir);
    PartSolution newsol{psol.vert, std::move(l_dir), std::move(new_set)};
    return newsol;
  };

  auto span_part_solution=[&](PartSolution const& psol) -> std::vector<PartSolution> {
    std::vector<PartSolution> list_sol;
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      if (psol.full_set[i_vect] == 0) {
        std::option<PartSolution> opt = span_new_solution(psol, i_vect);
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
      PartSolution esol{i_vect, {}, {}};
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
    int choice = l_levels[i_level].choice;
    PartSolution const& psol = l_levels[i_level].l_sol[choice];
    if (i_level == p) {
      f_insert(psol);
      return GoUpNextInTree();
    } else {
      std::vector<PartSolution> new_sols = span_part_solution(psol);
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
      return true;
    }
  };
  while(true) {
    bool test = NextInTree();
    if (!test) {
      break;
    }
  }
  return l_face;
}

int pow_two(int dim) {
  int pow = 1;
  for (int u=0; u<dim; u++) {
    pow *= 2;
  }
  return pow;
}

template<typename Tint>
std::vector<Face> enumerate_parallelepiped(MyMatrix<Tint> const& M) {
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
    Tint det = Determinant(Mdet);
    if (T_abs(det) == 1) {
      l_face.push_back(psol.full_set);
    }
  };
  kernel_enumerate_parallelepiped(dv, dim, f_insert);
  return l_face;
}

template<typename T, typename Tint>
struct ResultRobustClosest {
  T robust_minimum;
  std::vector<MyMatrix<Tint>> list_parallelepipeds;
};



template<typename T, typename Tint>
ResultRobustClosest<T,Tint> compute_robust_closest(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV) {
  T min(0);
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
  while(true) {
    if (n_iter == 0) {
      resultCVP<T, Tint> res_cvp = solve.SingleSolver(eV);
      min = res_cvp.TheNorm;
      std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect);
      if (l_face.size() > 0) {
        std::vector<MyMatrix<Tint>> list_parallelepipeds;
        for (auto & eFace: l_face) {
          list_parallelepipeds.push_back(get_msol(res_cvp.ListVect, eFace));
        }
        return {min, list_parallelepipeds};
      }
    } else {
      min = (min * T(3)) / T(2);
      std::vector<MyVector<Tint>> elist = soler.AtMostNormVectors(eV, min);
      std::vector<T> l_norm;
      for (auto & fV: elist) {
        MyVector<T> diff = UniversalVectorConversion<T,Tint>(fV) - eV;
        T norm = EvaluationQuadForm(solver.GramMat, diff);
        l_norm.push_back(norm);
      }
      MyMatrix<T> M = MatrixFromVectorFamily(elist);
      std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect);
      if (l_face.size() > 0) {
        T eff_min = min + T(1);
        std::vector<MyMatrix<Tint>> list_parallelepipeds;
        for (auto & eFace: l_face) {
          T local_max_norm(0);
          for (int& vert: FaceToVector<int>(eFace)) {
            if (l_norm[vert] > local_max_norm) {
              local_max_norm = l_norm[vert];
            }
          }
          MyMatrix<Tint> Mparall = get_msol(M, eFace);
          if (local_max_norm < eff_min) {
            list_parallelepipeds.clear();
            eff_min = local_max_norm;
            list_parallelepipeds.push_back(Mparall);
          } else {
            if (local_max_norm == eff_min) {
              list_parallelepipeds.push_back(Mparall);
            }
          }
        }
        return {eff_min, list_parallelepipeds};
      }
    }
  }
}



template<typename T, typename Tint>
T random_estimation_robust_covering(MyMatrix<T> const& GramMat, size_t n_iter, std::ostream & os) {
  CVPSolver<T,Tint> solver(GramMat, os);
  int dim = GramMat.rows();
  T max_cov(0);
  MyMatrix<T> eV(dim);
  for (size_t iter=0; iter<n_iter; iter++) {
    int denom = random() % 1000000000000000;
    T denom_T(denom);
    for (int i=0; i<dim; i++) {
      int val = random() % denom;
      T val_T(val);
      eV(i) = val_T / denom_T;
    }
    ResultRobustClosest<T,Tint> rrc = compute_robust_closest<T,Tint>(solver, eV);
    if (rrc.robust_minimum > max_cov) {
      max_cov = rrc.robust_minimum;
    }
  }
  return max_cov;
}







// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
