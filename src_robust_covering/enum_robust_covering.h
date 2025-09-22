// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_


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


template<typename Tint>
std::vector<Face> enumerate_parallelepiped(MyMatrix<Tint> const& M) {
  DataVect<Tint> dv = get_data_vect(M);
  std::vector<Face> l_face;
  int dim = M.cols();
  int pow = 1;
  for (int u=0; u<dim; u++) {
    pow *= 2;
  }
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







// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
