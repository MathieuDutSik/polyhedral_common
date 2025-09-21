// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_


/*
  We want to enumerate the possible parallelepipeds
  of fixed dimension.
  
 */
template<typename Tint>
std::vector<Face> enumerate_parallelepiped(MyMatrix<Tint> const& M, int const& p) {
  int miss_val = std::numeric_limits<int>::max();
  int n_vect = M.rows();
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(M, i_vect);
    map[V] = i_vect;
    ListV.push_back(V);
  }
  std::vector<MyMatrix<int>> ListMset;
  for (int i=0; i<=p; i++) {
    MyMatrix<int> Mset = BuildSet(i, 2);
    ListMset.push(Mset);
  }
  struct PartSolution {
    int vert;
    std::vector<int> l_dir;
    Face full_set;
  };


  auto span_new_solution=[&](PartSolution const& psol, int const& newdir) -> std::optional<PartSolution> {
    Face new_set = psol.full_set;
    MyVector<Tint> trans = ListV[newdir] - ListV[psol.vert];
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      if (psol.full_set[i_vect] == 1) {
        MyVector<Tint> newV = ListV[i_vect] + trans;
        auto iter = map.find(newV);
        if (iter == map.end()) {
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
  std::vector<Face> l_face;
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
    if (i_level == p) {
      l_face.push_back(l_levels[p].full_set);
      return GoUpNextInTree();
    } else {
      int choice = l_levels[i_level].choice;
      PartSolution const& psol = l_levels[i_level].l_sol[choice];
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






// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
