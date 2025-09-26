// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
#define SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_

// clang-format off
#include "ShortestUniversal.h"
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
    take into account 

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



template<typename T, typename Tint>
ResultRobustClosest<T,Tint> compute_robust_closest(CVPSolver<T,Tint> const& solver, MyVector<T> const& eV, [[maybe_unused]] std::ostream& os) {
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
      os << "ROBUST:   Before solver.SingleSolver\n";
#endif
      resultCVP<T, Tint> res_cvp = solver.SingleSolver(eV);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   After solver.SingleSolver\n";
#endif
      min_search = res_cvp.TheNorm;
      std::vector<Face> l_face = enumerate_parallelepiped(res_cvp.ListVect, os);
      if (l_face.size() > 0) {
        std::vector<MyMatrix<Tint>> list_parallelepipeds;
        for (auto & eFace: l_face) {
          list_parallelepipeds.push_back(get_msol(res_cvp.ListVect, eFace));
        }
        return {min_search, std::move(list_parallelepipeds)};
      }
    } else {
      min_search = (min_search * T(3)) / T(2);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   Before solver.AtMostNormVectors\n";
#endif
      std::vector<MyVector<Tint>> elist = solver.AtMostNormVectors(eV, min_search);
#ifdef DEBUG_ENUM_ROBUST_COVERING
      os << "ROBUST:   After solver.AtMostNormVectors |elist|=" << elist.size() << "\n";
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
        std::vector<MyMatrix<Tint>> list_parallelepipeds;
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
#ifdef DEBUG_ENUM_ROBUST_COVERING
        os << "ROBUST:   eff_min=" << eff_min << "\n";
#endif
        return {eff_min, std::move(list_parallelepipeds)};
      }
    }
    n_iter += 1;
  }
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







// clang-format off
#endif  // SRC_ROBUST_COVERING_ENUM_ROBUST_COVERING_H_
// clang-format on
