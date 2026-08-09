// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_PARALL_SEARCH_H_
#define SRC_ROBUST_COVERING_PARALL_SEARCH_H_

// clang-format off
#include "Shvec_exact.h"
#include "LatticeDelaunay.h"
#include <algorithm>
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
  // Indices of the points of the partial parallelepiped, sorted and without
  // repetition. Keeping the points as a sorted list of indices rather than as
  // a subset of the whole family lets a level cost the size of the
  // parallelepiped instead of the size of the family.
  std::vector<int> pts;
};

/*
  The enumeration asks a single question about the family of points: given a
  point of the family and a translation, which point of the family is its
  image, if any. Answering it on the multiprecision coordinates, through a
  hash of the whole vector, was the dominant cost of the whole program, so the
  points are relocated on a bounded integer grid whenever their coordinates
  allow it and the query becomes a few integer operations and one array
  access. MapPointLocator keeps the former behaviour for the families that do
  not fit on such a grid.
 */

// Above that many cells the bounding box is too sparse for the direct
// addressing to be worth its memory and the map based locator is used.
inline constexpr int64_t max_ncell_grid_locator = 1 << 22;

struct GridPointLocator {
  // A translation is a difference of grid coordinates.
  using Trans = std::vector<int64_t>;
  int n_vect;
  int dim;
  // Coordinates relative to the corner of the bounding box, n_vect blocks of
  // dim entries.
  std::vector<int64_t> coord;
  std::vector<int64_t> span;
  std::vector<int64_t> stride;
  // Index of the point sitting on a cell, -1 when the cell is empty.
  std::vector<int> grid;
  Trans get_trans() const { return Trans(dim); }
  void set_trans(int const &i_add, int const &i_sub, Trans &t) const {
    int64_t const *c_add = coord.data() + static_cast<size_t>(i_add) * dim;
    int64_t const *c_sub = coord.data() + static_cast<size_t>(i_sub) * dim;
    for (int k = 0; k < dim; k++) {
      t[k] = c_add[k] - c_sub[k];
    }
  }
  int translate(int const &i_pt, Trans const &t) const {
    int64_t const *c_pt = coord.data() + static_cast<size_t>(i_pt) * dim;
    int64_t flat = 0;
    for (int k = 0; k < dim; k++) {
      int64_t val = c_pt[k] + t[k];
      if (val < 0 || val >= span[k]) {
        return -1;
      }
      flat += val * stride[k];
    }
    return grid[flat];
  }
};

// Returns nothing when the coordinates do not fit in int64_t or when the
// bounding box needs too many cells.
template <typename Tint>
std::optional<GridPointLocator>
get_grid_point_locator(MyMatrix<Tint> const &M) {
  size_t n_vect = M.rows();
  size_t dim = M.cols();
  if (n_vect == 0) {
    return {};
  }
  std::vector<int64_t> coord(n_vect * dim);
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    for (size_t k = 0; k < dim; k++) {
      Tint const &val = M(i_vect, k);
      std::optional<int64_t> opt =
          UniversalScalarConversionCheck<int64_t, Tint>(val);
      if (!opt) {
        return {};
      }
      // Converting back and comparing keeps the test independent of whether
      // the conversion of Tint reports its overflows.
      if (UniversalScalarConversion<Tint, int64_t>(*opt) != val) {
        return {};
      }
      coord[i_vect * dim + k] = *opt;
    }
  }
  std::vector<int64_t> span(dim), stride(dim);
  int64_t n_cell = 1;
  for (size_t k = 0; k < dim; k++) {
    int64_t min_val = coord[k];
    int64_t max_val = coord[k];
    for (size_t i_vect = 1; i_vect < n_vect; i_vect++) {
      int64_t val = coord[i_vect * dim + k];
      min_val = std::min(min_val, val);
      max_val = std::max(max_val, val);
    }
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      coord[i_vect * dim + k] -= min_val;
    }
    span[k] = max_val - min_val + 1;
    stride[k] = n_cell;
    if (span[k] > max_ncell_grid_locator / n_cell) {
      return {};
    }
    n_cell *= span[k];
  }
  std::vector<int> grid(n_cell, -1);
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    int64_t flat = 0;
    for (size_t k = 0; k < dim; k++) {
      flat += coord[i_vect * dim + k] * stride[k];
    }
    grid[flat] = static_cast<int>(i_vect);
  }
  return GridPointLocator{static_cast<int>(n_vect), static_cast<int>(dim),
                          std::move(coord),         std::move(span),
                          std::move(stride),        std::move(grid)};
}

template <typename Tint> struct MapPointLocator {
  using Trans = MyVector<Tint>;
  int n_vect;
  int dim;
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
  Trans get_trans() const { return Trans(dim); }
  void set_trans(int const &i_add, int const &i_sub, Trans &t) const {
    t = ListV[i_add] - ListV[i_sub];
  }
  int translate(int const &i_pt, Trans const &t) const {
    MyVector<Tint> V = ListV[i_pt] + t;
    auto iter = map.find(V);
    if (iter == map.end()) {
      return -1;
    }
    return iter->second;
  }
};

template <typename Tint>
MapPointLocator<Tint> get_map_point_locator(MyMatrix<Tint> const &M) {
  int n_vect = M.rows();
  int dim = M.cols();
  std::vector<MyVector<Tint>> ListV;
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(M, i_vect);
    map[V] = i_vect;
    ListV.push_back(V);
  }
  return {n_vect, dim, std::move(ListV), std::move(map)};
}

/*
  We want to enumerate the possible parallelepipeds
  of fixed dimension.
  It is a simple tree search
 */
template <typename Tlocator, typename Finsert>
void kernel_enumerate_parallelepiped(Tlocator const &loc, int const &p,
                                     Finsert f_insert,
                                     [[maybe_unused]] std::ostream &os) {
  int n_vect = loc.n_vect;
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   kernel_enumerate_parallelepiped, n_vect=" << n_vect << "\n";
#endif
  // Reused across the whole search so that no allocation happens per node.
  typename Tlocator::Trans trans = loc.get_trans();
  std::vector<int> new_pts;

  auto span_new_solution =
      [&](PartSolution const &psol,
          int const &newdir) -> std::optional<PartSolution> {
    loc.set_trans(newdir, psol.vert, trans);
    size_t n_old = psol.pts.size();
    new_pts = psol.pts;
    for (auto &i_pt : psol.pts) {
      int pos = loc.translate(i_pt, trans);
      if (pos == -1) {
        return {};
      }
      new_pts.push_back(pos);
    }
    // The first half is already sorted, so only the image needs to be.
    std::sort(new_pts.begin() + n_old, new_pts.end());
    std::inplace_merge(new_pts.begin(), new_pts.begin() + n_old,
                       new_pts.end());
    new_pts.erase(std::unique(new_pts.begin(), new_pts.end()), new_pts.end());
    std::vector<int> l_dir = psol.l_dir;
    l_dir.push_back(newdir);
    return PartSolution{psol.vert, std::move(l_dir), new_pts};
  };

  /*
    The box spanned from psol.vert by the directions of psol.l_dir depends only
    on the set of directions, not on the order in which they were added. Every
    sub-box of a valid box is itself a valid partial solution, so restricting
    the search to strictly increasing sequences of directions still reaches
    every (vert, set of directions) pair, and reaches it exactly once instead of
    once per permutation. That divides the size of the tree by dim!.
   */
  auto span_part_solution =
      [&](PartSolution const &psol) -> std::vector<PartSolution> {
#ifdef DEBUG_ENUM_PARALL_SEARCH_DISABLE
    os << "PARALL:   span_part_solution |pts|=" << psol.pts.size() << "\n";
#endif
    std::vector<PartSolution> list_sol;
    int i_start = 0;
    if (!psol.l_dir.empty()) {
      i_start = psol.l_dir.back() + 1;
    }
    // psol.pts is sorted, so the points already in the parallelepiped are
    // skipped with a cursor instead of a membership test.
    size_t i_pts = 0;
    while (i_pts < psol.pts.size() && psol.pts[i_pts] < i_start) {
      i_pts++;
    }
    for (int i_vect = i_start; i_vect < n_vect; i_vect++) {
      if (i_pts < psol.pts.size() && psol.pts[i_pts] == i_vect) {
        i_pts++;
        continue;
      }
      std::optional<PartSolution> opt = span_new_solution(psol, i_vect);
      if (opt) {
        list_sol.push_back(std::move(*opt));
      }
    }
    return list_sol;
  };
  auto get_all_starts = [&]() -> std::vector<PartSolution> {
    std::vector<PartSolution> l_sol;
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      PartSolution esol{i_vect, {}, {i_vect}};
      l_sol.push_back(std::move(esol));
    }
    return l_sol;
  };
  struct OneLevel {
    std::vector<PartSolution> l_sol;
    size_t choice;
  };
  auto get_initial = [&]() -> OneLevel {
    std::vector<PartSolution> l_sol = get_all_starts();
    size_t choice = 0;
    return {std::move(l_sol), choice};
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
      OneLevel new_level{std::move(new_sols), new_choice};
      int new_i_level = i_level + 1;
      if (l_levels.size() >= static_cast<size_t>(new_i_level + 1)) {
        l_levels[new_i_level] = std::move(new_level);
      } else {
        l_levels.push_back(std::move(new_level));
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
  int n_vect = M.rows();
  int dim = M.cols();
  int pow = pow_two(dim);
  if (n_vect < pow) {
    // No point trying to enumerate when there are no solutions.
    return {};
  }
  std::unordered_set<Face> set_face;
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
      // The subset is materialized only for the parallelepipeds that are
      // kept, which are a small part of the leaves of the search.
      Face full_set(n_vect);
      for (auto &i_pt : psol.pts) {
        full_set[i_pt] = 1;
      }
      set_face.insert(std::move(full_set));
    }
  };
#ifdef DEBUG_ENUM_PARALL_SEARCH
  os << "PARALL:   Before kernel_enumerate_parallelepiped\n";
#endif
  std::optional<GridPointLocator> opt_grid = get_grid_point_locator(M);
  if (opt_grid) {
    kernel_enumerate_parallelepiped(*opt_grid, dim, f_insert, os);
  } else {
#ifdef DEBUG_ENUM_PARALL_SEARCH
    os << "PARALL:   The points do not fit on a grid, using the map\n";
#endif
    kernel_enumerate_parallelepiped(get_map_point_locator(M), dim, f_insert,
                                    os);
  }
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

/*
  The largest squared distance between two vertices of a parallelepiped.

  The quadratic form is evaluated for every pair of vertices, so the cost grows
  with the square of the number of vertices and, the vertices being the 2^dim
  corners of a box, quickly dominates. The evaluation is therefore done over
  the underlying ring: the Gram matrix is scaled once to an integral matrix,
  every norm is an integer numerator over the same positive denominator, and
  only the returned value goes back to the field. That replaces the canonical-
  izing field arithmetic of the fractions by plain multiply-accumulate.
 */
template <typename T> struct UpperBoundEvaluator {
  using Tring = typename underlying_ring<T>::ring_type;
  // ScaledGram = TheMult * GramMat with TheMult > 0, so the ordering of the
  // norms is the ordering of the integral values.
  MyMatrix<Tring> ScaledGram;
  T TheMult;
  int dim;
  template <typename Tint> T eval(MyMatrix<Tint> const &M) const {
    int n_ent = M.rows();
    MyMatrix<Tring> Mring = UniversalMatrixConversion<Tring, Tint>(M);
    MyVector<Tring> diff(dim);
    Tring upper_value(0);
    Tring norm, prod;
    for (int i_ent = 0; i_ent < n_ent; i_ent++) {
      for (int j_ent = i_ent + 1; j_ent < n_ent; j_ent++) {
        for (int i = 0; i < dim; i++) {
          diff(i) = Mring(i_ent, i) - Mring(j_ent, i);
        }
        norm = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            prod = diff(i) * diff(j);
            AddMul(norm, prod, ScaledGram(i, j));
          }
        }
        if (norm > upper_value) {
          upper_value = norm;
        }
      }
    }
    return UniversalScalarConversion<T, Tring>(upper_value) / TheMult;
  }
};

template <typename T>
UpperBoundEvaluator<T> get_upper_bound_evaluator(MyMatrix<T> const &GramMat) {
  FractionMatrixRing<T> frr = RemoveFractionMatrixPlusCoeffRing(GramMat);
  int dim = GramMat.rows();
  return UpperBoundEvaluator<T>{std::move(frr.TheMat), std::move(frr.TheMult),
                                dim};
}

template <typename T, typename Tint>
T compute_upper_bound_mat(MyMatrix<T> const &GramMat, MyMatrix<Tint> const &M) {
  return get_upper_bound_evaluator(GramMat).eval(M);
}

template <typename T, typename Tint>
T compute_upper_bound_rrc(MyMatrix<T> const &GramMat,
                          ResultRobustClosest<T, Tint> const &rrc) {
  // The scaling of the Gram matrix is shared by all the parallelepipeds.
  UpperBoundEvaluator<T> evaluator = get_upper_bound_evaluator(GramMat);
  T upper_value(0);
  for (auto &M : rrc.list_parallelepipeds) {
    T value = evaluator.eval(M);
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
      bool test = f_insert(rde, res_cvp.TheNorm);
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
  auto f_insert = [&](ResultDirectEnumeration<T, Tint> const &rde, [[maybe_unused]] T const& TheNorm) {
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


