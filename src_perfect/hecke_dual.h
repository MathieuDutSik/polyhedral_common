// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_DUAL_H_
#define SRC_PERFECT_HECKE_DUAL_H_

// clang-format off
#include "hecke_chain_map.h"
#include "SparseChainSolver.h"
// clang-format on

/*
  The "Dual" method of step A: induction from the perfect cones in the
  well rounded complex.
 */

//
// The dual construction: induction from the perfect cones
//

/*
  The induction from the vertices goes through the cells at infinity,
  whose stabilizers are infinite and whose stars are infinite; this is
  what makes the coverings explode for GL_4(Z) and GL_3(O_K). The dual
  construction runs the induction in the other direction, on the well
  rounded cells only: the images of the perfect cones (the vertices of
  the well rounded retract W) are found by the walk to the image of the
  barycenter; a facet (an edge of W, between two cones) is mapped to a
  path of adjacent cones; a cell of codimension k is mapped to a chain
  of cells of codimension k whose "coboundary" (the signed sum of the
  cells of codimension k-1 containing it) equals the image of its own
  coboundary. The chain complex is the transposed one, the cellular chain
  complex of W, and everything in it is finite. The harmonic
  representatives of the homology are the same as for the original
  complex, so step B applies unchanged.
 */

// The matrix A of the action Q -> M Q M^T on the T-space: the image of
// ListMat[k] is sum_l A(k,l) ListMat[l].
template <typename T>
MyMatrix<T> tspace_action_matrix(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &M) {
  int n_mat = LinSpa.ListMat.size();
  MyMatrix<T> A(n_mat, n_mat);
  for (int k = 0; k < n_mat; k++) {
    MyMatrix<T> Qimg = M * LinSpa.ListMat[k] * M.transpose();
    MyVector<T> V = SymmetricMatrixToVector(Qimg);
    std::optional<MyVector<T>> opt = SolutionMat(LinSpa.ListMatAsBigMat, V);
    if (!opt) {
      std::cerr << "HECKE: The matrix does not preserve the T-space\n";
      throw TerminalException{1};
    }
    for (int l = 0; l < n_mat; l++) {
      A(k, l) = (*opt)(l);
    }
  }
  return A;
}

// The sign of the determinant of the action Q -> M Q M^T on the T-space.
template <typename T>
int tspace_orientation_sign(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &M) {
  MyMatrix<T> A = tspace_action_matrix(LinSpa, M);
  T det = DeterminantMat(A);
  if (det > 0) {
    return 1;
  }
  if (det < 0) {
    return -1;
  }
  std::cerr << "HECKE: The action on the T-space is degenerate\n";
  throw TerminalException{1};
}

/*
  The cofaces of the representative cells, with the sign of the cell in
  the boundary of the coface, cached per orbit. Only for well rounded
  cells (finite stabilizers).
 */
template <typename T, typename Tint, typename Tgroup> struct UpperFaceCache {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  std::ostream &os;
  // cofaces[index][iOrb] : (coface (jOrb, P) at level index-1, sign)
  std::map<std::pair<int, int>, std::vector<std::pair<PerfectFace<Tint>, int>>> cofaces;
  // stars[index][iOrb] : the perfect cones containing the representative
  std::map<std::pair<int, int>, std::vector<PerfectFace<Tint>>> stars;
  UpperFaceCache(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, std::ostream &_os)
      : fce(_fce), os(_os) {}
  std::vector<std::pair<PerfectFace<Tint>, int>> const &get_cofaces(int const &index, int const &iOrb) {
    std::pair<int, int> key{index, iOrb};
    auto iter = cofaces.find(key);
    if (iter != cofaces.end()) {
      return iter->second;
    }
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    MyMatrix<Tint> key_cell = tot_set(face.EXT);
    std::pair<std::vector<MyMatrix<Tint>>, std::vector<PerfectFace<Tint>>> pair =
        get_all_upper_faces(fce, index, iOrb, os);
    std::vector<std::pair<PerfectFace<Tint>, int>> l_ret;
    for (auto &pf : pair.second) {
      // The sign of F_iOrb in d(F_jOrb P)
      int sign_tot = 0;
      int n_match = 0;
      for (auto &e : fce.boundaries[index - 1].ll_bound[pf.iOrb].l_bound) {
        if (e.iOrb == iOrb) {
          MyMatrix<Tint> Mcell = e.M * pf.M;
          MyMatrix<Tint> EXT_prod = face.EXT * Mcell;
          MyMatrix<Tint> EXT2 = tot_set(EXT_prod);
          if (EXT2 == key_cell) {
            MyMatrix<Tint> t = Inverse(Mcell);
            int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
            sign_tot += e.sign * sign;
            n_match += 1;
          }
        }
      }
      if (n_match != 1) {
        std::cerr << "HECKE: The cell should appear exactly once in the boundary of its coface, n_match=" << n_match << "\n";
        throw TerminalException{1};
      }
      l_ret.emplace_back(pf, sign_tot);
    }
    auto ret = cofaces.emplace(key, std::move(l_ret));
    return ret.first->second;
  }
  std::vector<PerfectFace<Tint>> const &get_star(int const &index, int const &iOrb) {
    std::pair<int, int> key{index, iOrb};
    auto iter = stars.find(key);
    if (iter != stars.end()) {
      return iter->second;
    }
    std::vector<PerfectFace<Tint>> l_star;
    if (index == 0) {
      int n = fce.pctdi.LinSpa.n;
      l_star.push_back({iOrb, IdentityMat<Tint>(n)});
    } else {
      std::unordered_set<MyMatrix<Tint>> set_key;
      for (auto &pair : get_cofaces(index, iOrb)) {
        std::vector<PerfectFace<Tint>> const &l_up = get_star(index - 1, pair.first.iOrb);
        // A copy: l_up is a reference into the map that the recursive
        // calls may not invalidate (std::map), but we keep it simple.
        std::vector<PerfectFace<Tint>> l_up_copy = l_up;
        for (auto &pf : l_up_copy) {
          MyMatrix<Tint> M = pf.M * pair.first.M;
          MyMatrix<Tint> EXT_cone = fce.levels[0].l_faces[pf.iOrb].EXT * M;
          MyMatrix<Tint> key_cone = tot_set(EXT_cone);
          if (set_key.insert(key_cone).second) {
            l_star.push_back({pf.iOrb, M});
          }
        }
      }
    }
    auto ret = stars.emplace(key, std::move(l_star));
    return ret.first->second;
  }
};

// The coboundary of a chain: the signed sum of the cofaces.
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
chain_coboundary(int const &index, std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                 UpperFaceCache<T, Tint, Tgroup> &ufc, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  ChainBuilder<T, Tint, Tgroup> cb(index - 1, fce, os);
  for (auto &fe : chain) {
    for (auto &pair : ufc.get_cofaces(index, fe.iOrb)) {
      T value = fe.value * T(pair.second);
      MyMatrix<Tint> M = pair.first.M * fe.M;
      cb.f_insert(value, pair.first.iOrb, M);
    }
  }
  return cb.get_faces();
}

/*
  The generic potential used by the peeling: a positive definite form of
  the T-space, the Gram matrix of the first perfect form slightly
  perturbed (so that the vectors of a cell have distinct heights).
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> get_peeling_potential(FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  MyMatrix<T> Q0 = fce.pctdi.l_perfect[0].gram;
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  uint64_t seed = 987654321;
  auto next_random = [&]() -> int {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((seed >> 33) % 11) - 5;
  };
  for (int i_try = 0; i_try < 20; i_try++) {
    MyMatrix<T> Q = Q0;
    T scale = T(1) / T(1000 * (i_try + 1));
    for (auto &M : ListMat) {
      int r = next_random();
      Q += (scale * T(r)) * M;
    }
    if (IsPositiveDefinite(Q, os)) {
      return Q;
    }
  }
  return Q0;
}

/*
  The peeling solver: the chain c (level k-1) is reduced by repeatedly
  taking its highest cell G (for the potential, the largest height of a
  vector of G) and subtracting the multiple of the coboundary of a level
  k face E of G that cancels G; the face is the one that gives the
  lexicographically smallest (maximal height, number of cells at that
  height, number of cells). Each accepted move decreases this objective,
  so the chain is pushed towards the low cells until it vanishes; the
  accumulated faces form the filling. The move fails cleanly when no
  face improves the objective or when a budget is exceeded.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
peel_partial_contracting_homotopy(int const &k, std::vector<PerfectFaceEntry<T, Tint>> const &rhs,
                          UpperFaceCache<T, Tint, Tgroup> &ufc, MyMatrix<T> const &Q0,
                          HeckeBudget const &budget, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  auto cell_height = [&](int const &index, int const &iOrb, MyMatrix<Tint> const &M) -> T {
    MyMatrix<Tint> EXT = fce.levels[index].l_faces[iOrb].EXT * M;
    T h(0);
    for (int i_row = 0; i_row < EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(EXT, i_row);
      T val = EvaluationQuadForm<T, Tint>(Q0, V);
      if (i_row == 0 || val > h) {
        h = val;
      }
    }
    return h;
  };
  struct Objective {
    T max_height;
    size_t n_top;
    size_t n_cell;
    bool operator<(Objective const &o) const {
      if (max_height != o.max_height) {
        return max_height < o.max_height;
      }
      if (n_top != o.n_top) {
        return n_top < o.n_top;
      }
      return n_cell < o.n_cell;
    }
  };
  // The coefficient of the cell (iOrb, M) in a chain.
  auto coefficient_of = [&](std::vector<PerfectFaceEntry<T, Tint>> const &c, int const &index,
                            int const &iOrb, MyMatrix<Tint> const &M) -> T {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    MyMatrix<Tint> EXT_M = face.EXT * M;
    MyMatrix<Tint> key = tot_set(EXT_M);
    for (auto &fe : c) {
      if (fe.iOrb == iOrb) {
        MyMatrix<Tint> EXT_fe = face.EXT * fe.M;
        MyMatrix<Tint> key_fe = tot_set(EXT_fe);
        if (key_fe == key) {
          MyMatrix<Tint> t = fe.M * Inverse(M);
          int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
          return fe.value * T(sign);
        }
      }
    }
    return T(0);
  };
  std::vector<PerfectFaceEntry<T, Tint>> c = rhs;
  std::vector<PerfectFaceEntry<T, Tint>> filling;
  // The heights of the cells of c, the multiset of heights and the map
  // from the cells (as sets of vectors) to their position, refreshed
  // after each accepted move so that the candidates are evaluated
  // incrementally.
  std::vector<T> h_c;
  std::multiset<T> heights;
  std::unordered_map<MyMatrix<Tint>, size_t> map_c;
  auto refresh = [&]() -> void {
    h_c.clear();
    heights.clear();
    map_c.clear();
    for (size_t i = 0; i < c.size(); i++) {
      T h = cell_height(k - 1, c[i].iOrb, c[i].M);
      h_c.push_back(h);
      heights.insert(h);
      MyMatrix<Tint> EXT_c = fce.levels[k - 1].l_faces[c[i].iOrb].EXT * c[i].M;
      map_c[tot_set(EXT_c)] = i;
    }
  };
  auto get_top = [&]() -> int {
    int i_top = -1;
    for (size_t i = 0; i < c.size(); i++) {
      if (i_top == -1 || h_c[i] > h_c[i_top]) {
        i_top = i;
      }
    }
    return i_top;
  };
  auto objective_of_heights = [&](std::multiset<T> const &hs) -> Objective {
    if (hs.empty()) {
      return Objective{T(0), 0, 0};
    }
    T top = *hs.rbegin();
    return Objective{top, hs.count(top), hs.size()};
  };
  refresh();
  Objective obj = objective_of_heights(heights);
  size_t n_move = 0;
  while (!c.empty()) {
    if (n_move >= budget.max_peel_moves) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, budget of moves exhausted, |c|=" << c.size() << "\n";
#endif
      return {};
    }
#ifdef DEBUG_HECKE_OPERATORS
    if (n_move % 100 == 0 && n_move > 0) {
      os << "HECKE: peel, move " << n_move << " |c|=" << c.size() << " max_height=" << obj.max_height << "\n";
    }
#endif
    int i_top = get_top();
    PerfectFaceEntry<T, Tint> const G = c[i_top];
    std::optional<Objective> best_obj;
    std::vector<PerfectFaceEntry<T, Tint>> best_cob;
    PerfectFaceEntry<T, Tint> best_E{0, T(0), G.M};
    for (auto &e : fce.boundaries[k - 1].ll_bound[G.iOrb].l_bound) {
      if (!fce.levels[k].l_faces[e.iOrb].is_well_rounded) {
        continue;
      }
      MyMatrix<Tint> Mface = e.M * G.M;
      std::vector<PerfectFaceEntry<T, Tint>> chain_E{{e.iOrb, T(1), Mface}};
      std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(k, chain_E, ufc, os);
      T s_G = coefficient_of(cob, k - 1, G.iOrb, G.M);
      if (s_G == 0) {
        continue;
      }
      T eps = G.value / s_G;
      // The heights after the move: G disappears, the touched cells of c
      // may vanish, the other cofaces of E appear.
      std::multiset<T> hs = heights;
      hs.erase(hs.find(h_c[i_top]));
      for (auto &fe : cob) {
        FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[k - 1].l_faces[fe.iOrb];
        MyMatrix<Tint> EXT_fe = face.EXT * fe.M;
        MyMatrix<Tint> key = tot_set(EXT_fe);
        auto iter = map_c.find(key);
        if (iter != map_c.end()) {
          size_t idx = iter->second;
          if (static_cast<int>(idx) == i_top) {
            continue;
          }
          MyMatrix<Tint> t = c[idx].M * Inverse(fe.M);
          int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
          T new_coef = c[idx].value - eps * fe.value * T(sign);
          if (new_coef == 0) {
            hs.erase(hs.find(h_c[idx]));
          }
        } else {
          T h = cell_height(k - 1, fe.iOrb, fe.M);
          hs.insert(h);
        }
      }
      Objective obj_new = objective_of_heights(hs);
      if (!best_obj || obj_new < *best_obj) {
        best_obj = obj_new;
        best_cob = cob;
        best_E = PerfectFaceEntry<T, Tint>{e.iOrb, eps, Mface};
      }
    }
    if (!best_obj || !(*best_obj < obj)) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, no improving move at move " << n_move << " |c|=" << c.size()
         << " max_height=" << obj.max_height << " n_top=" << obj.n_top << "\n";
#endif
      return {};
    }
    ChainBuilder<T, Tint, Tgroup> cb(k - 1, fce, os);
    for (auto &fe : c) {
      cb.f_insert(fe.value, fe.iOrb, fe.M);
    }
    for (auto &fe : best_cob) {
      cb.f_insert(-best_E.value * fe.value, fe.iOrb, fe.M);
    }
    c = cb.get_faces();
    if (c.size() > budget.max_chain) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, chain larger than the budget\n";
#endif
      return {};
    }
    filling.push_back(best_E);
    refresh();
    obj = objective_of_heights(heights);
    n_move += 1;
  }
  // The filling, with the duplicates merged
  ChainBuilder<T, Tint, Tgroup> cb(k, fce, os);
  for (auto &fe : filling) {
    cb.f_insert(fe.value, fe.iOrb, fe.M);
  }
  std::vector<PerfectFaceEntry<T, Tint>> sol = cb.get_faces();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: peel, solved in " << n_move << " moves, |sol|=" << sol.size() << "\n";
#endif
  return sol;
}

/*
  The dual partial contracting homotopy: given a chain b of level k-1 (which must
  be a coboundary), find a chain c of level k with coboundary(c) = b.
  The methods are tried in turn: the peeling, then the linear systems on
  the rings of unknowns around b (the faces of its cells first, then the
  cells sharing a coface with a current unknown, and so on) with all the
  cofaces of the unknowns as equations. Each method has its budget and
  the failure is reported cleanly (empty optional).
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
dual_partial_contracting_homotopy(int const &k, std::vector<PerfectFaceEntry<T, Tint>> const &rhs,
                          UpperFaceCache<T, Tint, Tgroup> &ufc, MyMatrix<T> const &Q0,
                          HeckeBudget const &budget, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  if (k > 1) {
    std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_peel =
        peel_partial_contracting_homotopy(k, rhs, ufc, Q0, budget, os);
    if (opt_peel) {
      return opt_peel;
    }
  }
  ComplexBuilder<T, Tint, Tgroup> cb_k(k, fce, os);
  auto insert_faces_of = [&](int const &iOrb_up, MyMatrix<Tint> const &M_up) -> void {
    for (auto &e : fce.boundaries[k - 1].ll_bound[iOrb_up].l_bound) {
      if (fce.levels[k].l_faces[e.iOrb].is_well_rounded) {
        MyMatrix<Tint> M = e.M * M_up;
        cb_k.f_insert(e.iOrb, M);
      }
    }
  };
  for (auto &fe : rhs) {
    insert_faces_of(fe.iOrb, fe.M);
  }
  size_t start = 0;
  for (size_t ring = 0; ring < budget.max_rings; ring++) {
    std::vector<PerfectFace<Tint>> faces_k = cb_k.get_faces();
    if (faces_k.size() > budget.max_unknowns) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: dual_partial_contracting_homotopy k=" << k << " ring=" << ring
         << " unknowns=" << faces_k.size() << " beyond the budget\n";
#endif
      return {};
    }
    ComplexBuilder<T, Tint, Tgroup> cb_km1(k - 1, fce, os);
    using T2 = Eigen::Triplet<T>;
    std::vector<T2> tripletList;
    int nbRow = faces_k.size();
    for (int u = 0; u < nbRow; u++) {
      PerfectFace<Tint> const &F = faces_k[u];
      for (auto &pair : ufc.get_cofaces(k, F.iOrb)) {
        MyMatrix<Tint> M = pair.first.M * F.M;
        cb_km1.f_insert(pair.first.iOrb, M);
        std::pair<size_t, int> p = cb_km1.get_position(pair.first.iOrb, M);
        T value = T(pair.second * p.second);
        tripletList.emplace_back(u, p.first, value);
      }
    }
    int nbCol = cb_km1.get_faces().size();
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: dual_partial_contracting_homotopy k=" << k << " ring=" << ring
       << " unknowns=" << nbRow << " equations=" << nbCol << "\n";
#endif
    MyVector<T> Bvect = ZeroVector<T>(nbCol);
    bool rhs_present = true;
    for (auto &fe : rhs) {
      std::optional<std::pair<size_t, int>> opt_p = cb_km1.find_position(fe.iOrb, fe.M);
      if (!opt_p) {
        rhs_present = false;
        break;
      }
      Bvect(opt_p->first) += T(opt_p->second) * fe.value;
    }
    bool is_graph = (k == 1);
    if (rhs_present) {
      SparseSolveResult<T> res = SparseSolutionMat_L1Exact(tripletList, nbRow, nbCol, Bvect, is_graph, budget.solve, os);
      if (res.status == SparseSolveStatus::Solved) {
        std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
        for (int u = 0; u < nbRow; u++) {
          T val = res.sol(u);
          if (val != 0) {
            PerfectFace<Tint> const &face = faces_k[u];
            chain_ret.push_back({face.iOrb, val, face.M});
          }
        }
        return chain_ret;
      }
    }
    // The next ring: the cells sharing a coface with a current unknown.
    size_t len = faces_k.size();
    for (size_t u = start; u < len; u++) {
      PerfectFace<Tint> const &F = faces_k[u];
      for (auto &pair : ufc.get_cofaces(k, F.iOrb)) {
        MyMatrix<Tint> M = pair.first.M * F.M;
        insert_faces_of(pair.first.iOrb, M);
      }
      if (cb_k.get_faces().size() > budget.max_unknowns) {
        break;
      }
    }
    start = len;
    if (cb_k.get_faces().size() == len) {
      break;
    }
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: dual_partial_contracting_homotopy k=" << k << " no filling found within the budgets\n";
#endif
  return {};
}

/*
  The image of a perfect cone C under y: a cone whose closure contains
  the image of the barycenter of C, found by the walk. When the point
  lies on a lower face the choice is not canonical; the equivariance is
  restored afterwards by the averaging over the coset stabilizer and the
  propagation along its orbit, as for the other levels. The sign is the
  orientation sign of y and of the matrix of the cone.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
get_cone_image(int const &iOrb, MyMatrix<T> const &y,
               FullComplexEnumeration<T, Tint, Tgroup> const &fce,
               ConeWalker<T, Tint, Tgroup> const &walker, int const &eps_y,
               [[maybe_unused]] std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  MyMatrix<Tint> const &EXT = fce.levels[0].l_faces[iOrb].EXT;
  MyMatrix<T> W = UniversalMatrixConversion<T, Tint>(EXT) * y;
  TopPerfectCone<Tint> tpc = walker.find_containing_cone_T(W);
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
  int eps_M = tspace_orientation_sign(LinSpa, M_T);
  PerfectFaceEntry<T, Tint> fe{tpc.i_perfect, T(eps_y * eps_M), tpc.M};
  return {fe};
}

// The right hand side coboundary(F) mapped by y_j, from the images of
// the level index-1.
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> get_hecke_coboundary_image(
    int const &index, int const &iOrb, int const &j,
    FullComplexEnumeration<T, Tint, Tgroup> const &fce, HeckeCosets<T> const &cosets,
    std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>> const &images_prev,
    UpperFaceCache<T, Tint, Tgroup> &ufc, std::ostream &os) {
  ChainBuilder<T, Tint, Tgroup> cb(index - 1, fce, os);
  MyMatrix<T> const &y = cosets.l_coset[j];
  for (auto &pair : ufc.get_cofaces(index, iOrb)) {
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(pair.first.M);
    MyMatrix<T> z = M_T * y;
    std::pair<int, MyMatrix<T>> dec = cosets.decompose(z);
    MyMatrix<Tint> M_conj = UniversalMatrixConversion<Tint, T>(dec.second);
    T e_sign = T(pair.second);
    for (auto &fe : images_prev[pair.first.iOrb][dec.first]) {
      T value = fe.value * e_sign;
      MyMatrix<Tint> M = fe.M * M_conj;
      cb.f_insert(value, fe.iOrb, M);
    }
  }
  return cb.get_faces();
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map_dual(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                             MyMatrix<T> const &x, std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  PerfectComplexOptions const &pco = fce.pctdi.pco;
  if (!pco.compute_boundary || !pco.compute_partial_contracting_homotopy) {
    std::cerr << "HECKE: The dual method requires ComputeBoundary = T and ComputePartialContractingHomotopy = T\n";
    throw TerminalException{1};
  }
  if (!is_stab_space(x, LinSpa)) {
    std::cerr << "HECKE: The matrix x should preserve the T-space\n";
    throw TerminalException{1};
  }
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> l_gens_G = get_ambient_group_generators(fce, os);
  std::function<bool(MyMatrix<T> const &)> f_member =
      [&LinSpa](MyMatrix<T> const &g) -> bool {
    return is_in_ambient_group(LinSpa, g);
  };
  HeckeCosets<T> cosets = enumerate_hecke_cosets<T, Tint>(x, l_gens_G, f_member, os);
  size_t n_coset = cosets.size();
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: enumerate_hecke_cosets|=" << time << "\n";
#endif
  int n_levels = fce.levels.size();
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images(n_levels);
  UpperFaceCache<T, Tint, Tgroup> ufc(fce, os);
  ConeWalker<T, Tint, Tgroup> walker(fce, os);
  HeckeBudget budget;
  MyMatrix<T> Q0 = get_peeling_potential(fce, os);
  auto is_wr = [&](int const &index, int const &iOrb) -> bool {
    return fce.levels[index].l_faces[iOrb].is_well_rounded;
  };
  for (int index = 0; index < n_levels; index++) {
    int n_orb = fce.levels[index].l_faces.size();
    images[index] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      if (!is_wr(index, iOrb)) {
        continue;
      }
      FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
      std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
          get_stab_coset_action(face.l_gens, cosets);
      std::vector<int> status(n_coset, 0);
      for (size_t j = 0; j < n_coset; j++) {
        if (status[j] == 1) {
          continue;
        }
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: dual index=" << index << " iOrb=" << iOrb << "/" << n_orb
           << " j=" << j << "/" << n_coset << "\n";
#endif
        std::vector<PerfectFaceEntry<T, Tint>> rhs;
        std::vector<PerfectFaceEntry<T, Tint>> sol;
        if (index == 0) {
          int eps_y = tspace_orientation_sign(LinSpa, cosets.l_coset[j]);
          sol = get_cone_image(iOrb, cosets.l_coset[j], fce, walker, eps_y, os);
        } else {
          rhs = get_hecke_coboundary_image(index, iOrb, j, fce, cosets, images[index - 1], ufc, os);
          if (!rhs.empty()) {
            std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_sol =
                dual_partial_contracting_homotopy(index, rhs, ufc, Q0, budget, os);
            if (!opt_sol) {
              throw HeckeFailure{"no filling found for the cell of level " + std::to_string(index) +
                                 " orbit " + std::to_string(iOrb) + " coset " + std::to_string(j) +
                                 " (right hand side of " + std::to_string(rhs.size()) + " cells)"};
            }
            sol = *opt_sol;
          }
        }
        CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
        std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
            get_coset_stabilizer_schreier(co, face.l_gens, act, n);
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (index > 0 && !is_chain_coset_stabilizer_invariant(index, iOrb, index - 1, rhs, l_schreier, fce, os)) {
          std::cerr << "HECKE: dual, the right hand side should be invariant under the coset stabilizer\n";
          throw TerminalException{1};
        }
#endif
        std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
            average_chain_coset_stabilizer(index, iOrb, sol, l_schreier, fce, os);
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: dual index=" << index << " iOrb=" << iOrb << " j=" << j
           << " |rhs|=" << rhs.size() << " |sol|=" << sol.size()
           << " |sol_avg|=" << sol_avg.size() << "\n";
#endif
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (index == 0) {
          // The coset stabilizer fixes the image of the barycenter, so the
          // averaged image is made of cones all containing it.
          MyMatrix<T> P = UniversalMatrixConversion<T, Tint>(face.EXT) * cosets.l_coset[j];
          for (auto &fe : sol_avg) {
            std::optional<int> opt = walker.n_facet_containing(fe.iOrb, fe.M, P);
            if (!opt) {
              std::cerr << "HECKE: dual, a cone of the averaged image of iOrb=" << iOrb
                        << " j=" << j << " does not contain the image of the barycenter\n";
              throw TerminalException{1};
            }
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: dual index=0 iOrb=" << iOrb << " j=" << j
               << " n_facet_containing=" << *opt << "\n";
#endif
          }
#ifdef DEBUG_HECKE_OPERATORS
          // The dimension of the space of forms fixed by the coset stabilizer.
          // If it is 1 then the image of the barycenter is the only fixed point
          // and the averaged image cannot be smaller.
          int n_mat = LinSpa.ListMat.size();
          MyMatrix<T> Equa = ZeroMatrix<T>(n_mat, n_mat * l_schreier.size());
          for (size_t i_s = 0; i_s < l_schreier.size(); i_s++) {
            MyMatrix<T> h = UniversalMatrixConversion<T, Tint>(l_schreier[i_s].second);
            MyMatrix<T> A = tspace_action_matrix(LinSpa, h);
            for (int k = 0; k < n_mat; k++) {
              for (int l = 0; l < n_mat; l++) {
                T val = A(k, l);
                if (k == l) {
                  val -= 1;
                }
                Equa(k, i_s * n_mat + l) = val;
              }
            }
          }
          os << "HECKE: dual index=0 iOrb=" << iOrb << " j=" << j
             << " dim_fixed_forms=" << NullspaceMat(Equa).rows() << "\n";
#endif
        }
        {
          if (index > 0) {
            std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(index, sol_avg, ufc, os);
            if (!is_equal_chain(cob, rhs, index - 1, fce, os)) {
              std::cerr << "HECKE: dual, the averaged chain is not a solution\n";
              throw TerminalException{1};
            }
          }
          if (!is_chain_coset_stabilizer_invariant(index, iOrb, index, sol_avg, l_schreier, fce, os)) {
            std::cerr << "HECKE: dual, the averaged chain should be invariant\n";
            throw TerminalException{1};
          }
        }
#endif
        size_t n_orbit = co.l_orbit.size();
        for (size_t pos = 0; pos < n_orbit; pos++) {
          int j2 = co.l_orbit[pos];
          int sign = get_face_orientation(face.EXT, LinSpa.ListMat, face.or_info, co.transversal[pos]);
          MyMatrix<Tint> Minv = Inverse(co.transversal_conj[pos]);
          images[index][iOrb][j2] = chain_right_multiply(sol_avg, T(sign), Minv);
          status[j2] = 1;
        }
      }
#ifdef SANITY_CHECK_HECKE_OPERATORS
      for (size_t j = 0; j < n_coset && index > 0; j++) {
        std::vector<PerfectFaceEntry<T, Tint>> rhs =
            get_hecke_coboundary_image(index, iOrb, j, fce, cosets, images[index - 1], ufc, os);
        std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(index, images[index][iOrb][j], ufc, os);
        if (!is_equal_chain(cob, rhs, index - 1, fce, os)) {
          std::cerr << "HECKE: dual, the chain map property fails at index=" << index
                    << " iOrb=" << iOrb << " j=" << j << "\n";
          throw TerminalException{1};
        }
      }
#endif
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: dual level " << index << " n_orb=" << n_orb << " done\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: dual level " << index << "|=" << time << "\n";
#endif
  }
  return {x, std::move(cosets), std::move(images), {}, "Dual", true};
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_DUAL_H_
// clang-format on
