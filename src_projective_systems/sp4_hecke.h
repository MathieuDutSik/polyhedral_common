// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PROJECTIVE_SYSTEMS_SP4_HECKE_H_
#define SRC_PROJECTIVE_SYSTEMS_SP4_HECKE_H_

// clang-format off
#include "sp4_complex.h"
#include "SparseChainSolver.h"
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SP4_HECKE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SP4_HECKE
#endif

#ifdef TIMINGS
#define TIMINGS_SP4_HECKE
#endif

/*
  Step A of the Hecke operators for Sp_4(Z) (SP4Z_ComputeHeckeOperators):
  the chain map of a double coset G x G on the complex of Sp_4(Z), as the
  images Phi(F, y_j) of the representative cells F of each level under the
  cosets y_j of G x G, in the format of HeckeChainMap of src_perfect, so
  that step B (compute_hecke_homology) gives the Hecke operators on the
  homology of any finite index subgroup.

  - A single vector v (the last index) goes to the cell of v y_j, found by
    a walk through the top cells (SP4Z_GetSequenceToRay): a top cell
    containing v y_j up to the positive multiples.
  - A cell F of a lower index goes to a filling of the image of its
    boundary, found on a covering of top cells containing the vectors of
    that image (SP4Z_SolutionInComplex), averaged over the stabilizer of
    the coset in Stab(F) and propagated to the orbit of the coset, as in
    the Vertex method of src_perfect.
 */

namespace sp4 {

// The vector of the upper triangle of v^T v (the rank one form of v).
inline MyVector<T> rank_one_form(MyVector<Tint> const &v) {
  MyVector<T> w(10);
  int pos = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i; j < 4; j++) {
      w(pos) = T(v(i) * v(j));
      pos += 1;
    }
  }
  return w;
}

// A top cell placed in the space: l_top[iTopCell].EXT * eMat.
struct PlacedCell {
  int iTopCell;
  MyMatrix<Tint> eMat;
  MyMatrix<Tint> EXTcell;
  ConfigKey key;
};

inline PlacedCell get_placed_cell(std::vector<TopCell> const &l_top, int const &iTopCell, MyMatrix<Tint> const &eMat) {
  MyMatrix<Tint> EXTcell = l_top[iTopCell].EXT * eMat;
  ConfigKey key = config_key(EXTcell);
  return {iTopCell, eMat, std::move(EXTcell), std::move(key)};
}

/*
  The value at the rank one form of v of the facet inequality of the cone
  of the rank one forms of the cell for the ridge eSet: positive in the
  cell, zero on the ridge.
 */
inline T facet_value(MyMatrix<Tint> const &EXTcell, std::vector<int> const &eSet, MyVector<Tint> const &v) {
  int n_set = eSet.size();
  MyMatrix<T> M(n_set, 10);
  for (int i = 0; i < n_set; i++) {
    MyVector<T> w = rank_one_form(GetMatrixRow(EXTcell, eSet[i]));
    for (int j = 0; j < 10; j++) {
      M(i, j) = w(j);
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(M);
  if (NSP.rows() != 1) {
    std::cerr << "SP4: facet_value, the ridge should span a hyperplane\n";
    throw TerminalException{1};
  }
  MyVector<T> fac = GetMatrixRow(NSP, 0);
  for (int i = 0; i < EXTcell.rows(); i++) {
    if (!std::binary_search(eSet.begin(), eSet.end(), i)) {
      T val = fac.dot(rank_one_form(GetMatrixRow(EXTcell, i)));
      if (val < 0) {
        fac = -fac;
      }
      break;
    }
  }
  return fac.dot(rank_one_form(v));
}

/*
  A walk through the top cells to a top cell containing v (up to sign) and
  the vectors l_cont (SP4Z_GetSequenceToRay_General): from the current
  cell, cross a facet whose inequality is negative at v while the new cell
  contains l_cont; when no such move exists, the cells containing l_cont
  around the visited ones are explored ("flowering"). Returns the sequence
  of the visited cells, reduced to a shortest path of adjacent cells.
 */
inline std::vector<PlacedCell> get_sequence_to_ray(std::vector<TopCell> const &l_top, MyVector<Tint> const &v,
                                                   MyMatrix<Tint> const &l_cont, PlacedCell const &start) {
  std::vector<PlacedCell> l_cell;
  std::map<ConfigKey, int> map_cell;
  auto f_contains_all = [&](MyMatrix<Tint> const &EXTcell) -> bool {
    for (int i = 0; i < l_cont.rows(); i++) {
      if (!position_antipodal(EXTcell, GetMatrixRow(l_cont, i))) {
        return false;
      }
    }
    return true;
  };
  auto f_insert = [&](PlacedCell const &c) -> std::pair<int, bool> {
    auto iter = map_cell.find(c.key);
    if (iter != map_cell.end()) {
      return {iter->second, false};
    }
    map_cell[c.key] = l_cell.size();
    l_cell.push_back(c);
    return {static_cast<int>(l_cell.size()) - 1, true};
  };
  auto f_next = [&](PlacedCell const &c) -> std::optional<PlacedCell> {
    for (auto &adj : l_top[c.iTopCell].l_adj) {
      if (facet_value(c.EXTcell, adj.eSet, v) < 0) {
        PlacedCell c_new = get_placed_cell(l_top, adj.iCell, adj.eEquiv * c.eMat);
        if (f_contains_all(c_new.EXTcell)) {
          return c_new;
        }
      }
    }
    return {};
  };
  f_insert(start);
  HeckeBudget budget;
  size_t n_step = 0;
  while (true) {
    // The greedy walk from the last cell.
    bool found = false;
    while (true) {
      n_step += 1;
      if (n_step > budget.max_walk_steps) {
        throw HeckeFailure{"SP4: get_sequence_to_ray, too many steps"};
      }
      PlacedCell const &c = l_cell.back();
      if (position_antipodal(c.EXTcell, v)) {
        found = true;
        break;
      }
      std::optional<PlacedCell> opt = f_next(c);
      if (!opt || !f_insert(*opt).second) {
        break;
      }
    }
    if (found) {
      break;
    }
    // The flowering from the last cell: the cells containing l_cont around
    // it, until one contains v or allows a move to a new cell, from which
    // the greedy walk resumes (it is the last cell).
    std::vector<int> l_queue{static_cast<int>(l_cell.size()) - 1};
    size_t pos = 0;
    bool resume = false;
    while (pos < l_queue.size() && !found && !resume) {
      PlacedCell c = l_cell[l_queue[pos]];
      for (auto &adj : l_top[c.iTopCell].l_adj) {
        PlacedCell c_new = get_placed_cell(l_top, adj.iCell, adj.eEquiv * c.eMat);
        if (!f_contains_all(c_new.EXTcell)) {
          continue;
        }
        std::pair<int, bool> ins = f_insert(c_new);
        if (!ins.second) {
          continue;
        }
        if (position_antipodal(c_new.EXTcell, v)) {
          found = true;
          break;
        }
        std::optional<PlacedCell> opt = f_next(c_new);
        if (opt && map_cell.count(opt->key) == 0) {
          resume = true;
          break;
        }
        l_queue.push_back(ins.first);
      }
      pos += 1;
      if (l_cell.size() > budget.max_walk_steps) {
        throw HeckeFailure{"SP4: get_sequence_to_ray, too many cells"};
      }
    }
    if (found) {
      break;
    }
    if (!resume) {
      throw HeckeFailure{"SP4: get_sequence_to_ray, the walk is stuck"};
    }
  }
  // The shortest path of adjacent cells (sharing a ridge) from the start
  // to the last cell.
  int n_cell = l_cell.size();
  int i_end = n_cell - 1;
  auto f_adjacent = [&](int i, int j) -> bool {
    int n_common = 0;
    for (int a = 0; a < l_cell[i].EXTcell.rows(); a++) {
      if (position_antipodal(l_cell[j].EXTcell, GetMatrixRow(l_cell[i].EXTcell, a))) {
        n_common += 1;
      }
    }
    return n_common == 9;
  };
  std::vector<int> prev(n_cell, -1);
  std::vector<int> l_queue{0};
  prev[0] = 0;
  size_t pos = 0;
  while (pos < l_queue.size() && prev[i_end] == -1) {
    int i = l_queue[pos];
    for (int j = 0; j < n_cell; j++) {
      if (prev[j] == -1 && f_adjacent(i, j)) {
        prev[j] = i;
        l_queue.push_back(j);
      }
    }
    pos += 1;
  }
  std::vector<PlacedCell> l_path;
  if (prev[i_end] == -1) {
    l_path.push_back(l_cell[i_end]);
    return l_path;
  }
  int i = i_end;
  while (i != 0) {
    l_path.push_back(l_cell[i]);
    i = prev[i];
  }
  l_path.push_back(l_cell[0]);
  std::reverse(l_path.begin(), l_path.end());
  return l_path;
}

/*
  Chains of a level: lists of entries value * (F_{iOrb} M). Two entries
  are the same cell iff their vectors are the same up to sign, and then
  M2 M1^{-1} is in the stabilizer of F_{iOrb}, whose orientation gives the
  relative sign.
 */
using Chain = std::vector<PerfectFaceEntry<T, Tint>>;

struct ChainReducer {
  Sp4CellComplex const &cx;
  int index;
  std::map<std::pair<int, ConfigKey>, size_t> map_cell;
  std::vector<PerfectFaceEntry<T, Tint>> l_entry;
  ChainReducer(Sp4CellComplex const &_cx, int const &_index) : cx(_cx), index(_index) {}
  // The position of the cell and the sign of M relative to the stored one.
  std::pair<size_t, int> locate(int const &iOrb, MyMatrix<Tint> const &M) {
    CellOrbit const &orbit = cx.levels[index][iOrb];
    std::pair<int, ConfigKey> key{iOrb, config_key(orbit.EXT * M)};
    auto iter = map_cell.find(key);
    if (iter == map_cell.end()) {
      size_t pos = l_entry.size();
      map_cell[key] = pos;
      l_entry.push_back({iOrb, T(0), M});
      return {pos, 1};
    }
    size_t pos = iter->second;
    MyMatrix<Tint> quot = M * Inverse(l_entry[pos].M);
    return {pos, orbit_orientation(orbit, quot)};
  }
  void insert(int const &iOrb, T const &value, MyMatrix<Tint> const &M) {
    std::pair<size_t, int> pair = locate(iOrb, M);
    l_entry[pair.first].value += value * T(pair.second);
  }
  Chain get_chain() const {
    Chain chain;
    for (auto &fe : l_entry) {
      if (fe.value != 0) {
        chain.push_back(fe);
      }
    }
    return chain;
  }
};

inline Chain reduce_chain(Sp4CellComplex const &cx, int const &index, Chain const &chain) {
  ChainReducer cr(cx, index);
  for (auto &fe : chain) {
    cr.insert(fe.iOrb, fe.value, fe.M);
  }
  return cr.get_chain();
}

// The boundary of a chain of the level index, a chain of the level index + 1.
inline Chain chain_boundary_sp4(Sp4CellComplex const &cx, int const &index, Chain const &chain) {
  ChainReducer cr(cx, index + 1);
  for (auto &fe : chain) {
    for (auto &eb : cx.levels[index][fe.iOrb].l_bound) {
      cr.insert(eb.iOrb, fe.value * T(eb.sign), eb.M * fe.M);
    }
  }
  return cr.get_chain();
}

inline bool is_equal_chain_sp4(Sp4CellComplex const &cx, int const &index, Chain const &c1, Chain const &c2) {
  Chain diff = c1;
  for (auto &fe : c2) {
    diff.push_back({fe.iOrb, -fe.value, fe.M});
  }
  return reduce_chain(cx, index, diff).empty();
}

/*
  The image of the boundary of F_{iOrb} (level index) under the coset y_j,
  from the images of the level index + 1:
    Phi(dF, y_j) = sum_faces sign Phi(F', y_k) g'   with M y_j = y_k g'.
 */
inline Chain get_boundary_image_sp4(Sp4CellComplex const &cx, int const &index, int const &iOrb, size_t const &j,
                                    HeckeCosets<T> const &cosets,
                                    std::vector<std::vector<Chain>> const &images_next) {
  ChainReducer cr(cx, index + 1);
  MyMatrix<T> const &y = cosets.l_coset[j];
  for (auto &eb : cx.levels[index][iOrb].l_bound) {
    MyMatrix<T> z = UniversalMatrixConversion<T, Tint>(eb.M) * y;
    std::pair<int, MyMatrix<T>> pair = cosets.decompose(z);
    MyMatrix<Tint> g_conj = UniversalMatrixConversion<Tint, T>(pair.second);
    for (auto &fe : images_next[eb.iOrb][pair.first]) {
      cr.insert(fe.iOrb, fe.value * T(eb.sign), fe.M * g_conj);
    }
  }
  return cr.get_chain();
}

/*
  The image of a single vector cell under y: the cell of the primitive
  vector of v y, found in a top cell by the walk.
 */
inline Chain get_vertex_image_sp4(Sp4CellComplex const &cx, int const &iOrb, MyMatrix<T> const &y) {
  int index = cx.levels.size() - 1;
  CellOrbit const &orbit = cx.levels[index][iOrb];
  MyVector<T> v0 = UniversalVectorConversion<T, Tint>(GetMatrixRow(orbit.EXT, 0));
  MyVector<T> w = y.transpose() * v0;
  w = RemoveFractionVector(w);
  MyVector<Tint> v = UniversalVectorConversion<Tint, T>(w);
  PlacedCell start = get_placed_cell(cx.l_top, 0, IdentityMat<Tint>(4));
  std::vector<PlacedCell> l_seq = get_sequence_to_ray(cx.l_top, v, MyMatrix<Tint>(0, 4), start);
  PlacedCell const &c = l_seq.back();
  int pos = *position_antipodal(c.EXTcell, v);
  Triple t{c.iTopCell, {pos}, c.eMat};
  for (size_t jOrb = 0; jOrb < cx.levels[index].size(); jOrb++) {
    for (auto &t_orb : cx.levels[index][jOrb].l_triple) {
      std::optional<MyMatrix<Tint>> opt = test_equivalence_triple(cx.l_top, t_orb, t);
      if (opt) {
        return {{static_cast<int>(jOrb), T(1), *opt}};
      }
    }
  }
  throw HeckeFailure{"SP4: get_vertex_image_sp4, the cell of the vector is not found"};
}

/*
  A filling of rhs (a cycle of the level index + 1, a boundary in the
  complex): a chain X of the level index with dX = rhs, found on a covering
  of top cells containing the vectors of rhs (SP4Z_SolutionInComplex).
  The cells are restricted to the span of the vectors of rhs when they do
  not span Z^4 (the singularities of the cells at infinity are kept). The
  covering grows by the neighbors when the system has no solution.
 */
inline Chain solve_in_complex(Sp4CellComplex const &cx, int const &index, Chain const &rhs, std::ostream &os) {
  std::vector<TopCell> const &l_top = cx.l_top;
  // The vectors of the right hand side.
  std::vector<MyVector<Tint>> l_vert;
  for (auto &fe : rhs) {
    MyMatrix<Tint> EXT = cx.levels[index + 1][fe.iOrb].EXT * fe.M;
    for (int i = 0; i < EXT.rows(); i++) {
      MyVector<Tint> v = GetMatrixRow(EXT, i);
      bool is_new = true;
      for (auto &u : l_vert) {
        if (u == v || u == MyVector<Tint>(-v)) {
          is_new = false;
        }
      }
      if (is_new) {
        l_vert.push_back(v);
      }
    }
  }
  MyMatrix<Tint> M_vert(l_vert.size(), 4);
  for (size_t i = 0; i < l_vert.size(); i++) {
    for (int j = 0; j < 4; j++) {
      M_vert(i, j) = l_vert[i](j);
    }
  }
  MyMatrix<T> M_vert_T = UniversalMatrixConversion<T, Tint>(M_vert);
  MyMatrix<T> l_equa = (RankMat(M_vert_T) == 4) ? MyMatrix<T>(0, 4) : NullspaceTrMat(M_vert_T);
  auto f_correct_face = [&](MyMatrix<Tint> const &EXT) -> bool {
    for (int i = 0; i < EXT.rows(); i++) {
      MyVector<T> v = UniversalVectorConversion<T, Tint>(GetMatrixRow(EXT, i));
      for (int k = 0; k < l_equa.rows(); k++) {
        if (GetMatrixRow(l_equa, k).dot(v) != 0) {
          return false;
        }
      }
    }
    return true;
  };
  // The first covering: the walks to the vectors.
  std::vector<PlacedCell> l_cov;
  std::set<ConfigKey> set_cov;
  auto f_insert_cov = [&](PlacedCell const &c) -> void {
    if (set_cov.insert(c.key).second) {
      l_cov.push_back(c);
    }
  };
  PlacedCell start = get_placed_cell(l_top, 0, IdentityMat<Tint>(4));
  std::vector<PlacedCell> l_seq = get_sequence_to_ray(l_top, l_vert[0], MyMatrix<Tint>(0, 4), start);
  PlacedCell good = l_seq.back();
  f_insert_cov(good);
  for (size_t i = 1; i < l_vert.size(); i++) {
    bool covered = false;
    for (auto &c : l_cov) {
      if (position_antipodal(c.EXTcell, l_vert[i])) {
        covered = true;
      }
    }
    if (!covered) {
      for (auto &c : get_sequence_to_ray(l_top, l_vert[i], MyMatrix<Tint>(0, 4), good)) {
        f_insert_cov(c);
      }
    }
  }
  // Each cell of rhs must be in a cell of the covering: otherwise the
  // walk to its vectors one by one, keeping the previous ones.
  for (auto &fe : rhs) {
    MyMatrix<Tint> EXT = cx.levels[index + 1][fe.iOrb].EXT * fe.M;
    bool contained = false;
    for (auto &c : l_cov) {
      bool all = true;
      for (int i = 0; i < EXT.rows(); i++) {
        if (!position_antipodal(c.EXTcell, GetMatrixRow(EXT, i))) {
          all = false;
        }
      }
      if (all) {
        contained = true;
        break;
      }
    }
    if (!contained) {
      PlacedCell c0 = l_cov[0];
      for (auto &c : l_cov) {
        if (position_antipodal(c.EXTcell, GetMatrixRow(EXT, 0))) {
          c0 = c;
          break;
        }
      }
      for (int iv = 0; iv < EXT.rows(); iv++) {
        MyMatrix<Tint> l_cont(iv, 4);
        for (int a = 0; a < iv; a++) {
          for (int b = 0; b < 4; b++) {
            l_cont(a, b) = EXT(a, b);
          }
        }
        std::vector<PlacedCell> l_seq2 = get_sequence_to_ray(l_top, GetMatrixRow(EXT, iv), l_cont, c0);
        for (auto &c : l_seq2) {
          f_insert_cov(c);
        }
        c0 = l_seq2.back();
      }
    }
  }
  SparseSolveBudget budget;
  bool is_graph = (index == static_cast<int>(cx.levels.size()) - 2);
  for (int n_flower = 0; n_flower < 100; n_flower++) {
    // The cells of the levels index and index + 1 in the covering.
    ChainReducer cr_row(cx, index);
    ChainReducer cr_col(cx, index + 1);
    for (int lev = index; lev <= index + 1; lev++) {
      ChainReducer &cr = (lev == index) ? cr_row : cr_col;
      for (auto &c : l_cov) {
        for (auto &fi : cx.face_index[lev][c.iTopCell]) {
          MyMatrix<Tint> uMat = fi.eEquiv * c.eMat;
          MyMatrix<Tint> EXT = cx.levels[lev][fi.iOrbit].EXT * uMat;
          if (f_correct_face(EXT)) {
            cr.locate(fi.iOrbit, uMat);
          }
        }
      }
    }
    int nbRow = cr_row.l_entry.size();
    std::vector<Eigen::Triplet<T>> tripletList;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      PerfectFaceEntry<T, Tint> const &fe = cr_row.l_entry[iRow];
      std::map<int, T> map_row;
      for (auto &eb : cx.levels[index][fe.iOrb].l_bound) {
        std::pair<size_t, int> pair = cr_col.locate(eb.iOrb, eb.M * fe.M);
        map_row[pair.first] += T(eb.sign * pair.second);
      }
      for (auto &kv : map_row) {
        if (kv.second != 0) {
          tripletList.emplace_back(iRow, kv.first, kv.second);
        }
      }
    }
    // A cell of rhs outside the covering gets a new column, which no row
    // reaches: the system has no solution and the covering grows.
    std::vector<std::pair<size_t, T>> l_rhs;
    for (auto &fe : rhs) {
      std::pair<size_t, int> pair = cr_col.locate(fe.iOrb, fe.M);
      l_rhs.push_back({pair.first, fe.value * T(pair.second)});
    }
    int nbCol = cr_col.l_entry.size();
    MyVector<T> Bvect = ZeroVector<T>(nbCol);
    for (auto &pair : l_rhs) {
      Bvect(pair.first) += pair.second;
    }
    SparseSolveResult<T> res = SparseSolutionMat_L1Exact(tripletList, nbRow, nbCol, Bvect, is_graph, budget, os);
    if (res.status == SparseSolveStatus::Solved) {
      Chain sol;
      for (int iRow = 0; iRow < nbRow; iRow++) {
        if (res.sol(iRow) != 0) {
          sol.push_back({cr_row.l_entry[iRow].iOrb, res.sol(iRow), cr_row.l_entry[iRow].M});
        }
      }
#ifdef SANITY_CHECK_SP4_HECKE
      if (!is_equal_chain_sp4(cx, index + 1, chain_boundary_sp4(cx, index, sol), rhs)) {
        std::cerr << "SP4: solve_in_complex, the solution is not a filling\n";
        throw TerminalException{1};
      }
#endif
      return sol;
    }
#ifdef DEBUG_SP4_HECKE
    os << "SP4: solve_in_complex index=" << index << " flowering " << n_flower << " |l_cov|=" << l_cov.size()
       << "\n";
#endif
    // The flowering: the neighbors of the covering.
    std::vector<PlacedCell> l_cov_old = l_cov;
    for (auto &c : l_cov_old) {
      for (auto &adj : l_top[c.iTopCell].l_adj) {
        f_insert_cov(get_placed_cell(l_top, adj.iCell, adj.eEquiv * c.eMat));
      }
    }
  }
  throw HeckeFailure{"SP4: solve_in_complex, no filling found after 100 flowerings"};
}

/*
  The average of a chain of the level index over the finite action of the
  stabilizer of the coset (as in average_chain_coset_stabilizer): the
  pairs (h, h') with h in Stab(F) and h' = y_j^{-1} h y_j act on the chain
  by chain -> sign_F(h) chain h'. The closure is over the distinct images
  of the chain (the action on the chain is finite even for a cell at
  infinity), and the average is their mean.
 */
inline Chain average_chain_sp4(Sp4CellComplex const &cx, int const &index, int const &iOrb, Chain const &chain,
                               std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> const &l_schreier) {
  if (l_schreier.empty()) {
    return chain;
  }
  CellOrbit const &orbit = cx.levels[index][iOrb];
  std::vector<int> l_sign;
  for (auto &pair : l_schreier) {
    l_sign.push_back(orbit_orientation(orbit, pair.first));
  }
  // The chains, as vectors on a common list of cells.
  ChainReducer cr(cx, index);
  auto f_vect = [&](Chain const &c) -> std::map<size_t, T> {
    std::map<size_t, T> V;
    for (auto &fe : c) {
      std::pair<size_t, int> pair = cr.locate(fe.iOrb, fe.M);
      V[pair.first] += fe.value * T(pair.second);
    }
    std::map<size_t, T> W;
    for (auto &kv : V) {
      if (kv.second != 0) {
        W[kv.first] = kv.second;
      }
    }
    return W;
  };
  std::vector<Chain> l_chain{chain};
  std::set<std::map<size_t, T>> set_vect{f_vect(chain)};
  std::map<size_t, T> sum = f_vect(chain);
  size_t pos = 0;
  HeckeBudget budget;
  while (pos < l_chain.size()) {
    if (l_chain.size() > budget.max_closure) {
      throw HeckeFailure{"SP4: average_chain_sp4, too many images"};
    }
    for (size_t i = 0; i < l_schreier.size(); i++) {
      Chain c_new = chain_right_multiply(l_chain[pos], T(l_sign[i]), l_schreier[i].second);
      std::map<size_t, T> V = f_vect(c_new);
      if (set_vect.insert(V).second) {
        l_chain.push_back(c_new);
        for (auto &kv : V) {
          sum[kv.first] += kv.second;
        }
      }
    }
    pos += 1;
  }
  T n_img(static_cast<long>(l_chain.size()));
  Chain avg;
  for (auto &kv : sum) {
    if (kv.second != 0) {
      avg.push_back({cr.l_entry[kv.first].iOrb, kv.second / n_img, cr.l_entry[kv.first].M});
    }
  }
  return avg;
}

/*
  The chain map of x on the levels from the single vectors (the last
  index) up to the index index_min (SP4Z_ComputeHeckeOperators with
  eDimMax = n_levels - index_min).
 */
inline HeckeChainMap<T, Tint> compute_sp4_hecke_chain_map(Sp4CellComplex const &cx, MyMatrix<T> const &x,
                                                          int const &index_min, std::ostream &os) {
#ifdef TIMINGS_SP4_HECKE
  MicrosecondTime time;
#endif
  int n_levels = cx.levels.size();
  std::function<bool(MyMatrix<T> const &)> f_member = [](MyMatrix<T> const &g) -> bool { return is_in_sp4(g); };
  HeckeCosets<T> cosets = enumerate_hecke_cosets<T, Tint>(x, cx.l_gens_G, f_member, os);
  size_t n_coset = cosets.size();
#ifdef TIMINGS_SP4_HECKE
  os << "|SP4: Hecke cosets, n_coset=" << n_coset << "|=" << time << "\n";
#endif
  std::vector<std::vector<std::vector<Chain>>> images(n_levels);
  std::vector<uint8_t> l_level_computed(n_levels, 0);
  for (int index = 0; index < n_levels; index++) {
    images[index] = std::vector<std::vector<Chain>>(cx.levels[index].size(), std::vector<Chain>(n_coset));
  }
  // The single vectors.
  int vlevel = n_levels - 1;
  for (size_t iOrb = 0; iOrb < cx.levels[vlevel].size(); iOrb++) {
    for (size_t j = 0; j < n_coset; j++) {
      images[vlevel][iOrb][j] = get_vertex_image_sp4(cx, iOrb, cosets.l_coset[j]);
    }
  }
  l_level_computed[vlevel] = 1;
#ifdef TIMINGS_SP4_HECKE
  os << "|SP4: vertex images|=" << time << "\n";
#endif
  int n = 4;
  for (int index = vlevel - 1; index >= index_min; index--) {
    for (size_t iOrb = 0; iOrb < cx.levels[index].size(); iOrb++) {
      CellOrbit const &orbit = cx.levels[index][iOrb];
      std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act = get_stab_coset_action(orbit.l_stab_gens, cosets);
      std::vector<int> status(n_coset, 0);
      for (size_t j = 0; j < n_coset; j++) {
        if (status[j] == 1) {
          continue;
        }
        Chain rhs = get_boundary_image_sp4(cx, index, iOrb, j, cosets, images[index + 1]);
        Chain sol;
        if (!rhs.empty()) {
          sol = solve_in_complex(cx, index, rhs, os);
        }
        CosetOrbit<Tint> co = get_coset_orbit(j, orbit.l_stab_gens, act, n_coset, n);
        std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
            get_coset_stabilizer_schreier(co, orbit.l_stab_gens, act, n);
        Chain sol_avg = average_chain_sp4(cx, index, iOrb, sol, l_schreier);
#ifdef SANITY_CHECK_SP4_HECKE
        if (!is_equal_chain_sp4(cx, index + 1, chain_boundary_sp4(cx, index, sol_avg), rhs)) {
          std::cerr << "SP4: the averaged chain is not a filling, index=" << index << " iOrb=" << iOrb << "\n";
          throw TerminalException{1};
        }
#endif
        // Phi(F, y_{j'}) = sign(z) Phi(F, y_j) z'^{-1} for z y_j = y_{j'} z'.
        for (size_t p = 0; p < co.l_orbit.size(); p++) {
          int j2 = co.l_orbit[p];
          int sign = orbit_orientation(orbit, co.transversal[p]);
          MyMatrix<Tint> Minv = Inverse(co.transversal_conj[p]);
          images[index][iOrb][j2] = chain_right_multiply(sol_avg, T(sign), Minv);
          status[j2] = 1;
        }
      }
#ifdef SANITY_CHECK_SP4_HECKE
      for (size_t j = 0; j < n_coset; j++) {
        Chain rhs = get_boundary_image_sp4(cx, index, iOrb, j, cosets, images[index + 1]);
        if (!is_equal_chain_sp4(cx, index + 1, chain_boundary_sp4(cx, index, images[index][iOrb][j]), rhs)) {
          std::cerr << "SP4: the chain map property fails, index=" << index << " iOrb=" << iOrb << " j=" << j
                    << "\n";
          throw TerminalException{1};
        }
      }
#endif
    }
    l_level_computed[index] = 1;
#ifdef TIMINGS_SP4_HECKE
    os << "|SP4: images of the level " << index << "|=" << time << "\n";
#endif
  }
  HeckeChainMap<T, Tint> hcm;
  hcm.x = x;
  hcm.cosets = std::move(cosets);
  hcm.images = std::move(images);
  hcm.l_level_computed = std::move(l_level_computed);
  hcm.method = "Sp4";
  hcm.transposed = false;
  return hcm;
}

}  // namespace sp4

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_SP4_HECKE_H_
// clang-format on
