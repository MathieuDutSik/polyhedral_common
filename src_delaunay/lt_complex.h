// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_LT_COMPLEX_H_
#define SRC_DELAUNAY_LT_COMPLEX_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "MAT_MatrixInt.h"
#include "SignatureSymmetric.h"
#include "Tspace_StabEquiInv.h"
#include "boost_serialization.h"
#include <boost/dynamic_bitset/serialization.hpp>
#include <fstream>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_LT_COMPLEX
#endif

#ifdef DEBUG
#define DEBUG_LT_COMPLEX
#endif

/*
  Enumeration of the cells of the L-type complex.

  The cone of positive semidefinite quadratic forms is subdivided into
  cells (closures of L-type / iso-Delaunay domains) of every dimension
  1 <= idx <= n*(n+1)/2. Each top-dimensional iso-Delaunay domain comes
  from the standard iso-Delaunay enumeration; the lower-dim cells are
  faces of the top-dim cones in the t-space.

  Per cell we expose:
  * The space of quadratic forms (LinSpa.ListMat, shared)
  * The inequalities (and equalities) defining the cell in the t-space
  * An inner ray (interior point of the cell, in t-space coordinates)

  With compute_full_dimensional = T, we keep only the cells whose inner
  ray represents a positive-definite Gram matrix (the "full-dimensional"
  cells in the user's terminology). Cells with degenerate inner forms
  (i.e. inner ray on the boundary of the PSD cone) are dropped.

  For each top-dim iso-Delaunay cell we enumerate its faces locally by
  walking the face lattice of its L-type cone in the t-space, quotienting
  by the lattice stabiliser acting on the inequalities. Cross-cell
  equivalences for the codim-1 facets are inherited from the iso-Delaunay
  adjacency structure (each facet appears as exactly one orbit per
  adjacency pair).
 */

struct LtypeComplexOptions {
  bool compute_full_dimensional;
};

inline bool LtypeComplexOptionsEqual(LtypeComplexOptions const &a,
                                     LtypeComplexOptions const &b) {
  return a.compute_full_dimensional == b.compute_full_dimensional;
}

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, LtypeComplexOptions &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("compute_full_dimensional", val.compute_full_dimensional);
}
}  // namespace boost::serialization

//
// A cell of the L-type complex.
//

template <typename T>
struct LtypeCell {
  // Subset of the parent top-dim cone's irredundant facet inequalities
  // that are tight on this cell (= equalities of the cell in t-space).
  Face f_inc;
  // The parent top-dim iso-Delaunay cell index in l_tot.
  int i_top;
  // Inequalities (rows of length dimSpace) defining the cone in t-space
  // (the non-tight irredundant inequalities of the parent restricted to
  // the affine span of this cell).
  MyMatrix<T> ListIneq;
  // Equalities (rows of length dimSpace) carving out the affine span of
  // this cell inside t-space (= the tight inequalities of the parent).
  MyMatrix<T> ListEqua;
  // Interior point of the cell, in t-space coordinates.
  MyVector<T> inner_ray;
  // Gram matrix corresponding to inner_ray, for the pos-def filter.
  MyMatrix<T> inner_gram;
  // Whether inner_gram is positive definite (i.e. the cell is "full
  // dimensional" in the user's sense).
  bool inner_is_pos_def;
};

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, LtypeCell<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("f_inc", val.f_inc);
  ar &make_nvp("i_top", val.i_top);
  ar &make_nvp("ListIneq", val.ListIneq);
  ar &make_nvp("ListEqua", val.ListEqua);
  ar &make_nvp("inner_ray", val.inner_ray);
  ar &make_nvp("inner_gram", val.inner_gram);
  ar &make_nvp("inner_is_pos_def", val.inner_is_pos_def);
}
}  // namespace boost::serialization

template <typename T>
struct LtypeCellsAtLevel {
  // idx of the cells at this level (= dimension in t-space).
  int idx;
  std::vector<LtypeCell<T>> l_cells;
};

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, LtypeCellsAtLevel<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("idx", val.idx);
  ar &make_nvp("l_cells", val.l_cells);
}
}  // namespace boost::serialization

template <typename T>
struct FullLtypeComplexEnumeration {
  LinSpaceMatrix<T> LinSpa;
  LtypeComplexOptions opts;
  // Grouped by dimension idx (= dim of cell in t-space).
  std::vector<LtypeCellsAtLevel<T>> levels;
};

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, FullLtypeComplexEnumeration<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("LinSpa", val.LinSpa);
  ar &make_nvp("opts", val.opts);
  ar &make_nvp("levels", val.levels);
}
}  // namespace boost::serialization

//
// Helpers.
//

// Build a Gram matrix from a t-space coordinate vector.
template <typename T>
MyMatrix<T> gram_from_tspace_vec(MyVector<T> const &v,
                                 std::vector<MyMatrix<T>> const &ListMat) {
  int dimSpace = ListMat.size();
  int n = ListMat[0].rows();
  MyMatrix<T> M = ZeroMatrix<T>(n, n);
  for (int u = 0; u < dimSpace; u++) {
    M += v(u) * ListMat[u];
  }
  return M;
}

// Action of g (lattice element in GL_n(Z)) on FAC indices. The
// inequality W transforms as W -> matrix_in_t_space(g) * W, exactly
// matching the convention used inside get_result_delaunay_adj.
template <typename T, typename Telt>
Telt permutation_on_fac(MyMatrix<T> const &MatSpace, MyMatrix<T> const &FAC) {
  using Tidx = typename Telt::Tidx;
  int n_row = FAC.rows();
  std::unordered_map<MyVector<T>, Tidx> MapV;
  for (int i = 0; i < n_row; i++) {
    MyVector<T> W = GetMatrixRow(FAC, i);
    MyVector<T> W_red = RemoveFractionVector(W);
    MapV[W_red] = static_cast<Tidx>(i);
  }
  std::vector<Tidx> img(n_row);
  for (int i = 0; i < n_row; i++) {
    MyVector<T> W = GetMatrixRow(FAC, i);
    MyVector<T> W_img = MatSpace * W;
    MyVector<T> W_img_red = RemoveFractionVector(W_img);
    auto it = MapV.find(W_img_red);
    if (it == MapV.end()) {
      std::cerr << "LTCOMP: permutation_on_fac: image inequality not in FAC\n";
      throw TerminalException{1};
    }
    img[i] = it->second;
  }
  return Telt(std::move(img));
}

// Group acting on FAC indices, derived from the lattice stabiliser of
// the cell's interior Gram matrix.
template <typename T, typename Tint, typename Tgroup>
Tgroup get_fac_permutation_group(LinSpaceMatrix<T> const &LinSpa,
                                  MyMatrix<T> const &GramMat,
                                  MyMatrix<T> const &FAC, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  std::vector<MyMatrix<T>> ListGenTot =
      LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, GramMat, os);
  std::vector<Telt> l_perm;
  l_perm.reserve(ListGenTot.size());
  for (auto &g : ListGenTot) {
    MyMatrix<T> MatSpace = matrix_in_t_space(g, LinSpa);
    l_perm.push_back(permutation_on_fac<T, Telt>(MatSpace, FAC));
  }
  return Tgroup(l_perm, FAC.rows());
}

// Compute (ListIneq, ListEqua) by splitting FAC according to which
// inequalities are marked as tight in f_inc.
template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
split_fac(MyMatrix<T> const &FAC, Face const &f_inc) {
  int nbIneq = FAC.rows();
  int dim = FAC.cols();
  int n_eq = f_inc.count();
  int n_in = nbIneq - n_eq;
  MyMatrix<T> ListIneq(n_in, dim);
  MyMatrix<T> ListEqua(n_eq, dim);
  int pos_in = 0;
  int pos_eq = 0;
  for (int j = 0; j < nbIneq; j++) {
    MyVector<T> W = GetMatrixRow(FAC, j);
    if (f_inc[j] == 1) {
      AssignMatrixRow(ListEqua, pos_eq, W);
      pos_eq += 1;
    } else {
      AssignMatrixRow(ListIneq, pos_in, W);
      pos_in += 1;
    }
  }
  return {std::move(ListIneq), std::move(ListEqua)};
}

// Compute an interior point of the face defined by f_inc, in t-space
// coordinates. Returns std::nullopt if the face is empty (e.g. the
// equalities are inconsistent with the strict inequalities).
template <typename T>
std::optional<MyVector<T>> get_face_interior(MyMatrix<T> const &FAC,
                                              Face const &f_inc,
                                              std::ostream &os) {
  try {
    int nbIneq = FAC.rows();
    int n_eq = f_inc.count();
    if (nbIneq - n_eq == 0) {
      // Pure linear subspace; pick the origin (does not count as a
      // real cell of the L-type complex).
      return {};
    }
    if (n_eq == 0) {
      MyVector<T> v = GetSpaceInteriorPoint_Basic(FAC, os);
      return v;
    }
    auto split = split_fac(FAC, f_inc);
    MyVector<T> v = GetSpaceInteriorPoint(split.first, split.second, os);
    return v;
  } catch (TerminalException const &) {
    return {};
  }
}

//
// Enumerate the cell orbits of one top-dimensional L-type cone modulo
// the lattice stabiliser GRPperm. We walk the face lattice top-down
// using a BFS by f_inc bit count, quotienting via OptCanonicalImage.
//
// The face lattice is bounded by 2^|FAC| faces; for the small cones
// arising in low-dim iso-Delaunay enumeration this is manageable.
//

template <typename T, typename Tint, typename Tgroup>
std::vector<LtypeCell<T>>
enumerate_face_orbits_for_top_cell(int i_top, LinSpaceMatrix<T> const &LinSpa,
                                    MyMatrix<T> const &FAC,
                                    Tgroup const &GRPperm, std::ostream &os) {
  int nbIneq = FAC.rows();
  int dimSpace = LinSpa.ListMat.size();
  std::vector<LtypeCell<T>> out;
  std::unordered_set<Face> seen_canonical;
  std::vector<Face> frontier;
  // Level 0: the full top-dim cone (no tight inequalities).
  Face empty_f(nbIneq);
  Face empty_can = GRPperm.OptCanonicalImage(empty_f);
  seen_canonical.insert(empty_can);
  frontier.push_back(empty_f);
  // Build the level-0 cell.
  {
    std::optional<MyVector<T>> opt_inner = get_face_interior(FAC, empty_f, os);
    if (!opt_inner) {
      std::cerr << "LTCOMP: failed to find an interior point for the "
                   "top-dim cell of i_top="
                << i_top << "\n";
      throw TerminalException{1};
    }
    MyVector<T> inner_ray = *opt_inner;
    MyMatrix<T> inner_gram = gram_from_tspace_vec(inner_ray, LinSpa.ListMat);
    bool pd = IsPositiveDefinite(inner_gram, os);
    auto split = split_fac(FAC, empty_f);
    out.push_back(LtypeCell<T>{empty_f, i_top, std::move(split.first),
                                std::move(split.second),
                                std::move(inner_ray), std::move(inner_gram),
                                pd});
  }
  // Descend: at each step, for each face in the frontier, generate all
  // covers (one extra tight inequality) and keep the new ones modulo
  // GRPperm and modulo "same face = same affine span" (a cover is only
  // kept if its dim is strictly less than the parent's, otherwise it
  // describes the same cell).
  while (!frontier.empty()) {
    std::vector<Face> next_frontier;
    for (auto const &f_inc : frontier) {
      int dim_parent;
      {
        auto split = split_fac(FAC, f_inc);
        dim_parent = dimSpace - RankMat(split.second);
      }
      if (dim_parent <= 1) {
        // The cell is at most 1-dimensional; its only "facet" is the
        // origin (excluded from the L-type complex).
        continue;
      }
      for (int j = 0; j < nbIneq; j++) {
        if (f_inc[j] == 1)
          continue;
        Face f_new = f_inc;
        f_new[j] = 1;
        auto split_new = split_fac(FAC, f_new);
        int rank_new = RankMat(split_new.second);
        int dim_new = dimSpace - rank_new;
        if (dim_new >= dim_parent)
          continue;  // Tightening this inequality is implied by f_inc.
        if (dim_new < 1)
          continue;
        // Saturate f_new: add all FAC indices whose row lies in the
        // span of the current equalities, so that two equivalent
        // descriptions of the same face collapse.
        Face f_sat = f_new;
        for (int k = 0; k < nbIneq; k++) {
          if (f_sat[k] == 1)
            continue;
          MyVector<T> Wk = GetMatrixRow(FAC, k);
          MyMatrix<T> Augmented(split_new.second.rows() + 1, dimSpace);
          for (int r = 0; r < split_new.second.rows(); r++) {
            for (int c = 0; c < dimSpace; c++) {
              Augmented(r, c) = split_new.second(r, c);
            }
          }
          for (int c = 0; c < dimSpace; c++) {
            Augmented(split_new.second.rows(), c) = Wk(c);
          }
          if (RankMat(Augmented) == rank_new) {
            f_sat[k] = 1;
          }
        }
        Face f_can = GRPperm.OptCanonicalImage(f_sat);
        if (seen_canonical.contains(f_can))
          continue;
        seen_canonical.insert(f_can);
        std::optional<MyVector<T>> opt_inner =
            get_face_interior(FAC, f_sat, os);
        if (!opt_inner)
          continue;  // Face is empty (inconsistent strict + equalities).
        MyVector<T> inner_ray = *opt_inner;
        MyMatrix<T> inner_gram =
            gram_from_tspace_vec(inner_ray, LinSpa.ListMat);
        bool pd = IsPositiveDefinite(inner_gram, os);
        auto split_keep = split_fac(FAC, f_sat);
        out.push_back(LtypeCell<T>{f_sat, i_top, std::move(split_keep.first),
                                    std::move(split_keep.second),
                                    std::move(inner_ray),
                                    std::move(inner_gram), pd});
        next_frontier.push_back(f_sat);
      }
    }
    frontier = std::move(next_frontier);
  }
  return out;
}

//
// Full enumeration: build per-top-cell face orbits, deduplicate across
// top-cells using the iso-Delaunay adjacency structure, then group by
// dimension idx and (optionally) filter to positive-definite cells.
//

template <typename T, typename Tint, typename Tgroup>
FullLtypeComplexEnumeration<T>
get_full_ltype_complex_enumeration_from_ltot(
    std::vector<DatabaseEntry_Serial<IsoDelaunayDomain_Obj<T, Tint, Tgroup>,
                                     IsoDelaunayDomain_AdjO<T, Tint>>> const
        &l_tot,
    LinSpaceMatrix<T> const &LinSpa, LtypeComplexOptions const &opts,
    std::ostream &os) {
  int dimSpace = LinSpa.ListMat.size();
  os << "LTCOMP: |l_tot|=" << l_tot.size() << "\n";
  int n_cell = l_tot.size();
  std::vector<std::vector<LtypeCell<T>>> per_cell;
  per_cell.resize(n_cell);
  // Cache of FAC per top-cell, used for cross-cell facet identification.
  std::vector<MyMatrix<T>> cached_FAC(n_cell);
  for (int i_top = 0; i_top < n_cell; i_top++) {
    auto const &eCell = l_tot[i_top];
    std::vector<FullAdjInfo<T>> ListIneq_all =
        ComputeDefiningIneqIsoDelaunayDomain<T, Tgroup>(
            eCell.x.DT_gram.DT, LinSpa.ListLineMat, os);
    MyMatrix<T> FAC_full = GetFACineq(ListIneq_all);
    MyMatrix<T> FAC_extend = AddFirstZeroColumn(FAC_full);
    std::vector<int> ListIrred = get_non_redundant_indices(FAC_extend, os);
    int nbIrred = ListIrred.size();
    MyMatrix<T> FAC(nbIrred, dimSpace);
    for (int i = 0; i < nbIrred; i++) {
      AssignMatrixRow(FAC, i, ListIneq_all[ListIrred[i]].eIneq);
    }
    cached_FAC[i_top] = FAC;
    Tgroup GRPperm = get_fac_permutation_group<T, Tint, Tgroup>(
        LinSpa, eCell.x.DT_gram.GramMat, FAC, os);
    os << "LTCOMP: i_top=" << i_top << " nbIrred=" << nbIrred
       << " |GRPperm|=" << GRPperm.size() << "\n";
    per_cell[i_top] = enumerate_face_orbits_for_top_cell<T, Tint, Tgroup>(
        i_top, LinSpa, FAC, GRPperm, os);
    os << "LTCOMP: i_top=" << i_top
       << " local face orbits=" << per_cell[i_top].size() << "\n";
  }
  // Deduplicate codim-1 facets across top cells using the iso-Delaunay
  // adjacency structure: each adjacency identifies one facet of A with
  // one facet of B, so we drop B's copy of the shared facet.
  std::vector<std::vector<bool>> drop(n_cell);
  for (int i_top = 0; i_top < n_cell; i_top++) {
    drop[i_top].resize(per_cell[i_top].size(), false);
  }
  for (int i_top = 0; i_top < n_cell; i_top++) {
    auto const &eCell = l_tot[i_top];
    MyMatrix<T> const &FAC_A = cached_FAC[i_top];
    int nbIrred_A = FAC_A.rows();
    for (auto const &eAdj : eCell.ListAdj) {
      int j_top = eAdj.iOrb;
      if (j_top <= i_top)
        continue;  // Process each unordered pair once.
      MyVector<T> const &W_A = eAdj.x.V;
      // Find the index of this inequality among FAC_A's rows.
      int idx_A = -1;
      MyVector<T> W_A_red = RemoveFractionVector(W_A);
      for (int i = 0; i < nbIrred_A; i++) {
        MyVector<T> W_row = GetMatrixRow(FAC_A, i);
        if (RemoveFractionVector(W_row) == W_A_red) {
          idx_A = i;
          break;
        }
      }
      if (idx_A == -1)
        continue;  // Should not normally happen, but be safe.
      // Find the cell in per_cell[j_top] whose f_inc activates the
      // shared inequality (after mapping through eBigMat).
      MyMatrix<T> eBigMat_T =
          UniversalMatrixConversion<T, Tint>(eAdj.x.eBigMat);
      MyMatrix<T> MatSpace = matrix_in_t_space(eBigMat_T, LinSpa);
      MyVector<T> W_B = MatSpace * W_A;
      MyVector<T> W_B_red = RemoveFractionVector(W_B);
      MyMatrix<T> const &FAC_B = cached_FAC[j_top];
      int idx_B = -1;
      for (int i = 0; i < FAC_B.rows(); i++) {
        MyVector<T> W_row = GetMatrixRow(FAC_B, i);
        if (RemoveFractionVector(W_row) == W_B_red) {
          idx_B = i;
          break;
        }
      }
      if (idx_B == -1)
        continue;
      // Drop B's orbit that has f_inc = {idx_B} (and any descendant
      // that already has idx_B set, since those are identified with
      // descendants on A's side via the same eBigMat).
      for (size_t k = 0; k < per_cell[j_top].size(); k++) {
        if (per_cell[j_top][k].f_inc[idx_B] == 1) {
          drop[j_top][k] = true;
        }
      }
    }
  }
  // Assemble the final list, grouped by dim idx.
  std::map<int, std::vector<LtypeCell<T>>> grouped;
  for (int i_top = 0; i_top < n_cell; i_top++) {
    for (size_t k = 0; k < per_cell[i_top].size(); k++) {
      if (drop[i_top][k])
        continue;
      LtypeCell<T> const &cell = per_cell[i_top][k];
      if (opts.compute_full_dimensional && !cell.inner_is_pos_def)
        continue;
      int idx = dimSpace - cell.ListEqua.rows();
      grouped[idx].push_back(cell);
    }
  }
  FullLtypeComplexEnumeration<T> ret;
  ret.LinSpa = LinSpa;
  ret.opts = opts;
  for (auto &kv : grouped) {
    LtypeCellsAtLevel<T> lv;
    lv.idx = kv.first;
    lv.l_cells = std::move(kv.second);
    ret.levels.push_back(std::move(lv));
  }
  return ret;
}

//
// Output: write the enumeration in GAP form.
//

template <typename T>
void WriteLtypeCellGAP(std::ostream &os_out, LtypeCell<T> const &cell) {
  os_out << "rec(idx_in_tspace:=" << (cell.ListIneq.cols() - cell.ListEqua.rows())
         << ", i_top:=" << (cell.i_top + 1) << ", ListIneq:=";
  WriteMatrixGAP(os_out, cell.ListIneq);
  os_out << ", ListEqua:=";
  WriteMatrixGAP(os_out, cell.ListEqua);
  os_out << ", inner_ray:=" << StringVectorGAP(cell.inner_ray);
  os_out << ", inner_gram:=" << StringMatrixGAP(cell.inner_gram);
  os_out << ", inner_is_pos_def:=" << GAP_logical(cell.inner_is_pos_def);
  os_out << ")";
}

template <typename T>
void WriteFullLtypeComplexEnumerationGAP(
    std::ostream &os_out, FullLtypeComplexEnumeration<T> const &flce) {
  os_out << "return rec(ListMat:=";
  WriteListMatrixGAP(os_out, flce.LinSpa.ListMat);
  os_out << ", levels:=[";
  int n_level = flce.levels.size();
  for (int i_level = 0; i_level < n_level; i_level++) {
    if (i_level > 0)
      os_out << ",\n";
    LtypeCellsAtLevel<T> const &lv = flce.levels[i_level];
    os_out << "rec(idx:=" << lv.idx << ", cells:=[";
    for (size_t k = 0; k < lv.l_cells.size(); k++) {
      if (k > 0)
        os_out << ",\n";
      WriteLtypeCellGAP(os_out, lv.l_cells[k]);
    }
    os_out << "])";
  }
  os_out << "]);\n";
}

// Per-cell summary that matches the spec ("space of quadratic forms",
// "inequalities that define the cone", "inner ray").
template <typename T>
struct LtypeCellSummary {
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<T> ListIneq;
  MyMatrix<T> ListEqua;
  MyVector<T> inner_ray;
};

template <typename T>
LtypeCellSummary<T>
get_ltype_cell_summary(FullLtypeComplexEnumeration<T> const &flce, int i_level,
                       int i_cell) {
  LtypeCell<T> const &cell = flce.levels[i_level].l_cells[i_cell];
  return {flce.LinSpa.ListMat, cell.ListIneq, cell.ListEqua, cell.inner_ray};
}

// clang-format off
#endif  // SRC_DELAUNAY_LT_COMPLEX_H_
// clang-format on
