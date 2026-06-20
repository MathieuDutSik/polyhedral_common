// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_LT_COMPLEX_H_
#define SRC_DELAUNAY_LT_COMPLEX_H_

// clang-format off
#include "InvariantVectorFamily.h"
#include "IsoDelaunayDomains.h"
#include "MAT_MatrixInt.h"
#include "SignatureSymmetric.h"
#include "Tspace_StabEquiInv.h"
#include "boost_serialization.h"
#include <boost/dynamic_bitset/serialization.hpp>
#include <deque>
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
  Enumeration of the full-dimensional cells of the L-type complex.

  Algorithm (specified by the user):
  * For each top-dimensional iso-Delaunay cell, compute the extreme rays
    of its L-type cone in the t-space.
  * Describe every cell of the L-type complex by its set of extreme rays
    (no bookkeeping about which top-dim cell or face index it came from).
  * For a cell with extreme rays R:
    - Compute the dual description of cone(R). Each facet of cone(R) is
      itself a cell, with its extreme rays being the subset of R lying
      on that facet.
    - Compute the sum of the extreme rays of a facet. The corresponding
      Gram matrix is positive definite iff that facet is "full
      dimensional" in the user's sense, i.e. its relative interior
      contains positive-definite forms.
    - Use that sum (and its invariant vector family) to test equivalence
      between cells via LINSPA_TestEquivalenceGramMatrix_SHV.
  * BFS over cells; only positive-definite cells are kept. Each new
    positive-definite cell is checked for equivalence with all existing
    cells of the same t-space dimension and added if it is new.

  Per cell we expose:
  * The space of quadratic forms (LinSpa.ListMat, shared)
  * The inequalities (and equalities) defining the cone
  * The inner ray (= sum of extreme rays, in t-space coordinates)
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
// A cell of the L-type complex (only positive-definite cells are
// instantiated here).
//

template <typename T>
struct LtypeCell {
  // Extreme rays of this cell in t-space, normalised to primitive
  // integer vectors with the first non-zero entry positive.
  MyMatrix<T> EXT;
  // Sum of EXT rows (= the "inner ray" used by the equivalence test).
  MyVector<T> inner_ray;
  // Gram matrix corresponding to inner_ray.
  MyMatrix<T> inner_gram;
  // Invariant vector family of inner_gram, used for the GL_n(Z)
  // equivalence test against other cells of the same dimension.
  MyMatrix<T> SHV_T;
  // Inequalities of the cone within its affine span (rows of length
  // dimSpace).
  MyMatrix<T> ListIneq;
  // Linear equalities defining the affine span of this cell inside
  // t-space (rows of length dimSpace; empty for top-dim cells).
  MyMatrix<T> ListEqua;
  // The t-space dimension of this cell (= rank of EXT = dim of affine
  // span).
  int idx;
};

namespace boost::serialization {
template <class Archive, typename T>
inline void serialize(Archive &ar, LtypeCell<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("inner_ray", val.inner_ray);
  ar &make_nvp("inner_gram", val.inner_gram);
  ar &make_nvp("SHV_T", val.SHV_T);
  ar &make_nvp("ListIneq", val.ListIneq);
  ar &make_nvp("ListEqua", val.ListEqua);
  ar &make_nvp("idx", val.idx);
}
}  // namespace boost::serialization

template <typename T>
struct LtypeCellsAtLevel {
  int idx;  // t-space dimension of the cells at this level
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

// Normalise each row of EXT to a primitive integer vector with the
// first non-zero entry positive. Two rays of a cone differ by a
// positive scaling iff their canonical forms agree.
template <typename T>
MyMatrix<T> canonicalize_ext_rows(MyMatrix<T> const &EXT) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  MyMatrix<T> out(n_row, n_col);
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> v = GetMatrixRow(EXT, i_row);
    MyVector<T> v_can = CanonicalizeVectorToInvertible(v);
    AssignMatrixRow(out, i_row, v_can);
  }
  return out;
}

// Equalities defining the affine span of EXT in t-space (orthogonal
// complement of the row span of EXT).
template <typename T>
MyMatrix<T> equalities_from_ext(MyMatrix<T> const &EXT) {
  // NullspaceTrMat(EXT) returns rows W with EXT * W = 0, i.e. EXT
  // each row dotted with W gives 0 — so W is orthogonal to every row
  // of EXT, hence in the row-span's orthogonal complement.
  return NullspaceTrMat(EXT);
}

//
// Build an LtypeCell from a set of extreme rays.
//
// For the cell to be considered, its inner Gram matrix must be
// positive definite (otherwise the cell is on the boundary of the PSD
// cone and is not full-dimensional in the user's sense).
//

template <typename T, typename Tint>
std::optional<LtypeCell<T>>
build_cell_from_ext(MyMatrix<T> const &EXT_in, LinSpaceMatrix<T> const &LinSpa,
                    std::ostream &os) {
  MyMatrix<T> EXT = canonicalize_ext_rows(EXT_in);
  // Drop duplicate rays that may arise from canonicalisation.
  {
    std::set<MyVector<T>> seen;
    int n = EXT.rows();
    std::vector<int> keep;
    keep.reserve(n);
    for (int i = 0; i < n; i++) {
      MyVector<T> v = GetMatrixRow(EXT, i);
      if (seen.insert(v).second) {
        keep.push_back(i);
      }
    }
    if (static_cast<int>(keep.size()) != n) {
      MyMatrix<T> EXT2(keep.size(), EXT.cols());
      for (size_t i = 0; i < keep.size(); i++) {
        AssignMatrixRow(EXT2, i, GetMatrixRow(EXT, keep[i]));
      }
      EXT = std::move(EXT2);
    }
  }
  int n_ext = EXT.rows();
  int dimSpace = LinSpa.ListMat.size();
  if (n_ext == 0) {
    return {};
  }
  MyVector<T> inner_ray = ZeroVector<T>(dimSpace);
  for (int i = 0; i < n_ext; i++) {
    inner_ray += GetMatrixRow(EXT, i);
  }
  MyMatrix<T> inner_gram = gram_from_tspace_vec(inner_ray, LinSpa.ListMat);
  if (!IsPositiveDefinite(inner_gram, os)) {
    return {};
  }
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(inner_gram, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int rank_ext = RankMat(EXT);
  MyMatrix<T> ListEqua = equalities_from_ext(EXT);
  MyMatrix<T> ListIneq;
  if (rank_ext == dimSpace) {
    ListIneq = DirectDualDescription_mat(EXT, os);
  } else {
    // For a lower-dim cone, the dual description is most naturally
    // computed inside the cone's affine span. We pick a basis of the
    // span, project EXT into it, dual-describe there, and lift each
    // facet inequality back to t-space.
    MyMatrix<T> Basis = NullspaceTrMat(ListEqua);
    // Basis rows form a basis of the row span of EXT (= the affine
    // span). Express each row of EXT in basis coordinates: solve
    // ext_row = coeffs * Basis.
    MyMatrix<T> EXT_basis(n_ext, rank_ext);
    for (int i = 0; i < n_ext; i++) {
      MyVector<T> r = GetMatrixRow(EXT, i);
      std::optional<MyVector<T>> opt = SolutionMat(Basis, r);
      MyVector<T> coeffs =
          unfold_opt(opt, "build_cell_from_ext: row not in span");
      AssignMatrixRow(EXT_basis, i, coeffs);
    }
    MyMatrix<T> FAC_basis = DirectDualDescription_mat(EXT_basis, os);
    int n_fac = FAC_basis.rows();
    ListIneq.resize(n_fac, dimSpace);
    for (int i = 0; i < n_fac; i++) {
      MyVector<T> w_basis = GetMatrixRow(FAC_basis, i);
      // The inequality W in t-space coordinates satisfies, for any
      // affine-span vector x = c * Basis: W . x = w_basis . c. We can
      // take W = Basis^T * w_basis.
      MyVector<T> W = Basis.transpose() * w_basis;
      AssignMatrixRow(ListIneq, i, W);
    }
  }
  LtypeCell<T> cell{std::move(EXT),       std::move(inner_ray),
                    std::move(inner_gram), std::move(SHV_T),
                    std::move(ListIneq),  std::move(ListEqua),
                    rank_ext};
  return cell;
}

// GL_n(Z) equivalence between two positive-definite cells, using their
// inner Gram matrices and SHV families.
template <typename T, typename Tgroup>
bool cells_are_gram_equivalent(LtypeCell<T> const &c1, LtypeCell<T> const &c2,
                                LinSpaceMatrix<T> const &LinSpa,
                                std::ostream &os) {
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(LinSpa, c1.inner_gram,
                                                      c2.inner_gram, c1.SHV_T,
                                                      c2.SHV_T, os);
  return opt.has_value();
}

//
// Full enumeration: BFS over cells, starting from the top-dim
// iso-Delaunay cells, descending through facets, keeping only
// positive-definite cells, deduplicating by GL_n(Z) equivalence of
// inner Gram matrices.
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
  std::vector<LtypeCell<T>> cells;
  std::deque<int> frontier;
  auto try_add = [&](MyMatrix<T> const &cell_EXT_in) -> bool {
    std::optional<LtypeCell<T>> opt =
        build_cell_from_ext<T, Tint>(cell_EXT_in, LinSpa, os);
    if (!opt)
      return false;
    LtypeCell<T> cell = std::move(*opt);
    for (auto const &existing : cells) {
      if (existing.idx != cell.idx)
        continue;
      if (cells_are_gram_equivalent<T, Tgroup>(cell, existing, LinSpa, os))
        return false;
    }
    cells.push_back(std::move(cell));
    frontier.push_back(static_cast<int>(cells.size()) - 1);
    return true;
  };
  for (size_t i_top = 0; i_top < l_tot.size(); i_top++) {
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
    MyMatrix<T> EXT_top = DirectDualDescription_mat(FAC, os);
    bool added = try_add(EXT_top);
    os << "LTCOMP: i_top=" << i_top << " nbIrred=" << nbIrred
       << " n_ext=" << EXT_top.rows() << " added=" << added << "\n";
  }
  // BFS through facets.
  while (!frontier.empty()) {
    int idx_cur = frontier.front();
    frontier.pop_front();
    MyMatrix<T> EXT_cur = cells[idx_cur].EXT;
    int n_ext = EXT_cur.rows();
    if (n_ext <= 1) {
      continue;
    }
    MyMatrix<T> EXT_for_dd;
    bool used_projection = false;
    MyMatrix<T> Basis;
    int rank_ext = cells[idx_cur].idx;
    if (rank_ext == dimSpace) {
      EXT_for_dd = EXT_cur;
    } else {
      // Project EXT into a basis of its affine span so that dual
      // description sees a full-dimensional cone.
      Basis = NullspaceTrMat(cells[idx_cur].ListEqua);
      EXT_for_dd.resize(n_ext, rank_ext);
      for (int i = 0; i < n_ext; i++) {
        MyVector<T> r = GetMatrixRow(EXT_cur, i);
        std::optional<MyVector<T>> opt = SolutionMat(Basis, r);
        MyVector<T> coeffs =
            unfold_opt(opt, "BFS: row not in affine span");
        AssignMatrixRow(EXT_for_dd, i, coeffs);
      }
      used_projection = true;
    }
    vectface vf = DirectDualDescription_vf(EXT_for_dd, os);
    for (auto const &face : vf) {
      int n_face = face.count();
      if (n_face == 0)
        continue;
      MyMatrix<T> facet_ext(n_face, dimSpace);
      int pos = 0;
      for (int i = 0; i < n_ext; i++) {
        if (face[i] == 1) {
          AssignMatrixRow(facet_ext, pos, GetMatrixRow(EXT_cur, i));
          pos += 1;
        }
      }
      (void)used_projection;
      (void)Basis;
      try_add(facet_ext);
    }
  }
  // Group by idx.
  std::map<int, std::vector<LtypeCell<T>>> grouped;
  for (auto &cell : cells) {
    grouped[cell.idx].push_back(std::move(cell));
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
  os_out << "rec(idx:=" << cell.idx << ", EXT:=";
  WriteMatrixGAP(os_out, cell.EXT);
  os_out << ", ListIneq:=";
  WriteMatrixGAP(os_out, cell.ListIneq);
  os_out << ", ListEqua:=";
  WriteMatrixGAP(os_out, cell.ListEqua);
  os_out << ", inner_ray:=" << StringVectorGAP(cell.inner_ray);
  os_out << ", inner_gram:=" << StringMatrixGAP(cell.inner_gram);
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
