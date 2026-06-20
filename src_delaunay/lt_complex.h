// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_LT_COMPLEX_H_
#define SRC_DELAUNAY_LT_COMPLEX_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "Tspace_StabEquiInv.h"
#include "MatrixGroupNest.h"
#include "triples.h"
#include "boost_serialization.h"
#include <boost/dynamic_bitset/serialization.hpp>
#include <fstream>
#include <set>
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

  This mirrors the design of src_perfect/perfect_complex.h, but for
  iso-Delaunay (L-type) domains:

  * Top-dimensional cells come from the iso-Delaunay enumeration.
  * Lower-dimensional cells are obtained by walking down the faces of
    the top-dimensional cones and identifying equivalent faces using
    the triples machinery (src_polydecomp/triples.h).

  Per cell we expose:
  * The space of quadratic forms (LinSpa.ListMat, shared)
  * The inequalities (and equalities) defining the cell in the t-space
  * An inner ray (interior point of the cell, in t-space coordinates)

  The action used is the action of the lattice automorphisms on the
  t-space: a group element g acts on t-space coordinates c by
    c |-> matrix_in_t_space(g, LinSpa).transpose() * c
  which is the same convention as PermutationBuilder uses for the
  shortest-vector action in perfect_complex.h.
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
// Top-dimensional cell info
//

template <typename T, typename Tint_impl, typename Tgroup_impl>
struct LtypeInfoForComplex {
  using Tint = T;
  using TintLattice = Tint_impl;
  using Tgroup = Tgroup_impl;
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  MyMatrix<T> GramMat;
  MyMatrix<T> EXT;
  MyMatrix<T> ListIneq;
  std::optional<PreImagerElementContainer<T, Telt, TintGroup>> opt_pre_imager;
  Tgroup GRP_ext;
  std::vector<sing_adj<T>> l_sing_adj;
  MyMatrix<T> find_matrix(Telt const &x, [[maybe_unused]] std::ostream &os) const {
    if (opt_pre_imager) {
      PreImagerElementContainer<T, Telt, TintGroup> const &pre_imager =
          *opt_pre_imager;
      std::optional<MyMatrix<T>> opt = pre_imager.get_preimage(x);
      MyMatrix<T> M = unfold_opt(opt, "The element should belong to the group");
#ifdef SANITY_CHECK_LT_COMPLEX
      Telt x_img = get_elt_from_matrix<T, Telt>(M, EXT, os);
      if (x_img != x) {
        std::cerr << "LTCOMP: Inconsistency in pre_image computation\n";
        throw TerminalException{1};
      }
#endif
      return M;
    }
    MyMatrix<T> M = FindTransformation(EXT, EXT, x);
#ifdef SANITY_CHECK_LT_COMPLEX
    Telt x_img = get_elt_from_matrix<T, Telt>(M, EXT, os);
    if (x_img != x) {
      std::cerr << "LTCOMP: Inconsistency in pre_image computation\n";
      throw TerminalException{1};
    }
#endif
    return M;
  }
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void save(Archive &ar, LtypeInfoForComplex<T, Tint, Tgroup> const &val,
                 [[maybe_unused]] const unsigned int version) {
  bool has_pre_imager = val.opt_pre_imager.has_value();
  ar &make_nvp("GramMat", val.GramMat);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("ListIneq", val.ListIneq);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("l_sing_adj", val.l_sing_adj);
  ar &make_nvp("has_pre_imager", has_pre_imager);
  if (has_pre_imager) {
    std::vector<MyMatrix<T>> l_matr = val.opt_pre_imager->get_list_matr_gens();
    ar &make_nvp("l_matr", l_matr);
  }
}

template <class Archive, typename T, typename Tint, typename Tgroup>
inline void load(Archive &ar, LtypeInfoForComplex<T, Tint, Tgroup> &val,
                 [[maybe_unused]] const unsigned int version) {
  bool has_pre_imager = false;
  ar &make_nvp("GramMat", val.GramMat);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("ListIneq", val.ListIneq);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("l_sing_adj", val.l_sing_adj);
  ar &make_nvp("has_pre_imager", has_pre_imager);
  if (has_pre_imager) {
    std::vector<MyMatrix<T>> l_matr;
    ar &make_nvp("l_matr", l_matr);
    if (l_matr.empty()) {
      std::cerr << "LTCOMP: Missing l_matr for preimage reconstruction\n";
      throw TerminalException{1};
    }
    using Telt = typename Tgroup::Telt;
    using TintGroup = typename Tgroup::Tint;
    std::vector<Telt> l_perm =
        get_list_elt_from_list_matrices<T, Telt>(l_matr, val.EXT, std::cerr);
    MyMatrix<T> id = IdentityMat<T>(val.EXT.cols());
    val.opt_pre_imager =
        PreImagerElementContainer<T, Telt, TintGroup>(l_matr, l_perm, id);
  } else {
    val.opt_pre_imager = {};
  }
}

template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, LtypeInfoForComplex<T, Tint, Tgroup> &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}
}  // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct LtypeComplexTopDimInfo {
  std::vector<LtypeInfoForComplex<T, Tint, Tgroup>> l_ltype;
  LinSpaceMatrix<T> LinSpa;
  LtypeComplexOptions opts;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar,
                      LtypeComplexTopDimInfo<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_ltype", val.l_ltype);
  ar &make_nvp("LinSpa", val.LinSpa);
  ar &make_nvp("opts", val.opts);
}
}  // namespace boost::serialization

//
// Build the top-dimensional info from the iso-Delaunay enumeration.
//
// l_tot is the output of EnumerateAndStore_Serial for iso-Delaunay
// domains; we recover EXT (extreme rays of each L-type cone), the
// group acting on EXT, and the adjacency information needed by the
// triples machinery.
//

template <typename T, typename Tint, typename Tgroup>
LtypeComplexTopDimInfo<T, Tint, Tgroup> generate_ltype_complex_top_dim_info(
    std::vector<DatabaseEntry_Serial<IsoDelaunayDomain_Obj<T, Tint, Tgroup>,
                                     IsoDelaunayDomain_AdjO<T, Tint>>> const
        &l_tot,
    LinSpaceMatrix<T> const &LinSpa, LtypeComplexOptions const &opts,
    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  std::vector<LtypeInfoForComplex<T, Tint, Tgroup>> l_ltype;
  size_t n_cell = l_tot.size();
#ifdef DEBUG_LT_COMPLEX
  os << "LTCOMP: generate_ltype_complex_top_dim_info n_cell=" << n_cell << "\n";
#endif
  for (size_t i_cell = 0; i_cell < n_cell; i_cell++) {
    auto const &eCell = l_tot[i_cell];
    // Inequalities defining this cell.
    std::vector<FullAdjInfo<T>> const &ListIneqRed = eCell.x.ListIneqRed;
    size_t nbIneq = ListIneqRed.size();
    int dimSpace = LinSpa.ListMat.size();
    MyMatrix<T> FAC(nbIneq, dimSpace);
    for (size_t i = 0; i < nbIneq; i++) {
      AssignMatrixRow(FAC, i, ListIneqRed[i].eIneq);
    }
    // Extreme rays of the L-type cone (in t-space coordinates).
    MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
    int n_ext = EXT.rows();
#ifdef DEBUG_LT_COMPLEX
    os << "LTCOMP: i_cell=" << i_cell << " nbIneq=" << nbIneq
       << " n_ext=" << n_ext << "\n";
#endif
    // Group acting on EXT.
    // We obtain matrix generators of the lattice stabilizer from the
    // (already computed) GRPperm acting on inequalities — but we need
    // their action on rays. We re-derive them from the t-space action
    // matrices coming from LINSPA_ComputeStabilizer.
    std::vector<MyMatrix<T>> ListGenTot =
        LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, eCell.x.DT_gram.GramMat,
                                                  os);
    std::vector<MyMatrix<T>> ListMatSpace;
    ListMatSpace.reserve(ListGenTot.size());
    for (auto &eGenTot : ListGenTot) {
      ListMatSpace.push_back(matrix_in_t_space(eGenTot, LinSpa));
    }
    std::vector<Telt> l_perm =
        get_list_elt_from_list_matrices<T, Telt>(ListMatSpace, EXT, os);
    Tgroup GRP_ext(l_perm, n_ext);
    std::optional<PreImagerElementContainer<T, Telt, TintGroup>> opt_pre_imager;
    if (RankMat(EXT) < EXT.cols()) {
      MyMatrix<T> id = IdentityMat<T>(EXT.cols());
      opt_pre_imager = PreImagerElementContainer<T, Telt, TintGroup>(
          ListMatSpace, l_perm, id);
    }
    // Adjacency entries: convert each iso-Delaunay adjacency into a
    // sing_adj describing which face of EXT is the corresponding facet.
    std::vector<sing_adj<T>> l_sing_adj;
    l_sing_adj.reserve(eCell.ListAdj.size());
    for (auto &eAdj : eCell.ListAdj) {
      size_t jCone = static_cast<size_t>(eAdj.iOrb);
      MyVector<T> const &Vineq = eAdj.x.V;
      // Face of EXT lying on this facet hyperplane.
      Face f_ext(n_ext);
      for (int i_ext = 0; i_ext < n_ext; i_ext++) {
        MyVector<T> Vray = GetMatrixRow(EXT, i_ext);
        T scal = Vineq.dot(Vray);
        if (scal == 0) {
          f_ext[i_ext] = 1;
        }
      }
      // Action on the t-space induced by the lattice element eBigMat.
      MyMatrix<T> eBigMat_T = UniversalMatrixConversion<T, Tint>(eAdj.x.eBigMat);
      MyMatrix<T> MatSpace = matrix_in_t_space(eBigMat_T, LinSpa);
      // For the perfect-complex convention the eMat is the inverse of
      // the matrix coming out of the database (the database matrix maps
      // the source cell forward, the triples machinery uses the matrix
      // pulling the neighbour's rays back into the source's frame).
      MyMatrix<T> eMat = Inverse(MatSpace);
      sing_adj<T> adj{jCone, f_ext, std::move(eMat)};
      l_sing_adj.emplace_back(std::move(adj));
    }
    LtypeInfoForComplex<T, Tint, Tgroup> info{eCell.x.DT_gram.GramMat,
                                              std::move(EXT),
                                              std::move(FAC),
                                              std::move(opt_pre_imager),
                                              std::move(GRP_ext),
                                              std::move(l_sing_adj)};
    l_ltype.emplace_back(std::move(info));
  }
  return {std::move(l_ltype), LinSpa, opts};
}

//
// A cell of the L-type complex.
//

template <typename T, typename Tint, typename Tgroup>
struct LtypeFace {
  std::vector<triple<T>> l_triple;
  std::vector<MyMatrix<T>> l_gens;
  MyMatrix<T> EXT;
  Tgroup GRP_ext;
  MyMatrix<T> ListIneq;
  MyMatrix<T> ListEqua;
  MyVector<T> inner_ray;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, LtypeFace<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_triple", val.l_triple);
  ar &make_nvp("l_gens", val.l_gens);
  ar &make_nvp("EXT", val.EXT);
  ar &make_nvp("GRP_ext", val.GRP_ext);
  ar &make_nvp("ListIneq", val.ListIneq);
  ar &make_nvp("ListEqua", val.ListEqua);
  ar &make_nvp("inner_ray", val.inner_ray);
}
}  // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct LtypeFaces {
  std::vector<LtypeFace<T, Tint, Tgroup>> l_faces;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, LtypeFaces<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("l_faces", val.l_faces);
}
}  // namespace boost::serialization

//
// Helpers used to compute the inequalities / inner ray attached to a
// face of a top-dimensional cell.
//

template <typename T>
Face find_tight_facets(MyMatrix<T> const &FAC, MyMatrix<T> const &EXT,
                       Face const &f_ext) {
  int nbIneq = FAC.rows();
  Face f_ineq(nbIneq);
  std::vector<int> indices;
  indices.reserve(f_ext.count());
  for (int i = 0; i < static_cast<int>(f_ext.size()); i++) {
    if (f_ext[i] == 1) {
      indices.push_back(i);
    }
  }
  for (int j = 0; j < nbIneq; j++) {
    MyVector<T> Vineq = GetMatrixRow(FAC, j);
    bool all_zero = true;
    for (int i : indices) {
      MyVector<T> Vray = GetMatrixRow(EXT, i);
      if (Vineq.dot(Vray) != 0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      f_ineq[j] = 1;
    }
  }
  return f_ineq;
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
split_ineq_equa(MyMatrix<T> const &FAC, Face const &f_ineq) {
  int nbIneq = FAC.rows();
  int dim = FAC.cols();
  int n_eq = f_ineq.count();
  int n_in = nbIneq - n_eq;
  MyMatrix<T> ListIneq(n_in, dim);
  MyMatrix<T> ListEqua(n_eq, dim);
  int pos_in = 0;
  int pos_eq = 0;
  for (int j = 0; j < nbIneq; j++) {
    MyVector<T> Vineq = GetMatrixRow(FAC, j);
    if (f_ineq[j] == 1) {
      AssignMatrixRow(ListEqua, pos_eq, Vineq);
      pos_eq += 1;
    } else {
      AssignMatrixRow(ListIneq, pos_in, Vineq);
      pos_in += 1;
    }
  }
  return {std::move(ListIneq), std::move(ListEqua)};
}

template <typename T, typename Tint, typename Tgroup>
MyVector<T>
get_inner_ray_for_face(LtypeComplexTopDimInfo<T, Tint, Tgroup> const &lctdi,
                       triple<T> const &t, std::ostream &os) {
  LtypeInfoForComplex<T, Tint, Tgroup> const &cone =
      lctdi.l_ltype[t.iCone];
  MyMatrix<T> const &FAC = cone.ListIneq;
  Face f_ineq = find_tight_facets(FAC, cone.EXT, t.f_ext);
  int nbIneq = FAC.rows();
  int dim = FAC.cols();
  int n_eq = f_ineq.count();
  if (n_eq == 0) {
    return GetSpaceInteriorPoint_Basic(FAC, os);
  }
  int n_in = nbIneq - n_eq;
  MyMatrix<T> ListIneq(n_in, dim);
  MyMatrix<T> ListEqua(n_eq, dim);
  int pos_in = 0;
  int pos_eq = 0;
  for (int j = 0; j < nbIneq; j++) {
    MyVector<T> Vineq = GetMatrixRow(FAC, j);
    if (f_ineq[j] == 1) {
      AssignMatrixRow(ListEqua, pos_eq, Vineq);
      pos_eq += 1;
    } else {
      AssignMatrixRow(ListIneq, pos_in, Vineq);
      pos_in += 1;
    }
  }
  MyVector<T> inner = GetSpaceInteriorPoint(ListIneq, ListEqua, os);
  // Pull it back to the top-dim cell's own frame.
  MyMatrix<T> eMat_T = t.eMat;
  return eMat_T.transpose() * inner;
}

//
// First step: register each top-dim cell as a face cell.
//

template <typename T, typename Tint, typename Tgroup>
LtypeFaces<T, Tint, Tgroup> get_first_step_ltype_complex_enumeration(
    LtypeComplexTopDimInfo<T, Tint, Tgroup> const &lctdi, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  std::vector<LtypeFace<T, Tint, Tgroup>> l_faces;
  int dimSpace = lctdi.LinSpa.ListMat.size();
  size_t n_top = lctdi.l_ltype.size();
#ifdef DEBUG_LT_COMPLEX
  os << "LTCOMP: get_first_step_ltype_complex_enumeration n_top=" << n_top
     << "\n";
#endif
  for (size_t i_top = 0; i_top < n_top; i_top++) {
    LtypeInfoForComplex<T, Tint, Tgroup> const &top = lctdi.l_ltype[i_top];
    int n_ext = top.EXT.rows();
    Face f_ext(n_ext);
    for (int i_ext = 0; i_ext < n_ext; i_ext++) {
      f_ext[i_ext] = 1;
    }
    MyMatrix<T> eMatId = IdentityMat<T>(dimSpace);
    triple<T> t{i_top, f_ext, eMatId};
    std::vector<triple<T>> l_triple = {t};
    std::vector<MyMatrix<T>> l_gens;
    std::vector<Telt> l_elt = top.GRP_ext.SmallGeneratingSet();
    for (auto &ePermGen : l_elt) {
      MyMatrix<T> eMatrGen = top.find_matrix(ePermGen, os);
      l_gens.push_back(eMatrGen);
    }
    MyVector<T> inner_ray = GetSpaceInteriorPoint_Basic(top.ListIneq, os);
    MyMatrix<T> ListEqua(0, dimSpace);
    LtypeFace<T, Tint, Tgroup> face{
        std::move(l_triple), std::move(l_gens),  top.EXT,
        top.GRP_ext,         top.ListIneq,       std::move(ListEqua),
        std::move(inner_ray)};
    l_faces.emplace_back(std::move(face));
  }
  return {std::move(l_faces)};
}

//
// Step from a level to the next (lower) one.
//

template <typename T, typename Tint, typename Tgroup>
LtypeFaces<T, Tint, Tgroup>
compute_next_level_ltype(LtypeComplexTopDimInfo<T, Tint, Tgroup> const &lctdi,
                         LtypeFaces<T, Tint, Tgroup> const &level,
                         std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int dimSpace = lctdi.LinSpa.ListMat.size();
  std::vector<LtypeFace<T, Tint, Tgroup>> l_faces;
  auto find_matching_entry =
      [&](triple<T> const &t) -> std::optional<std::pair<int, MyMatrix<T>>> {
    int i_dom = 0;
    for (auto &face1 : l_faces) {
      std::optional<MyMatrix<T>> opt =
          test_triple_in_listtriple(lctdi.l_ltype, face1.l_triple, t, os);
      if (opt) {
        return std::pair<int, MyMatrix<T>>{i_dom, *opt};
      }
      i_dom += 1;
    }
    return {};
  };
  using Tfull_triple =
      std::pair<std::vector<triple<T>>, std::vector<MyMatrix<T>>>;
  auto get_initial_triple = [&](LtypeFace<T, Tint, Tgroup> const &face,
                                Face const &eIncd_sma) -> triple<T> {
    triple<T> const &t_big = face.l_triple[0];
    size_t iCone = t_big.iCone;
    int n_ext_big = lctdi.l_ltype[iCone].EXT.rows();
    Face f(n_ext_big);
    size_t index = 0;
    for (int i = 0; i < n_ext_big; i++) {
      if (t_big.f_ext[i] == 1) {
        if (eIncd_sma[index] == 1) {
          f[i] = 1;
        }
        index += 1;
      }
    }
    triple<T> t{iCone, f, t_big.eMat};
    return canonicalize_triple(lctdi.l_ltype, t, os);
  };
  auto f_insert = [&](triple<T> const &t,
                      std::optional<Tfull_triple> &opt_t) -> int {
    std::optional<std::pair<int, MyMatrix<T>>> opt = find_matching_entry(t);
    if (opt) {
      return opt->first;
    }
    if (!opt_t) {
      opt_t = get_spanning_list_triple(lctdi.l_ltype, t, os);
    }
    Tfull_triple const &pair = *opt_t;
    int iCone = t.iCone;
    int n_ext = t.f_ext.count();
    int n_ext_big = lctdi.l_ltype[iCone].EXT.rows();
    MyMatrix<T> EXT_face(n_ext, dimSpace);
    int pos = 0;
    for (int i_big = 0; i_big < n_ext_big; i_big++) {
      if (t.f_ext[i_big] == 1) {
        MyVector<T> V1 = GetMatrixRow(lctdi.l_ltype[iCone].EXT, i_big);
        MyVector<T> V2 = t.eMat.transpose() * V1;
        AssignMatrixRow(EXT_face, pos, V2);
        pos += 1;
      }
    }
    PermutationBuilder<T, Telt> builder(EXT_face);
    std::vector<Telt> l_gens_perm;
    for (auto &eMatrGen : pair.second) {
      Telt elt = builder.get_permutation(eMatrGen, os);
      l_gens_perm.push_back(std::move(elt));
    }
    Tgroup GRP_ext(l_gens_perm, n_ext);
    MyMatrix<T> const &FAC = lctdi.l_ltype[iCone].ListIneq;
    Face f_ineq = find_tight_facets(FAC, lctdi.l_ltype[iCone].EXT, t.f_ext);
    auto split = split_ineq_equa(FAC, f_ineq);
    MyMatrix<T> ListIneq_local = split.first;
    MyMatrix<T> ListEqua_local = split.second;
    MyVector<T> inner_local;
    int n_eq = f_ineq.count();
    if (n_eq == 0) {
      inner_local = GetSpaceInteriorPoint_Basic(ListIneq_local, os);
    } else {
      inner_local = GetSpaceInteriorPoint(ListIneq_local, ListEqua_local, os);
    }
    // Express inequalities, equalities and inner ray in the face's own
    // frame. Rays transform as v -> eMat^T * v, so to preserve the
    // pairing W^T v inequalities transform as W -> eMat^{-1} * W, i.e.
    // for the row layout, ListIneq -> ListIneq * eMat^{-T}.
    MyMatrix<T> eMatT = t.eMat.transpose();
    MyMatrix<T> eMat_inv = Inverse(t.eMat);
    MyMatrix<T> eMat_inv_T = eMat_inv.transpose();
    MyMatrix<T> ListIneq_pulled = ListIneq_local * eMat_inv_T;
    MyMatrix<T> ListEqua_pulled = ListEqua_local * eMat_inv_T;
    MyVector<T> inner_pulled = eMatT * inner_local;
    LtypeFace<T, Tint, Tgroup> face{
        pair.first,           pair.second,
        std::move(EXT_face),  std::move(GRP_ext),
        std::move(ListIneq_pulled), std::move(ListEqua_pulled),
        std::move(inner_pulled)};
    int i_dom = static_cast<int>(l_faces.size());
    l_faces.emplace_back(std::move(face));
    return i_dom;
  };
  size_t n_level = level.l_faces.size();
  for (size_t i = 0; i < n_level; i++) {
    LtypeFace<T, Tint, Tgroup> const &face = level.l_faces[i];
    // The facets of this face cell come from the inequalities of the
    // parent top-dim cone, restricted to the rays of the face. Each
    // inequality which is not already tight gives a candidate facet.
    triple<T> const &t_big = face.l_triple[0];
    size_t iCone = t_big.iCone;
    LtypeInfoForComplex<T, Tint, Tgroup> const &top = lctdi.l_ltype[iCone];
    MyMatrix<T> const &FAC_top = top.ListIneq;
    int nbIneq = FAC_top.rows();
    int n_ext_big = top.EXT.rows();
    Face f_ineq_active =
        find_tight_facets(FAC_top, top.EXT, t_big.f_ext);
    // List of relevant rays (the ones in the current face).
    std::vector<int> face_rays;
    face_rays.reserve(t_big.f_ext.count());
    for (int i_big = 0; i_big < n_ext_big; i_big++) {
      if (t_big.f_ext[i_big] == 1) {
        face_rays.push_back(i_big);
      }
    }
    std::set<Face> seen;
    for (int j = 0; j < nbIneq; j++) {
      if (f_ineq_active[j] == 1) {
        continue;
      }
      MyVector<T> Vineq = GetMatrixRow(FAC_top, j);
      Face eIncd_big(n_ext_big);
      int pos_in_face = 0;
      Face eIncd_sma(face_rays.size());
      bool any_nonzero = false;
      for (int i_big : face_rays) {
        MyVector<T> Vray = GetMatrixRow(top.EXT, i_big);
        if (Vineq.dot(Vray) == 0) {
          eIncd_big[i_big] = 1;
          eIncd_sma[pos_in_face] = 1;
        } else {
          any_nonzero = true;
        }
        pos_in_face += 1;
      }
      if (!any_nonzero || eIncd_sma.count() == 0) {
        continue;
      }
      if (seen.contains(eIncd_sma)) {
        continue;
      }
      seen.insert(eIncd_sma);
      triple<T> t = get_initial_triple(face, eIncd_sma);
      std::optional<Tfull_triple> opt_t;
      (void)f_insert(t, opt_t);
    }
  }
  return {std::move(l_faces)};
}

//
// Full enumeration of all cells of the L-type complex.
//

template <typename T, typename Tint, typename Tgroup>
struct FullLtypeComplexEnumeration {
  LtypeComplexTopDimInfo<T, Tint, Tgroup> lctdi;
  std::vector<LtypeFaces<T, Tint, Tgroup>> levels;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tint, typename Tgroup>
inline void serialize(Archive &ar,
                      FullLtypeComplexEnumeration<T, Tint, Tgroup> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("lctdi", val.lctdi);
  ar &make_nvp("levels", val.levels);
}
}  // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
FullLtypeComplexEnumeration<T, Tint, Tgroup>
full_ltype_complex_enumeration(LtypeComplexTopDimInfo<T, Tint, Tgroup> lctdi,
                               std::ostream &os) {
  std::vector<LtypeFaces<T, Tint, Tgroup>> levels;
  LtypeFaces<T, Tint, Tgroup> level0 =
      get_first_step_ltype_complex_enumeration(lctdi, os);
  levels.push_back(level0);
#ifdef DEBUG_LT_COMPLEX
  os << "LTCOMP: level 0 has " << level0.l_faces.size() << " faces\n";
#endif
  if (lctdi.opts.compute_full_dimensional) {
    return {std::move(lctdi), std::move(levels)};
  }
  while (true) {
    LtypeFaces<T, Tint, Tgroup> next =
        compute_next_level_ltype(lctdi, levels.back(), os);
#ifdef DEBUG_LT_COMPLEX
    os << "LTCOMP: next level has " << next.l_faces.size() << " faces\n";
#endif
    if (next.l_faces.empty()) {
      break;
    }
    levels.emplace_back(std::move(next));
  }
  return {std::move(lctdi), std::move(levels)};
}

template <typename T, typename Tint, typename Tgroup>
FullLtypeComplexEnumeration<T, Tint, Tgroup>
get_full_ltype_complex_enumeration_kernel(
    DataIsoDelaunayDomains<T, Tint, Tgroup> &data,
    LtypeComplexOptions const &opts, int max_runtime_second, std::ostream &os) {
  using Tdata = DataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::vector<Tout> l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  os << "LTCOMP: |l_tot|=" << l_tot.size() << "\n";
  LtypeComplexTopDimInfo<T, Tint, Tgroup> lctdi =
      generate_ltype_complex_top_dim_info<T, Tint, Tgroup>(
          l_tot, data_func.data.LinSpa, opts, os);
  return full_ltype_complex_enumeration(std::move(lctdi), os);
}

//
// Output: write the enumeration in GAP form.
//
// For each cell we emit a record with the space of quadratic forms,
// the inequalities and equalities that define the cone, and an inner
// ray.
//

template <typename T, typename Tint, typename Tgroup>
void WriteLtypeFaceGAP(std::ostream &os_out,
                       LinSpaceMatrix<T> const &LinSpa,
                       LtypeFace<T, Tint, Tgroup> const &face) {
  os_out << "rec(ListMat:=";
  WriteListMatrixGAP(os_out, LinSpa.ListMat);
  os_out << ", ListIneq:=";
  WriteMatrixGAP(os_out, face.ListIneq);
  os_out << ", ListEqua:=";
  WriteMatrixGAP(os_out, face.ListEqua);
  os_out << ", inner_ray:=" << StringVectorGAP(face.inner_ray);
  os_out << ", EXT:=";
  WriteMatrixGAP(os_out, face.EXT);
  os_out << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteFullLtypeComplexEnumerationGAP(
    std::ostream &os_out,
    FullLtypeComplexEnumeration<T, Tint, Tgroup> const &flce) {
  os_out << "return rec(ListMat:=";
  WriteListMatrixGAP(os_out, flce.lctdi.LinSpa.ListMat);
  os_out << ", levels:=[";
  int n_level = flce.levels.size();
  for (int i_level = 0; i_level < n_level; i_level++) {
    if (i_level > 0) {
      os_out << ",\n";
    }
    os_out << "[";
    int n_face = flce.levels[i_level].l_faces.size();
    for (int i_face = 0; i_face < n_face; i_face++) {
      if (i_face > 0) {
        os_out << ",\n";
      }
      WriteLtypeFaceGAP(os_out, flce.lctdi.LinSpa,
                        flce.levels[i_level].l_faces[i_face]);
    }
    os_out << "]";
  }
  os_out << "]);\n";
}

//
// Access: per-cell summary that mirrors the description in the spec.
//

template <typename T>
struct LtypeCellSummary {
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<T> ListIneq;
  MyMatrix<T> ListEqua;
  MyVector<T> inner_ray;
};

template <typename T, typename Tint, typename Tgroup>
LtypeCellSummary<T> get_ltype_cell_summary(
    FullLtypeComplexEnumeration<T, Tint, Tgroup> const &flce, int i_level,
    int i_face) {
  LtypeFace<T, Tint, Tgroup> const &face =
      flce.levels[i_level].l_faces[i_face];
  return {flce.lctdi.LinSpa.ListMat, face.ListIneq, face.ListEqua,
          face.inner_ray};
}

// clang-format off
#endif  // SRC_DELAUNAY_LT_COMPLEX_H_
// clang-format on
