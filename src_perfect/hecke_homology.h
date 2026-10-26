// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_HOMOLOGY_H_
#define SRC_PERFECT_HECKE_HOMOLOGY_H_

// clang-format off
#include "hecke_subgroups.h"
#include "hecke_chain_map.h"
#include "MAT_SparseMatrixMod.h"
#include <cstdio>
#include <fstream>
#include <string>
// clang-format on

/*
  Step B of the Hecke operators: the complex of the quotient by a finite
  index subgroup Gamma, its homology (rationally, or modulo a prime by
  sparse linear algebra) and the Hecke matrices on it.
 */

//
// Step B: the quotient by a finite index subgroup and its homology
//

/*
  For a cell F (level index, orbit iOrb), the orbits of Stab_G(F) on the
  cosets g_i Gamma. They correspond to the Gamma-orbits of cells in the
  G-orbit of F: the orbit of the coset i gives the cell F g_i. For each
  coset i we store the representative rep of its orbit and the sign of a
  transversal element t of Stab_G(F) with t g_rep Gamma = g_i Gamma, so
  that F g_i = sign (F g_rep) gamma for some gamma in Gamma.
  The cell F g_rep is orientable for Gamma if all Schreier generators of
  the stabilizer of the coset rep preserve the orientation of F.
 */
struct GammaOrbitSplit {
  std::vector<int> orbit_rep;
  std::vector<int> transversal_sign;
  std::vector<int> l_rep;
  std::vector<int> rep_position;
  std::vector<bool> l_orientable;
  // The size of the orbit of each representative.
  std::vector<size_t> l_orbit_size;
};

template <typename T, typename Tint, typename Tgroup>
GammaOrbitSplit get_gamma_orbit_split(
    FacePerfectComplex<T, Tint, Tgroup> const &face,
    std::vector<MyMatrix<T>> const &ListMat, CosetTableGamma<Tint> const &ct,
    [[maybe_unused]] std::ostream &os) {
  size_t n_coset = ct.size();
  size_t n_gen = face.l_gens.size();
  std::vector<std::vector<int>> perm(n_gen, std::vector<int>(n_coset));
  std::vector<int> sign_gen(n_gen);
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    MyMatrix<Tint> const &s = face.l_gens[i_gen];
    sign_gen[i_gen] = get_face_orientation(face.EXT, ListMat, face.or_info, s);
    GammaAction a = ct.gamma.get_action(s);
    for (size_t i = 0; i < n_coset; i++) {
      perm[i_gen][i] = ct.index_of_action(i, a);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      if (perm[i_gen][i] != ct.index_of(s * ct.l_coset[i])) {
        std::cerr << "HECKE: the action on the keys differs from the product\n";
        throw TerminalException{1};
      }
#endif
    }
  }
  GammaOrbitSplit gos;
  gos.orbit_rep = std::vector<int>(n_coset, -1);
  gos.transversal_sign = std::vector<int>(n_coset, 0);
  gos.rep_position = std::vector<int>(n_coset, -1);
  for (size_t i_start = 0; i_start < n_coset; i_start++) {
    if (gos.orbit_rep[i_start] != -1) {
      continue;
    }
    int rep = i_start;
    std::vector<int> l_orbit{rep};
    gos.orbit_rep[rep] = rep;
    gos.transversal_sign[rep] = 1;
    size_t pos = 0;
    while (pos < l_orbit.size()) {
      int i1 = l_orbit[pos];
      for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
        int i2 = perm[i_gen][i1];
        if (gos.orbit_rep[i2] == -1) {
          gos.orbit_rep[i2] = rep;
          gos.transversal_sign[i2] = sign_gen[i_gen] * gos.transversal_sign[i1];
          l_orbit.push_back(i2);
        }
      }
      pos += 1;
    }
    // The Schreier generators t_{i2}^{-1} s t_{i1} have sign
    // sign(t_{i2}) sign(s) sign(t_{i1}).
    bool is_orientable = true;
    for (auto &i1 : l_orbit) {
      for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
        int i2 = perm[i_gen][i1];
        int sign = gos.transversal_sign[i2] * sign_gen[i_gen] * gos.transversal_sign[i1];
        if (sign == -1) {
          is_orientable = false;
        }
      }
    }
    gos.rep_position[rep] = gos.l_rep.size();
    gos.l_rep.push_back(rep);
    gos.l_orientable.push_back(is_orientable);
    gos.l_orbit_size.push_back(l_orbit.size());
  }
  return gos;
}

/*
  The cells of the complex that form the basis of the quotient: all of
  them, the well rounded ones (the complex Pi / Pi_infinity, whose homology
  is the cohomology of Gamma), or the cells at infinity (their vectors do
  not span, the forms are degenerate: the subcomplex Pi_infinity, whose
  quotient by Gamma is the boundary of the Satake compactification).
 */
enum class CellSelection { All, WellRounded, AtInfinity };

template <typename T, typename Tint, typename Tgroup>
bool is_cell_selected(FacePerfectComplex<T, Tint, Tgroup> const &face, CellSelection const &selection) {
  if (selection == CellSelection::WellRounded) {
    return face.is_well_rounded;
  }
  if (selection == CellSelection::AtInfinity) {
    return !face.is_well_rounded;
  }
  return true;
}

/*
  The chain complex of the quotient by Gamma. The basis of each level is
  the list of orientable Gamma-cells among the selected ones.
  The boundary matrix of a level has one row per basis cell of the level
  and one column per basis cell of the next level.
 */
template <typename T, typename Tint> struct GammaComplex {
  CosetTableGamma<Tint> ct;
  CellSelection selection;
  std::vector<std::vector<GammaOrbitSplit>> splits;
  std::vector<std::vector<std::pair<int, int>>> basis;
  std::vector<std::vector<std::vector<int>>> basis_index;
  // The order of the stabilizer in Gamma of each basis cell F g_rep: the
  // order of Stab_G(F) divided by the size of the orbit of the coset rep.
  // It is 0 for a cell at infinity, whose stabilizer is infinite.
  std::vector<std::vector<int64_t>> stab_order;
  // The dense boundary matrices (for the rational homology) and the sparse
  // ones: boundary_sparse[index][i] is the list of the (column, value) of
  // the row i, sorted by column.
  std::vector<MyMatrix<T>> boundary;
  std::vector<std::vector<std::vector<std::pair<int, int>>>> boundary_sparse;
  // The position and sign of the cell F_{iOrb} M of level index in the
  // basis, if it is present.
  std::optional<std::pair<int, int>>
  cell_position(int const &index, int const &iOrb, MyMatrix<Tint> const &M) const {
    return cell_position_coset(index, iOrb, ct.index_of(M));
  }
  // The same for the cell F_{iOrb} g_{i0}, given by its coset index.
  std::optional<std::pair<int, int>>
  cell_position_coset(int const &index, int const &iOrb, int const &i0) const {
    GammaOrbitSplit const &gos = splits[index][iOrb];
    int rep = gos.orbit_rep[i0];
    int idx = basis_index[index][iOrb][rep];
    if (idx == -1) {
      return {};
    }
    std::pair<int, int> pair{idx, gos.transversal_sign[i0]};
    return pair;
  }
};

/*
  The boundary of the basis cells of the level index of gc_src, in the
  basis cells of the level index + 1 of gc_dst (the two complexes being
  for the same subgroup, with possibly different cell selections): one
  sparse row per cell of gc_src, sorted by column.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<std::vector<std::pair<int, int>>>
get_gamma_boundary_sparse(FullComplexEnumeration<T, Tint, Tgroup> const &fce, GammaComplex<T, Tint> const &gc_src,
                          GammaComplex<T, Tint> const &gc_dst, int const &index) {
  CosetTableGamma<Tint> const &ct = gc_src.ct;
  int n_row = gc_src.basis[index].size();
  std::vector<std::vector<std::pair<int, int>>> D_sparse(n_row);
  // The actions of the matrices of the boundary on the cosets, per orbit.
  int n_orb = fce.levels[index].l_faces.size();
  std::vector<std::vector<GammaAction>> ll_act(n_orb);
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
      ll_act[iOrb].push_back(ct.gamma.get_action(ebnd.M));
    }
  }
  for (int i_row = 0; i_row < n_row; i_row++) {
    int iOrb = gc_src.basis[index][i_row].first;
    int i_coset = gc_src.basis[index][i_row].second;
    std::map<int, int> map_row;
    std::vector<BoundEntry<Tint>> const &l_bound = fce.boundaries[index].ll_bound[iOrb].l_bound;
    for (size_t k = 0; k < l_bound.size(); k++) {
      BoundEntry<Tint> const &ebnd = l_bound[k];
      int i0 = ct.index_of_action(i_coset, ll_act[iOrb][k]);
      std::optional<std::pair<int, int>> opt = gc_dst.cell_position_coset(index + 1, ebnd.iOrb, i0);
      if (opt) {
        map_row[opt->first] += ebnd.sign * opt->second;
      }
    }
    for (auto &kv : map_row) {
      if (kv.second != 0) {
        D_sparse[i_row].push_back(kv);
      }
    }
  }
  return D_sparse;
}

template <typename T, typename Tint, typename Tgroup>
GammaComplex<T, Tint>
compute_gamma_complex(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                      CosetTableGamma<Tint> const &ct,
                      CellSelection const &selection,
                      std::vector<uint8_t> const &l_level_needed,
                      bool const &dense_boundary, std::ostream &os) {
  GammaComplex<T, Tint> gc;
  gc.ct = ct;
  gc.selection = selection;
  int n_levels = fce.levels.size();
  size_t n_coset = ct.size();
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  gc.splits = std::vector<std::vector<GammaOrbitSplit>>(n_levels);
  gc.basis = std::vector<std::vector<std::pair<int, int>>>(n_levels);
  gc.basis_index = std::vector<std::vector<std::vector<int>>>(n_levels);
  gc.stab_order = std::vector<std::vector<int64_t>>(n_levels);
  for (int index = 0; index < n_levels; index++) {
    if (l_level_needed[index] == 0) {
      continue;
    }
    int n_orb = fce.levels[index].l_faces.size();
    gc.basis_index[index] = std::vector<std::vector<int>>(n_orb, std::vector<int>(n_coset, -1));
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
      GammaOrbitSplit gos = get_gamma_orbit_split(face, ListMat, ct, os);
      if (is_cell_selected(face, selection)) {
        // The well rounded cells span the space, so Stab_G(F) acts
        // faithfully on their vectors and GRP_ext is the stabilizer. For a
        // cell at infinity, the stabilizer is infinite (l_gens generates
        // it, which is all the orbit split needs) and GRP_ext is only its
        // image on the vectors.
        int64_t order_stab = 0;
        if (face.is_well_rounded) {
          auto order_grp = face.GRP_ext.size();
          order_stab = UniversalScalarConversion<int64_t, decltype(order_grp)>(order_grp);
        }
        size_t n_rep = gos.l_rep.size();
        for (size_t i_rep = 0; i_rep < n_rep; i_rep++) {
          if (gos.l_orientable[i_rep]) {
            int rep = gos.l_rep[i_rep];
            int64_t orbit_size = gos.l_orbit_size[i_rep];
            if (face.is_well_rounded && order_stab % orbit_size != 0) {
              std::cerr << "HECKE: the orbit size " << orbit_size << " does not divide the order "
                        << order_stab << " of the stabilizer\n";
              throw TerminalException{1};
            }
            gc.basis_index[index][iOrb][rep] = gc.basis[index].size();
            gc.basis[index].push_back({iOrb, rep});
            gc.stab_order[index].push_back(face.is_well_rounded ? order_stab / orbit_size : 0);
          }
        }
      }
      gc.splits[index].push_back(std::move(gos));
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: compute_gamma_complex index=" << index << " n_orb=" << n_orb
       << " |basis|=" << gc.basis[index].size() << "\n";
#endif
  }
  //
  // The boundary matrices
  //
  for (int index = 0; index < n_levels - 1; index++) {
    if (l_level_needed[index] == 0 || l_level_needed[index + 1] == 0) {
      gc.boundary.push_back(MyMatrix<T>(0, 0));
      gc.boundary_sparse.push_back({});
      continue;
    }
    int n_row = gc.basis[index].size();
    int n_col = gc.basis[index + 1].size();
    std::vector<std::vector<std::pair<int, int>>> D_sparse = get_gamma_boundary_sparse(fce, gc, gc, index);
    if (dense_boundary) {
      MyMatrix<T> D = ZeroMatrix<T>(n_row, n_col);
      for (int i_row = 0; i_row < n_row; i_row++) {
        for (auto &kv : D_sparse[i_row]) {
          D(i_row, kv.first) = T(kv.second);
        }
      }
      gc.boundary.push_back(D);
    } else {
      gc.boundary.push_back(MyMatrix<T>(0, 0));
    }
    gc.boundary_sparse.emplace_back(std::move(D_sparse));
  }
#ifdef SANITY_CHECK_HECKE_OPERATORS
  for (int index = 0; index < n_levels - 2; index++) {
    if (l_level_needed[index] == 0 || l_level_needed[index + 1] == 0 || l_level_needed[index + 2] == 0) {
      continue;
    }
    bool is_zero = true;
    for (auto &row : gc.boundary_sparse[index]) {
      std::map<int, int> map_prod;
      for (auto &e1 : row) {
        for (auto &e2 : gc.boundary_sparse[index + 1][e1.first]) {
          map_prod[e2.first] += e1.second * e2.second;
        }
      }
      for (auto &kv : map_prod) {
        if (kv.second != 0) {
          is_zero = false;
        }
      }
    }
    if (!is_zero) {
      std::cerr << "HECKE: The product of the boundary matrices should be zero at index=" << index << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return gc;
}

/*
  A basis of harmonic representatives of the homology at a level:
  the vectors v with v D_index = 0 and v orthogonal to the rows of
  D_{index-1}. They span a complement of the boundaries in the cycles.
 */
template <typename T, typename Tint>
MyMatrix<T> get_harmonic_basis(GammaComplex<T, Tint> const &gc, int const &index) {
  int n_levels = gc.basis.size();
  int n_cell = gc.basis[index].size();
  int n_col = 0;
  if (index < n_levels - 1) {
    n_col += gc.boundary[index].cols();
  }
  if (index > 0) {
    n_col += gc.boundary[index - 1].rows();
  }
  if (n_cell == 0) {
    return MyMatrix<T>(0, 0);
  }
  if (n_col == 0) {
    return IdentityMat<T>(n_cell);
  }
  MyMatrix<T> E(n_cell, n_col);
  int pos = 0;
  if (index < n_levels - 1) {
    MyMatrix<T> const &D = gc.boundary[index];
    for (int i = 0; i < n_cell; i++) {
      for (int j = 0; j < D.cols(); j++) {
        E(i, pos + j) = D(i, j);
      }
    }
    pos += D.cols();
  }
  if (index > 0) {
    MyMatrix<T> const &D = gc.boundary[index - 1];
    for (int i = 0; i < n_cell; i++) {
      for (int j = 0; j < D.rows(); j++) {
        E(i, pos + j) = D(j, i);
      }
    }
  }
  return NullspaceMat(E);
}

template <typename T> struct HeckeHomologyLevel {
  HeckeHomologyLevel() = default;
  HeckeHomologyLevel(int const &_index, int const &_n_cell, MyMatrix<T> _harmonic, MyMatrix<T> _hecke_matrix,
                     std::string _hecke_failure)
      : index(_index), n_cell(_n_cell), harmonic(std::move(_harmonic)), hecke_matrix(std::move(_hecke_matrix)),
        hecke_failure(std::move(_hecke_failure)) {}
  int index = 0;
  int n_cell = 0;
  MyMatrix<T> harmonic;
  MyMatrix<T> hecke_matrix;
  // Nonempty when the Hecke matrix could not be computed (modulo p, when p
  // divides a denominator of the images); the homology is still given.
  std::string hecke_failure;
  // With the boundary map (see get_boundary_map_mod): the dimension of its
  // kernel (-1 if not computed), a basis of the kernel in the coordinates
  // of the homology, the Hecke matrix on it in that basis, and whether the
  // kernel is stable under the Hecke operator (it must be).
  int dim_ker_boundary = -1;
  MyMatrix<T> ker_boundary;
  MyMatrix<T> hecke_matrix_ker;
  bool ker_boundary_stable = true;
};

template <typename T> struct HeckeHomologyResult {
  std::string method; // The method of the chain map.
  int64_t prime = 0;  // If nonzero, the Hecke matrices are modulo this prime.
  std::string subgroup_type;
  bool only_well_rounded;
  size_t n_coset_gamma;
  size_t n_hecke_coset;
  std::vector<HeckeHomologyLevel<T>> l_level;
};

/*
  The actions on the coset keys of the matrices M of the entries of the
  images of a level: the coset of M g is then found from the key of g by
  an action, without computing the key of the product.
 */
using HeckeImageKeyActions = std::vector<std::vector<std::vector<GammaAction>>>;

template <typename T, typename Tint>
HeckeImageKeyActions get_hecke_image_key_actions(int const &index, GammaComplex<T, Tint> const &gc,
                                                 HeckeChainMap<T, Tint> const &hcm) {
  HeckeImageKeyActions actions;
  for (auto &l_chain : hcm.images[index]) {
    std::vector<std::vector<GammaAction>> l_act_orb;
    for (auto &chain : l_chain) {
      std::vector<GammaAction> l_act;
      for (auto &fe : chain) {
        l_act.push_back(gc.ct.gamma.get_action(fe.M));
      }
      l_act_orb.emplace_back(std::move(l_act));
    }
    actions.emplace_back(std::move(l_act_orb));
  }
  return actions;
}

/*
  The Hecke operator on the basis cell (iOrb, i) of level index:
    T(F g_i) = sum_k Phi(F, y_j) g'   with g_i x_k = y_j g'
  expressed in the basis of the level. With the actions of the entries
  (get_hecke_image_key_actions), the coset of M g' is the action of M on
  the key of g'.
 */
template <typename T, typename Tint>
std::map<int, T> get_hecke_image_cell_sparse(
    int const &index, int const &i_cell, GammaComplex<T, Tint> const &gc,
    HeckeChainMap<T, Tint> const &hcm, HeckeCosets<T> const &cosets_gamma,
    HeckeImageKeyActions const *actions = nullptr) {
  std::map<int, T> V;
  int iOrb = gc.basis[index][i_cell].first;
  int i_coset = gc.basis[index][i_cell].second;
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(gc.ct.l_coset[i_coset]);
  size_t n_hecke = cosets_gamma.size();
  for (size_t k = 0; k < n_hecke; k++) {
    MyMatrix<T> z = g_T * cosets_gamma.l_coset[k];
    std::pair<int, MyMatrix<T>> pair = hcm.cosets.decompose_ambient(z);
    MyMatrix<Tint> g_conj = UniversalMatrixConversion<Tint, T>(pair.second);
    std::vector<PerfectFaceEntry<T, Tint>> const &chain = hcm.images[index][iOrb][pair.first];
    std::optional<GammaKey> key_conj;
    if (actions) {
      key_conj = gc.ct.gamma.get_key(g_conj);
    }
    for (size_t e = 0; e < chain.size(); e++) {
      PerfectFaceEntry<T, Tint> const &fe = chain[e];
      std::optional<std::pair<int, int>> opt;
      if (actions) {
        auto iter = gc.ct.map_key.find(gc.ct.gamma.act(*key_conj, (*actions)[iOrb][pair.first][e]));
        if (iter == gc.ct.map_key.end()) {
          std::cerr << "HECKE: get_hecke_image_cell_sparse, the coset of M g was not found\n";
          throw TerminalException{1};
        }
        opt = gc.cell_position_coset(index, fe.iOrb, iter->second);
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (opt != gc.cell_position(index, fe.iOrb, fe.M * g_conj)) {
          std::cerr << "HECKE: get_hecke_image_cell_sparse, the action on the keys is wrong\n";
          throw TerminalException{1};
        }
#endif
      } else {
        opt = gc.cell_position(index, fe.iOrb, fe.M * g_conj);
      }
      if (opt) {
        V[opt->first] += fe.value * T(opt->second);
      }
    }
  }
  return V;
}

template <typename T, typename Tint, typename Tgroup>
MyVector<T> get_hecke_image_cell(
    int const &index, int const &i_cell, GammaComplex<T, Tint> const &gc,
    HeckeChainMap<T, Tint> const &hcm, HeckeCosets<T> const &cosets_gamma,
    [[maybe_unused]] FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeImageKeyActions const *actions = nullptr) {
  int n_cell = gc.basis[index].size();
  MyVector<T> V = ZeroVector<T>(n_cell);
  for (auto &kv : get_hecke_image_cell_sparse(index, i_cell, gc, hcm, cosets_gamma, actions)) {
    V(kv.first) += kv.second;
  }
  return V;
}

/*
  For a chain map commuting with the transposed differential (the Dual
  method), the images are chains of the dual complex W (the well rounded
  retract), and the complex of the quotient Gamma \ W is not the transpose
  of the one of Gamma \ Pi when the cells have stabilizers in Gamma: by
  counting the incident pairs over the orbits, its boundary is
    B(sigma, rho) = D(rho, sigma) |Gamma_sigma| / |Gamma_rho|.
  In the coordinates u = z S (S the diagonal of the stabilizer orders),
  in which the cycles and boundaries are those of the transposed complex
  of Gamma \ Pi, the operator z -> z T becomes u -> u S^{-1} T S: the entry
  (c, c') of the images is multiplied by |Gamma_c'| / |Gamma_c|.
 */
template <typename T, typename Tint>
T get_dual_stab_weight(GammaComplex<T, Tint> const &gc, int const &index, int const &c, int const &c_img) {
  return UniversalScalarConversion<T, int64_t>(gc.stab_order[index][c_img]) /
         UniversalScalarConversion<T, int64_t>(gc.stab_order[index][c]);
}

/*
  The cycle conditions of the cells of a level modulo p, one sparse row per
  cell c (a vector v is a cycle iff sum_c v_c row_c = 0), and the
  generators of the boundaries.
  For a chain map commuting with the boundary (Sharbly, Vertex) the cycles
  are the v with v D_index = 0 and the boundaries the rows of D_{index-1};
  for the transposed one (Dual), the v with D_{index-1} v = 0 and the
  columns of D_index.
 */
struct CycleConditionsMod {
  std::vector<SparseVectorMod> l_cycle_row;
  int n_next = 0;
  std::vector<SparseVectorMod> l_bnd;
};

template <typename T, typename Tint>
CycleConditionsMod get_cycle_conditions_mod(int const &index, GammaComplex<T, Tint> const &gc,
                                            bool const &transposed, int64_t const &p) {
  int n_levels = gc.basis.size();
  int n_cell = gc.basis[index].size();
  CycleConditionsMod cc;
  cc.l_cycle_row.resize(n_cell);
  auto f_push = [&](SparseVectorMod &v, int const &idx, int const &val) -> void {
    int64_t val_p = ResidueModPrime(val, p);
    if (val_p != 0) {
      v.push_back({idx, val_p});
    }
  };
  if (!transposed) {
    if (index < n_levels - 1) {
      cc.n_next = gc.basis[index + 1].size();
      for (int c = 0; c < n_cell; c++) {
        for (auto &e : gc.boundary_sparse[index][c]) {
          f_push(cc.l_cycle_row[c], e.first, e.second);
        }
      }
    }
    if (index > 0) {
      for (auto &row : gc.boundary_sparse[index - 1]) {
        SparseVectorMod v;
        for (auto &e : row) {
          f_push(v, e.first, e.second);
        }
        cc.l_bnd.emplace_back(std::move(v));
      }
    }
  } else {
    if (index > 0) {
      cc.n_next = gc.basis[index - 1].size();
      int n_prev = gc.boundary_sparse[index - 1].size();
      for (int r = 0; r < n_prev; r++) {
        for (auto &e : gc.boundary_sparse[index - 1][r]) {
          f_push(cc.l_cycle_row[e.first], r, e.second);
        }
      }
    }
    if (index < n_levels - 1) {
      int n_col_next = gc.basis[index + 1].size();
      std::vector<SparseVectorMod> l_col(n_col_next);
      for (int c = 0; c < n_cell; c++) {
        for (auto &e : gc.boundary_sparse[index][c]) {
          f_push(l_col[e.first], c, e.second);
        }
      }
      for (auto &v : l_col) {
        cc.l_bnd.emplace_back(std::move(v));
      }
    }
  }
  for (auto &v : cc.l_cycle_row) {
    std::sort(v.begin(), v.end());
  }
  for (auto &v : cc.l_bnd) {
    std::sort(v.begin(), v.end());
  }
  return cc;
}

#ifdef TRACK_INFO_HECKE_OPERATORS
/*
  Writes the cycle conditions of the level as an instance of the cycle
  lifting job of basic_common_cpp/src_matrix/Bench_SparseKernelMod: the
  line "p n_row n_col", then one line per cell with its number of entries
  and the pairs (column, value). The file is written under a temporary
  name and renamed when complete.
 */
inline void write_cycle_conditions_instance(int const &index, CycleConditionsMod const &cc, int64_t const &p,
                                            std::ostream &os) {
  std::string file = "HECKE_cycle_conditions_index" + std::to_string(index) + "_n" +
                     std::to_string(cc.l_cycle_row.size()) + "_m" + std::to_string(cc.n_next) + ".txt";
  std::string file_tmp = file + ".tmp";
  {
    std::ofstream os_file(file_tmp);
    os_file << p << " " << cc.l_cycle_row.size() << " " << cc.n_next << "\n";
    for (auto &v : cc.l_cycle_row) {
      os_file << v.size();
      for (auto &e : v) {
        os_file << " " << e.first << " " << e.second;
      }
      os_file << "\n";
    }
  }
  std::rename(file_tmp.c_str(), file.c_str());
  os << "HECKE: the cycle conditions of the level " << index << " are written to " << file << "\n";
}
#endif

/*
  The boundary map of the Satake compactification at a level k: the well
  rounded quotient complex is the quotient of the full complex by the
  subcomplex of the cells at infinity, and the connecting map
    delta : H_k(well rounded) -> H_{k+1}(at infinity)
  sends the class of a cycle z to the class of the part z D_wi of its
  boundary on the cells at infinity (a cycle there, since D^2 = 0). For
  GL_2 over Q at the symbol level, it is the boundary map {a, b} -> [b] - [a]
  of the modular symbols to the cusps (Cremona). The kernel of delta
  contains the interior cohomology (the image of the cohomology with
  compact support); it is Hecke stable, since delta is Hecke equivariant.
  The data: for each well rounded basis cell of the level k its boundary
  on the cells at infinity of the level k + 1 (D_wi), the boundaries of the
  cells at infinity of the level k in those of the level k + 1 (D_ii), and
  the number of cells at infinity of the level k + 1.
 */
struct BoundaryMapData {
  std::vector<SparseVectorMod> l_row_wi;
  std::vector<SparseVectorMod> l_row_ii;
  int n_inf_next = 0;
};

/*
  The kernel of delta on the homology classes of the basis cycles l_cycle,
  and the Hecke matrix on it. The class of delta(z) is given by its
  coordinates in the quotient of the chains at infinity of the level k + 1
  by the boundaries of the cells at infinity of the level k (a cycle has a
  zero class there iff it is a boundary).
 */
template <typename T>
void get_boundary_map_mod(BoundaryMapData const &bmd, std::vector<SparseVectorMod> const &l_cycle,
                          int64_t const &p, HeckeHomologyLevel<T> &level, [[maybe_unused]] std::ostream &os) {
  int dim_hom = l_cycle.size();
  SparseQuotientMod quot_inf(bmd.l_row_ii, bmd.n_inf_next, p);
  std::vector<SparseVectorMod> l_delta;
  for (auto &z : l_cycle) {
    std::map<int, int64_t> acc;
    for (auto &e : z) {
      for (auto &f : bmd.l_row_wi[e.first]) {
        int64_t &val = acc[f.first];
        val = ResidueModPrime(val + e.second * f.second, p);
      }
    }
    SparseVectorMod w;
    for (auto &kv : acc) {
      if (kv.second != 0) {
        w.push_back(kv);
      }
    }
    std::vector<int64_t> coord = quot_inf.coordinates(w);
    SparseVectorMod delta;
    for (size_t u = 0; u < coord.size(); u++) {
      if (coord[u] != 0) {
        delta.push_back({static_cast<int>(u), coord[u]});
      }
    }
    l_delta.emplace_back(std::move(delta));
  }
  // The kernel: the combinations x of the classes with sum_i x_i delta_i = 0.
  std::vector<SparseVectorMod> l_ker = SparseLeftKernelMod(l_delta, quot_inf.dim(), p);
  int dim_ker = l_ker.size();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: boundary map, dim_hom=" << dim_hom << " dim H(infinity)=" << quot_inf.dim()
     << " dim ker=" << dim_ker << "\n";
#endif
  level.dim_ker_boundary = dim_ker;
  level.ker_boundary = ZeroMatrix<T>(dim_ker, dim_hom);
  for (int i = 0; i < dim_ker; i++) {
    for (auto &e : l_ker[i]) {
      level.ker_boundary(i, e.first) = UniversalScalarConversion<T, int64_t>(e.second);
    }
  }
  int n_col_hecke = level.hecke_matrix.cols();
  if (dim_ker == 0 || n_col_hecke != dim_hom) {
    level.hecke_matrix_ker = MyMatrix<T>(0, 0);
    return;
  }
  // The rows of K extended by markers: the reduction of (y, 0), y = x T
  // for a row x of K, leaves (0, -c) with y = sum_i c_i K_i if y is in
  // the span of the rows of K.
  SparseEchelonMod echelon(p);
  for (int i = 0; i < dim_ker; i++) {
    SparseVectorMod v = l_ker[i];
    v.push_back({dim_hom + i, 1});
    echelon.insert(v);
  }
  level.hecke_matrix_ker = ZeroMatrix<T>(dim_ker, dim_ker);
  for (int i = 0; i < dim_ker; i++) {
    std::vector<int64_t> y(dim_hom, 0);
    for (auto &e : l_ker[i]) {
      for (int j = 0; j < dim_hom; j++) {
        int64_t t_ij = UniversalScalarConversion<int64_t, T>(level.hecke_matrix(e.first, j));
        y[j] = ResidueModPrime(y[j] + e.second * t_ij, p);
      }
    }
    SparseVectorMod v_y;
    for (int j = 0; j < dim_hom; j++) {
      if (y[j] != 0) {
        v_y.push_back({j, y[j]});
      }
    }
    for (auto &e : echelon.reduce(v_y)) {
      if (e.first < dim_hom) {
        level.ker_boundary_stable = false;
        level.hecke_matrix_ker = MyMatrix<T>(0, 0);
        return;
      }
      level.hecke_matrix_ker(i, e.first - dim_hom) = UniversalScalarConversion<T, int64_t>(ResidueModPrime(-e.second, p));
    }
  }
}

/*
  The Hecke operator on the homology of a level modulo a prime p, by
  sparse linear algebra: no matrix of the size of the complex is formed,
  which is what limits the rational computation. With B the space of the
  boundaries and Z the one of the cycles:
  - the cycles are the kernel of the constraints given by the columns of
    the matrix of the cycle conditions (SparseKernelMod, with the pivots
    chosen globally): a cycle is determined by its restriction to a set F
    of free cells, and every vector on F extends to a cycle,
  - H = Z / B is then the quotient of the vectors on F by the restrictions
    of the boundaries (SparseQuotientMod),
  - a basis of H is given by the basis cells of that quotient, lifted to
    cycles, and the image of a basis vector is restricted to F, its class
    in the quotient giving its coordinates.
  The entries of the Hecke matrix are residues mod p. Without with_hecke
  only the dimension is computed.
 */
template <typename T, typename Tint>
HeckeHomologyLevel<T> get_hecke_level_mod(int const &index, GammaComplex<T, Tint> const &gc,
                                 HeckeChainMap<T, Tint> const &hcm,
                                 HeckeCosets<T> const &cosets_gamma, int64_t const &p,
                                 bool const &with_hecke, BoundaryMapData const *bmd,
                                 [[maybe_unused]] std::ostream &os) {
  int n_cell = gc.basis[index].size();
  bool transposed = hcm.transposed;
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  CycleConditionsMod cc = get_cycle_conditions_mod(index, gc, transposed, p);
  // The constraints: the columns of the matrix of the cycle conditions.
  std::vector<SparseVectorMod> l_constraint(cc.n_next);
  for (int c = 0; c < n_cell; c++) {
    for (auto &e : cc.l_cycle_row[c]) {
      l_constraint[e.first].push_back({c, e.second});
    }
  }
  SparseKernelMod ker(l_constraint, n_cell, p);
  l_constraint.clear();
  std::vector<int> const &l_F = ker.free_columns();
  int n_F = l_F.size();
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: mod p index=" << index << " kernel of the cycle conditions, n_cell=" << n_cell
     << " n_next=" << cc.n_next << " rank=" << ker.rank() << " n_sparse_pivot=" << ker.n_sparse_pivot()
     << " n_dense_col=" << ker.n_dense_col() << " n_entry=" << ker.n_entry() << "|=" << time << "\n";
#endif
  // The restrictions of the boundaries to F, the cell l_F[u] at the
  // position u.
  std::vector<int> pos_F(n_cell, -1);
  for (int u = 0; u < n_F; u++) {
    pos_F[l_F[u]] = u;
  }
  std::vector<SparseVectorMod> l_bnd_F;
  for (auto &v : cc.l_bnd) {
    SparseVectorMod w;
    for (auto &e : v) {
      if (pos_F[e.first] != -1) {
        w.push_back({pos_F[e.first], e.second});
      }
    }
    if (!w.empty()) {
      std::sort(w.begin(), w.end());
      l_bnd_F.emplace_back(std::move(w));
    }
  }
  cc.l_bnd.clear();
  SparseQuotientMod quot_bnd(l_bnd_F, n_F, p);
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: mod p index=" << index << " quotient by the boundaries, |l_bnd_F|=" << l_bnd_F.size()
     << " n_F=" << n_F << " rank=" << quot_bnd.rank() << " n_sparse_pivot=" << quot_bnd.n_sparse_pivot()
     << " n_dense_col=" << quot_bnd.n_dense_col() << " n_entry=" << quot_bnd.n_entry() << "|=" << time << "\n";
#endif
  l_bnd_F.clear();
  // The basis columns of the quotient give the cells of the basis of H.
  std::vector<int> l_basis = quot_bnd.basis_columns();
  int dim_hom = l_basis.size();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: mod p index=" << index << " n_cell=" << n_cell << " rank(Z conditions)=" << ker.rank()
     << " rank(B)=" << quot_bnd.rank() << " dim_hom=" << dim_hom << "\n";
#endif
  // The cycles lifting the basis cells.
  std::vector<SparseVectorMod> l_cycle;
  for (int i = 0; i < dim_hom; i++) {
    l_cycle.emplace_back(ker.lift(l_basis[i]));
  }
  if (!with_hecke) {
    HeckeHomologyLevel<T> level{index, n_cell, MyMatrix<T>(dim_hom, 0), MyMatrix<T>(0, 0),
                                "no chain map at this level"};
    if (bmd) {
      get_boundary_map_mod(*bmd, l_cycle, p, level, os);
    }
    return level;
  }
  MyMatrix<T> hecke_matrix = ZeroMatrix<T>(dim_hom, dim_hom);
  HeckeImageKeyActions actions = get_hecke_image_key_actions(index, gc, hcm);
  for (int i = 0; i < dim_hom; i++) {
    SparseVectorMod const &v = l_cycle[i];
    // Its image, restricted to F, and its class modulo the boundaries.
    std::map<int, int64_t> acc;
    for (auto &e : v) {
      for (auto &kv : get_hecke_image_cell_sparse(index, e.first, gc, hcm, cosets_gamma, &actions)) {
        if (pos_F[kv.first] == -1) {
          continue;
        }
        T value = kv.second;
        if (transposed) {
          // The dual complex of the quotient (see get_dual_stab_weight).
          value *= get_dual_stab_weight(gc, index, e.first, kv.first);
        }
        std::optional<int64_t> opt = residue_mod_p(value, p);
        if (!opt) {
          // The homology is known, but not the Hecke matrix.
          HeckeHomologyLevel<T> level{index, n_cell, MyMatrix<T>(dim_hom, 0), MyMatrix<T>(0, 0),
                                      "a denominator of an image is divisible by p=" + std::to_string(p)};
          if (bmd) {
            get_boundary_map_mod(*bmd, l_cycle, p, level, os);
          }
          return level;
        }
        int64_t &val = acc[pos_F[kv.first]];
        val = ResidueModPrime(val + e.second * (*opt), p);
      }
    }
    SparseVectorMod w;
    for (auto &kv : acc) {
      if (kv.second != 0) {
        w.push_back(kv);
      }
    }
    std::vector<int64_t> coord = quot_bnd.coordinates(w);
    for (int j = 0; j < dim_hom; j++) {
      hecke_matrix(i, j) = UniversalScalarConversion<T, int64_t>(coord[j]);
    }
  }
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: mod p index=" << index << " Hecke matrix, dim_hom=" << dim_hom << "|=" << time << "\n";
#endif
  HeckeHomologyLevel<T> level{index, n_cell, MyMatrix<T>(dim_hom, 0), std::move(hecke_matrix), ""};
  if (bmd) {
    get_boundary_map_mod(*bmd, l_cycle, p, level, os);
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: mod p index=" << index << " boundary map, dim_ker=" << level.dim_ker_boundary << "|=" << time
       << "\n";
#endif
  }
  return level;
}

/*
  The cosets of Gamma \ Gamma x Gamma without any enumeration, when Gamma
  contains the principal congruence subgroup of level N and det(x) is
  coprime to N. The left cosets y_k G of G x G are then in bijection with
  the Gamma \ Gamma x Gamma, and the one in y_k G is y_k g for the coset
  g Gamma such that y_k g satisfies the congruences of Gamma x Gamma
  modulo N:
  - "Gamma0": e_n y ~ e_n (up to a unit), so the key e_n g^{-1} of g is
    the one of e_n y_k,
  - "Gamma1": e_n y = e_n x = c e_n, so the key of g is c^{-1} e_n y_k,
  - "Principal": y = x, so the key of g is y_k^{-1} x,
  - "Full": any g,
  and for the SL variants det(y_k g) = det(x). Nothing is returned when
  one of the conditions is not met.
 */
template <typename T, typename Tint>
std::optional<HeckeCosets<T>>
get_hecke_cosets_gamma_direct(HeckeChainMap<T, Tint> const &hcm, CosetTableGamma<Tint> const &ct,
                              std::function<bool(MyMatrix<T> const &)> const &f_member_gamma) {
  FiniteIndexSubgroup<Tint> const &gamma = ct.gamma;
  if (gamma.is_hermitian()) {
    // The congruences over O_K are not implemented here: the cosets are
    // enumerated.
    return {};
  }
  int n = gamma.n;
  int64_t m = gamma.get_n_mod();
  std::string type = gamma.base_type();
  MyMatrix<T> const &x = hcm.x;
  T det_x = DeterminantMat(x);
  if (!IsInteger(det_x)) {
    return {};
  }
  Tint det_x_int = UniversalScalarConversion<Tint, T>(det_x);
  if (GcdPair(det_x_int, gamma.N) != 1) {
    return {};
  }
  // The residue modulo m of a rational whose denominator is coprime to m.
  auto f_res = [&](T const &val) -> int64_t {
    Tint num = UniversalScalarConversion<Tint, T>(GetNumerator_z(val));
    Tint den = UniversalScalarConversion<Tint, T>(GetDenominator_z(val));
    int64_t num_r = UniversalScalarConversion<int64_t, Tint>(ResInt(num, gamma.N));
    int64_t den_r = UniversalScalarConversion<int64_t, Tint>(ResInt(den, gamma.N));
    if (m == 1) {
      return 0;
    }
    return (num_r * InverseModPrime(den_r, m)) % m;
  };
  int64_t c_inv = 1;
  if (type == "Gamma1") {
    for (int j = 0; j < n - 1; j++) {
      if (f_res(x(n - 1, j)) != 0) {
        return {};
      }
    }
    c_inv = InverseModPrime(f_res(x(n - 1, n - 1)), m);
  }
  HeckeCosets<T> hc;
  hc.x = x;
  hc.f_member = f_member_gamma;
  for (auto &y : hcm.cosets.l_coset) {
    GammaKey key;
    if (type == "Gamma0" || type == "Gamma1") {
      for (int j = 0; j < n; j++) {
        key.push_back((f_res(y(n - 1, j)) * c_inv) % m);
      }
      if (type == "Gamma0") {
        gamma.canonicalize_gamma0(key);
      }
    } else if (type == "Principal") {
      MyMatrix<T> z = Inverse(y) * x;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          key.push_back(f_res(z(i, j)));
        }
      }
    } else if (type != "Full") {
      return {};
    }
    if (gamma.is_sl()) {
      T ratio = det_x / DeterminantMat(y);
      key.push_back((ratio > 0) ? 1 : -1);
    }
    auto iter = ct.map_key.find(key);
    if (iter == ct.map_key.end()) {
      return {};
    }
    MyMatrix<T> g = UniversalMatrixConversion<T, Tint>(ct.l_coset[iter->second]);
    MyMatrix<T> y_new = y * g;
    hc.l_coset_inv.push_back(Inverse(y_new));
    hc.l_coset.push_back(y_new);
  }
  return hc;
}

template <typename T, typename Tint, typename Tgroup>
HeckeHomologyResult<T> compute_hecke_homology(
    FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeChainMap<T, Tint> const &hcm,
    FiniteIndexSubgroup<Tint> const &gamma, bool const &only_well_rounded,
    bool const &all_levels, bool const &boundary_map, int64_t const &prime, std::ostream &os) {
  [[maybe_unused]] bool const &dual = hcm.transposed;
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> l_gens_G = get_ambient_group_generators(fce, os);
  CosetTableGamma<Tint> ct = enumerate_cosets_gamma(gamma, l_gens_G, n, os);
#ifdef DEBUG_HECKE_OPERATORS
  {
    // The orbifold Euler characteristic of G on the well rounded retract:
    // sum over the levels of (-1)^k sum over the orbits 1 / |Stab_G(F)|.
    // It vanishes when the symmetric space has odd dimension (the Bianchi
    // groups), which checks the stabilizers and the orbits of the complex.
    T chi_orb(0);
    for (size_t index = 0; index < fce.levels.size(); index++) {
      for (auto &face : fce.levels[index].l_faces) {
        if (face.is_well_rounded) {
          auto order = face.GRP_ext.size();
          T val = T(1) / UniversalScalarConversion<T, decltype(order)>(order);
          chi_orb += (index % 2 == 0) ? val : -val;
        }
      }
    }
    os << "HECKE: orbifold Euler characteristic of G on the well rounded cells: " << chi_orb << "\n";
  }
#endif
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: enumerate_cosets_gamma|=" << time << "\n";
#endif
  std::function<bool(MyMatrix<T> const &)> f_member_gamma =
      [&LinSpa, &gamma](MyMatrix<T> const &g) -> bool {
    if (!is_in_ambient_group(LinSpa, g)) {
      return false;
    }
    MyMatrix<Tint> g_int = UniversalMatrixConversion<Tint, T>(g);
    return gamma.is_in_gamma(g_int);
  };
  // If Gamma contains the principal congruence subgroup of level N and
  // det(x) is coprime to N, the cosets of Gamma \ Gamma x Gamma are in
  // bijection with those of G \ G x G (strong approximation), so their
  // number is known and the enumeration with the many Schreier generators
  // of Gamma can stop early.
  size_t n_expected = 0;
  {
    T det_x = DeterminantMat(hcm.x);
    if (IsInteger(det_x)) {
      Tint det_int = UniversalScalarConversion<Tint, T>(det_x);
      if (GcdPair(det_int, gamma.N) == 1) {
        n_expected = hcm.cosets.size();
      }
    }
  }
  // The direct construction first; the enumeration with the Schreier
  // generators of Gamma otherwise (their number is about the index times
  // the number of generators of G, which is large for the large levels).
  std::optional<HeckeCosets<T>> opt_cosets = get_hecke_cosets_gamma_direct(hcm, ct, f_member_gamma);
  HeckeCosets<T> cosets_gamma;
  if (opt_cosets) {
    cosets_gamma = *opt_cosets;
  } else {
    std::vector<MyMatrix<Tint>> l_gens_gamma = get_schreier_generators(ct, l_gens_G, n, os);
    cosets_gamma = enumerate_hecke_cosets<T, Tint>(hcm.x, l_gens_gamma, f_member_gamma, os, n_expected);
  }
#ifdef SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE
  {
    // The enumeration without the early stop, with the many Schreier
    // generators of Gamma: its memory and time grow with the index.
    std::cerr << "HECKE: SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE, comparing the Hecke cosets of "
                 "Gamma with their full enumeration\n";
    std::vector<MyMatrix<Tint>> l_gens_gamma = get_schreier_generators(ct, l_gens_G, n, os);
    HeckeCosets<T> cosets_gamma_full =
        enumerate_hecke_cosets<T, Tint>(hcm.x, l_gens_gamma, f_member_gamma, os);
    if (n_expected > 0 && cosets_gamma_full.size() != n_expected) {
      std::cerr << "HECKE: the number of Hecke cosets of Gamma is " << cosets_gamma_full.size()
                << " and not " << n_expected << " as for G\n";
      throw TerminalException{1};
    }
    if (cosets_gamma.size() != cosets_gamma_full.size()) {
      std::cerr << "HECKE: the Hecke cosets of Gamma have " << cosets_gamma.size()
                << " elements instead of " << cosets_gamma_full.size() << "\n";
      throw TerminalException{1};
    }
    // Each coset found is one of those of the enumeration, and distinct.
    std::set<int> set_idx;
    for (auto &y : cosets_gamma.l_coset) {
      std::optional<std::pair<int, MyMatrix<T>>> opt = cosets_gamma_full.find_decompose(y);
      if (!opt) {
        std::cerr << "HECKE: a Hecke coset of Gamma is not in Gamma x Gamma\n";
        throw TerminalException{1};
      }
      set_idx.insert(opt->first);
    }
    if (set_idx.size() != cosets_gamma.size()) {
      std::cerr << "HECKE: two Hecke cosets of Gamma coincide\n";
      throw TerminalException{1};
    }
    std::cerr << "HECKE: SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE, done\n";
  }
#endif
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: enumerate_hecke_cosets (Gamma)|=" << time << "\n";
#endif
  // Only the levels of the chain map and their neighbors (for the
  // boundaries) are needed: the dense boundary matrices of the other
  // levels can be by far the largest objects.
  int n_levels_fce = fce.levels.size();
  std::vector<uint8_t> l_level_needed(n_levels_fce, 0);
  for (int index = 0; index < n_levels_fce; index++) {
    if (all_levels || hcm.is_level_computed(index)) {
      for (int k = std::max(0, index - 1); k <= std::min(n_levels_fce - 1, index + 1); k++) {
        l_level_needed[k] = 1;
      }
    }
  }
  bool dense_boundary = (prime == 0);
  CellSelection selection = only_well_rounded ? CellSelection::WellRounded : CellSelection::All;
  GammaComplex<T, Tint> gc = compute_gamma_complex(fce, ct, selection, l_level_needed, dense_boundary, os);
  // The boundary map of the Satake compactification: the cells at infinity
  // of the levels k and k + 1 for each level k computed.
  std::optional<GammaComplex<T, Tint>> opt_gc_inf;
  if (boundary_map) {
    if (prime == 0 || !only_well_rounded || hcm.transposed || fce.pctdi.pco.only_well_rounded) {
      std::cerr << "HECKE: the boundary map needs the homology modulo a prime of the well rounded complex, "
                   "a chain map commuting with the boundary (Sharbly or Vertex) and the full complex "
                   "(OnlyWellRounded = F)\n";
      throw TerminalException{1};
    }
    std::vector<uint8_t> l_level_needed_inf(n_levels_fce, 0);
    for (int index = 0; index < n_levels_fce; index++) {
      if (all_levels || hcm.is_level_computed(index)) {
        l_level_needed_inf[index] = 1;
        if (index + 1 < n_levels_fce) {
          l_level_needed_inf[index + 1] = 1;
        }
      }
    }
    opt_gc_inf = compute_gamma_complex(fce, ct, CellSelection::AtInfinity, l_level_needed_inf, false, os);
  }
#ifdef TRACK_INFO_HECKE_OPERATORS
  if (prime > 0) {
    for (int index = 0; index < n_levels_fce; index++) {
      if (all_levels || hcm.is_level_computed(index)) {
        write_cycle_conditions_instance(index, get_cycle_conditions_mod(index, gc, hcm.transposed, prime), prime, os);
      }
    }
  }
#endif
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: compute_gamma_complex|=" << time << "\n";
#endif
  HeckeHomologyResult<T> result;
  result.method = hcm.method;
  result.prime = prime;
  result.subgroup_type = gamma.type;
  result.only_well_rounded = only_well_rounded;
  result.n_coset_gamma = ct.size();
  result.n_hecke_coset = cosets_gamma.size();
  int n_levels = fce.levels.size();
  for (int index = 0; index < n_levels; index++) {
    // Without the chain map at this level, only the homology is given.
    bool with_hecke = hcm.is_level_computed(index);
    if (!with_hecke && !all_levels) {
      continue;
    }
    int n_cell = gc.basis[index].size();
    if (prime > 0) {
      std::optional<BoundaryMapData> opt_bmd;
      if (opt_gc_inf) {
        GammaComplex<T, Tint> const &gc_inf = *opt_gc_inf;
        BoundaryMapData bmd;
        auto f_mod = [&](std::vector<std::pair<int, int>> const &row) -> SparseVectorMod {
          SparseVectorMod v;
          for (auto &e : row) {
            int64_t val = ResidueModPrime(e.second, prime);
            if (val != 0) {
              v.push_back({e.first, val});
            }
          }
          return v;
        };
        if (index + 1 < n_levels) {
          bmd.n_inf_next = gc_inf.basis[index + 1].size();
          for (auto &row : get_gamma_boundary_sparse(fce, gc, gc_inf, index)) {
            bmd.l_row_wi.emplace_back(f_mod(row));
          }
          for (auto &row : gc_inf.boundary_sparse[index]) {
            bmd.l_row_ii.emplace_back(f_mod(row));
          }
        } else {
          bmd.l_row_wi.resize(gc.basis[index].size());
        }
        opt_bmd = std::move(bmd);
      }
      BoundaryMapData const *ptr_bmd = opt_bmd ? &(*opt_bmd) : nullptr;
      result.l_level.emplace_back(get_hecke_level_mod(index, gc, hcm, cosets_gamma, prime, with_hecke, ptr_bmd, os));
#ifdef TIMINGS_HECKE_OPERATORS
      os << "|HECKE: homology mod p level " << index << "|=" << time << "\n";
#endif
      continue;
    }
    MyMatrix<T> harmonic = get_harmonic_basis(gc, index);
    int dim_hom = harmonic.rows();
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: index=" << index << " n_cell=" << n_cell << " dim_hom=" << dim_hom << "\n";
#endif
    if (!with_hecke) {
      HeckeHomologyLevel<T> level{index, n_cell, std::move(harmonic), MyMatrix<T>(0, 0),
                                  "no chain map at this level"};
      result.l_level.emplace_back(std::move(level));
      continue;
    }
    MyMatrix<T> hecke_matrix = ZeroMatrix<T>(dim_hom, dim_hom);
    if (dim_hom > 0) {
      // The images of the basis cells, computed once.
      std::vector<MyVector<T>> l_img(n_cell);
      HeckeImageKeyActions actions = get_hecke_image_key_actions(index, gc, hcm);
      for (int i_cell = 0; i_cell < n_cell; i_cell++) {
        l_img[i_cell] = get_hecke_image_cell(index, i_cell, gc, hcm, cosets_gamma, fce, &actions);
        if (dual) {
          // The dual complex of the quotient (see get_dual_stab_weight).
          for (int j_cell = 0; j_cell < n_cell; j_cell++) {
            if (l_img[i_cell](j_cell) != 0) {
              l_img[i_cell](j_cell) *= get_dual_stab_weight(gc, index, i_cell, j_cell);
            }
          }
        }
      }
      MyMatrix<T> Gram = harmonic * harmonic.transpose();
      MyMatrix<T> GramInv = Inverse(Gram);
      for (int i_hom = 0; i_hom < dim_hom; i_hom++) {
        MyVector<T> w = ZeroVector<T>(n_cell);
        for (int i_cell = 0; i_cell < n_cell; i_cell++) {
          T coeff = harmonic(i_hom, i_cell);
          if (coeff != 0) {
            w += coeff * l_img[i_cell];
          }
        }
        // w is a cycle, decomposed as harmonic part plus boundary. The
        // boundaries are orthogonal to the harmonic vectors, so the
        // coordinates are the orthogonal projection.
        MyVector<T> proj = GramInv * (harmonic * w);
        for (int j_hom = 0; j_hom < dim_hom; j_hom++) {
          hecke_matrix(i_hom, j_hom) = proj(j_hom);
        }
#ifdef SANITY_CHECK_HECKE_OPERATORS
        // With the dual construction the chain map commutes with the
        // transposed differential: the cycles and boundaries are those
        // of the transposed complex.
        if (!dual && index < n_levels - 1) {
          MyVector<T> wD = gc.boundary[index].transpose() * w;
          if (!IsZeroVector(wD)) {
            std::cerr << "HECKE: The image of a cycle should be a cycle, index=" << index << "\n";
            throw TerminalException{1};
          }
        }
        if (dual && index > 0) {
          MyVector<T> wD = gc.boundary[index - 1] * w;
          if (!IsZeroVector(wD)) {
            std::cerr << "HECKE: The image of a dual cycle should be a dual cycle, index=" << index << "\n";
            throw TerminalException{1};
          }
        }
        MyVector<T> diff = w - harmonic.transpose() * proj;
        if (!IsZeroVector(diff)) {
          std::optional<MyVector<T>> opt;
          if (!dual) {
            if (index > 0) {
              opt = SolutionMat(gc.boundary[index - 1], diff);
            }
          } else {
            if (index < n_levels - 1) {
              MyMatrix<T> DT = gc.boundary[index].transpose();
              opt = SolutionMat(DT, diff);
            }
          }
          if (!opt) {
            std::cerr << "HECKE: The difference should be a boundary, index=" << index << "\n";
            throw TerminalException{1};
          }
        }
#endif
      }
    }
    HeckeHomologyLevel<T> level{index, n_cell, std::move(harmonic), std::move(hecke_matrix), ""};
    result.l_level.emplace_back(std::move(level));
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: homology level " << index << "|=" << time << "\n";
#endif
  }
  return result;
}

template <typename T>
void WriteHeckeHomologyResultGAP(std::ostream &os_out, HeckeHomologyResult<T> const &result) {
  os_out << "rec(method:=\"" << result.method << "\"";
  if (result.prime > 0) {
    os_out << ", prime:=" << result.prime;
  }
  os_out << ", subgroup_type:=\"" << result.subgroup_type << "\"";
  os_out << ", only_well_rounded:=" << GAP_logical(result.only_well_rounded);
  os_out << ", n_coset_gamma:=" << result.n_coset_gamma;
  os_out << ", n_hecke_coset:=" << result.n_hecke_coset;
  os_out << ",\nListLevel:=[";
  size_t n_level = result.l_level.size();
  for (size_t i = 0; i < n_level; i++) {
    if (i > 0) {
      os_out << ",\n";
    }
    HeckeHomologyLevel<T> const &level = result.l_level[i];
    os_out << "rec(index:=" << level.index;
    os_out << ", n_cell:=" << level.n_cell;
    os_out << ", dim_homology:=" << level.harmonic.rows();
    os_out << ", harmonic:=";
    WriteMatrixGAP(os_out, level.harmonic);
    os_out << ", HeckeMatrix:=";
    WriteMatrixGAP(os_out, level.hecke_matrix);
    if (!level.hecke_failure.empty()) {
      os_out << ", hecke_failure:=\"" << level.hecke_failure << "\"";
    }
    if (level.dim_ker_boundary >= 0) {
      os_out << ", dim_ker_boundary:=" << level.dim_ker_boundary;
      os_out << ", ker_boundary:=";
      WriteMatrixGAP(os_out, level.ker_boundary);
      os_out << ", HeckeMatrixKerBoundary:=";
      WriteMatrixGAP(os_out, level.hecke_matrix_ker);
      os_out << ", ker_boundary_stable:=" << GAP_logical(level.ker_boundary_stable);
    }
    os_out << ")";
  }
  os_out << "])";
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_HOMOLOGY_H_
// clang-format on
