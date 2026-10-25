// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_HOMOLOGY_H_
#define SRC_PERFECT_HECKE_HOMOLOGY_H_

// clang-format off
#include "hecke_subgroups.h"
#include "hecke_chain_map.h"
#include "MAT_SparseMatrixMod.h"
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
  The chain complex of the quotient by Gamma. The basis of each level is
  the list of orientable Gamma-cells (restricted to the well rounded ones
  if only_well_rounded is set, which gives the complex Pi / Pi_infinity).
  The boundary matrix of a level has one row per basis cell of the level
  and one column per basis cell of the next level.
 */
template <typename T, typename Tint> struct GammaComplex {
  CosetTableGamma<Tint> ct;
  bool only_well_rounded;
  std::vector<std::vector<GammaOrbitSplit>> splits;
  std::vector<std::vector<std::pair<int, int>>> basis;
  std::vector<std::vector<std::vector<int>>> basis_index;
  // The order of the stabilizer in Gamma of each basis cell F g_rep: the
  // order of Stab_G(F) divided by the size of the orbit of the coset rep.
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

template <typename T, typename Tint, typename Tgroup>
GammaComplex<T, Tint>
compute_gamma_complex(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                      CosetTableGamma<Tint> const &ct,
                      bool const &only_well_rounded,
                      std::vector<uint8_t> const &l_level_needed,
                      bool const &dense_boundary, std::ostream &os) {
  GammaComplex<T, Tint> gc;
  gc.ct = ct;
  gc.only_well_rounded = only_well_rounded;
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
      bool is_ok = true;
      if (only_well_rounded && !face.is_well_rounded) {
        is_ok = false;
      }
      if (is_ok) {
        // The well rounded cells span the space, so Stab_G(F) acts
        // faithfully on their vectors.
        auto order_grp = face.GRP_ext.size();
        int64_t order_stab = UniversalScalarConversion<int64_t, decltype(order_grp)>(order_grp);
        size_t n_rep = gos.l_rep.size();
        for (size_t i_rep = 0; i_rep < n_rep; i_rep++) {
          if (gos.l_orientable[i_rep]) {
            int rep = gos.l_rep[i_rep];
            int64_t orbit_size = gos.l_orbit_size[i_rep];
            if (order_stab % orbit_size != 0) {
              std::cerr << "HECKE: the orbit size " << orbit_size << " does not divide the order "
                        << order_stab << " of the stabilizer\n";
              throw TerminalException{1};
            }
            gc.basis_index[index][iOrb][rep] = gc.basis[index].size();
            gc.basis[index].push_back({iOrb, rep});
            gc.stab_order[index].push_back(order_stab / orbit_size);
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
      int iOrb = gc.basis[index][i_row].first;
      int i_coset = gc.basis[index][i_row].second;
      std::map<int, int> map_row;
      std::vector<BoundEntry<Tint>> const &l_bound = fce.boundaries[index].ll_bound[iOrb].l_bound;
      for (size_t k = 0; k < l_bound.size(); k++) {
        BoundEntry<Tint> const &ebnd = l_bound[k];
        int i0 = ct.index_of_action(i_coset, ll_act[iOrb][k]);
        std::optional<std::pair<int, int>> opt = gc.cell_position_coset(index + 1, ebnd.iOrb, i0);
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
  int index;
  int n_cell;
  MyMatrix<T> harmonic;
  MyMatrix<T> hecke_matrix;
  // Nonempty when the Hecke matrix could not be computed (modulo p, when p
  // divides a denominator of the images); the homology is still given.
  std::string hecke_failure;
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
  The Hecke operator on the homology of a level modulo a prime p, by
  sparse linear algebra: no matrix of the size of the complex is formed,
  which is what limits the rational computation. With B the space of the
  boundaries and Z the one of the cycles, H = Z / B is computed in the
  coordinates of the cells outside a set S whose cycle condition rows are
  a basis of the cycle conditions (see below); a basis of H is given by
  the basis cells of the quotient by the boundaries there
  (SparseQuotientMod), lifted to cycles, and the image of a basis vector
  is restricted to the complement of S, its class in the quotient giving
  its coordinates.
  For a chain map commuting with the boundary (Sharbly, Vertex) the cycles
  are the v with v D_index = 0 and the boundaries the rows of D_{index-1};
  for the transposed one (Dual), the v with D_{index-1} v = 0 and the
  columns of D_index. The entries of the Hecke matrix are residues mod p.
 */
template <typename T, typename Tint>
HeckeHomologyLevel<T> get_hecke_level_mod(int const &index, GammaComplex<T, Tint> const &gc,
                                 HeckeChainMap<T, Tint> const &hcm,
                                 HeckeCosets<T> const &cosets_gamma, int64_t const &p,
                                 bool const &with_hecke, [[maybe_unused]] std::ostream &os) {
  int n_levels = gc.basis.size();
  int n_cell = gc.basis[index].size();
  bool transposed = hcm.transposed;
  // The cycle condition sum_c v_c l_cycle_row[c] = 0 and the generators of
  // the boundaries, as sparse vectors modulo p.
  std::vector<SparseVectorMod> l_cycle_row(n_cell);
  std::vector<SparseVectorMod> l_bnd;
  int n_next = 0;
  auto f_push = [&](SparseVectorMod &v, int const &idx, int const &val) -> void {
    int64_t val_p = ResidueModPrime(val, p);
    if (val_p != 0) {
      v.push_back({idx, val_p});
    }
  };
  if (!transposed) {
    if (index < n_levels - 1) {
      n_next = gc.basis[index + 1].size();
      for (int c = 0; c < n_cell; c++) {
        for (auto &e : gc.boundary_sparse[index][c]) {
          f_push(l_cycle_row[c], e.first, e.second);
        }
      }
    }
    if (index > 0) {
      for (auto &row : gc.boundary_sparse[index - 1]) {
        SparseVectorMod v;
        for (auto &e : row) {
          f_push(v, e.first, e.second);
        }
        l_bnd.emplace_back(std::move(v));
      }
    }
  } else {
    if (index > 0) {
      n_next = gc.basis[index - 1].size();
      int n_prev = gc.boundary_sparse[index - 1].size();
      for (int r = 0; r < n_prev; r++) {
        for (auto &e : gc.boundary_sparse[index - 1][r]) {
          f_push(l_cycle_row[e.first], r, e.second);
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
        l_bnd.emplace_back(std::move(v));
      }
    }
  }
  // The column order of an echelon form: by increasing number of vectors
  // containing the column, so that the pivots (the first entries) are
  // rarely used columns and the fill-in stays small. Returns the label of
  // each column.
  auto f_label = [](std::vector<SparseVectorMod> const &l_vect, int const &n_col) -> std::vector<int> {
    std::vector<int> l_count(n_col, 0);
    for (auto &v : l_vect) {
      for (auto &e : v) {
        l_count[e.first] += 1;
      }
    }
    std::vector<int> l_order(n_col);
    for (int c = 0; c < n_col; c++) {
      l_order[c] = c;
    }
    std::stable_sort(l_order.begin(), l_order.end(),
                     [&](int const &a, int const &b) -> bool { return l_count[a] < l_count[b]; });
    std::vector<int> l_label(n_col);
    for (int c = 0; c < n_col; c++) {
      l_label[l_order[c]] = c;
    }
    return l_label;
  };
  auto f_relabel = [](SparseVectorMod const &v, std::vector<int> const &l_label) -> SparseVectorMod {
    SparseVectorMod w;
    w.reserve(v.size());
    for (auto &e : v) {
      w.push_back({l_label[e.first], e.second});
    }
    std::sort(w.begin(), w.end());
    return w;
  };
  /*
    The approach of Ash, Gunnells and McConnell (their matrix eta): first
    a set S of cells whose cycle condition rows are a basis of the row
    space of the cycle condition. A cycle v is then determined by its
    restriction to the complement S^c, and every vector on S^c extends to
    a cycle, so H = Z / B is the quotient of the vectors on S^c by the
    restrictions of the boundaries. That quotient has a small dimension,
    so the elimination of the boundaries has a low corank and a small
    fill-in, unlike the elimination of the boundaries in all the cells.
   */
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  std::vector<int> l_label_next = f_label(l_cycle_row, n_next);
  SparseEchelonMod echelon_cyc(p);
  std::vector<uint8_t> in_S(n_cell, 0);
  for (int c = 0; c < n_cell; c++) {
    if (!l_cycle_row[c].empty() && echelon_cyc.insert(f_relabel(l_cycle_row[c], l_label_next))) {
      in_S[c] = 1;
    }
  }
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: mod p index=" << index << " echelon of the cycle conditions, n_cell=" << n_cell
     << " n_next=" << n_next << " |S|=" << echelon_cyc.rank() << "|=" << time << "\n";
#endif
  // The same rows with markers at n_next + c: the reduction of the cycle
  // condition row of a cell gives its expression in the rows of S.
  SparseEchelonMod echelon_cyc_mark(p);
  for (int c = 0; c < n_cell; c++) {
    if (with_hecke && in_S[c] == 1) {
      SparseVectorMod v = f_relabel(l_cycle_row[c], l_label_next);
      v.push_back({n_next + c, 1});
      echelon_cyc_mark.insert(v);
    }
  }
#ifdef TIMINGS_HECKE_OPERATORS
  {
    size_t n_entry = 0;
    for (size_t i = 0; i < echelon_cyc_mark.rank(); i++) {
      n_entry += echelon_cyc_mark.get_row(i).size();
    }
    os << "|HECKE: mod p index=" << index << " marked echelon, n_entry=" << n_entry << "|=" << time << "\n";
  }
#endif
  // The restrictions of the boundaries to S^c, with the cells of S^c
  // relabeled by the count.
  std::vector<SparseVectorMod> l_bnd_res;
  for (auto &v : l_bnd) {
    SparseVectorMod w;
    for (auto &e : v) {
      if (in_S[e.first] == 0) {
        w.push_back(e);
      }
    }
    if (!w.empty()) {
      l_bnd_res.emplace_back(std::move(w));
    }
  }
  // The quotient of the vectors on S^c (the cell c at the position
  // pos_Sc[c]) by the restrictions of the boundaries. The pivots are chosen
  // globally and the last columns are treated densely: the incremental
  // reduced echelon form fills in on the large levels.
  std::vector<int> pos_Sc(n_cell, -1);
  std::vector<int> l_Sc;
  for (int c = 0; c < n_cell; c++) {
    if (in_S[c] == 0) {
      pos_Sc[c] = l_Sc.size();
      l_Sc.push_back(c);
    }
  }
  int n_Sc = l_Sc.size();
  std::vector<SparseVectorMod> l_bnd_Sc;
  for (auto &v : l_bnd_res) {
    SparseVectorMod w;
    for (auto &e : v) {
      w.push_back({pos_Sc[e.first], e.second});
    }
    l_bnd_Sc.emplace_back(std::move(w));
  }
  SparseQuotientMod quot_bnd(l_bnd_Sc, n_Sc, p);
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: mod p index=" << index << " quotient by the boundaries, |l_bnd_res|=" << l_bnd_res.size()
     << " n_Sc=" << n_Sc << " rank=" << quot_bnd.rank() << " n_sparse_pivot=" << quot_bnd.n_sparse_pivot()
     << " n_dense_col=" << quot_bnd.n_dense_col() << " n_entry=" << quot_bnd.n_entry() << "|=" << time << "\n";
#endif
  // The basis columns of the quotient give the cells of the basis of H.
  std::vector<int> l_free;
  for (auto &u : quot_bnd.basis_columns()) {
    l_free.push_back(l_Sc[u]);
  }
  int dim_hom = l_free.size();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: mod p index=" << index << " n_cell=" << n_cell << " |S|=" << echelon_cyc.rank()
     << " rank(B)=" << quot_bnd.rank() << " dim_hom=" << dim_hom << "\n";
#endif
  if (!with_hecke) {
    return {index, n_cell, MyMatrix<T>(dim_hom, 0), MyMatrix<T>(0, 0), "no chain map at this level"};
  }
  MyMatrix<T> hecke_matrix = ZeroMatrix<T>(dim_hom, dim_hom);
  HeckeImageKeyActions actions = get_hecke_image_key_actions(index, gc, hcm);
  for (int i = 0; i < dim_hom; i++) {
    int c_free = l_free[i];
    // The cycle lifting the free cell: e_c plus its part on S.
    SparseVectorMod v{{c_free, 1}};
    if (!l_cycle_row[c_free].empty()) {
      SparseVectorMod red = echelon_cyc_mark.reduce(f_relabel(l_cycle_row[c_free], l_label_next));
      for (auto &e : red) {
        if (e.first < n_next) {
          throw HeckeFailure{"Hecke matrix mod p: a cycle condition row is not in the span of S"};
        }
        v.push_back({e.first - n_next, e.second});
      }
      std::sort(v.begin(), v.end());
    }
    // Its image, restricted to S^c, and its class modulo the boundaries.
    std::map<int, int64_t> acc;
    for (auto &e : v) {
      for (auto &kv : get_hecke_image_cell_sparse(index, e.first, gc, hcm, cosets_gamma, &actions)) {
        if (in_S[kv.first] == 1) {
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
          return {index, n_cell, MyMatrix<T>(dim_hom, 0), MyMatrix<T>(0, 0),
                  "a denominator of an image is divisible by p=" + std::to_string(p)};
        }
        int64_t &val = acc[pos_Sc[kv.first]];
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
  return {index, n_cell, MyMatrix<T>(dim_hom, 0), std::move(hecke_matrix), ""};
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
    bool const &all_levels, int64_t const &prime, std::ostream &os) {
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
  GammaComplex<T, Tint> gc = compute_gamma_complex(fce, ct, only_well_rounded, l_level_needed, dense_boundary, os);
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
      result.l_level.emplace_back(get_hecke_level_mod(index, gc, hcm, cosets_gamma, prime, with_hecke, os));
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
    os_out << ")";
  }
  os_out << "])";
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_HOMOLOGY_H_
// clang-format on
