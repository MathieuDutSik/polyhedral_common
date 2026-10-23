// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_CHAIN_MAP_H_
#define SRC_PERFECT_HECKE_CHAIN_MAP_H_

// clang-format off
#include "hecke_chains.h"
// clang-format on

/*
  Step A of the Hecke operators: the chain map Phi(F, y_j) on the complex
  (HeckeChainMap), the averaging over the coset stabilizers and the
  construction from the vertices (the "Vertex" method).
 */

//
// Step A: the Hecke chain map on the complex
//

/*
  The action of the Hecke double coset G x G on the cell complex (step A).
  It is built once, by any of the methods (compute_hecke_chain_map_dual,
  compute_hecke_chain_map, compute_hecke_chain_map_sharbly, ...), and is
  then consumed by compute_hecke_homology for any finite index subgroup.
  Nothing in it depends on the subgroup.
 */
template <typename T, typename Tint> struct HeckeChainMap {
  MyMatrix<T> x;
  HeckeCosets<T> cosets;
  // images[index][iOrb][j] is the chain of level index representing the
  // image of the cell F_{iOrb} under the coset y_j.
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images;
  // The levels at which the images were computed. Empty means all of
  // them (the chain map methods); the sharbly method computes only the
  // symbol level and the level below it, and the homology is reported
  // only there.
  std::vector<uint8_t> l_level_computed;
  // The method that built it ("Dual", "Vertex", "Sharbly").
  std::string method;
  // Whether the images commute with the transposed differential (the
  // coboundary from a cell to the cells containing it, as for "Dual")
  // rather than with the boundary. The cycles and boundaries of step B
  // are taken accordingly.
  bool transposed;
  bool is_level_computed(int const &index) const {
    if (l_level_computed.empty()) {
      return true;
    }
    return l_level_computed[index] == 1;
  }
};

// The f_member of the cosets is not stored: it is the membership in the
// ambient group and is set again on reading.
namespace boost::serialization {
template <class Archive, typename T, typename Tint>
inline void serialize(Archive &ar, HeckeChainMap<T, Tint> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", val.x);
  ar &make_nvp("cosets_x", val.cosets.x);
  ar &make_nvp("l_coset", val.cosets.l_coset);
  ar &make_nvp("l_coset_inv", val.cosets.l_coset_inv);
  ar &make_nvp("images", val.images);
  ar &make_nvp("l_level_computed", val.l_level_computed);
  ar &make_nvp("method", val.method);
  ar &make_nvp("transposed", val.transposed);
}
}  // namespace boost::serialization

/*
  The action of the generators of the stabilizer of a cell on the Hecke
  cosets: act[i_gen][j] = (j', s') with s y_j = y_{j'} s'.
 */
template <typename T, typename Tint>
std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>>
get_stab_coset_action(std::vector<MyMatrix<Tint>> const &l_gens,
                      HeckeCosets<T> const &hc) {
  size_t n_gen = l_gens.size();
  size_t n_coset = hc.size();
  std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act(n_gen);
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    MyMatrix<T> s_T = UniversalMatrixConversion<T, Tint>(l_gens[i_gen]);
    std::vector<std::pair<int, MyMatrix<Tint>>> l_act;
    for (size_t j = 0; j < n_coset; j++) {
      MyMatrix<T> y = s_T * hc.l_coset[j];
      std::pair<int, MyMatrix<T>> pair = hc.decompose(y);
      MyMatrix<Tint> s_conj = UniversalMatrixConversion<Tint, T>(pair.second);
      l_act.emplace_back(pair.first, std::move(s_conj));
    }
    act[i_gen] = std::move(l_act);
  }
  return act;
}

/*
  The orbit of a coset j under the stabilizer with the transversal:
  for the coset j' in the orbit, transversal[pos] = z is an element of the
  stabilizer with z y_j = y_{j'} z' and transversal_conj[pos] = z'.
 */
template <typename Tint> struct CosetOrbit {
  std::vector<int> l_orbit;
  std::vector<int> position;
  std::vector<MyMatrix<Tint>> transversal;
  std::vector<MyMatrix<Tint>> transversal_conj;
};

template <typename Tint>
CosetOrbit<Tint> get_coset_orbit(
    int const &j, std::vector<MyMatrix<Tint>> const &l_gens,
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> const &act,
    size_t const &n_coset, int const &n) {
  CosetOrbit<Tint> co;
  co.position = std::vector<int>(n_coset, -1);
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  co.l_orbit.push_back(j);
  co.position[j] = 0;
  co.transversal.push_back(id);
  co.transversal_conj.push_back(id);
  size_t pos = 0;
  size_t n_gen = l_gens.size();
  while (pos < co.l_orbit.size()) {
    int j1 = co.l_orbit[pos];
    for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
      std::pair<int, MyMatrix<Tint>> const &pair = act[i_gen][j1];
      int j2 = pair.first;
      if (co.position[j2] == -1) {
        co.position[j2] = co.l_orbit.size();
        co.l_orbit.push_back(j2);
        MyMatrix<Tint> z = l_gens[i_gen] * co.transversal[pos];
        MyMatrix<Tint> z_conj = pair.second * co.transversal_conj[pos];
        co.transversal.emplace_back(std::move(z));
        co.transversal_conj.emplace_back(std::move(z_conj));
      }
    }
    pos += 1;
  }
  return co;
}

/*
  The Schreier generators of the stabilizer S_j of the coset j in the
  stabilizer of the cell. Each generator h comes with h' = y_j^{-1} h y_j.
 */
template <typename Tint>
std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>>
get_coset_stabilizer_schreier(
    CosetOrbit<Tint> const &co, std::vector<MyMatrix<Tint>> const &l_gens,
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> const &act,
    int const &n) {
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  std::unordered_set<MyMatrix<Tint>> set_h;
  std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier;
  size_t n_gen = l_gens.size();
  size_t n_orbit = co.l_orbit.size();
  for (size_t pos1 = 0; pos1 < n_orbit; pos1++) {
    int j1 = co.l_orbit[pos1];
    for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
      std::pair<int, MyMatrix<Tint>> const &pair = act[i_gen][j1];
      int j2 = pair.first;
      int pos2 = co.position[j2];
      MyMatrix<Tint> h = Inverse(co.transversal[pos2]) * l_gens[i_gen] * co.transversal[pos1];
      if (h != id) {
        if (set_h.insert(h).second) {
          MyMatrix<Tint> h_conj = Inverse(co.transversal_conj[pos2]) * pair.second * co.transversal_conj[pos1];
          l_schreier.emplace_back(std::move(h), std::move(h_conj));
        }
      }
    }
  }
  return l_schreier;
}

/*
  Averaging of a chain of level index, representing the image of the
  cell F_{iOrb} under a coset y_j, over the stabilizer S_j of the coset
  in Stab_G(F). The averaged chain satisfies
    chain h' = sign_F(h) chain     for h in S_j, h' = y_j^{-1} h y_j.
  S_j may be infinite for a cell at infinity. But its elements act on the
  chain only through their restriction to the span P of the vectors of F
  (the chain lives in P y_j) and the restriction is determined by the
  signed permutation of the vectors of F, which forms a finite group. So
  the closure is computed with EXT * h as key and the average is taken
  over that finite group.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> average_chain_coset_stabilizer(
    int const &index, int const &iOrb,
    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
    std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> const &l_schreier,
    FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  if (l_schreier.empty()) {
    return chain;
  }
  FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
  MyMatrix<Tint> const &EXT = face.EXT;
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  std::unordered_map<MyMatrix<Tint>, size_t> map_key;
  std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_elt;
  std::vector<int> l_sign;
  auto f_insert = [&](MyMatrix<Tint> const &h, MyMatrix<Tint> const &h_conj) -> void {
    MyMatrix<Tint> key = EXT * h;
    if (map_key.contains(key)) {
#ifdef SANITY_CHECK_HECKE_OPERATORS
      // The average is taken over the restrictions to the span of the
      // vectors of F, keyed by EXT * h. Two elements with the same key
      // must act identically on the chain: the chain must live in the
      // subspace on which the restriction is faithful. If this fails the
      // solution was not forced into the right subspace before the
      // averaging and the averaged chain is representative-dependent.
      size_t idx = map_key.at(key);
      if (l_elt[idx].first == h) {
        // The same element: nothing to check. Different elements with the
        // same key only occur when the restriction is not faithful, i.e.
        // for the cells at infinity.
        return;
      }
      int sign0 = l_sign[idx];
      int sign1 = get_face_orientation(EXT, ListMat, face.or_info, h);
      MyMatrix<Tint> prod = Inverse(h_conj) * l_elt[idx].second;
      std::vector<PerfectFaceEntry<T, Tint>> chain_img =
          chain_right_multiply(chain, T(sign0 * sign1), prod);
      if (!is_equal_chain(chain_img, chain, index, fce, os)) {
        std::cerr << "HECKE: two stabilizer elements with the same restriction to the\n";
        std::cerr << "HECKE: span of the cell act differently on the chain. The chain is\n";
        std::cerr << "HECKE: not supported on the subspace where the restriction is\n";
        std::cerr << "HECKE: faithful, so the average over the finite restriction group\n";
        std::cerr << "HECKE: is not an average over the stabilizer. The solution must be\n";
        std::cerr << "HECKE: forced into the right subspace before the averaging.\n";
        std::cerr << "HECKE: index=" << index << " iOrb=" << iOrb << "\n";
        throw TerminalException{1};
      }
#endif
      return;
    }
    map_key[key] = l_elt.size();
    int sign = get_face_orientation(EXT, ListMat, face.or_info, h);
    l_elt.emplace_back(h, h_conj);
    l_sign.push_back(sign);
  };
  f_insert(id, id);
  size_t pos = 0;
  HeckeBudget budget;
  while (pos < l_elt.size()) {
    if (l_elt.size() > budget.max_closure) {
      throw HeckeFailure{"average_chain_coset_stabilizer: closure larger than " + std::to_string(budget.max_closure)};
    }
    for (auto &pair : l_schreier) {
      MyMatrix<Tint> h = l_elt[pos].first * pair.first;
      MyMatrix<Tint> h_conj = l_elt[pos].second * pair.second;
      f_insert(h, h_conj);
    }
    pos += 1;
  }
  size_t n_elt = l_elt.size();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: average_chain_coset_stabilizer index=" << index << " iOrb=" << iOrb
     << " |l_schreier|=" << l_schreier.size() << " |group|=" << n_elt << "\n";
#endif
  if (n_elt == 1) {
    return chain;
  }
  int n_elt_i = n_elt;
  T inv_size = T(1) / T(n_elt_i);
  ChainBuilder<T, Tint, Tgroup> cb(index, fce, os);
  for (size_t i_elt = 0; i_elt < n_elt; i_elt++) {
    MyMatrix<Tint> Minv = Inverse(l_elt[i_elt].second);
    T scal = inv_size * T(l_sign[i_elt]);
    for (auto &fe : chain) {
      T value = fe.value * scal;
      MyMatrix<Tint> M = fe.M * Minv;
      cb.f_insert(value, fe.iOrb, M);
    }
  }
  return cb.get_faces();
}

/*
  Whether the chain of level index_chain satisfies
  chain h' = sign_F(h) chain for the Schreier generators, with F the
  cell (index_face, iOrb) whose coset stabilizer is considered.
 */
template <typename T, typename Tint, typename Tgroup>
bool is_chain_coset_stabilizer_invariant(
    int const &index_face, int const &iOrb, int const &index_chain,
    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
    std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> const &l_schreier,
    FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index_face].l_faces[iOrb];
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  for (auto &pair : l_schreier) {
    int sign = get_face_orientation(face.EXT, ListMat, face.or_info, pair.first);
    std::vector<PerfectFaceEntry<T, Tint>> chain_img =
        chain_right_multiply(chain, T(sign), pair.second);
    if (!is_equal_chain(chain_img, chain, index_chain, fce, os)) {
      return false;
    }
  }
  return true;
}

/*
  The image of a vertex under a coset y: the vertex of the primitive
  vectors of EXT * y.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
get_vertex_image(int const &vlevel, int const &iOrb, MyMatrix<T> const &y,
                 FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                 ConeWalker<T, Tint, Tgroup> const &walker,
                 std::unordered_map<MyMatrix<Tint>, std::vector<PerfectFaceEntry<T, Tint>>> &cache,
                 [[maybe_unused]] std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> const &EXT = fce.levels[vlevel].l_faces[iOrb].EXT;
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  MyMatrix<T> prod = EXT_T * y;
  int n_row = prod.rows();
  MyMatrix<Tint> EXT_img(n_row, n);
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(prod, i_row);
    MyVector<T> Vred = RemoveFractionVector(V);
    MyVector<Tint> Vint = UniversalVectorConversion<Tint, T>(Vred);
    AssignMatrixRow(EXT_img, i_row, Vint);
  }
  MyMatrix<Tint> key = tot_set(EXT_img);
  auto iter = cache.find(key);
  if (iter != cache.end()) {
    return iter->second;
  }
  FceFaceSearch<Tint> ffs = walker.search_vertex(prod);
  PerfectFaceEntry<T, Tint> fe{ffs.iOrb, T(1), ffs.M};
  std::vector<PerfectFaceEntry<T, Tint>> chain{fe};
  cache[key] = chain;
  return chain;
}

template <typename T, typename Tint, typename Tgroup>
void check_fce_for_hecke(FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  PerfectComplexOptions const &pco = fce.pctdi.pco;
  if (pco.only_well_rounded) {
    std::cerr << "HECKE: The Hecke operators require the full complex with the vertices\n";
    std::cerr << "HECKE: Please set OnlyWellRounded = F\n";
    throw TerminalException{1};
  }
  if (!pco.compute_boundary) {
    std::cerr << "HECKE: The Hecke operators require the boundary, please set ComputeBoundary = T\n";
    throw TerminalException{1};
  }
  if (!pco.compute_partial_contracting_homotopy) {
    std::cerr << "HECKE: The Hecke operators require ComputePartialContractingHomotopy = T\n";
    throw TerminalException{1};
  }
  int n_levels = fce.levels.size();
  for (auto &face : fce.levels[n_levels - 1].l_faces) {
    if (face.or_info.ListRowSelect.size() != 1) {
      std::cerr << "HECKE: The last level should be made of vertices\n";
      throw TerminalException{1};
    }
  }
}

/*
  The right hand side Phi(dF, y_j) for the cell F = F_{iOrb} of level
  index, computed from the images of the level index+1.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> get_hecke_boundary_image(
    int const &index, int const &iOrb, int const &j,
    FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeCosets<T> const &cosets,
    std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>> const &images_next,
    std::ostream &os) {
  ChainBuilder<T, Tint, Tgroup> cb(index + 1, fce, os);
  MyMatrix<T> const &y = cosets.l_coset[j];
  for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(ebnd.M);
    MyMatrix<T> z = M_T * y;
    std::pair<int, MyMatrix<T>> pair = cosets.decompose(z);
    MyMatrix<Tint> M_conj = UniversalMatrixConversion<Tint, T>(pair.second);
    T e_sign = T(ebnd.sign);
    for (auto &fe : images_next[ebnd.iOrb][pair.first]) {
      T value = fe.value * e_sign;
      MyMatrix<Tint> M = fe.M * M_conj;
      cb.f_insert(value, fe.iOrb, M);
    }
  }
  return cb.get_faces();
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                        MyMatrix<T> const &x, std::ostream &os) {
  check_fce_for_hecke(fce);
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
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
  int vlevel = n_levels - 1;
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images(n_levels);
  //
  // The vertices
  //
  {
    int n_orb = fce.levels[vlevel].l_faces.size();
    std::unordered_map<MyMatrix<Tint>, std::vector<PerfectFaceEntry<T, Tint>>> cache;
    ConeWalker<T, Tint, Tgroup> walker(fce, os);
    images[vlevel] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      for (size_t j = 0; j < n_coset; j++) {
        images[vlevel][iOrb][j] =
            get_vertex_image(vlevel, iOrb, cosets.l_coset[j], fce, walker, cache, os);
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: vertex level " << vlevel << " n_orb=" << n_orb << " done\n";
#endif
  }
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: vertex images|=" << time << "\n";
#endif
  //
  // The other levels, from the bottom to the top
  //
  for (int index = vlevel - 1; index >= 0; index--) {
    int n_orb = fce.levels[index].l_faces.size();
    images[index] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
      std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
          get_stab_coset_action(face.l_gens, cosets);
      std::vector<int> status(n_coset, 0);
      for (size_t j = 0; j < n_coset; j++) {
        if (status[j] == 1) {
          continue;
        }
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: index=" << index << " iOrb=" << iOrb << "/" << n_orb
           << " j=" << j << "/" << n_coset << "\n";
#endif
        std::vector<PerfectFaceEntry<T, Tint>> rhs =
            get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
        std::vector<PerfectFaceEntry<T, Tint>> sol;
        if (!rhs.empty()) {
          std::pair<std::function<bool(MyMatrix<Tint> const &)>, bool> pair_filter =
              get_subspace_filter(index + 1, rhs, fce, os);
          sol = partial_contracting_homotopy(index + 1, rhs, fce, pair_filter.first, os);
        }
        CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
        std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
            get_coset_stabilizer_schreier(co, face.l_gens, act, n);
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (!is_chain_coset_stabilizer_invariant(index, iOrb, index + 1, rhs, l_schreier, fce, os)) {
          std::cerr << "HECKE: The right hand side should be invariant under the coset stabilizer\n";
          throw TerminalException{1};
        }
#endif
        std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
            average_chain_coset_stabilizer(index, iOrb, sol, l_schreier, fce, os);
#ifdef SANITY_CHECK_HECKE_OPERATORS
        {
          std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary(index, sol_avg, fce, os);
          if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
            std::cerr << "HECKE: The averaged chain is not a solution\n";
            throw TerminalException{1};
          }
          if (!is_chain_coset_stabilizer_invariant(index, iOrb, index, sol_avg, l_schreier, fce, os)) {
            std::cerr << "HECKE: The averaged chain should be invariant under the coset stabilizer\n";
            throw TerminalException{1};
          }
        }
#endif
        // Propagation to the orbit of the coset under the stabilizer:
        // Phi(F, y_{j'}) = sign(z) Phi(F, y_j) z'^{-1} for z y_j = y_{j'} z'.
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
      for (size_t j = 0; j < n_coset; j++) {
        std::vector<PerfectFaceEntry<T, Tint>> rhs =
            get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
        std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary(index, images[index][iOrb][j], fce, os);
        if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
          std::cerr << "HECKE: The chain map property fails at index=" << index
                    << " iOrb=" << iOrb << " j=" << j << "\n";
          throw TerminalException{1};
        }
      }
#endif
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: level " << index << " n_orb=" << n_orb << " done\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: level " << index << "|=" << time << "\n";
#endif
  }
  return {x, std::move(cosets), std::move(images), {}, "Vertex", false};
}

template <typename T, typename Tint>
void WriteHeckeChainMapGAP(std::ostream &os_out, HeckeChainMap<T, Tint> const &hcm) {
  os_out << "rec(x:=";
  WriteMatrixGAP(os_out, hcm.x);
  os_out << ",\nListCoset:=";
  WriteListMatrixGAP(os_out, hcm.cosets.l_coset);
  os_out << ",\nListImages:=[";
  size_t n_levels = hcm.images.size();
  for (size_t index = 0; index < n_levels; index++) {
    if (index > 0) {
      os_out << ",\n";
    }
    os_out << "[";
    size_t n_orb = hcm.images[index].size();
    for (size_t iOrb = 0; iOrb < n_orb; iOrb++) {
      if (iOrb > 0) {
        os_out << ",\n";
      }
      os_out << "[";
      size_t n_coset = hcm.images[index][iOrb].size();
      for (size_t j = 0; j < n_coset; j++) {
        if (j > 0) {
          os_out << ",\n";
        }
        WriteListPerfectFaceEntryGAP(os_out, hcm.images[index][iOrb][j]);
      }
      os_out << "]";
    }
    os_out << "]";
  }
  os_out << "])";
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_CHAIN_MAP_H_
// clang-format on
