// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_OPERATORS_H_
#define SRC_PERFECT_HECKE_OPERATORS_H_

// clang-format off
#include "hecke_dual.h"
#include "hecke_sharbly.h"
#include "hecke_homology.h"
// clang-format on

/*
  Hecke operators on the homology of the perfect form complex.

  Notation. G is the arithmetic group of the T-space (GL_n(Z) for the
  classic T-space), acting on the right on vector configurations: a cell
  is given by a matrix EXT of vectors (each vector together with its
  opposite) and g in G maps it to EXT * g. The perfect form complex is
  the full Voronoi complex Pi: all faces of the perfect cones, including
  the faces "at infinity" whose vectors do not span R^n. Its cells are
  organized by levels: level 0 holds the perfect cones, level k the faces
  of codimension k, and the last level holds the vertices (single
  vectors). The differential goes from a level to the next one (from a
  cone to its facets). Restricting to the well rounded cells (those whose
  vectors span R^n) gives the quotient complex Pi / Pi_infinity, whose
  homology modulo a finite index subgroup Gamma of G is (rationally) the
  cohomology H^k(Gamma, Q) with k the level, by the standard duality.

  A rational matrix x in the commensurator (it must preserve the T-space)
  defines the Hecke operator of the double coset Gamma x Gamma. The
  computation is organized in two steps, following ProjectiveSystem.g of
  the MyPolyhedral GAP package.

  Step A (the action on the cell complex, independent of Gamma).
  Write G x G = union of the cosets y_j G. For every orbit representative
  F of a cell under G and every coset y_j we compute a chain Phi(F, y_j)
  of the complex which represents "F x y_j". The chains satisfy
    (1) d Phi(F, y_j) = Phi(dF, y_j),
    (2) Phi(F g, y_j) = Phi(F, g y_j) = Phi(F, y_{j'}) g'  for g y_j = y_{j'} g',
  where (2) is used to define Phi on non-representative cells and, for g in
  the stabilizer of F, is a consistency requirement. Vertices are mapped
  to vertices (the primitive vector of v x). For higher cells the right
  hand side of (1) is computed from the lower level and Phi(F, y_j) is
  obtained by the partial contracting homotopy. The solution is searched among
  the cells whose vectors lie in the span of the vectors of the right
  hand side (this is what makes the construction respect the cells at
  infinity, and hence pass to the quotient Pi / Pi_infinity). Then the
  chain is averaged over the stabilizer of the coset in Stab_G(F) so that
  (2) holds for the stabilizer, and (2) propagates the chain to the other
  cosets of the orbit of Stab_G(F) on the cosets.

  Step B (the action on the homology of the quotient by Gamma).
  Gamma is a finite index subgroup, G = union of the cosets g_i Gamma.
  The Gamma-orbits of cells in the G-orbit of F correspond to the orbits
  of Stab_G(F) on the cosets g_i Gamma. This gives the cells of the
  quotient complex by Gamma with their orientability, hence the boundary
  matrices and the homology. The Hecke double coset Gamma x Gamma is
  split into cosets x_k Gamma and the operator on the cell F g_i is
    T(F g_i) = sum_k Phi(F, y_j) g'   with g_i x_k = y_j g',
  the cells of the result being then identified modulo Gamma. The Hecke
  matrices are expressed in a basis of harmonic representatives of the
  homology.
 */

/*
  The files:
  - hecke_common.h: the failures, the budgets, the ambient group G and the
    cosets of a Hecke double coset,
  - hecke_subgroups.h: the finite index subgroups Gamma and their cosets,
  - hecke_chains.h: the chains and the walk through the perfect cones,
  - hecke_chain_map.h: step A, the chain map and the "Vertex" method,
  - hecke_dual.h: the "Dual" method of step A,
  - hecke_symbols.h, hecke_sharbly.h: the "Sharbly" method of step A,
  - hecke_homology.h: step B,
  - this file: the choice of the method of step A and its cache.
 */

//
// Step A, all methods: building the chain map once and caching it
//

/*
  The chain map of x by the given method. If cache_file is not "none" and
  holds the chain map of the same x on the same complex, it is read
  instead of computed; otherwise it is computed and written there. The
  stored images refer to the cells by their orbit index, so the complex
  is identified by the vectors of the representatives of the orbits.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<std::vector<MyMatrix<Tint>>>
get_complex_fingerprint(FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  std::vector<std::vector<MyMatrix<Tint>>> l_ext;
  for (auto &level : fce.levels) {
    std::vector<MyMatrix<Tint>> l_ext_level;
    for (auto &face : level.l_faces) {
      l_ext_level.push_back(face.EXT);
    }
    l_ext.emplace_back(std::move(l_ext_level));
  }
  return l_ext;
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map_method(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                               MyMatrix<T> const &x, std::string const &method,
                               std::ostream &os) {
  if (method == "Dual") {
    return compute_hecke_chain_map_dual(fce, x, os);
  }
  if (method == "Vertex") {
    return compute_hecke_chain_map(fce, x, os);
  }
  if (method == "Sharbly") {
    return compute_hecke_chain_map_sharbly(fce, x, os);
  }
  std::cerr << "HECKE: HeckeMethod should be Dual, Vertex or Sharbly, method=" << method << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
get_hecke_chain_map(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                    MyMatrix<T> const &x, std::string const &method,
                    std::string const &cache_file, std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  std::vector<std::vector<MyMatrix<Tint>>> l_ext = get_complex_fingerprint(fce);
  if (cache_file != "none" && FILE_IsExistingFile(cache_file)) {
    std::vector<std::vector<MyMatrix<Tint>>> l_ext_read;
    HeckeChainMap<T, Tint> hcm;
    {
      std::ifstream ifs(cache_file);
      boost::archive::text_iarchive ia(ifs);
      ia >> l_ext_read;
      ia >> hcm;
    }
    if (l_ext_read == l_ext && hcm.x == x) {
      hcm.cosets.f_member = [&LinSpa](MyMatrix<T> const &g) -> bool {
        return is_in_ambient_group(LinSpa, g);
      };
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: chain map read from " << cache_file << " method=" << hcm.method << "\n";
#endif
      return hcm;
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: the chain map of " << cache_file << " is for another x or complex, recomputing\n";
#endif
  }
  HeckeChainMap<T, Tint> hcm = compute_hecke_chain_map_method(fce, x, method, os);
  if (cache_file != "none") {
    std::ofstream ofs(cache_file);
    boost::archive::text_oarchive oa(ofs);
    oa << l_ext;
    oa << hcm;
  }
  return hcm;
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_OPERATORS_H_
// clang-format on
