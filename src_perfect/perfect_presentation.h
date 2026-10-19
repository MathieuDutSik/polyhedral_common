// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_PRESENTATION_H_
#define SRC_PERFECT_PERFECT_PRESENTATION_H_

// clang-format off
#include "perfect_complex.h"
#include "Presentation.h"
#include "Resolution.h"
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERFECT_PRESENTATION
#endif

#ifdef DEBUG
#define DEBUG_PERFECT_PRESENTATION
#endif

/*
  Presentation by generators and relators of the group Gamma acting on the
  perfect form complex.

  The group acts on the well rounded retract W, which is contractible
  (A. Ash, Small-dimensional classifying spaces for arithmetic subgroups of
  general linear groups, Duke Math. J. 51 (1984)). Its cells are the faces
  of the perfect domains whose relative interior meets the positive definite
  cone:
  --- The vertices are the perfect domains.
  --- The edges are the facets of the perfect domains, each one shared by
      two perfect domains.
  --- The 2-cells are the well rounded faces of codimension 2. The perfect
      domains containing such a face form a cycle, each one meeting the next
      along a facet containing the face.
  K. S. Brown (Presentations for groups acting on simply-connected
  complexes, J. Pure Appl. Algebra 32 (1984)) gives a presentation from the
  stabilizers of the cells and the side transformations. It is the one used
  in R. Coulangeon, G. Nebe, O. Braun, S. Schoennenbeck, Computing in
  arithmetic groups with Voronoi's algorithm, J. Algebra 435 (2015).

  The group acts on the right: the domain v_i g has the vector configuration
  EXT_i g. The perfect domain v_j s_{i,k} is the neighbor of v_i across the
  facet f_k of l_sing_adj[k], and s_{i,k} is a side transformation. The
  neighbor across the facet f = f_k h, h in Gamma_i, is v_j s_{i,k} h.

  Generators:
  --- For each perfect domain v_i, the strong generators of its stabilizer
      Gamma_i, which is finite.
  --- For each perfect domain v_i and each orbit of facets f_k, the side
      transformation s_{i,k}.
  Relators:
  --- R0: The relators of each Gamma_i, from its stabilizer chain.
  --- R1: For g in the stabilizer of the facet f_k, which fixes v_i and its
      neighbor v_j s_{i,k}: s_{i,k} g s_{i,k}^{-1} = c with c in Gamma_j.
  --- R2: The reversal of the edge. The domain v_i s_{i,k}^{-1} is the
      neighbor of v_j across the facet f_k s_{i,k}^{-1} = f_{k'} h' of v_j.
      So it is v_i s_{j,k'} h' and s_{j,k'} h' s_{i,k} = c' with c' in
      Gamma_i.
  --- R3: For each orbit of 2-cells, the walk around the cycle of the
      perfect domains containing it. Starting from v_{i_0}, every step
      crosses a facet f = f_k h of the current domain v_i g containing the
      2-cell: the next domain is v_j s_{i,k} h g. The walk ends back at
      v_{i_0} with an element g of Gamma_{i_0}: the word of g is equal to
      the element of Gamma_{i_0}.
  --- R4: The theorem requires the representatives of the perfect domains
      to form a tree: s_{i,k} = 1 for the edges of a spanning tree of the
      quotient graph. The representatives are changed to
      v'_i = v_i t_i, with t_0 = 1 and t_j = s_{i,k} t_i for the edge
      (i, k) of the tree from v_i to v_j, so that v'_j is the neighbor of
      v'_i. The generators of the stabilizers become t_i^{-1} g t_i and the
      side transformations t_j^{-1} s_{i,k} t_i, which is 1 on the edges of
      the tree. The relators R0, R1, R2, R3 hold for the new generators
      since they are obtained by the same conjugations.
  Then the generating set is simplified by Tietze transformations.

  The elements of Gamma_i are written as words in its strong generators by
  sifting them through the stabilizer chain. The action of Gamma_i on the
  vector configuration EXT_i of the perfect domain is faithful when EXT_i is
  of full rank. Otherwise, the action is on the orbit of EXT_i and of the
  standard basis, which is of full rank.
 */

/*
  The words of the elements of a finite permutation group in the strong
  generators of its stabilizer chain, and the relators of the presentation
  on those generators. For the small groups, the relators are the ones of
  the free resolution, which are fewer and shorter than the Schreier
  relators of the stabilizer chain.
 */
static constexpr size_t PERFECT_PRESENTATION_MAX_ORDER_RESOLUTION = 25000;

template <typename Tgroup> struct FiniteGroupWordSolver {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using TintGroup = typename Tgroup::Tint;
  Tgroup GRP;
  permutalib::Presentation<Telt> pres;
  std::unordered_map<Telt, int> map_letter;
  FiniteGroupWordSolver(Tgroup const &_GRP)
      : GRP(_GRP), pres(GRP.PresentationOnStrongGenerators()) {
    int n_gen = pres.generators.size();
    for (int i_gen = 0; i_gen < n_gen; i_gen++) {
      map_letter.try_emplace(pres.generators[i_gen], i_gen + 1);
    }
    if (n_gen > 0 &&
        GRP.size() <= TintGroup(PERFECT_PRESENTATION_MAX_ORDER_RESOLUTION)) {
      Telt id(static_cast<Tidx>(GRP.n_act()));
      permutalib::Presentation<Telt> pres_res =
          permutalib::PresentationByResolution<Telt, TintGroup>(
              pres.generators, id);
      // The strong generators are distinct and not the identity, so they
      // are the generators of pres_res, in the same order.
      if (pres_res.generators != pres.generators) {
        std::cerr << "PERFPRES: The generators of the resolution are not the "
                     "strong generators\n";
        throw TerminalException{1};
      }
      pres.relators = pres_res.relators;
    }
  }
  // The sifting gives g l_1 ... l_m = 1 with l_i the labels of the
  // transversals, so g = (l_1 ... l_m)^{-1}.
  permutalib::PresentationWord get_word(Telt const &g) const {
    auto S = GRP.stab_chain();
    using Tidx_label = std::decay_t<decltype(S->transversal[0])>;
    Tidx_label miss = std::numeric_limits<Tidx_label>::max();
    permutalib::PresentationWord w;
    Telt gW = g;
    auto Sptr = S;
    while (Sptr != nullptr && Sptr->orbit.size() > 0) {
      auto bpt = Sptr->orbit[0];
      auto img = PowAct(bpt, gW);
      if (Sptr->transversal[img] == miss) {
        std::cerr << "PERFPRES: The element does not belong to the group\n";
        throw TerminalException{1};
      }
      while (img != bpt) {
        Tidx_label lbl = Sptr->transversal[img];
        Telt const &x = S->comm->labels[lbl];
        gW *= x;
        if (!x.is_identity()) {
          w.push_back(map_letter.at(x));
        }
        img = PowAct(bpt, gW);
      }
      Sptr = Sptr->stabilizer;
    }
    if (!gW.is_identity()) {
      std::cerr << "PERFPRES: The sifting does not end at the identity\n";
      throw TerminalException{1};
    }
    return permutalib::InverseWord(w);
  }
};

template <typename Tint>
MyMatrix<Tint> EvaluateWordMatrix(std::vector<MyMatrix<Tint>> const &l_gen,
                                  int n,
                                  permutalib::PresentationWord const &w) {
  MyMatrix<Tint> M = IdentityMat<Tint>(n);
  for (auto &x : w) {
    if (x > 0) {
      M = M * l_gen[x - 1];
    } else {
      M = M * Inverse(l_gen[-x - 1]);
    }
  }
  return M;
}

// The data of the stabilizer of a perfect domain.
template <typename Tint, typename Tgroup> struct PerfectDomainGroup {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  // The vectors on which the action is faithful, the first ones are EXT.
  MyMatrix<Tint> V;
  std::unordered_map<MyVector<Tint>, Tidx> map_V;
  FiniteGroupWordSolver<Tgroup> solver;
  // The matrices of the strong generators of the solver.
  std::vector<MyMatrix<Tint>> l_matr;
  Telt get_elt(MyMatrix<Tint> const &M) const {
    int n_V = V.rows();
    std::vector<Tidx> v(n_V);
    for (int u = 0; u < n_V; u++) {
      MyVector<Tint> V1 = GetMatrixRow(V, u);
      MyVector<Tint> V2 = M.transpose() * V1;
      auto iter = map_V.find(V2);
      if (iter == map_V.end()) {
        std::cerr << "PERFPRES: The matrix does not stabilize the perfect "
                     "domain\n";
        throw TerminalException{1};
      }
      v[u] = iter->second;
    }
    return Telt(v);
  }
  permutalib::PresentationWord get_word(MyMatrix<Tint> const &M) const {
    return solver.get_word(get_elt(M));
  }
};

template <typename Tint, typename Tidx>
std::unordered_map<MyVector<Tint>, Tidx>
get_map_vectors(MyMatrix<Tint> const &V) {
  std::unordered_map<MyVector<Tint>, Tidx> map;
  int n_V = V.rows();
  for (int u = 0; u < n_V; u++) {
    MyVector<Tint> eV = GetMatrixRow(V, u);
    map[eV] = u;
  }
  return map;
}

template <typename Tint, typename Tgroup>
PerfectDomainGroup<Tint, Tgroup>
get_perfect_domain_group(MyMatrix<Tint> const &EXT, Tgroup const &GRP_ext,
                         std::vector<MyMatrix<Tint>> const &l_matr_gens) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = EXT.cols();
  auto get_group = [&]() -> std::pair<MyMatrix<Tint>, Tgroup> {
    if (RankMat(EXT) == n) {
      return {EXT, GRP_ext};
    }
    // The orbit of the vectors of EXT and of the standard basis
    std::vector<MyVector<Tint>> l_vect;
    std::unordered_set<MyVector<Tint>> set_vect;
    auto f_insert = [&](MyVector<Tint> const &eV) -> void {
      if (set_vect.insert(eV).second) {
        l_vect.push_back(eV);
      }
    };
    for (int u = 0; u < EXT.rows(); u++) {
      f_insert(GetMatrixRow(EXT, u));
    }
    for (int i = 0; i < n; i++) {
      MyVector<Tint> eV = ZeroVector<Tint>(n);
      eV(i) = 1;
      f_insert(eV);
    }
    size_t pos = 0;
    while (pos < l_vect.size()) {
      MyVector<Tint> eV = l_vect[pos];
      for (auto &eMatr : l_matr_gens) {
        MyVector<Tint> eV_img = eMatr.transpose() * eV;
        f_insert(eV_img);
      }
      pos++;
    }
    MyMatrix<Tint> V = MatrixFromVectorFamily(l_vect);
    std::unordered_map<MyVector<Tint>, Tidx> map_V =
        get_map_vectors<Tint, Tidx>(V);
    int n_V = V.rows();
    std::vector<Telt> l_gens;
    for (auto &eMatr : l_matr_gens) {
      std::vector<Tidx> v(n_V);
      for (int u = 0; u < n_V; u++) {
        MyVector<Tint> eV_img = eMatr.transpose() * l_vect[u];
        v[u] = map_V.at(eV_img);
      }
      l_gens.push_back(Telt(v));
    }
    return {V, Tgroup(l_gens, n_V)};
  };
  std::pair<MyMatrix<Tint>, Tgroup> pair = get_group();
  MyMatrix<Tint> const &V = pair.first;
  std::unordered_map<MyVector<Tint>, Tidx> map_V =
      get_map_vectors<Tint, Tidx>(V);
  FiniteGroupWordSolver<Tgroup> solver(pair.second);
  std::vector<MyMatrix<Tint>> l_matr;
  for (auto &ePerm : solver.pres.generators) {
    l_matr.push_back(FindTransformation(V, V, ePerm));
  }
  return {V, std::move(map_V), std::move(solver), std::move(l_matr)};
}

/*
  The presentation of the group acting on the perfect form complex. The
  levels 0, 1, 2 of the complex have to be computed.
 */
template <typename T, typename Tint, typename Tgroup>
permutalib::Presentation<MyMatrix<Tint>>
get_presentation_unsimplified_fce(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                  std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using PresentationWord = permutalib::PresentationWord;
  std::vector<PerfectFormInfoForComplex<T, Tint, Tgroup>> const &l_perfect =
      fce.pctdi.l_perfect;
  int n = fce.pctdi.LinSpa.n;
  size_t n_perfect = l_perfect.size();
  std::vector<MyMatrix<Tint>> l_gen;
  std::vector<PresentationWord> l_rel;
  //
  // The generators
  //
  std::vector<PerfectDomainGroup<Tint, Tgroup>> l_dom;
  std::vector<int> l_offset_stab;
  for (auto &perfect : l_perfect) {
    std::vector<MyMatrix<Tint>> l_matr_gens;
    if (perfect.opt_pre_imager) {
      l_matr_gens = perfect.opt_pre_imager->get_list_matr_gens();
    }
    l_dom.push_back(get_perfect_domain_group(perfect.EXT, perfect.GRP_ext,
                                             l_matr_gens));
    l_offset_stab.push_back(l_gen.size());
    for (auto &eMatr : l_dom.back().l_matr) {
      l_gen.push_back(eMatr);
    }
  }
  std::vector<int> l_offset_side;
  for (auto &perfect : l_perfect) {
    l_offset_side.push_back(l_gen.size());
    for (auto &adj : perfect.l_sing_adj) {
      l_gen.push_back(adj.eMat);
    }
  }
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: n_perfect=" << n_perfect << " |l_gen|=" << l_gen.size()
     << "\n";
#endif
  auto get_stab_word = [&](size_t i, MyMatrix<Tint> const &M)
      -> PresentationWord {
    PresentationWord w = l_dom[i].get_word(M);
    for (auto &x : w) {
      x += (x > 0) ? l_offset_stab[i] : -l_offset_stab[i];
    }
    return w;
  };
  auto get_side_letter = [&](size_t i, size_t k) -> int {
    return l_offset_side[i] + static_cast<int>(k) + 1;
  };
  auto append = [](PresentationWord &w1, PresentationWord const &w2) -> void {
    w1.insert(w1.end(), w2.begin(), w2.end());
  };
  auto get_face = [&](size_t i, MyMatrix<Tint> const &M) -> Face {
    std::unordered_map<MyVector<Tint>, typename Telt::Tidx> const &map_V =
        l_dom[i].map_V;
    int n_ext = l_perfect[i].EXT.rows();
    Face f(n_ext);
    for (int u = 0; u < M.rows(); u++) {
      MyVector<Tint> eV = GetMatrixRow(M, u);
      auto iter = map_V.find(eV);
      if (iter == map_V.end() || static_cast<int>(iter->second) >= n_ext) {
        std::cerr << "PERFPRES: The vector does not belong to the perfect "
                     "domain\n";
        throw TerminalException{1};
      }
      f[iter->second] = 1;
    }
    return f;
  };
  // Whether g maps the vector configuration of v_i to itself.
  auto is_stabilizing = [&](size_t i, MyMatrix<Tint> const &g) -> bool {
    MyMatrix<Tint> const &EXT = l_perfect[i].EXT;
    int n_ext = EXT.rows();
    for (int u = 0; u < n_ext; u++) {
      MyVector<Tint> eV = GetMatrixRow(EXT, u);
      MyVector<Tint> eV_img = g.transpose() * eV;
      auto iter = l_dom[i].map_V.find(eV_img);
      if (iter == l_dom[i].map_V.end() ||
          static_cast<int>(iter->second) >= n_ext) {
        return false;
      }
    }
    return true;
  };
  auto get_face_vectors = [&](size_t i, Face const &f) -> MyMatrix<Tint> {
    return SelectRow(l_perfect[i].EXT, f);
  };
  // The facet f of v_i written as f_k h with h in Gamma_i.
  auto get_facet_decomposition = [&](size_t i, Face const &f)
      -> std::pair<size_t, MyMatrix<Tint>> {
    PerfectFormInfoForComplex<T, Tint, Tgroup> const &perfect = l_perfect[i];
    size_t n_adj = perfect.l_sing_adj.size();
    for (size_t k = 0; k < n_adj; k++) {
      Face const &f_k = perfect.l_sing_adj[k].f_ext;
      if (f_k.count() != f.count()) {
        continue;
      }
      for (auto &pair : perfect.orbit_representatives(f_k)) {
        if (pair.first == f) {
          return {k, perfect.find_matrix(pair.second, os)};
        }
      }
    }
    std::cerr << "PERFPRES: Failed to find the facet in the list of facets\n";
    throw TerminalException{1};
  };
  //
  // R0: The relators of the stabilizers
  //
  for (size_t i = 0; i < n_perfect; i++) {
    for (auto &rel : l_dom[i].solver.pres.relators) {
      PresentationWord w = rel;
      for (auto &x : w) {
        x += (x > 0) ? l_offset_stab[i] : -l_offset_stab[i];
      }
      l_rel.push_back(w);
    }
  }
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: R0, |l_rel|=" << l_rel.size() << "\n";
#endif
  //
  // R1 and R2: The relators of the edges
  //
  for (size_t i = 0; i < n_perfect; i++) {
    PerfectFormInfoForComplex<T, Tint, Tgroup> const &perfect = l_perfect[i];
    size_t n_adj = perfect.l_sing_adj.size();
    for (size_t k = 0; k < n_adj; k++) {
      sing_adj<Tint> const &adj = perfect.l_sing_adj[k];
      size_t j = adj.jCone;
      MyMatrix<Tint> const &s = adj.eMat;
      MyMatrix<Tint> s_inv = Inverse(s);
      int letter_s = get_side_letter(i, k);
      // R1
      int n_V = l_dom[i].V.rows();
      Face f_V(n_V);
      for (int u = 0; u < static_cast<int>(adj.f_ext.size()); u++) {
        if (adj.f_ext[u] == 1) {
          f_V[u] = 1;
        }
      }
      Tgroup GRP_edge = l_dom[i].solver.GRP.Stabilizer_OnSets(f_V);
      for (auto &ePerm : GRP_edge.GeneratorsOfGroup()) {
        MyMatrix<Tint> g = FindTransformation(l_dom[i].V, l_dom[i].V, ePerm);
        MyMatrix<Tint> c = s * g * s_inv;
        PresentationWord w{letter_s};
        append(w, get_stab_word(i, g));
        w.push_back(-letter_s);
        append(w, permutalib::InverseWord(get_stab_word(j, c)));
        l_rel.push_back(w);
      }
      // R2
      MyMatrix<Tint> EXT_facet = get_face_vectors(i, adj.f_ext) * s_inv;
      Face f_rev = get_face(j, EXT_facet);
      std::pair<size_t, MyMatrix<Tint>> pair = get_facet_decomposition(j, f_rev);
      size_t k_rev = pair.first;
      MyMatrix<Tint> const &h_rev = pair.second;
      sing_adj<Tint> const &adj_rev = l_perfect[j].l_sing_adj[k_rev];
      if (adj_rev.jCone != i) {
        std::cerr << "PERFPRES: The reversed edge does not go back to the "
                     "same perfect domain\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> c_rev = adj_rev.eMat * h_rev * s;
      PresentationWord w{get_side_letter(j, k_rev)};
      append(w, get_stab_word(j, h_rev));
      w.push_back(letter_s);
      append(w, permutalib::InverseWord(get_stab_word(i, c_rev)));
      l_rel.push_back(w);
    }
  }
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: R1/R2, |l_rel|=" << l_rel.size() << "\n";
#endif
  //
  // R3: The relators of the 2-cells
  //
  if (fce.levels.size() > 2) {
    for (auto &face : fce.levels[2].l_faces) {
      if (!face.is_well_rounded) {
        continue;
      }
      triple<Tint> const &t = face.l_triple[0];
      size_t i_start = t.iCone;
      MyMatrix<Tint> EXT_cell = get_face_vectors(i_start, t.f_ext);
      size_t i = i_start;
      MyMatrix<Tint> g = IdentityMat<Tint>(n);
      PresentationWord w_g;
      std::optional<Face> f_prev;
      size_t n_step = 0;
      while (true) {
        Face f_cell = get_face(i, EXT_cell * Inverse(g));
        // The facets of v_i containing the 2-cell
        std::vector<Face> l_facet;
        for (auto &adj : l_perfect[i].full_adjacencies(os)) {
          if (f_cell.is_subset_of(adj.f_ext)) {
            l_facet.push_back(adj.f_ext);
          }
        }
        if (l_facet.size() != 2) {
          std::cerr << "PERFPRES: The 2-cell should be contained in two "
                       "facets, |l_facet|="
                    << l_facet.size() << "\n";
          throw TerminalException{1};
        }
        Face f_next = l_facet[0];
        if (f_prev && *f_prev == f_next) {
          f_next = l_facet[1];
        }
        std::pair<size_t, MyMatrix<Tint>> pair = get_facet_decomposition(i, f_next);
        size_t k = pair.first;
        MyMatrix<Tint> const &h = pair.second;
        sing_adj<Tint> const &adj = l_perfect[i].l_sing_adj[k];
        size_t j = adj.jCone;
        MyMatrix<Tint> step = adj.eMat * h;
        PresentationWord w{get_side_letter(i, k)};
        append(w, get_stab_word(i, h));
        append(w, w_g);
        w_g = w;
        g = step * g;
        f_prev = get_face(j, get_face_vectors(i, f_next) * Inverse(step));
        i = j;
        n_step++;
        if (i == i_start && is_stabilizing(i, g)) {
          break;
        }
        if (n_step > 1000000) {
          std::cerr << "PERFPRES: The cycle around the 2-cell is too long\n";
          throw TerminalException{1};
        }
      }
#ifdef DEBUG_PERFECT_PRESENTATION
      os << "PERFPRES: R3, cycle of length " << n_step << "\n";
#endif
      append(w_g, permutalib::InverseWord(get_stab_word(i_start, g)));
      l_rel.push_back(w_g);
    }
  }
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: R3, |l_rel|=" << l_rel.size() << "\n";
#endif
  //
  // R4: The tree of representatives
  //
  std::vector<std::optional<MyMatrix<Tint>>> l_t(n_perfect);
  l_t[0] = IdentityMat<Tint>(n);
  std::vector<size_t> l_queue{0};
  size_t pos_queue = 0;
  while (pos_queue < l_queue.size()) {
    size_t i = l_queue[pos_queue];
    pos_queue++;
    size_t n_adj = l_perfect[i].l_sing_adj.size();
    for (size_t k = 0; k < n_adj; k++) {
      sing_adj<Tint> const &adj = l_perfect[i].l_sing_adj[k];
      size_t j = adj.jCone;
      if (!l_t[j]) {
        l_t[j] = adj.eMat * *l_t[i];
        l_queue.push_back(j);
        l_rel.push_back({get_side_letter(i, k)});
      }
    }
  }
  if (l_queue.size() != n_perfect) {
    std::cerr << "PERFPRES: The perfect domains are not connected\n";
    throw TerminalException{1};
  }
  for (size_t i = 0; i < n_perfect; i++) {
    MyMatrix<Tint> const &t_i = *l_t[i];
    MyMatrix<Tint> t_i_inv = Inverse(t_i);
    size_t n_stab = l_dom[i].l_matr.size();
    for (size_t u = 0; u < n_stab; u++) {
      MyMatrix<Tint> &eMatr = l_gen[l_offset_stab[i] + u];
      eMatr = t_i_inv * eMatr * t_i;
    }
    size_t n_adj = l_perfect[i].l_sing_adj.size();
    for (size_t k = 0; k < n_adj; k++) {
      size_t j = l_perfect[i].l_sing_adj[k].jCone;
      MyMatrix<Tint> &eMatr = l_gen[l_offset_side[i] + k];
      eMatr = Inverse(*l_t[j]) * eMatr * t_i;
    }
  }
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: R4, |l_rel|=" << l_rel.size() << "\n";
#endif
#ifdef SANITY_CHECK_PERFECT_PRESENTATION
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  for (auto &rel : l_rel) {
    if (EvaluateWordMatrix(l_gen, n, rel) != id) {
      std::cerr << "PERFPRES: A relator is not the identity\n";
      throw TerminalException{1};
    }
  }
#endif
  return {std::move(l_gen), permutalib::ReducedListRelators(l_rel)};
}

template <typename T, typename Tint, typename Tgroup>
permutalib::Presentation<MyMatrix<Tint>>
get_presentation_fce(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                     std::ostream &os) {
  permutalib::Presentation<MyMatrix<Tint>> P =
      get_presentation_unsimplified_fce(fce, os);
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: Before simplification |gens|=" << P.generators.size()
     << " |rels|=" << P.relators.size()
     << " length=" << permutalib::TotalLength(P.relators) << "\n";
#endif
  // The elimination of the generators is stopped when the total length of
  // the relators would be more than twice the one of the unsimplified
  // presentation: the last eliminations make the relators much longer.
  size_t n_keep = 0;
  size_t max_total_length = 2 * permutalib::TotalLength(P.relators);
  permutalib::Presentation<MyMatrix<Tint>> P_simp =
      permutalib::SimplifyPresentation(P, n_keep, max_total_length);
#ifdef DEBUG_PERFECT_PRESENTATION
  os << "PERFPRES: After simplification |gens|=" << P_simp.generators.size()
     << " |rels|=" << P_simp.relators.size()
     << " length=" << permutalib::TotalLength(P_simp.relators) << "\n";
#endif
#ifdef SANITY_CHECK_PERFECT_PRESENTATION
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  for (auto &rel : P_simp.relators) {
    if (EvaluateWordMatrix(P_simp.generators, n, rel) != id) {
      std::cerr << "PERFPRES: A simplified relator is not the identity\n";
      throw TerminalException{1};
    }
  }
#endif
  return P_simp;
}

// The GAP code of a record with the finitely presented group and the
// matrices of its generators.
template <typename Tint>
void WritePresentationGAP(std::ostream &os_out,
                          permutalib::Presentation<MyMatrix<Tint>> const &P) {
  std::string F = "F";
  os_out << "local F, G;\n";
  os_out << "F:=FreeGroup(" << P.generators.size() << ");\n";
  os_out << "G:=F/[";
  bool is_first = true;
  for (auto &rel : P.relators) {
    if (!is_first) {
      os_out << ",";
    }
    is_first = false;
    os_out << "\n" << permutalib::GapStringWord(rel, F);
  }
  os_out << "];\n";
  os_out << "return rec(fp:=G, ListGen:=";
  WriteListMatrixGAP(os_out, P.generators);
  os_out << ");\n";
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_PRESENTATION_H_
// clang-format on
