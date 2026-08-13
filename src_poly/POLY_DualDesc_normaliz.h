// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESC_NORMALIZ_H_
#define SRC_POLY_POLY_DUALDESC_NORMALIZ_H_

// clang-format off
#include "POLY_DualDesc_beneath_and_beyond.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixInverse.h"
#include <algorithm>
#include <deque>
#include <list>
#include <map>
#include <memory>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_NORMALIZ_DUAL_DESC
#endif

#ifdef DISABLE_DEBUG_NORMALIZ_DUAL_DESC
#undef DEBUG_NORMALIZ_DUAL_DESC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_NORMALIZ_DUAL_DESC
#endif

#ifdef TIMINGS
#define TIMINGS_NORMALIZ_DUAL_DESC
#endif

// Native port of the Normaliz primal support-hyperplane algorithm
// (Fourier-Motzkin with pyramid decomposition), following
//  * W. Bruns, B. Ichim, "Normaliz: Algorithms for affine monoids and rational
//    cones", J. Algebra 324 (2010), https://arxiv.org/abs/0910.2845
//  * W. Bruns, B. Ichim, C. Soeger, "The power of pyramid decomposition in
//    Normaliz", https://arxiv.org/abs/1206.1916
// and the reference implementation in normaliz-3.11.1 full_cone.cpp
// (build_cone / find_new_facets / process_pyramids / process_pyramid /
// select_supphyps_from / match_neg_hyp_with_pos_hyps).
//
// Same contract as the beneath-and-beyond backend: EXT is the set of extreme
// rays of a *full-dimensional pointed cone*, one ray per row. A facet is a
// linear form f with f . x >= 0 on every row, f . x == 0 exactly on the
// incident rows. Low-dimensional or redundancy-carrying inputs must be
// preprocessed (or fed to cdd/lrs) as elsewhere in src_poly.
//
// What is ported:
//  * The incremental build: generators sorted lexicographically, start
//    simplex = lex-first rank-d row subset, then one Fourier-Motzkin step per
//    generator with the three-way facet classification (positive / negative /
//    neutral) and removal of the negative facets.
//  * find_new_facets: the separate treatment of simplicial facets (map-based
//    subfacet matching) and non-simplicial ones (mother/daughter shortcut,
//    birth-epoch "extension" test, missing-generator counting, then the rank
//    test vs comparison test decision).
//  * Pyramid decomposition: when |Pos|*|Neg| crosses the recursion bound the
//    FM matching is replaced by pyramids (apex = new generator, base = a
//    visible facet). Simplicial pyramids are solved by one kernel-vector
//    computation per facet; small ones recurse into a child cone that hands
//    its support hyperplanes back through select_supphyps_from (whose
//    scalar-product filter is at once the validity check and the
//    deduplication); large ones (largePyramidFactor * Comparisons[k] >
//    #hyperplanes) are deferred and matched one negative facet against all
//    positive facets in match_neg_hyp_with_pos_hyps.
//  * All arithmetic is ring arithmetic (FM combination, gcd content
//    reduction, division-free rank tests via SelectIndependentRows, all the
//    facet normals of a simplicial cone at once from one fraction-free
//    adjugate), so the kernel runs on the underlying integer ring of a field
//    input -- and on TryInt64 with exact-ring fallback, as in the
//    beneath-and-beyond wrapper.
//
// What is deliberately not ported: OpenMP parallelization and its buffers,
// triangulations / Hilbert series / multiplicity, the floating-point rank
// predictor, and the runtime-measurement heuristics (their a-priori formulas
// are kept).

namespace normaliz_dd {

// The recursion bound for pyramid decomposition: pyramids are built when
// nr_neg*nr_pos - nr_neg_simp*nr_pos_simp >= dim * SuppHypRecursionFactor
// (times arith_cost_factor for non-machine arithmetic). Constants from
// full_cone.cpp.
constexpr size_t nmz_SuppHypRecursionFactor = 320000;
constexpr size_t nmz_largePyramidFactor = 20;
// Factor by which multiprecision arithmetic is assumed slower than machine
// integers (GMP_time_factor in full_cone.cpp). Machine types get 1.
template <typename Tint>
constexpr size_t nmz_arith_cost_factor() {
  if constexpr (std::is_fundamental<Tint>::value)
    return 1;
  else
    return 10;
}

// FACETDATA of full_cone.cpp. GenInHyp is indexed by kernel-local generator
// indices; Hyp lives in the ambient space.
template <typename Tint> struct NmzFacet {
  MyVector<Tint> Hyp;
  Face GenInHyp;
  Tint ValNewGen;
  size_t BornAt;
  size_t Ident;
  size_t Mother;
  bool simplicial;
  bool neutral;
  bool positive;
  bool negative;
};

// One cone in the recursion (the top cone or a pyramid). Generators are not
// copied: every level references rows of the one top matrix through RowIdx,
// so the rank tests and simplex computations all run on the same matrix.
template <typename Tint> struct NmzKernel {
  MyMatrix<Tint> const &TopGen;
  std::ostream &os;
  // local generator i  <->  row RowIdx[i] of TopGen
  std::vector<int> RowIdx;
  size_t dim;
  size_t nr_gen;
  bool is_pyramid;
  NmzKernel *Mother;
  // local indices (in Mother) of this pyramid's generators; entry 0 is the
  // apex, the rest ascending
  std::vector<size_t> Mother_Key;
  size_t apex;

  std::vector<uint8_t> in_triang;
  std::vector<size_t> GensInCone;
  size_t nrGensInCone;
  std::list<NmzFacet<Tint>> Facets;
  size_t old_nr_supp_hyps;
  size_t HypCounter;
  std::vector<size_t> Comparisons;
  size_t nrTotalComparisons;
  std::list<NmzFacet<Tint>> LargeRecPyrs;

  // top cone
  NmzKernel(MyMatrix<Tint> const &_TopGen, std::ostream &_os)
      : TopGen(_TopGen), os(_os),
        dim(_TopGen.cols()), nr_gen(_TopGen.rows()), is_pyramid(false),
        Mother(nullptr), apex(0), in_triang(nr_gen, 0), nrGensInCone(0),
        old_nr_supp_hyps(0), HypCounter(1), nrTotalComparisons(0) {
    RowIdx.resize(nr_gen);
    for (size_t i = 0; i < nr_gen; i++)
      RowIdx[i] = static_cast<int>(i);
  }

  // pyramid: generators are Mother's gens selected by Key (apex first)
  NmzKernel(NmzKernel &C, std::vector<size_t> const &Key)
      : TopGen(C.TopGen), os(C.os), dim(C.dim),
        nr_gen(Key.size()), is_pyramid(true), Mother(&C), Mother_Key(Key),
        apex(0), in_triang(nr_gen, 0), nrGensInCone(0), old_nr_supp_hyps(0),
        HypCounter(1), nrTotalComparisons(0) {
    RowIdx.resize(nr_gen);
    for (size_t i = 0; i < nr_gen; i++)
      RowIdx[i] = C.RowIdx[Key[i]];
  }

  Tint v_scal(MyVector<Tint> const &Hyp, size_t i_gen) const {
    Tint sum(0);
    int row = RowIdx[i_gen];
    for (size_t k = 0; k < dim; k++)
      AddMul(sum, Hyp(k), TopGen(row, k));
    return sum;
  }

  void number_hyperplane(NmzFacet<Tint> &hyp, size_t born_at, size_t mother) {
    hyp.Mother = mother;
    hyp.BornAt = born_at;
    hyp.Ident = HypCounter;
    HypCounter++;
  }

  void set_simplicial(NmzFacet<Tint> &hyp) {
    size_t nr_gen_in_hyp = 0;
    for (size_t i = 0; i < nr_gen; i++)
      if (in_triang[i] && hyp.GenInHyp.test(i))
        nr_gen_in_hyp++;
    hyp.simplicial = (nr_gen_in_hyp == dim - 2);
  }

  // Fourier-Motzkin combination of a positive / negative pair:
  // NewHyp = pos.ValNewGen * neg.Hyp - neg.ValNewGen * pos.Hyp, made
  // primitive. Ring arithmetic only.
  void add_hyperplane(size_t new_generator, NmzFacet<Tint> const &positive,
                      NmzFacet<Tint> const &negative,
                      std::list<NmzFacet<Tint>> &NewHyps,
                      bool known_to_be_simplicial) {
    NmzFacet<Tint> NewFacet;
    NewFacet.Hyp = MyVector<Tint>(dim);
    for (size_t k = 0; k < dim; k++) {
      NewFacet.Hyp(k) = positive.ValNewGen * negative.Hyp(k) -
                        negative.ValNewGen * positive.Hyp(k);
    }
    NormalizeVectorContent(NewFacet.Hyp);
    NewFacet.ValNewGen = 0;
    NewFacet.GenInHyp = positive.GenInHyp & negative.GenInHyp;
    if (known_to_be_simplicial)
      NewFacet.simplicial = true;
    else
      set_simplicial(NewFacet);
    NewFacet.GenInHyp.set(new_generator);
    NewFacet.neutral = false;
    NewFacet.positive = false;
    NewFacet.negative = false;
    number_hyperplane(NewFacet, nrGensInCone, positive.Ident);
    NewHyps.emplace_back(std::move(NewFacet));
  }

  // All dim facet normals of the simplicial cone spanned by the local
  // generators key (|key| == dim, linearly independent) at once: with S the
  // generator submatrix, S * adj(S) = det(S) * Id, so column j of the
  // adjugate vanishes on every key generator except key[j] -- one
  // fraction-free elimination instead of dim kernel solves (simplex_data in
  // the original). Each normal is content-reduced and oriented positive on
  // its opposite generator.
  std::vector<MyVector<Tint>>
  simplex_facet_normals(std::vector<size_t> const &key) {
    MyMatrix<Tint> S(dim, dim);
    for (size_t k = 0; k < dim; k++) {
      int row = RowIdx[key[k]];
      for (size_t u = 0; u < dim; u++)
        S(k, u) = TopGen(row, u);
    }
    std::pair<MyMatrix<Tint>, Tint> pair = AdjugateDeterminant(S);
    std::vector<MyVector<Tint>> normals;
    normals.reserve(dim);
    for (size_t j = 0; j < dim; j++) {
      MyVector<Tint> raw(dim);
      for (size_t u = 0; u < dim; u++)
        raw(u) = pair.first(u, j);
      MyVector<Tint> normal = beneath_beyond::canonicalize_normal(raw);
      if (v_scal(normal, key[j]) < 0)
        normal = -normal;
      normals.emplace_back(std::move(normal));
    }
    return normals;
  }

  // find_and_evaluate_start_simplex: the lex-first rank-dim subset of the
  // local generators, its facets, and the bookkeeping start values.
  void find_and_evaluate_start_simplex() {
    std::vector<int> candidates(RowIdx.begin(), RowIdx.end());
    std::vector<int> sel = SelectIndependentRows(TopGen, candidates, dim);
#ifdef SANITY_CHECK_NORMALIZ_DUAL_DESC
    if (sel.size() != dim) {
      std::cerr << "NMZ: start simplex has rank " << sel.size()
                << " but dim=" << dim << " (input not full-dimensional?)\n";
      throw TerminalException{1};
    }
#endif
    // map top rows back to local indices (RowIdx is increasing in local
    // index order per construction, and sel preserves candidate order)
    std::vector<size_t> key;
    key.reserve(dim);
    size_t pos = 0;
    for (auto &row : sel) {
      while (RowIdx[pos] != row)
        pos++;
      key.push_back(pos);
    }
    std::vector<MyVector<Tint>> normals = simplex_facet_normals(key);
    for (size_t j = 0; j < dim; j++) {
      NmzFacet<Tint> NewFacet;
      NewFacet.Hyp = std::move(normals[j]);
      NewFacet.GenInHyp = Face(nr_gen);
      for (size_t k = 0; k < dim; k++)
        if (k != j)
          NewFacet.GenInHyp.set(key[k]);
      NewFacet.ValNewGen = 0;
      NewFacet.simplicial = true;
      NewFacet.neutral = false;
      NewFacet.positive = false;
      NewFacet.negative = false;
      number_hyperplane(NewFacet, 0, 0);
      Facets.emplace_back(std::move(NewFacet));
    }
    for (size_t j = 0; j < dim; j++) {
      in_triang[key[j]] = 1;
      GensInCone.push_back(key[j]);
    }
    nrGensInCone = dim;
    nrTotalComparisons = dim * dim / 2;
    if constexpr (nmz_arith_cost_factor<Tint>() != 1)
      nrTotalComparisons *= (nmz_arith_cost_factor<Tint>() / 4);
  }

  // Division-free rank test: do the top rows of the listed local generators
  // have rank >= target?
  bool rank_at_least(std::vector<size_t> const &local_key, size_t target) {
    std::vector<int> rows;
    rows.reserve(local_key.size());
    for (auto &i : local_key)
      rows.push_back(RowIdx[i]);
    return SelectIndependentRows(TopGen, rows, target).size() == target;
  }

  // find_new_facets: the Fourier-Motzkin step, one insertion. The facets are
  // already classified (positive / neutral / negative, ValNewGen) by
  // build_cone.
  void find_new_facets(size_t new_generator) {
    if (dim <= 1)
      return;
    size_t i, j, k;
    size_t subfacet_dim = dim - 2;
    size_t facet_dim = dim - 1;

    std::deque<NmzFacet<Tint> *> Pos_Simp, Pos_Non_Simp;
    std::deque<NmzFacet<Tint> *> Neg_Simp, Neg_Non_Simp;
    std::deque<NmzFacet<Tint> *> Neutral_Simp, Neutral_Non_Simp;

    Face GenInPosHyp(nr_gen), GenInNegHyp(nr_gen);
    for (auto &facet : Facets) {
      if (facet.positive)
        GenInPosHyp |= facet.GenInHyp;
      if (facet.negative)
        GenInNegHyp |= facet.GenInHyp;
    }
    Face Gen_BothSides = GenInPosHyp & GenInNegHyp;
    std::vector<size_t> Gen_BothSides_key;
    for (i = 0; i < nr_gen; ++i)
      if (Gen_BothSides[i])
        Gen_BothSides_key.push_back(i);

    for (auto &facet : Facets) {
      bool simplex = facet.simplicial;
      if (facet.neutral) {
        facet.GenInHyp.set(new_generator);
        facet.simplicial = false;
        if (simplex)
          Neutral_Simp.push_back(&facet);
        else
          Neutral_Non_Simp.push_back(&facet);
        continue;
      }
      size_t nr_relevant_gens = 0;
      for (auto &g : Gen_BothSides_key)
        if (facet.GenInHyp[g])
          nr_relevant_gens++;
      if (nr_relevant_gens < subfacet_dim)
        continue;
      if (facet.positive) {
        if (simplex)
          Pos_Simp.push_back(&facet);
        else
          Pos_Non_Simp.push_back(&facet);
      } else if (facet.negative) {
        if (simplex)
          Neg_Simp.push_back(&facet);
        else
          Neg_Non_Simp.push_back(&facet);
      }
    }
    size_t nr_PosSimp = Pos_Simp.size();
    size_t nr_PosNonSimp = Pos_Non_Simp.size();
    size_t nr_NegSimp = Neg_Simp.size();
    size_t nr_NegNonSimp = Neg_Non_Simp.size();
    size_t nr_NeuSimp = Neutral_Simp.size();
    size_t nr_NeuNonSimp = Neutral_Non_Simp.size();

    // subfacets of the negative simplicial facets
    std::list<std::pair<Face, int>> Neg_Subfacet_Multi;
    {
      for (i = 0; i < nr_NegSimp; i++) {
        Face RelGen_NegHyp = Gen_BothSides & Neg_Simp[i]->GenInHyp;
        size_t nr_RelGen_NegHyp = 0;
        for (j = 0; j < nr_gen; j++) {
          if (RelGen_NegHyp.test(j))
            nr_RelGen_NegHyp++;
          if (nr_RelGen_NegHyp > subfacet_dim)
            break;
        }
        if (nr_RelGen_NegHyp == subfacet_dim)
          Neg_Subfacet_Multi.emplace_back(RelGen_NegHyp,
                                          static_cast<int>(i));
        if (nr_RelGen_NegHyp == facet_dim) {
          for (k = 0; k < nr_gen; k++) {
            if (RelGen_NegHyp.test(k)) {
              Face subfacet = RelGen_NegHyp;
              subfacet.reset(k);
              Neg_Subfacet_Multi.emplace_back(subfacet, static_cast<int>(i));
            }
          }
        }
      }
    }
    Neg_Subfacet_Multi.sort();
    // a subfacet shared by two negative simplicial facets is interior to the
    // negative side: remove both copies
    for (auto jj = Neg_Subfacet_Multi.begin();
         jj != Neg_Subfacet_Multi.end();) {
      auto del = jj++;
      if (jj != Neg_Subfacet_Multi.end() && (*jj).first == (*del).first) {
        Neg_Subfacet_Multi.erase(del);
        del = jj++;
        Neg_Subfacet_Multi.erase(del);
      }
    }
    size_t nr_NegSubfMult = Neg_Subfacet_Multi.size();

    // remove those that lie in a neutral facet or a negative non-simplicial
    // facet (the size guards against a quadratic disaster are the original's)
    if (nr_NegSubfMult * (nr_NeuSimp + nr_NeuNonSimp + nr_NegNonSimp) <=
        100000000) {
      for (auto &entry : Neg_Subfacet_Multi) {
        Face const &subfacet = entry.first;
        bool found = false;
        if (nr_NeuSimp < 100000) {
          for (i = 0; i < nr_NeuSimp; i++) {
            found = subfacet.is_subset_of(Neutral_Simp[i]->GenInHyp);
            if (found)
              break;
          }
        }
        if (!found && nr_NeuNonSimp < 100000) {
          for (i = 0; i < nr_NeuNonSimp; i++) {
            found = subfacet.is_subset_of(Neutral_Non_Simp[i]->GenInHyp);
            if (found)
              break;
          }
        }
        if (!found && nr_NegNonSimp < 100000) {
          for (i = 0; i < nr_NegNonSimp; i++) {
            found = subfacet.is_subset_of(Neg_Non_Simp[i]->GenInHyp);
            if (found)
              break;
          }
        }
        if (found)
          entry.second = -1;
      }
    }
    std::map<Face, int> Neg_Subfacet;
    {
      auto last_inserted = Neg_Subfacet.begin();
      for (auto &entry : Neg_Subfacet_Multi)
        if (entry.second != -1)
          last_inserted = Neg_Subfacet.insert(last_inserted, entry);
    }
    size_t nr_NegSubf = Neg_Subfacet.size();
    Neg_Subfacet_Multi.clear();

    std::vector<std::list<NmzFacet<Tint>>> NewHypsSimp(nr_PosSimp);
    std::vector<std::list<NmzFacet<Tint>>> NewHypsNonSimp(nr_PosNonSimp);

    nrTotalComparisons += nr_NegNonSimp * nr_PosNonSimp;

    //
    // Positive simplicial vs negative simplicial and non-simplicial
    //
    std::vector<size_t> key(nr_gen);
    size_t nr_missing;
    bool common_subfacet;
    for (i = 0; i < nr_PosSimp; i++) {
      Face RelGen_PosHyp = Gen_BothSides & Pos_Simp[i]->GenInHyp;
      size_t nr_RelGen_PosHyp = 0;
      for (j = 0; j < nr_gen && nr_RelGen_PosHyp <= facet_dim; j++)
        if (RelGen_PosHyp.test(j)) {
          key[nr_RelGen_PosHyp] = j;
          nr_RelGen_PosHyp++;
        }
      if (nr_RelGen_PosHyp < subfacet_dim)
        continue;

      if (nr_RelGen_PosHyp == subfacet_dim) {
        auto jj_map = Neg_Subfacet.find(RelGen_PosHyp);
        if (jj_map != Neg_Subfacet.end()) {
          add_hyperplane(new_generator, *Pos_Simp[i],
                         *Neg_Simp[(*jj_map).second], NewHypsSimp[i], true);
          (*jj_map).second = -1;
        }
      }
      if (nr_RelGen_PosHyp == facet_dim) {
        for (k = 0; k < nr_gen; k++) {
          if (RelGen_PosHyp.test(k)) {
            Face subfacet = RelGen_PosHyp;
            subfacet.reset(k);
            auto jj_map = Neg_Subfacet.find(subfacet);
            if (jj_map != Neg_Subfacet.end()) {
              add_hyperplane(new_generator, *Pos_Simp[i],
                             *Neg_Simp[(*jj_map).second], NewHypsSimp[i],
                             true);
              (*jj_map).second = -1;
            }
          }
        }
      }

      for (j = 0; j < nr_NegNonSimp; j++) {
        nr_missing = 0;
        common_subfacet = true;
        for (k = 0; k < nr_RelGen_PosHyp; k++) {
          if (!Neg_Non_Simp[j]->GenInHyp.test(key[k])) {
            nr_missing++;
            if (nr_missing == 2 || nr_RelGen_PosHyp == subfacet_dim) {
              common_subfacet = false;
              break;
            }
          }
        }
        if (common_subfacet) {
          add_hyperplane(new_generator, *Pos_Simp[i], *Neg_Non_Simp[j],
                         NewHypsSimp[i], true);
          if (nr_RelGen_PosHyp == subfacet_dim)
            break;
        }
      }
    }

    //
    // Positive non-simplicial vs negative simplicial and non-simplicial
    //
    {
      // GenInHyp of every non-simplicial facet, for the comparison test;
      // successful reducers move to the front ("darwinistic" order)
      std::list<Face> Facets_0_1;
      for (i = 0; i < nr_PosNonSimp; ++i)
        Facets_0_1.push_back(Pos_Non_Simp[i]->GenInHyp);
      for (i = 0; i < nr_NegNonSimp; ++i)
        Facets_0_1.push_back(Neg_Non_Simp[i]->GenInHyp);
      for (i = 0; i < nr_NeuNonSimp; ++i)
        Facets_0_1.push_back(Neutral_Non_Simp[i]->GenInHyp);
      size_t nr_NonSimp = nr_PosNonSimp + nr_NegNonSimp + nr_NeuNonSimp;

      size_t missing_bound, nr_CommonGens;
      std::vector<size_t> common_key;
      common_key.reserve(nr_gen);
      std::vector<int> key_start(nrGensInCone);

      for (i = 0; i < nr_PosNonSimp; i++) {
        auto jj_map = Neg_Subfacet.begin();
        for (j = 0; j < nr_NegSubf; ++j, ++jj_map) {
          if ((*jj_map).second != -1) {
            if (jj_map->first.is_subset_of(Pos_Non_Simp[i]->GenInHyp)) {
              add_hyperplane(new_generator, *Pos_Non_Simp[i],
                             *Neg_Simp[(*jj_map).second], NewHypsNonSimp[i],
                             true);
              (*jj_map).second = -1;
            }
          }
        }

        NmzFacet<Tint> *PosHyp_Pointer = Pos_Non_Simp[i];
        Face RelGen_PosHyp = Gen_BothSides & PosHyp_Pointer->GenInHyp;
        size_t nr_RelGen_PosHyp = 0;
        int last_existing = -1;
        for (size_t jj = 0; jj < nrGensInCone; jj++) {
          j = GensInCone[jj];
          if (RelGen_PosHyp.test(j)) {
            key[nr_RelGen_PosHyp] = j;
            for (size_t kk = last_existing + 1; kk <= jj; kk++)
              key_start[kk] = static_cast<int>(nr_RelGen_PosHyp);
            nr_RelGen_PosHyp++;
            last_existing = static_cast<int>(jj);
          }
        }
        if (last_existing < static_cast<int>(nrGensInCone) - 1)
          for (size_t kk = last_existing + 1; kk < nrGensInCone; kk++)
            key_start[kk] = static_cast<int>(nr_RelGen_PosHyp);
        if (nr_RelGen_PosHyp < subfacet_dim)
          continue;

        missing_bound = nr_RelGen_PosHyp - subfacet_dim;

        for (j = 0; j < nr_NegNonSimp; j++) {
          NmzFacet<Tint> *NegHyp_Pointer = Neg_Non_Simp[j];

          if (PosHyp_Pointer->Ident == NegHyp_Pointer->Mother ||
              NegHyp_Pointer->Ident == PosHyp_Pointer->Mother) {
            // mother and daughter: their intersection is a subfacet
            add_hyperplane(new_generator, *PosHyp_Pointer, *NegHyp_Pointer,
                           NewHypsNonSimp[i], false);
            continue;
          }

          bool extension_test =
              PosHyp_Pointer->BornAt == NegHyp_Pointer->BornAt ||
              (PosHyp_Pointer->BornAt < NegHyp_Pointer->BornAt &&
               NegHyp_Pointer->Mother != 0) ||
              (NegHyp_Pointer->BornAt < PosHyp_Pointer->BornAt &&
               PosHyp_Pointer->Mother != 0);
          size_t both_existing_from =
              key_start[std::max(PosHyp_Pointer->BornAt,
                                 NegHyp_Pointer->BornAt)];

          nr_missing = 0;
          nr_CommonGens = 0;
          common_key.clear();
          size_t second_loop_bound = nr_RelGen_PosHyp;
          common_subfacet = true;

          // Two facets that are not mother and daughter cannot have shared a
          // subfacet at the birth of the younger one: they intersect in a
          // subfacet now only if a common generator arrived afterwards.
          if (extension_test) {
            bool extended = false;
            second_loop_bound = both_existing_from;
            for (k = both_existing_from; k < nr_RelGen_PosHyp; k++) {
              if (!NegHyp_Pointer->GenInHyp.test(key[k])) {
                nr_missing++;
                if (nr_missing > missing_bound) {
                  common_subfacet = false;
                  break;
                }
              } else {
                extended = true;
                common_key.push_back(key[k]);
                nr_CommonGens++;
              }
            }
            if (!extended || !common_subfacet)
              continue;
          }

          for (k = 0; k < second_loop_bound; k++) {
            if (!NegHyp_Pointer->GenInHyp.test(key[k])) {
              nr_missing++;
              if (nr_missing > missing_bound) {
                common_subfacet = false;
                break;
              }
            } else {
              common_key.push_back(key[k]);
              nr_CommonGens++;
            }
          }
          if (!common_subfacet)
            continue;

          if (subfacet_dim <= 2) {
            add_hyperplane(new_generator, *PosHyp_Pointer, *NegHyp_Pointer,
                           NewHypsNonSimp[i], false);
            continue;
          }

          // rank test vs comparison test, by the a-priori cost estimate
          bool ranktest =
              (nr_NonSimp > nmz_arith_cost_factor<Tint>() * dim * dim *
                                nr_CommonGens / 3);
          if (ranktest) {
            if (!rank_at_least(common_key, subfacet_dim))
              common_subfacet = false;
          } else {
            Face CommonGens = RelGen_PosHyp & NegHyp_Pointer->GenInHyp;
            for (auto a = Facets_0_1.begin(); a != Facets_0_1.end(); ++a) {
              if (CommonGens.is_subset_of(*a) &&
                  (*a != PosHyp_Pointer->GenInHyp) &&
                  (*a != NegHyp_Pointer->GenInHyp)) {
                common_subfacet = false;
                Facets_0_1.splice(Facets_0_1.begin(), Facets_0_1, a);
                break;
              }
            }
          }
          if (common_subfacet)
            add_hyperplane(new_generator, *PosHyp_Pointer, *NegHyp_Pointer,
                           NewHypsNonSimp[i], false);
        }
      }
    }

    for (i = 0; i < nr_PosSimp; i++)
      Facets.splice(Facets.end(), NewHypsSimp[i]);
    for (i = 0; i < nr_PosNonSimp; i++)
      Facets.splice(Facets.end(), NewHypsNonSimp[i]);
  }

  // select_supphyps_from: the mother cone (=this) selects the new global
  // support hyperplanes from the facet list of a daughter pyramid. A pyramid
  // facet containing the apex is a new global facet iff it is strictly
  // positive on every inserted generator outside the pyramid -- a value 0
  // would mean the same hyperplane is delivered by another pyramid, so this
  // one scalar-product filter is both the validity check and the
  // deduplication.
  void select_supphyps_from(std::list<NmzFacet<Tint>> &NewFacets,
                            size_t new_generator,
                            std::vector<size_t> const &Pyramid_key) {
    Face in_Pyr(nr_gen);
    for (size_t i = 0; i < Pyramid_key.size(); i++)
      in_Pyr.set(Pyramid_key[i]);
    for (auto &pyr_hyp : NewFacets) {
      if (!pyr_hyp.GenInHyp.test(0))
        continue;
      bool new_global_hyp = true;
      for (size_t i = 0; i < nr_gen; ++i) {
        if (in_Pyr.test(i) || !in_triang[i])
          continue;
        if (v_scal(pyr_hyp.Hyp, i) <= 0) {
          new_global_hyp = false;
          break;
        }
      }
      if (!new_global_hyp)
        continue;
      NmzFacet<Tint> NewFacet;
      NewFacet.Hyp = std::move(pyr_hyp.Hyp);
      NewFacet.GenInHyp = Face(nr_gen);
      for (size_t i = 0; i < Pyramid_key.size(); ++i)
        if (pyr_hyp.GenInHyp.test(i) && in_triang[Pyramid_key[i]])
          NewFacet.GenInHyp.set(Pyramid_key[i]);
      NewFacet.GenInHyp.set(new_generator);
      NewFacet.ValNewGen = 0;
      NewFacet.simplicial = pyr_hyp.simplicial;
      NewFacet.neutral = false;
      NewFacet.positive = false;
      NewFacet.negative = false;
      number_hyperplane(NewFacet, nrGensInCone, 0);
      Facets.emplace_back(std::move(NewFacet));
    }
  }

  // match_neg_hyp_with_pos_hyps: the Fourier-Motzkin step for one *large*
  // recursive pyramid, i.e. one negative facet matched directly against all
  // positive facets (more efficient than building the pyramid when the
  // pyramid is large).
  void match_neg_hyp_with_pos_hyps(NmzFacet<Tint> const &Neg,
                                   size_t new_generator,
                                   std::vector<NmzFacet<Tint> *> const &PosHyps,
                                   Face const &GenIn_PosHyp,
                                   std::list<Face> &Facets_0_1) {
    size_t missing_bound, nr_common_gens;
    std::vector<size_t> common_key;
    common_key.reserve(nr_gen);
    std::vector<size_t> key(nr_gen);
    bool common_subfacet;
    size_t subfacet_dim = dim - 2;
    size_t nr_missing;
    std::list<NmzFacet<Tint>> NewHyps;

    Face RelGens_InNegHyp = Neg.GenInHyp & GenIn_PosHyp;

    std::vector<int> key_start(nrGensInCone);
    size_t nr_RelGens_InNegHyp = 0;
    size_t j;
    int last_existing = -1;
    for (size_t jj = 0; jj < nrGensInCone; jj++) {
      j = GensInCone[jj];
      if (RelGens_InNegHyp.test(j)) {
        key[nr_RelGens_InNegHyp] = j;
        for (size_t kk = last_existing + 1; kk <= jj; kk++)
          key_start[kk] = static_cast<int>(nr_RelGens_InNegHyp);
        nr_RelGens_InNegHyp++;
        last_existing = static_cast<int>(jj);
      }
    }
    if (last_existing < static_cast<int>(nrGensInCone) - 1)
      for (size_t kk = last_existing + 1; kk < nrGensInCone; kk++)
        key_start[kk] = static_cast<int>(nr_RelGens_InNegHyp);

    if (nr_RelGens_InNegHyp < dim - 2)
      return;
    missing_bound = nr_RelGens_InNegHyp - subfacet_dim;

    for (auto const &Pos : PosHyps) {
      if (Neg.Ident == Pos->Mother || Pos->Ident == Neg.Mother) {
        add_hyperplane(new_generator, *Pos, Neg, NewHyps, false);
        continue;
      }
      bool extension_test =
          Neg.BornAt == Pos->BornAt ||
          (Neg.BornAt < Pos->BornAt && Pos->Mother != 0) ||
          (Pos->BornAt < Neg.BornAt && Neg.Mother != 0);
      size_t both_existing_from = key_start[std::max(Neg.BornAt, Pos->BornAt)];

      nr_missing = 0;
      nr_common_gens = 0;
      common_key.clear();
      size_t second_loop_bound = nr_RelGens_InNegHyp;
      common_subfacet = true;
      Face common_gens(nr_gen);

      if (extension_test) {
        bool extended = false;
        second_loop_bound = both_existing_from;
        for (size_t k = both_existing_from; k < nr_RelGens_InNegHyp; k++) {
          if (!Pos->GenInHyp.test(key[k])) {
            nr_missing++;
            if (nr_missing > missing_bound) {
              common_subfacet = false;
              break;
            }
          } else {
            extended = true;
            common_key.push_back(key[k]);
            common_gens.set(key[k]);
            nr_common_gens++;
          }
        }
        if (!extended || !common_subfacet)
          continue;
      }
      for (size_t k = 0; k < second_loop_bound; k++) {
        if (!Pos->GenInHyp.test(key[k])) {
          nr_missing++;
          if (nr_missing > missing_bound) {
            common_subfacet = false;
            break;
          }
        } else {
          common_key.push_back(key[k]);
          common_gens.set(key[k]);
          nr_common_gens++;
        }
      }
      if (!common_subfacet)
        continue;

      if (!Pos->simplicial) {
        bool ranktest = true;
        if constexpr (nmz_arith_cost_factor<Tint>() != 1)
          ranktest = (old_nr_supp_hyps > nmz_arith_cost_factor<Tint>() * dim *
                                             dim * nr_common_gens / 3);
        if (ranktest) {
          if (!rank_at_least(common_key, subfacet_dim))
            common_subfacet = false;
        } else {
          for (auto hp_t = Facets_0_1.begin(); hp_t != Facets_0_1.end();
               ++hp_t) {
            if (common_gens.is_subset_of(*hp_t) && (*hp_t != Neg.GenInHyp) &&
                (*hp_t != Pos->GenInHyp)) {
              Facets_0_1.splice(Facets_0_1.begin(), Facets_0_1, hp_t);
              common_subfacet = false;
              break;
            }
          }
        }
      }
      if (common_subfacet)
        add_hyperplane(new_generator, *Pos, Neg, NewHyps, false);
    }
    Facets.splice(Facets.end(), NewHyps);
  }

  void evaluate_large_rec_pyramids(size_t new_generator) {
    size_t nrLargeRecPyrs = LargeRecPyrs.size();
    if (nrLargeRecPyrs == 0)
      return;
    std::list<Face> Facets_0_1;
    {
      auto Fac = Facets.begin();
      for (size_t i = 0; i < old_nr_supp_hyps; ++i, ++Fac) {
        if (Fac->simplicial)
          continue;
        Facets_0_1.push_back(Fac->GenInHyp);
      }
    }
    std::vector<NmzFacet<Tint> *> PosHyps;
    Face GenIn_PosHyp(nr_gen);
    {
      auto ii = Facets.begin();
      for (size_t ij = 0; ij < old_nr_supp_hyps; ++ij, ++ii)
        if (ii->ValNewGen > 0) {
          GenIn_PosHyp |= ii->GenInHyp;
          PosHyps.push_back(&(*ii));
        }
    }
    nrTotalComparisons += PosHyps.size() * nrLargeRecPyrs;
#ifdef DEBUG_NORMALIZ_DUAL_DESC
    os << "NMZ: large pyramids " << nrLargeRecPyrs << "\n";
#endif
    for (auto &Neg : LargeRecPyrs)
      match_neg_hyp_with_pos_hyps(Neg, new_generator, PosHyps, GenIn_PosHyp,
                                  Facets_0_1);
    LargeRecPyrs.clear();
  }

  // process_pyramid: a simplicial pyramid is solved on the spot; a small one
  // recurses into a child cone (which hands its facets back at the end of its
  // build_cone); a large one is deferred to evaluate_large_rec_pyramids.
  void process_pyramid(std::vector<size_t> const &Pyramid_key,
                       size_t new_generator, NmzFacet<Tint> const &hyp) {
    if (Pyramid_key.size() == dim) {
      // simplicial pyramid: facet j of the pyramid contains all its
      // generators but the j-th
      std::list<NmzFacet<Tint>> NewFacets;
      std::vector<MyVector<Tint>> normals = simplex_facet_normals(Pyramid_key);
      for (size_t j = 0; j < dim; j++) {
        NmzFacet<Tint> NewFacet;
        NewFacet.Hyp = std::move(normals[j]);
        NewFacet.GenInHyp = Face(Pyramid_key.size());
        NewFacet.GenInHyp.set();
        NewFacet.GenInHyp.reset(j);
        NewFacet.ValNewGen = 0;
        NewFacet.simplicial = true;
        NewFacet.neutral = false;
        NewFacet.positive = false;
        NewFacet.negative = false;
        NewFacet.BornAt = 0;
        NewFacet.Ident = 0;
        NewFacet.Mother = 0;
        NewFacets.emplace_back(std::move(NewFacet));
      }
      select_supphyps_from(NewFacets, new_generator, Pyramid_key);
      return;
    }
    // large vs small, by the accumulated-comparisons estimate of the
    // original (with a bounds clamp the original leaves implicit)
    bool large = false;
    if (!Comparisons.empty()) {
      size_t idx = Pyramid_key.size() - dim;
      if (idx >= Comparisons.size())
        idx = Comparisons.size() - 1;
      large = (nmz_largePyramidFactor * Comparisons[idx] > old_nr_supp_hyps);
    }
    if (large) {
      LargeRecPyrs.push_back(hyp);
      return;
    }
    NmzKernel<Tint> Pyramid(*this, Pyramid_key);
    Pyramid.apex = new_generator;
    Pyramid.build_cone();
    nrTotalComparisons += Pyramid.nrTotalComparisons;
  }

  // process_pyramids (recursive variant): one pyramid per visible facet.
  void process_pyramids(size_t new_generator) {
    std::vector<size_t> Pyramid_key;
    Pyramid_key.reserve(nr_gen);
    auto hyp = Facets.begin();
    for (size_t kk = 0; kk < old_nr_supp_hyps; ++kk, ++hyp) {
      if (hyp->ValNewGen == 0) {
        hyp->GenInHyp.set(new_generator);
        hyp->simplicial = false;
        continue;
      }
      if (hyp->ValNewGen > 0)
        continue;
      Pyramid_key.clear();
      Pyramid_key.push_back(new_generator);
      for (size_t i = 0; i < nr_gen; i++)
        if (in_triang[i] && hyp->GenInHyp.test(i))
          Pyramid_key.push_back(i);
      process_pyramid(Pyramid_key, new_generator, *hyp);
    }
    evaluate_large_rec_pyramids(new_generator);
  }

  // build_cone: the incremental main loop.
  void build_cone() {
    size_t RecBoundSuppHyp =
        dim * nmz_SuppHypRecursionFactor * nmz_arith_cost_factor<Tint>();
    find_and_evaluate_start_simplex();
    for (size_t i = 0; i < nr_gen; ++i) {
      if (in_triang[i])
        continue;
      terminate_in_arithmetic_error<Tint>();
      old_nr_supp_hyps = Facets.size();
      bool is_new_generator = false;
      size_t nr_pos = 0, nr_neg = 0, nr_pos_simp = 0, nr_neg_simp = 0;
      for (auto &facet : Facets) {
        Tint scalar_product = v_scal(facet.Hyp, i);
        facet.ValNewGen = scalar_product;
        facet.negative = false;
        facet.positive = false;
        facet.neutral = false;
        if (scalar_product < 0) {
          is_new_generator = true;
          facet.negative = true;
          nr_neg++;
          if (facet.simplicial)
            nr_neg_simp++;
          continue;
        }
        if (scalar_product == 0) {
          facet.neutral = true;
          continue;
        }
        nr_pos++;
        facet.positive = true;
        if (facet.simplicial)
          nr_pos_simp++;
      }
      if (!is_new_generator)
        continue;

      if (nr_neg * nr_pos - nr_neg_simp * nr_pos_simp >= RecBoundSuppHyp) {
#ifdef DEBUG_NORMALIZ_DUAL_DESC
        os << "NMZ: gen=" << i + 1 << " pyramid decomposition nr_pos="
           << nr_pos << " nr_neg=" << nr_neg << "\n";
#endif
        process_pyramids(i);
      } else {
        find_new_facets(i);
      }

      // remove the negative (visible) facets
      {
        auto l = Facets.begin();
        for (size_t j = 0; j < old_nr_supp_hyps; j++) {
          if (l->negative)
            l = Facets.erase(l);
          else
            ++l;
        }
      }
      GensInCone.push_back(i);
      nrGensInCone++;
      Comparisons.push_back(nrTotalComparisons);
      in_triang[i] = 1;
#ifdef DEBUG_NORMALIZ_DUAL_DESC
      if (!is_pyramid)
        os << "NMZ: gen=" << i + 1 << ", " << Facets.size() << " hyp\n";
#endif
    }
    if (is_pyramid)
      Mother->select_supphyps_from(Facets, apex, Mother_Key);
  }
};

// The facet output type: primitive inward normal and full incidence over the
// original rows.
template <typename T> struct NormalizDDFacet {
  MyVector<T> normal;
  Face incd;
};

// Kernel driver on a ring matrix: sort the rows lexicographically (the
// insertion order the original uses for support-hyperplane computations),
// run build_cone, then translate the facets back to the original row order.
// The incidence bits over inserted generators are exact by construction; the
// rows never inserted (only possible when the input carries redundant rays)
// are completed by scalar products.
template <typename Tint>
std::vector<NormalizDDFacet<Tint>>
NormalizDualDesc_Kernel(MyMatrix<Tint> const &EXT, std::ostream &os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
#ifdef TIMINGS_NORMALIZ_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef SANITY_CHECK_NORMALIZ_DUAL_DESC
  if (nbCol < 2) {
    std::cerr << "NMZ: requires EXT.cols() >= 2, nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
#endif
  // lexicographic sort of the rows (perm_by_weights with no weights)
  std::vector<int> perm(nbRow);
  for (int i = 0; i < nbRow; i++)
    perm[i] = i;
  std::sort(perm.begin(), perm.end(), [&](int a, int b) -> bool {
    for (int k = 0; k < nbCol; k++) {
      if (EXT(a, k) < EXT(b, k))
        return true;
      if (EXT(b, k) < EXT(a, k))
        return false;
    }
    return a < b;
  });
  MyMatrix<Tint> EXTsort(nbRow, nbCol);
  for (int i = 0; i < nbRow; i++)
    for (int k = 0; k < nbCol; k++)
      EXTsort(i, k) = EXT(perm[i], k);

  NmzKernel<Tint> kernel(EXTsort, os);
  kernel.build_cone();
  terminate_in_arithmetic_error<Tint>();

  // rows never inserted: their incidences must be completed by hand
  std::vector<size_t> missing;
  for (size_t i = 0; i < kernel.nr_gen; i++)
    if (!kernel.in_triang[i])
      missing.push_back(i);

  std::vector<NormalizDDFacet<Tint>> result;
  result.reserve(kernel.Facets.size());
  for (auto &facet : kernel.Facets) {
    Face incd(nbRow);
    for (size_t i = 0; i < kernel.nr_gen; i++)
      if (kernel.in_triang[i] && facet.GenInHyp.test(i))
        incd.set(perm[i]);
    for (auto &i : missing)
      if (kernel.v_scal(facet.Hyp, i) == 0)
        incd.set(perm[i]);
    result.push_back({std::move(facet.Hyp), std::move(incd)});
  }
#ifdef TIMINGS_NORMALIZ_DUAL_DESC
  os << "NMZ: |EXT|=" << nbRow << "/" << nbCol
     << " |facets|=" << result.size() << " time=" << time << "\n";
#endif
  return result;
}

} // namespace normaliz_dd

// Runs the kernel and returns facets with T-typed normals. As in the
// beneath-and-beyond wrapper: the arithmetic is division-free, so for a field
// input the whole computation runs on the underlying integer ring, attempting
// machine integers (TryInt64) first and falling back to the exact ring on
// overflow. Facet incidences are combinatorial and type independent.
template <typename T>
std::vector<normaliz_dd::NormalizDDFacet<T>>
NormalizDualDesc_run(MyMatrix<T> const &EXT, std::ostream &os) {
  if constexpr (is_ring_field<T>::value) {
    using Tring = typename underlying_ring<T>::ring_type;
    int nbRow = EXT.rows();
    int nbCol = EXT.cols();
    MyMatrix<Tring> EXTring(nbRow, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      MyVector<T> eRow = NonUniqueScaleToIntegerVector(GetMatrixRow(EXT, iRow));
      AssignMatrixRow(EXTring, iRow, UniversalVectorConversion<Tring, T>(eRow));
    }
    if constexpr (use_try_int64<Tring>::value) {
      try {
        MyMatrix<TryInt64> EXTtry = ConvertMatrixToTryInt64<TryInt64>(EXTring);
        std::vector<normaliz_dd::NormalizDDFacet<TryInt64>> try_facets =
            normaliz_dd::NormalizDualDesc_Kernel(EXTtry, os);
        terminate_in_arithmetic_error<TryInt64>();
        std::vector<normaliz_dd::NormalizDDFacet<T>> result;
        result.reserve(try_facets.size());
        for (auto &f : try_facets) {
          int len = f.normal.size();
          MyVector<T> normal(len);
          for (int u = 0; u < len; u++)
            normal(u) = ConvertFromTryInt64<T>(f.normal(u));
          result.push_back({std::move(normal), std::move(f.incd)});
        }
        return result;
      } catch (TryIntException const &) {
      }
    }
    std::vector<normaliz_dd::NormalizDDFacet<Tring>> ring_facets =
        normaliz_dd::NormalizDualDesc_Kernel(EXTring, os);
    std::vector<normaliz_dd::NormalizDDFacet<T>> result;
    result.reserve(ring_facets.size());
    for (auto &f : ring_facets)
      result.push_back(
          {UniversalVectorConversion<T, Tring>(f.normal), std::move(f.incd)});
    return result;
  } else {
    return normaliz_dd::NormalizDualDesc_Kernel(EXT, os);
  }
}

// Facet incidences, matching DirectFacetComputationIncidence.
template <typename T>
vectface POLY_DualDescription_NormalizIncidence(MyMatrix<T> const &EXT,
                                                std::ostream &os) {
  std::vector<normaliz_dd::NormalizDDFacet<T>> facets =
      NormalizDualDesc_run(EXT, os);
  vectface vf(EXT.rows());
  for (auto &facet : facets)
    vf.push_back(facet.incd);
  return vf;
}

// Facet normals, matching DirectFacetComputationInequalities.
template <typename T>
MyMatrix<T> POLY_DualDescription_NormalizInequalities(MyMatrix<T> const &EXT,
                                                      std::ostream &os) {
  std::vector<normaliz_dd::NormalizDDFacet<T>> facets =
      NormalizDualDesc_run(EXT, os);
  int n_facet = facets.size();
  int nbCol = EXT.cols();
  MyMatrix<T> FAC(n_facet, nbCol);
  for (int i = 0; i < n_facet; i++)
    AssignMatrixRow(FAC, i, facets[i].normal);
  return FAC;
}

// (Face, inequality) pairs, matching DirectFacetComputationFaceIneq.
template <typename T, typename Fprocess>
void POLY_DualDescription_NormalizFaceIneq(MyMatrix<T> const &EXT,
                                           Fprocess f_process,
                                           std::ostream &os) {
  std::vector<normaliz_dd::NormalizDDFacet<T>> facets =
      NormalizDualDesc_run(EXT, os);
  for (auto &facet : facets) {
    std::pair<Face, MyVector<T>> pair{facet.incd, facet.normal};
    f_process(pair);
  }
}

// clang-format off
#endif  // SRC_POLY_POLY_DUALDESC_NORMALIZ_H_
// clang-format on
