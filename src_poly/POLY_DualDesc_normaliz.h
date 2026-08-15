// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESC_NORMALIZ_H_
#define SRC_POLY_POLY_DUALDESC_NORMALIZ_H_

// clang-format off
#include "POLY_DualDesc_beneath_and_beyond.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_MatrixInverse.h"
#include <algorithm>
#include <list>
#include <map>
#include <optional>
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

// The machine-integer type used for the fast attempt over a field input, as
// in the lrs backend: it defaults to the exact-detection TryCarryInt64;
// define POLY_NORMALIZ_TRY_SIMD to use the conservative, SIMD-friendly
// TrySimdInt64 instead (for A/B testing the two deferred-overflow flavours).
#ifdef POLY_NORMALIZ_TRY_SIMD
using NmzTryInt = TrySimdInt64;
#else
using NmzTryInt = TryCarryInt64;
#endif

#ifdef NMZ_COEFF_STATS
// Maximum absolute facet coefficient seen across the whole run (all kernel
// levels), for calibrating the bound-discipline arithmetic. Reset by the
// top-level driver.
inline int64_t nmz_stat_max_abs_coeff = 0;
#endif

// The recursion bound for pyramid decomposition: pyramids are built when
// nr_neg*nr_pos - nr_neg_simp*nr_pos_simp >= dim * SuppHypRecursionFactor
// (times arith_cost_factor for non-machine arithmetic). Constants from
// full_cone.cpp.
constexpr size_t nmz_SuppHypRecursionFactor = 320000;
constexpr size_t nmz_largePyramidFactor = 20;
// Factor by which multiprecision arithmetic is assumed slower than machine
// integers (GMP_time_factor in full_cone.cpp). Machine types get 1 -- and so
// do the deferred-overflow TryInt types, which are machine arithmetic with a
// cheap check: classifying them as heavy would flip the rank-test /
// comparison-test decisions and the pyramid recursion bound away from what
// the reference uses for long long, which measurably loses on cones with a
// large intermediate facet count (the comparison test scales with the facet
// count, the rank test does not).
template <typename Tint>
constexpr size_t nmz_arith_cost_factor() {
  if constexpr (std::is_fundamental<Tint>::value || is_try_int<Tint>::value)
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

  std::vector<uint8_t> in_triang;
  std::vector<size_t> GensInCone;
  size_t nrGensInCone;
  // InsertedUpTo[s] = the set of generators inserted in the first s steps
  // (GensInCone[0..s)), used for the word-parallel extension test: a facet
  // pair born at steps (a, b) has gained a common generator since the birth
  // of the younger iff its common-generator set is NOT contained in
  // InsertedUpTo[max(a, b)].
  std::vector<Face> InsertedUpTo;
  std::vector<NmzFacet<Tint>> Facets;
  size_t old_nr_supp_hyps;
  size_t HypCounter;
  std::vector<size_t> Comparisons;
  size_t nrTotalComparisons;
  std::vector<NmzFacet<Tint>> LargeRecPyrs;

  // top cone
  NmzKernel(MyMatrix<Tint> const &_TopGen, std::ostream &_os)
      : TopGen(_TopGen), os(_os),
        dim(_TopGen.cols()), nr_gen(_TopGen.rows()), is_pyramid(false),
        in_triang(nr_gen, 0), nrGensInCone(0),
        old_nr_supp_hyps(0), HypCounter(1), nrTotalComparisons(0) {
    RowIdx.resize(nr_gen);
    for (size_t i = 0; i < nr_gen; i++)
      RowIdx[i] = static_cast<int>(i);
  }

  // pyramid: generators are C's gens selected by Key (apex first). The
  // pyramid holds no back-reference to its creator: the give-back of its
  // support hyperplanes is done by the caller (process_pyramid) once
  // build_cone has finished.
  NmzKernel(NmzKernel const &C, std::vector<size_t> const &Key)
      : TopGen(C.TopGen), os(C.os), dim(C.dim),
        nr_gen(Key.size()), is_pyramid(true),
        in_triang(nr_gen, 0), nrGensInCone(0), old_nr_supp_hyps(0),
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
#ifdef NMZ_COEFF_STATS
    // Measurement hook for the bound-discipline design: every persisted
    // facet passes through here at creation, right after gcd reduction, so
    // this maximum is exactly the H of the safety conditions
    // dim*G*H < 2^63 (scalar product) and 2*(dim*G*H)*H < 2^63 (FM combo).
    if constexpr (uses_deferred_overflow<Tint>::value ||
                  std::is_fundamental<Tint>::value) {
      for (size_t k = 0; k < dim; k++) {
        int64_t v;
        if constexpr (uses_deferred_overflow<Tint>::value)
          v = hyp.Hyp(k).get_const_val();
        else
          v = static_cast<int64_t>(hyp.Hyp(k));
        if (v < 0)
          v = -v;
        if (v > nmz_stat_max_abs_coeff)
          nmz_stat_max_abs_coeff = v;
      }
    }
#endif
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
                      std::vector<NmzFacet<Tint>> &NewHyps,
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
    InsertedUpTo.push_back(Face(nr_gen));
    for (size_t j = 0; j < dim; j++) {
      in_triang[key[j]] = 1;
      GensInCone.push_back(key[j]);
      Face nxt = InsertedUpTo.back();
      nxt.set(key[j]);
      InsertedUpTo.push_back(std::move(nxt));
    }
    nrGensInCone = dim;
    nrTotalComparisons = dim * dim / 2;
    if constexpr (nmz_arith_cost_factor<Tint>() != 1)
      nrTotalComparisons *= (nmz_arith_cost_factor<Tint>() / 4);
  }

  // workspace of the Bareiss rank test, sized on first use and reused
  // across calls (the test runs once per surviving candidate pair)
  MyMatrix<Tint> rank_ws;

  // Fraction-free (Bareiss) rank test: do the top rows of the generators
  // selected by the face have rank >= target? Exact divisions by the
  // previous pivot, no gcd reductions, so the cost matches the reference's
  // long-long row echelon (the gcd-reducing SelectIndependentRows is an
  // order of magnitude more expensive here and made the rank test a worse
  // deal than the comparison test it replaces). For a deferred-overflow
  // type an overflow OF THIS ELIMINATION is absorbed and reported as
  // nullopt, so the caller falls back to the comparison test for that one
  // candidate instead of aborting the whole run; a flag raised by earlier
  // arithmetic is left untouched for the outer terminate check.
  std::optional<bool> rank_face_at_least(Face const &f, size_t target) {
    size_t nb_row = f.count();
    if (nb_row < target)
      return false;
    if (rank_ws.rows() == 0)
      rank_ws = MyMatrix<Tint>(nr_gen, dim);
    size_t pos = 0;
    boost::dynamic_bitset<>::size_type i_bit = f.find_first();
    while (i_bit != boost::dynamic_bitset<>::npos) {
      int row = RowIdx[i_bit];
      for (size_t u = 0; u < dim; u++)
        rank_ws(pos, u) = TopGen(row, u);
      pos++;
      i_bit = f.find_next(i_bit);
    }
    bool pre_correct = true;
    if constexpr (uses_deferred_overflow<Tint>::value)
      pre_correct = is_correct;
    Tint prev(1);
    size_t rank = 0;
    size_t col = 0;
    while (rank < target && col < dim) {
      // pivot search in the current column, below the found pivots
      size_t piv = rank;
      while (piv < nb_row && rank_ws(piv, col) == 0)
        piv++;
      if (piv == nb_row) {
        col++;
        continue;
      }
      if (piv != rank)
        for (size_t u = col; u < dim; u++)
          std::swap(rank_ws(rank, u), rank_ws(piv, u));
      Tint pivot = rank_ws(rank, col);
      for (size_t i = rank + 1; i < nb_row; i++) {
        Tint coef = rank_ws(i, col);
        for (size_t u = col + 1; u < dim; u++) {
          Tint val = pivot * rank_ws(i, u) - coef * rank_ws(rank, u);
          Tint quot = val / prev; // exact division (Bareiss)
#ifdef SANITY_CHECK_NORMALIZ_DUAL_DESC
          if constexpr (!uses_deferred_overflow<Tint>::value) {
            if (quot * prev != val) {
              std::cerr << "NMZ: non-exact division in the Bareiss rank\n";
              throw TerminalException{1};
            }
          }
#endif
          rank_ws(i, u) = quot;
        }
        rank_ws(i, col) = Tint(0);
      }
      prev = pivot;
      rank++;
      col++;
    }
    if constexpr (uses_deferred_overflow<Tint>::value) {
      if (pre_correct && !is_correct) {
        // the overflow came from this elimination only: absorb it
        is_correct = true;
        return {};
      }
      if (!is_correct)
        // earlier arithmetic already overflowed; the result is unusable and
        // the outer terminate check will abort the run
        return {};
    }
    return rank == target;
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

    std::vector<NmzFacet<Tint> *> Pos_Simp, Pos_Non_Simp;
    std::vector<NmzFacet<Tint> *> Neg_Simp, Neg_Non_Simp;
    std::vector<NmzFacet<Tint> *> Neutral_Simp, Neutral_Non_Simp;

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
    std::vector<std::pair<Face, int>> Neg_Subfacet_Multi;
    for (i = 0; i < nr_NegSimp; i++) {
      Face RelGen_NegHyp = Gen_BothSides & Neg_Simp[i]->GenInHyp;
      size_t nr_RelGen_NegHyp = 0;
      for (j = 0; j < nr_gen; j++) {
        if (RelGen_NegHyp.test(j))
          nr_RelGen_NegHyp++;
        if (nr_RelGen_NegHyp > subfacet_dim)
          break;
      }
      if (nr_RelGen_NegHyp == facet_dim) {
        for (k = 0; k < nr_gen; k++) {
          if (RelGen_NegHyp.test(k)) {
            Face subfacet = RelGen_NegHyp;
            subfacet.reset(k);
            Neg_Subfacet_Multi.emplace_back(std::move(subfacet),
                                            static_cast<int>(i));
          }
        }
      } else if (nr_RelGen_NegHyp == subfacet_dim) {
        Neg_Subfacet_Multi.emplace_back(std::move(RelGen_NegHyp),
                                        static_cast<int>(i));
      }
    }
    std::sort(Neg_Subfacet_Multi.begin(), Neg_Subfacet_Multi.end());
    // a subfacet shared by two negative simplicial facets is interior to the
    // negative side: mark both copies (-2), pair by pair as the original's
    // list erasure does
    {
      size_t p = 0;
      while (p < Neg_Subfacet_Multi.size()) {
        if (p + 1 < Neg_Subfacet_Multi.size() &&
            Neg_Subfacet_Multi[p + 1].first == Neg_Subfacet_Multi[p].first) {
          Neg_Subfacet_Multi[p].second = -2;
          Neg_Subfacet_Multi[p + 1].second = -2;
          p += 2;
        } else {
          p += 1;
        }
      }
    }
    size_t nr_NegSubfMult = 0;
    for (auto &entry : Neg_Subfacet_Multi)
      if (entry.second != -2)
        nr_NegSubfMult++;

    // remove those that lie in a neutral facet or a negative non-simplicial
    // facet (the size guards against a quadratic disaster are the original's)
    if (nr_NegSubfMult * (nr_NeuSimp + nr_NeuNonSimp + nr_NegNonSimp) <=
        100000000) {
      for (auto &entry : Neg_Subfacet_Multi) {
        if (entry.second == -2)
          continue;
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
        if (entry.second >= 0)
          last_inserted = Neg_Subfacet.insert(last_inserted, std::move(entry));
    }
    size_t nr_NegSubf = Neg_Subfacet.size();
    Neg_Subfacet_Multi.clear();

    std::vector<std::vector<NmzFacet<Tint>>> NewHypsSimp(nr_PosSimp);
    std::vector<std::vector<NmzFacet<Tint>>> NewHypsNonSimp(nr_PosNonSimp);

    nrTotalComparisons += nr_NegNonSimp * nr_PosNonSimp;

    //
    // Positive simplicial vs negative simplicial and non-simplicial
    //
    // The per-pair generator tests below run word-parallel on the incidence
    // bitsets (in-place AND + popcount + subset test against InsertedUpTo)
    // instead of the original's one-generator-at-a-time key loops: same
    // accept/reject decisions, a fraction of the work, and no per-pair
    // allocation.
    Face RelGen_PosHyp(nr_gen);
    Face CommonGens(nr_gen);
    for (i = 0; i < nr_PosSimp; i++) {
      RelGen_PosHyp = Gen_BothSides;
      RelGen_PosHyp &= Pos_Simp[i]->GenInHyp;
      size_t nr_RelGen_PosHyp = RelGen_PosHyp.count();
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

      // reject on the first missing generator when only one subfacet fits
      // (nr_RelGen == subfacet_dim), on the second one otherwise
      size_t allowed_missing = (nr_RelGen_PosHyp == facet_dim) ? 1 : 0;
      for (j = 0; j < nr_NegNonSimp; j++) {
        CommonGens = RelGen_PosHyp;
        CommonGens &= Neg_Non_Simp[j]->GenInHyp;
        if (nr_RelGen_PosHyp - CommonGens.count() <= allowed_missing) {
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

      size_t missing_bound;

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
        RelGen_PosHyp = Gen_BothSides;
        RelGen_PosHyp &= PosHyp_Pointer->GenInHyp;
        size_t nr_RelGen_PosHyp = RelGen_PosHyp.count();
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

          CommonGens = RelGen_PosHyp;
          CommonGens &= NegHyp_Pointer->GenInHyp;
          size_t nr_CommonGens = CommonGens.count();
          if (nr_RelGen_PosHyp - nr_CommonGens > missing_bound)
            continue;

          // Two facets that are not mother and daughter cannot have shared a
          // subfacet at the birth of the younger one: they intersect in a
          // subfacet now only if a common generator arrived afterwards.
          if (extension_test) {
            size_t both_born =
                std::max(PosHyp_Pointer->BornAt, NegHyp_Pointer->BornAt);
            if (CommonGens.is_subset_of(InsertedUpTo[both_born]))
              continue;
          }

          if (subfacet_dim <= 2) {
            add_hyperplane(new_generator, *PosHyp_Pointer, *NegHyp_Pointer,
                           NewHypsNonSimp[i], false);
            continue;
          }

          // rank test vs comparison test, by the a-priori cost estimate; a
          // rank test that overflowed sends the candidate to the comparison
          // test instead
          bool common_subfacet = true;
          bool ranktest =
              (nr_NonSimp > nmz_arith_cost_factor<Tint>() * dim * dim *
                                nr_CommonGens / 3);
          if (ranktest) {
            std::optional<bool> opt =
                rank_face_at_least(CommonGens, subfacet_dim);
            if (opt)
              common_subfacet = *opt;
            else
              ranktest = false;
          }
          if (!ranktest) {
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
      for (auto &hyp : NewHypsSimp[i])
        Facets.emplace_back(std::move(hyp));
    for (i = 0; i < nr_PosNonSimp; i++)
      for (auto &hyp : NewHypsNonSimp[i])
        Facets.emplace_back(std::move(hyp));
  }

  // select_supphyps_from: the mother cone (=this) selects the new global
  // support hyperplanes from the facet list of a daughter pyramid. A pyramid
  // facet containing the apex is a new global facet iff it is strictly
  // positive on every inserted generator outside the pyramid -- a value 0
  // would mean the same hyperplane is delivered by another pyramid, so this
  // one scalar-product filter is both the validity check and the
  // deduplication.
  void select_supphyps_from(std::vector<NmzFacet<Tint>> &NewFacets,
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
  // New hyperplanes go into the caller-owned NewHyps accumulator, NOT into
  // Facets: the PosHyps pointers into Facets must stay valid across all the
  // match calls of one evaluate_large_rec_pyramids round.
  void match_neg_hyp_with_pos_hyps(NmzFacet<Tint> const &Neg,
                                   size_t new_generator,
                                   std::vector<NmzFacet<Tint> *> const &PosHyps,
                                   Face const &GenIn_PosHyp,
                                   std::list<Face> &Facets_0_1,
                                   std::vector<NmzFacet<Tint>> &NewHyps) {
    size_t subfacet_dim = dim - 2;

    Face RelGens_InNegHyp = Neg.GenInHyp & GenIn_PosHyp;
    size_t nr_RelGens_InNegHyp = RelGens_InNegHyp.count();
    if (nr_RelGens_InNegHyp < subfacet_dim)
      return;
    size_t missing_bound = nr_RelGens_InNegHyp - subfacet_dim;

    // Word-parallel pair tests, as in find_new_facets.
    Face CommonGens(nr_gen);
    for (auto const &Pos : PosHyps) {
      if (Neg.Ident == Pos->Mother || Pos->Ident == Neg.Mother) {
        add_hyperplane(new_generator, *Pos, Neg, NewHyps, false);
        continue;
      }
      bool extension_test =
          Neg.BornAt == Pos->BornAt ||
          (Neg.BornAt < Pos->BornAt && Pos->Mother != 0) ||
          (Pos->BornAt < Neg.BornAt && Neg.Mother != 0);

      CommonGens = RelGens_InNegHyp;
      CommonGens &= Pos->GenInHyp;
      size_t nr_common_gens = CommonGens.count();
      if (nr_RelGens_InNegHyp - nr_common_gens > missing_bound)
        continue;
      if (extension_test) {
        size_t both_born = std::max(Neg.BornAt, Pos->BornAt);
        if (CommonGens.is_subset_of(InsertedUpTo[both_born]))
          continue;
      }

      bool common_subfacet = true;
      if (!Pos->simplicial) {
        bool ranktest = true;
        if constexpr (nmz_arith_cost_factor<Tint>() != 1)
          ranktest = (old_nr_supp_hyps > nmz_arith_cost_factor<Tint>() * dim *
                                             dim * nr_common_gens / 3);
        if (ranktest) {
          std::optional<bool> opt =
              rank_face_at_least(CommonGens, subfacet_dim);
          if (opt)
            common_subfacet = *opt;
          else
            ranktest = false;
        }
        if (!ranktest) {
          for (auto hp_t = Facets_0_1.begin(); hp_t != Facets_0_1.end();
               ++hp_t) {
            if (CommonGens.is_subset_of(*hp_t) && (*hp_t != Neg.GenInHyp) &&
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
  }

  void evaluate_large_rec_pyramids(size_t new_generator) {
    size_t nrLargeRecPyrs = LargeRecPyrs.size();
    if (nrLargeRecPyrs == 0)
      return;
    // a std::list: the comparison test moves successful reducers to the front
    std::list<Face> Facets_0_1;
    for (size_t i = 0; i < old_nr_supp_hyps; ++i) {
      if (Facets[i].simplicial)
        continue;
      Facets_0_1.push_back(Facets[i].GenInHyp);
    }
    std::vector<NmzFacet<Tint> *> PosHyps;
    Face GenIn_PosHyp(nr_gen);
    for (size_t i = 0; i < old_nr_supp_hyps; ++i)
      if (Facets[i].ValNewGen > 0) {
        GenIn_PosHyp |= Facets[i].GenInHyp;
        PosHyps.push_back(&Facets[i]);
      }
    nrTotalComparisons += PosHyps.size() * nrLargeRecPyrs;
#ifdef DEBUG_NORMALIZ_DUAL_DESC
    os << "NMZ: large pyramids " << nrLargeRecPyrs << "\n";
#endif
    std::vector<NmzFacet<Tint>> NewHyps;
    for (auto &Neg : LargeRecPyrs)
      match_neg_hyp_with_pos_hyps(Neg, new_generator, PosHyps, GenIn_PosHyp,
                                  Facets_0_1, NewHyps);
    for (auto &hyp : NewHyps)
      Facets.emplace_back(std::move(hyp));
    LargeRecPyrs.clear();
  }

  // process_pyramid: a simplicial pyramid is solved on the spot; a small one
  // recurses into a child cone whose facets are handed back through
  // select_supphyps_from; a large one is deferred to
  // evaluate_large_rec_pyramids. The give-back appends to Facets, so the hyp
  // reference (a Facets element) is only read before any such append: the
  // large-pyramid copy is taken first, and the other branches do not use it.
  void process_pyramid(std::vector<size_t> const &Pyramid_key,
                       size_t new_generator, NmzFacet<Tint> const &hyp) {
    if (Pyramid_key.size() == dim) {
      // simplicial pyramid: facet j of the pyramid contains all its
      // generators but the j-th
      std::vector<NmzFacet<Tint>> NewFacets;
      NewFacets.reserve(dim);
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
    Pyramid.build_cone();
    select_supphyps_from(Pyramid.Facets, new_generator, Pyramid_key);
    nrTotalComparisons += Pyramid.nrTotalComparisons;
  }

  // process_pyramids (recursive variant): one pyramid per visible facet.
  void process_pyramids(size_t new_generator) {
    std::vector<size_t> Pyramid_key;
    Pyramid_key.reserve(nr_gen);
    for (size_t kk = 0; kk < old_nr_supp_hyps; ++kk) {
      // process_pyramid may append to Facets (small-pyramid give-back), so
      // the reference is re-taken each iteration and not used afterwards
      NmzFacet<Tint> &hyp = Facets[kk];
      if (hyp.ValNewGen == 0) {
        hyp.GenInHyp.set(new_generator);
        hyp.simplicial = false;
        continue;
      }
      if (hyp.ValNewGen > 0)
        continue;
      Pyramid_key.clear();
      Pyramid_key.push_back(new_generator);
      for (size_t i = 0; i < nr_gen; i++)
        if (in_triang[i] && hyp.GenInHyp.test(i))
          Pyramid_key.push_back(i);
      process_pyramid(Pyramid_key, new_generator, hyp);
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

      // remove the negative (visible) facets; the facets appended during
      // this step are never negative, so the erasure can scan everything
      std::erase_if(Facets,
                    [](NmzFacet<Tint> const &f) -> bool { return f.negative; });
      GensInCone.push_back(i);
      nrGensInCone++;
      {
        Face nxt = InsertedUpTo.back();
        nxt.set(i);
        InsertedUpTo.push_back(std::move(nxt));
      }
      Comparisons.push_back(nrTotalComparisons);
      in_triang[i] = 1;
#ifdef DEBUG_NORMALIZ_DUAL_DESC
      if (!is_pyramid)
        os << "NMZ: gen=" << i + 1 << ", " << Facets.size() << " hyp\n";
#endif
    }
  }
};

// Kernel driver on a ring matrix: sort the rows lexicographically (the
// insertion order the original uses for support-hyperplane computations),
// run build_cone, then hand each facet to the caller as
// f_facet(MyVector<Tint>&& normal, Face&& incd) with the incidence over the
// ORIGINAL row order -- no facet vector is materialized here. The incidence
// bits over inserted generators are exact by construction; the rows never
// inserted (only possible when the input carries redundant rays) are
// completed by scalar products.
template <typename Tint, typename Ffacet>
void NormalizDualDesc_Kernel_f(MyMatrix<Tint> const &EXT, std::ostream &os,
                               Ffacet f_facet) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
#ifdef NMZ_COEFF_STATS
  nmz_stat_max_abs_coeff = 0;
#endif
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

#ifdef TIMINGS_NORMALIZ_DUAL_DESC
  size_t n_facet = kernel.Facets.size();
#endif
  for (auto &facet : kernel.Facets) {
    Face incd(nbRow);
    for (size_t i = 0; i < kernel.nr_gen; i++)
      if (kernel.in_triang[i] && facet.GenInHyp.test(i))
        incd.set(perm[i]);
    for (auto &i : missing)
      if (kernel.v_scal(facet.Hyp, i) == 0)
        incd.set(perm[i]);
    f_facet(std::move(facet.Hyp), std::move(incd));
  }
#ifdef TIMINGS_NORMALIZ_DUAL_DESC
  os << "NMZ: |EXT|=" << nbRow << "/" << nbCol
     << " |facets|=" << n_facet << " time=" << time << "\n";
#endif
#ifdef NMZ_COEFF_STATS
  os << "NMZ: max_abs_coeff=" << nmz_stat_max_abs_coeff << "\n";
#endif
}

} // namespace normaliz_dd

// Convert a kernel-arithmetic normal back to T: identity for T itself,
// unpacking for the TryInt64 machine integers, the universal conversion for
// the underlying ring.
template <typename T, typename Tw>
MyVector<T> NormalizNormalConvert(MyVector<Tw> &&normal) {
  if constexpr (std::is_same_v<T, Tw>) {
    return std::move(normal);
  } else if constexpr (is_try_int<Tw>::value) {
    int len = normal.size();
    MyVector<T> ret(len);
    for (int u = 0; u < len; u++)
      ret(u) = ConvertFromTryInt64<T>(normal(u));
    return ret;
  } else {
    return UniversalVectorConversion<T, Tw>(normal);
  }
}

// Attempts the kernel on machine integers (TryInt64). On success the
// buffered facets are replayed to f_facet and true is returned; on
// overflow nothing has been emitted and false is returned, so the caller
// falls back to the exact ring without partial results.
template <typename Tring, typename Ffacet>
bool NormalizDualDesc_try_int64(MyMatrix<Tring> const &EXTring,
                                std::ostream &os, Ffacet f_facet) {
  try {
    MyMatrix<normaliz_dd::NmzTryInt> EXTtry = ConvertMatrixToTryInt64<normaliz_dd::NmzTryInt>(EXTring);
    std::vector<std::pair<MyVector<normaliz_dd::NmzTryInt>, Face>> buffer;
    normaliz_dd::NormalizDualDesc_Kernel_f(
        EXTtry, os, [&buffer](MyVector<normaliz_dd::NmzTryInt> &&normal, Face &&incd) {
          buffer.emplace_back(std::move(normal), std::move(incd));
        });
    terminate_in_arithmetic_error<normaliz_dd::NmzTryInt>();
    for (auto &entry : buffer)
      f_facet(std::move(entry.first), std::move(entry.second));
    return true;
  } catch (TryIntException const &) {
    return false;
  }
}

// Runs the kernel and hands every facet to the caller as
// f_facet(MyVector<Tw>&& normal, Face&& incd), where Tw is the arithmetic
// the kernel ran on (T itself, its underlying integer ring, or TryInt64):
// the callers convert the normal with NormalizNormalConvert -- or ignore it,
// so the incidence-only entry point never converts a normal at all. As in
// the beneath-and-beyond wrapper the arithmetic is division-free, so for a
// field input the whole computation runs on the underlying integer ring.
// Both the field case (after scaling to the ring) and the plain ring case
// attempt machine integers (TryInt64) first and fall back to the exact
// ring on overflow.
template <typename T, typename Ffacet>
void NormalizDualDesc_process(MyMatrix<T> const &EXT, std::ostream &os,
                              Ffacet f_facet) {
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
      if (NormalizDualDesc_try_int64(EXTring, os, f_facet))
        return;
    }
    normaliz_dd::NormalizDualDesc_Kernel_f(EXTring, os, f_facet);
  } else {
    if constexpr (use_try_int64<T>::value) {
      if (NormalizDualDesc_try_int64(EXT, os, f_facet))
        return;
    }
    normaliz_dd::NormalizDualDesc_Kernel_f(EXT, os, f_facet);
  }
}

// Facet incidences, matching DirectFacetComputationIncidence. The normals
// are dropped without ever being converted.
template <typename T>
vectface POLY_DualDescription_NormalizIncidence(MyMatrix<T> const &EXT,
                                                std::ostream &os) {
  vectface vf(EXT.rows());
  NormalizDualDesc_process(
      EXT, os, [&vf]([[maybe_unused]] auto &&normal, Face &&incd) -> void {
        vf.push_back(incd);
      });
  return vf;
}

// Facet normals, matching DirectFacetComputationInequalities.
template <typename T>
MyMatrix<T> POLY_DualDescription_NormalizInequalities(MyMatrix<T> const &EXT,
                                                      std::ostream &os) {
  std::vector<MyVector<T>> ListVect;
  NormalizDualDesc_process(
      EXT, os,
      [&ListVect](auto &&normal, [[maybe_unused]] Face &&incd) -> void {
        ListVect.push_back(NormalizNormalConvert<T>(std::move(normal)));
      });
  return MatrixFromVectorFamily(ListVect);
}

// (Face, inequality) pairs, matching DirectFacetComputationFaceIneq: each
// facet is streamed straight into f_process.
template <typename T, typename Fprocess>
void POLY_DualDescription_NormalizFaceIneq(MyMatrix<T> const &EXT,
                                           Fprocess f_process,
                                           std::ostream &os) {
  NormalizDualDesc_process(
      EXT, os, [&f_process](auto &&normal, Face &&incd) -> void {
        std::pair<Face, MyVector<T>> pair{
            std::move(incd), NormalizNormalConvert<T>(std::move(normal))};
        f_process(pair);
      });
}

// clang-format off
#endif  // SRC_POLY_POLY_DUALDESC_NORMALIZ_H_
// clang-format on
