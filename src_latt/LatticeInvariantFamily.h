// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEINVARIANTFAMILY_H_
#define SRC_LATT_LATTICEINVARIANTFAMILY_H_

// clang-format off
#include "LatticeMethods.h"
#include "LatticeRootSystem.h"
#include "LatticeInvariantBV.h"
#include "InvariantVectorFamily.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_INVARIANT_FAMILY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_INVARIANT_FAMILY
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_INVARIANT_FAMILY
#endif

/*
  The recognition of a lattice among a known family of pairwise
  non-isometric lattices, the classes of a genus for instance, as needed
  once the genus has been enumerated: the Hecke operators of the algebraic
  modular forms take the neighbors of the known classes, and each neighbor
  is one of them.

  The family is turned into a decision tree on isometry invariants, built
  from the family itself: at a node with several candidate classes, the
  first invariant of the list below whose values split them is used, and
  the tree goes on below each value. The invariants, from the cheapest:
  * the root system (the vectors of norm 1 and 2), as its Dynkin type;
  * the number of vectors of norm k, k = 1, 2, ..., as long as the vectors
    of norm at most k stay few (max_vectors);
  * the BV invariant (LatticeInvariantBV.h), where it is predicted cheaper
    than the canonical form (InvariantFamilyCostBV);
  * the norms of the iterated shortest vector family;
  * the canonical form (heuristic), which separates any two classes.
  So the tree only uses an expensive invariant where the cheap ones do not
  separate: for the Niemeier lattices the root system decides everything.
  The recognition of a lattice follows the tree, computing the invariants
  of the nodes on its path only, once each. A lattice isometric to a class
  of the family has its invariants, so it reaches the leaf of that class; a
  value that the tree does not know means that the lattice is not in the
  family.
 */

enum class InvariantFamilyStep {
  RootSystem,
  ThetaCount,
  BV,
  ShortestFamily,
  Canonical
};

inline std::string NameInvariantFamilyStep(InvariantFamilyStep step) {
  switch (step) {
  case InvariantFamilyStep::RootSystem:
    return "root_system";
  case InvariantFamilyStep::ThetaCount:
    return "theta_count";
  case InvariantFamilyStep::BV:
    return "bv";
  case InvariantFamilyStep::ShortestFamily:
    return "shortest_family";
  case InvariantFamilyStep::Canonical:
    return "canonical";
  }
  return "unknown";
}

// A node: the invariant used (with its norm for ThetaCount) and, for each
// of its values, the next node (index >= 0) or the class (-1 - index).
struct InvariantFamilyNode {
  InvariantFamilyStep step;
  int param;
  std::map<std::string, int> children;
};

template <typename Tint> struct InvariantFamily {
  std::vector<MyMatrix<Tint>> ListGram;
  // A second form for each lattice, or empty: the classes are then those of
  // the pairs (G, S) (a Hermitian lattice, S = G Omega^T), the invariants of
  // G being only filters and the canonical form that of the pair.
  std::vector<MyMatrix<Tint>> ListSecond;
  std::vector<InvariantFamilyNode> nodes;
  // the first node, or -1 - index for a family of one class
  int root;
  // the bound on the number of vectors (up to sign) of the theta counts
  size_t max_vectors;
};

/*
  The invariants of one lattice, computed on demand and kept: the vectors
  up to sign of norm at most the largest k asked so far, the others once.
 */
// The integral Gram matrix G in gmp, for the computations done there whatever
// Tint is (through int64_t for the other integer types, whose values fit).
template <typename Tint>
MyMatrix<mpq_class> InvariantFamilyGmpMatrix(MyMatrix<Tint> const &G) {
  if constexpr (std::is_same_v<Tint, mpz_class>) {
    return UniversalMatrixConversion<mpq_class, mpz_class>(G);
  } else {
    return UniversalMatrixConversion<mpq_class, int64_t>(
        UniversalMatrixConversion<int64_t, Tint>(G));
  }
}

// The LLL-reduced Gram matrix of G, computed in gmp.
template <typename Tint>
MyMatrix<Tint> InvariantFamilyReduced(MyMatrix<Tint> const &G,
                                      std::ostream &os) {
  MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(G);
  MyMatrix<mpq_class> G_red =
      LLLreducedBasis<mpq_class, mpz_class>(G_q, os).GramMatRed;
  if constexpr (std::is_same_v<Tint, mpz_class>) {
    return UniversalMatrixConversion<mpz_class, mpq_class>(G_red);
  } else {
    return UniversalMatrixConversion<Tint, int64_t>(
        UniversalMatrixConversion<int64_t, mpq_class>(G_red));
  }
}

// The LLL reduction of the Gram matrix G, in gmp, applied to the pair (G,
// S): P G P^T and P S P^T.
template <typename Tint>
std::pair<MyMatrix<Tint>, MyMatrix<Tint>>
InvariantFamilyReducedPair(MyMatrix<Tint> const &G, MyMatrix<Tint> const &S,
                           std::ostream &os) {
  MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(G);
  MyMatrix<mpq_class> S_q = InvariantFamilyGmpMatrix(S);
  LLLreduction<mpq_class, mpz_class> rec =
      LLLreducedBasis<mpq_class, mpz_class>(G_q, os);
  MyMatrix<mpq_class> P = UniversalMatrixConversion<mpq_class, mpz_class>(rec.Pmat);
  MyMatrix<mpq_class> S_red = P * S_q * P.transpose();
  auto back = [&](MyMatrix<mpq_class> const &M) -> MyMatrix<Tint> {
    if constexpr (std::is_same_v<Tint, mpz_class>) {
      return UniversalMatrixConversion<mpz_class, mpq_class>(M);
    } else {
      return UniversalMatrixConversion<Tint, int64_t>(
          UniversalMatrixConversion<int64_t, mpq_class>(M));
    }
  };
  return {back(rec.GramMatRed), back(S_red)};
}

template <typename Tint> struct InvariantFamilyCache {
  using T = typename overlying_field<Tint>::field_type;
  MyMatrix<Tint> G_int;
  MyMatrix<T> G;
  std::optional<std::string> root_system;
  int theta_bound = 0;
  // the number of vectors up to sign of each norm at most theta_bound
  std::map<T, size_t> theta;
  // set when the vectors of norm at most theta_over_bound exceed the bound
  bool theta_over = false;
  int theta_over_bound = 0;
  std::optional<std::string> bv;
  std::optional<std::string> shortest_family;
  std::optional<std::string> canonical;
  // The lattice is first LLL-reduced in gmp: the invariants do not depend
  // on the basis, and a skewed basis (a neighbor, for instance) would make
  // the entries of the computations overflow a fixed size Tint.
  InvariantFamilyCache(MyMatrix<Tint> const &G_in, std::ostream &os)
      : G_int(InvariantFamilyReduced(G_in, os)),
        G(UniversalMatrixConversion<T, Tint>(G_int)) {}
  // the second form, on the same reduced basis (empty without one)
  MyMatrix<Tint> S_int;
  InvariantFamilyCache(MyMatrix<Tint> const &G_in, MyMatrix<Tint> const &S_in,
                       std::ostream &os) {
    std::pair<MyMatrix<Tint>, MyMatrix<Tint>> pr =
        InvariantFamilyReducedPair(G_in, S_in, os);
    G_int = pr.first;
    S_int = pr.second;
    G = UniversalMatrixConversion<T, Tint>(G_int);
  }
};

/*
  The number of vectors (up to sign) of norm at most k, or nothing if it
  exceeds max_vectors. The counts by norm come from the enumeration (no
  vector is kept, no norm recomputed), and when they must go further the
  bound is doubled, so that the enumerations of k = 1, 2, 3, ... are not
  all redone. If the doubled bound has too many vectors, the largest bound
  within max_vectors is found by bisection and kept, with the next one as
  the bound known to exceed it: every later k is then answered without an
  enumeration.
 */
template <typename Tint>
std::optional<size_t> InvariantFamilyCountAtMost(InvariantFamilyCache<Tint> &cache,
                                                 int k, size_t max_vectors,
                                                 std::ostream &os) {
  using T = typename InvariantFamilyCache<Tint>::T;
  if (k > cache.theta_bound) {
    if (cache.theta_over && cache.theta_over_bound <= k) {
      return {};
    }
    // G is reduced already (in gmp): a second reduction, in Tint, could
    // overflow a fixed size Tint.
    CVPSolver<T, Tint> solver =
        GetCVPSolver<T, Tint>(cache.G, GramReduction::AlreadyReduced, os);
    auto counts_at = [&](int bound) -> std::optional<std::map<T, size_t>> {
      return solver.norm_counts_at_most(T(bound), max_vectors);
    };
    int k_try = std::max(k, 2 * cache.theta_bound);
    std::optional<std::map<T, size_t>> counts = counts_at(k_try);
    if (counts) {
      cache.theta = *counts;
      cache.theta_bound = k_try;
    } else {
      // within the bound at lo (or lo = theta_bound), beyond it at hi
      int lo = cache.theta_bound;
      int hi = k_try;
      if (k < k_try) {
        std::optional<std::map<T, size_t>> counts_k = counts_at(k);
        if (counts_k) {
          lo = k;
          cache.theta = *counts_k;
        } else {
          hi = k;
        }
      }
      while (hi - lo > 1 && lo >= k) {
        int mid = lo + (hi - lo) / 2;
        std::optional<std::map<T, size_t>> counts_mid = counts_at(mid);
        if (counts_mid) {
          lo = mid;
          cache.theta = *counts_mid;
        } else {
          hi = mid;
        }
      }
      cache.theta_bound = lo;
      cache.theta_over = true;
      cache.theta_over_bound = hi;
      if (lo < k) {
        return {};
      }
    }
  }
  size_t count = 0;
  for (auto const &[norm, c] : cache.theta) {
    if (norm <= T(k)) {
      count += c;
    }
  }
  return count;
}

/*
  The number of vectors (up to sign) of the BV invariant of the lattice: those
  of norm at most M, M the least integer with at least twice the rank of
  them (LatticeBVOptions defaults, see BVShortVectors), from the theta counts
  of the cache. Nothing if they exceed max_vectors.
 */
template <typename Tint>
std::optional<size_t> InvariantFamilyBVSize(InvariantFamilyCache<Tint> &cache,
                                            size_t max_vectors,
                                            std::ostream &os) {
  size_t min_vectors = 2 * static_cast<size_t>(cache.G.rows());
  for (int k = 1;; k++) {
    std::optional<size_t> count =
        InvariantFamilyCountAtMost(cache, k, max_vectors, os);
    if (!count) {
      return {};
    }
    if (*count >= min_vectors) {
      return count;
    }
  }
}

/*
  The predicted costs in seconds of the BV invariant on N vectors and of the
  canonical form (heuristic) in rank n, fitted on 414 lattices of ranks 3 to
  24 (the largest DefLatDB genus of each rank, c02216, the Niemeier
  lattices) with Bench_InvariantCost: t_BV = 8.9e-7 N^1.65 and
  t_can = 7.9e-6 n^2.83. Choosing the predicted cheaper one costs 3.9 s on
  that sample, against 3.7 s for the best choice, 8.1 s for the canonical
  form always and 248 s for BV always. The canonical form is cheaper beyond
  N = 195 vectors in rank 10, 586 in rank 19, 875 in rank 24. The sample has
  none of the lattices of very large automorphism groups whose canonical
  form is much more expensive than the model.
 */
inline double InvariantFamilyCostBV(size_t N) {
  return 8.9e-7 * std::pow(static_cast<double>(N), 1.65);
}

inline double InvariantFamilyCostCanonical(int n) {
  return 7.9e-6 * std::pow(static_cast<double>(n), 2.83);
}

/*
  The value of one invariant on a lattice, as a string; nothing for a
  ThetaCount beyond the bound on the number of vectors.
 */
template <typename Tint, typename Tgroup>
std::optional<std::string>
InvariantFamilyValue(InvariantFamilyCache<Tint> &cache,
                     InvariantFamilyStep step, int param, size_t max_vectors,
                     std::ostream &os) {
  using T = typename InvariantFamilyCache<Tint>::T;
  switch (step) {
  case InvariantFamilyStep::RootSystem: {
    if (!cache.root_system) {
      RootSystemData<Tint> rs =
          ComputeRootSystem<T, Tint>(cache.G, GramReduction::AlreadyReduced, os);
      std::vector<std::pair<char, int>> comps;
      for (auto const &c : rs.components) {
        comps.push_back({c.letter, c.rank});
      }
      std::sort(comps.begin(), comps.end());
      std::ostringstream str;
      str << rs.n_roots;
      for (auto const &[letter, rank] : comps) {
        str << " " << letter << rank;
      }
      cache.root_system = str.str();
    }
    return cache.root_system;
  }
  case InvariantFamilyStep::ThetaCount: {
    if (!InvariantFamilyCountAtMost(cache, param, max_vectors, os)) {
      return {};
    }
    auto iter = cache.theta.find(T(param));
    size_t c = (iter == cache.theta.end()) ? 0 : iter->second;
    return std::to_string(c);
  }
  case InvariantFamilyStep::BV: {
    if (!cache.bv) {
      // too many vectors: the lattice is not one that reached a BV node
      if (!InvariantFamilyBVSize(cache, max_vectors, os)) {
        return {};
      }
      LatticeBVOptions opt;
      std::ostringstream str;
      str << LatticeBV<T, Tint>(cache.G, opt, os);
      cache.bv = str.str();
    }
    return cache.bv;
  }
  case InvariantFamilyStep::ShortestFamily: {
    if (!cache.shortest_family) {
      // In gmp whatever Tint is: the recursion rescales projected Gram
      // matrices whose entries can overflow 64 bits.
      MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(cache.G_int);
      MyMatrix<mpz_class> SHV =
          IteratedShortestVectorFamilyHalf<mpq_class, mpz_class>(G_q, os);
      std::map<mpq_class, int> norms;
      for (int i = 0; i < SHV.rows(); i++) {
        MyVector<mpz_class> v = GetMatrixRow(SHV, i);
        norms[EvaluationQuadForm<mpq_class, mpz_class>(G_q, v)]++;
      }
      std::ostringstream str;
      for (auto const &[norm, mult] : norms) {
        str << norm << ":" << mult << " ";
      }
      cache.shortest_family = str.str();
    }
    return cache.shortest_family;
  }
  case InvariantFamilyStep::Canonical: {
    if (!cache.canonical) {
      // In gmp whatever Tint is: the canonical forms are not available for
      // every arithmetic, and the fallback is rarely needed.
      MyMatrix<mpq_class> G_q = InvariantFamilyGmpMatrix(cache.G_int);
      std::vector<MyMatrix<mpq_class>> l_mat{G_q};
      if (cache.S_int.rows() > 0) {
        l_mat.push_back(InvariantFamilyGmpMatrix(cache.S_int));
      }
      MyMatrix<mpz_class> B =
          LatticeCanonicalBasisMultiple<mpq_class, mpz_class, Tgroup>(
              l_mat, LatticeStabCanIsomMethod::Heuristic, os);
      MyMatrix<mpq_class> B_T = UniversalMatrixConversion<mpq_class, mpz_class>(B);
      std::ostringstream str;
      for (auto const &M : l_mat) {
        MyMatrix<mpq_class> Can = B_T * M * B_T.transpose();
        for (int i = 0; i < Can.rows(); i++) {
          for (int j = 0; j < Can.cols(); j++) {
            str << Can(i, j) << " ";
          }
        }
        str << "| ";
      }
      cache.canonical = str.str();
    }
    return cache.canonical;
  }
  }
  return {};
}

/*
  The decision tree of the family ListGram, pairwise non-isometric
  lattices. max_vectors bounds the number of vectors up to sign enumerated
  for the theta counts, on every lattice of the family.
 */
template <typename Tint, typename Tgroup>
InvariantFamily<Tint>
get_invariant_family(std::vector<MyMatrix<Tint>> const &ListGram,
                     std::ostream &os, size_t max_vectors = 20000,
                     std::vector<MyMatrix<Tint>> const &ListSecond = {}) {
  InvariantFamily<Tint> family;
  family.ListGram = ListGram;
  family.ListSecond = ListSecond;
  family.max_vectors = max_vectors;
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
  // per kind of work: the evaluations and their total time in microseconds
  std::map<std::string, std::pair<size_t, int64_t>> timings;
  MicrosecondTime time_total;
  auto print_timings = [&](std::string const &when) -> void {
    os << "|INVFAM: " << when << " nodes=" << family.nodes.size()
       << " elapsed_s=" << time_total.const_eval_int64() / 1000000.0;
    for (auto const &[name, pair] : timings) {
      os << " " << name << "=" << pair.first << "/" << pair.second / 1000000.0
         << "s";
    }
    os << "\n";
  };
  MicrosecondTime time_reduce;
#endif
  std::vector<InvariantFamilyCache<Tint>> caches;
  for (auto const &G : ListGram) {
    if (ListSecond.empty()) {
      caches.emplace_back(G, os);
    } else {
      caches.emplace_back(G, ListSecond[caches.size()], os);
    }
  }
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
  timings["reduce"] = {ListGram.size(), time_reduce.eval_int64()};
#endif
  // The values of one step on the candidates, if it splits them.
  auto try_split = [&](std::vector<size_t> const &cands,
                       InvariantFamilyStep step, int param)
      -> std::optional<std::map<std::string, std::vector<size_t>>> {
    std::map<std::string, std::vector<size_t>> parts;
    for (size_t i : cands) {
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
      MicrosecondTime time;
#endif
      std::optional<std::string> val = InvariantFamilyValue<Tint, Tgroup>(
          caches[i], step, param, max_vectors, os);
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
      std::pair<size_t, int64_t> &entry =
          timings[NameInvariantFamilyStep(step)];
      entry.first++;
      entry.second += time.eval_int64();
#endif
      if (!val) {
        return {};
      }
      parts[*val].push_back(i);
    }
    if (parts.size() < 2) {
      return {};
    }
    return parts;
  };
  // The node of a set of candidates, built recursively: its index, or
  // -1 - index for a single class.
  std::function<int(std::vector<size_t> const &)> build =
      [&](std::vector<size_t> const &cands) -> int {
    if (cands.size() == 1) {
      return -1 - static_cast<int>(cands[0]);
    }
    std::optional<std::map<std::string, std::vector<size_t>>> parts;
    InvariantFamilyStep step = InvariantFamilyStep::RootSystem;
    int param = 0;
    parts = try_split(cands, step, param);
    // the theta counts, norm by norm, while the vectors stay few
    for (int k = 1; !parts; k++) {
      bool within = true;
      for (size_t i : cands) {
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
        MicrosecondTime time;
#endif
        bool ok = InvariantFamilyCountAtMost(caches[i], k, max_vectors, os)
                      .has_value();
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
        std::pair<size_t, int64_t> &entry = timings["theta_enumeration"];
        entry.first++;
        entry.second += time.eval_int64();
#endif
        if (!ok) {
          within = false;
          break;
        }
      }
      if (!within) {
        break;
      }
      step = InvariantFamilyStep::ThetaCount;
      param = k;
      parts = try_split(cands, step, param);
    }
    // BV only where it is predicted cheaper than the canonical form, on
    // the largest BV vector set of the candidates
    bool use_bv = !parts;
    if (use_bv) {
      size_t N_max = 0;
      for (size_t i : cands) {
        std::optional<size_t> N = InvariantFamilyBVSize(caches[i], max_vectors, os);
        if (!N) {
          use_bv = false;
          break;
        }
        N_max = std::max(N_max, *N);
      }
      int n = static_cast<int>(caches[cands[0]].G.rows());
      if (use_bv &&
          InvariantFamilyCostBV(N_max) >= InvariantFamilyCostCanonical(n)) {
        use_bv = false;
      }
    }
    for (InvariantFamilyStep s :
         {InvariantFamilyStep::BV, InvariantFamilyStep::ShortestFamily,
          InvariantFamilyStep::Canonical}) {
      if (parts) {
        break;
      }
      if (s == InvariantFamilyStep::BV && !use_bv) {
        continue;
      }
      step = s;
      param = 0;
      parts = try_split(cands, step, param);
    }
    if (!parts) {
      std::cerr << "INVFAM: " << cands.size()
                << " lattices of the family have the same canonical form, "
                   "they are not pairwise non-isometric\n";
      throw TerminalException{1};
    }
    int idx = static_cast<int>(family.nodes.size());
    family.nodes.push_back({step, param, {}});
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
    if (family.nodes.size() % 500 == 0) {
      print_timings("progress");
    }
#endif
    for (auto const &[val, sub] : *parts) {
      int child = build(sub);
      family.nodes[idx].children[val] = child;
    }
    return idx;
  };
  std::vector<size_t> all(ListGram.size());
  for (size_t i = 0; i < all.size(); i++) {
    all[i] = i;
  }
  if (all.empty()) {
    std::cerr << "INVFAM: an empty family\n";
    throw TerminalException{1};
  }
  family.root = build(all);
#ifdef SANITY_CHECK_LATTICE_INVARIANT_FAMILY
  // every class of the family reaches its own leaf
  for (size_t i = 0; i < ListGram.size(); i++) {
    int node = family.root;
    while (node >= 0) {
      InvariantFamilyNode const &nd = family.nodes[node];
      std::optional<std::string> val = InvariantFamilyValue<Tint, Tgroup>(
          caches[i], nd.step, nd.param, max_vectors, os);
      auto iter = val ? nd.children.find(*val) : nd.children.end();
      if (iter == nd.children.end()) {
        std::cerr << "INVFAM: the class " << i << " of the family has no "
                  << "branch at a node of the tree\n";
        throw TerminalException{1};
      }
      node = iter->second;
    }
    if (static_cast<size_t>(-1 - node) != i) {
      std::cerr << "INVFAM: the class " << i << " of the family reaches the "
                << "leaf of the class " << (-1 - node) << "\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef TIMINGS_LATTICE_INVARIANT_FAMILY
  print_timings("done");
#endif
#ifdef DEBUG_LATTICE_INVARIANT_FAMILY
  std::map<std::string, size_t> count_step;
  for (auto const &node : family.nodes) {
    count_step[NameInvariantFamilyStep(node.step)]++;
  }
  os << "INVFAM: " << ListGram.size() << " lattices, " << family.nodes.size()
     << " nodes:";
  for (auto const &[name, c] : count_step) {
    os << " " << name << "=" << c;
  }
  os << "\n";
#endif
  return family;
}

/*
  The index in the family of the class of the lattice G, by the decision
  tree. An error if G is not isometric to a lattice of the family.
 */
template <typename Tint, typename Tgroup>
size_t family_invariant_index(InvariantFamily<Tint> const &family,
                              MyMatrix<Tint> const &G, std::ostream &os,
                              MyMatrix<Tint> const &S = MyMatrix<Tint>(0, 0)) {
  InvariantFamilyCache<Tint> cache =
      (S.rows() > 0) ? InvariantFamilyCache<Tint>(G, S, os)
                     : InvariantFamilyCache<Tint>(G, os);
  int node = family.root;
  while (node >= 0) {
    InvariantFamilyNode const &nd = family.nodes[node];
    std::optional<std::string> val = InvariantFamilyValue<Tint, Tgroup>(
        cache, nd.step, nd.param, family.max_vectors, os);
    if (!val) {
      // a theta count beyond the bound: G has more short vectors than every
      // lattice that reached this node
      std::cerr << "INVFAM: the lattice is not in the family (too many "
                   "short vectors)\n";
      throw TerminalException{1};
    }
    auto iter = nd.children.find(*val);
    if (iter == nd.children.end()) {
      std::cerr << "INVFAM: the lattice is not in the family ("
                << NameInvariantFamilyStep(nd.step) << " " << *val
                << " unknown)\n";
      throw TerminalException{1};
    }
    node = iter->second;
  }
  return static_cast<size_t>(-1 - node);
}

// clang-format off
#endif  // SRC_LATT_LATTICEINVARIANTFAMILY_H_
// clang-format on
