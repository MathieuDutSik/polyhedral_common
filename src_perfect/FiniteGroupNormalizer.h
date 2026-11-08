// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_FINITEGROUPNORMALIZER_H_
#define SRC_PERFECT_FINITEGROUPNORMALIZER_H_

// clang-format off
#include "perfect_tspace.h"
#include "Tspace_ListMatTransform.h"
#include "LatticeStabEquiCan.h"
#include "OnlineExhaustiveReduction.h"
#include <algorithm>
#include <functional>
#include <string>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_FINITE_GROUP_NORMALIZER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_FINITE_GROUP_NORMALIZER
#endif

#ifdef TIMINGS
#define TIMINGS_FINITE_GROUP_NORMALIZER
#endif

/*
  The normalizer in GL_n(Z) of a finite group G of GL_n(Z) and the
  Z-equivalence (conjugacy in GL_n(Z)) of two finite groups, by the method
  of Opgenorth, Plesken and Schulz (Crystallographic algorithms and
  tables, Acta Cryst. A54 (1998) 517-531, Remark 5 and Sections 3.3.2,
  3.3.4 and 4).

  Convention: a matrix g acts on the forms by M -> g M g^T, the space of
  invariant forms is F(G) = {M : g M g^T = M for g in G}, and conjugation
  is H -> x H x^{-1}. If x G1 x^{-1} = G2 then x F(G1) x^T = F(G2).

  B = B(G), the Bravais group, is the pointwise stabilizer of F(G) in
  GL_n(Z); it is finite, contains G and F(B) = F(G). Then:
  - N(B) is the global stabilizer of the T-space F(G), and Voronoi's
    algorithm on the G-perfect forms of F(G) gives generators of it from
    the stabilizers of the perfect forms and the equivalences between
    them (Opgenorth's algorithm, done by DataPerfectTspaceFunc);
  - N(G) is contained in N(B) (it preserves F(G)), and it is the
    stabilizer of G in the conjugation action of N(B) on the subgroups of
    the finite group B: the orbit is finite and the Schreier generators
    give N(G).
  - G1 and G2 are Z-equivalent iff some x maps B1 to B2 by conjugation and
    x G1 x^{-1} is conjugate to G2 under N(B2). If x exists, x P1 x^T is a
    perfect form of F(G2) for a perfect form P1 of F(G1); so x is found by
    comparing P1 with the perfect forms of F(G2) up to N(B2): x0 P1 x0^T =
    P2 gives all the isometries x0 O(P1), and among them the ones mapping
    B1 to B2 are found by the conjugation orbit of B1 in the finite group
    O(P1). The last step is the orbit of x G1 x^{-1} under N(B2).

  The finite groups are handled by the list of their elements, so the
  method is meant for the groups of moderate order of small dimension.
 */

struct FiniteGroupNormalizerBounds {
  size_t max_group = 2000000;
  size_t max_orbit = 200000;
};

template <typename Tint>
std::vector<MyMatrix<Tint>>
FiniteMatrixGroupElements(std::vector<MyMatrix<Tint>> const &gens, int n,
                          size_t max_size) {
  std::vector<MyMatrix<Tint>> elts{IdentityMat<Tint>(n)};
  std::unordered_set<MyMatrix<Tint>> seen{elts[0]};
  for (size_t pos = 0; pos < elts.size(); pos++) {
    for (auto const &g : gens) {
      MyMatrix<Tint> x = elts[pos] * g;
      if (seen.insert(x).second) {
        elts.push_back(x);
        if (elts.size() > max_size) {
          std::cerr << "FGN: the finite group has more than " << max_size
                    << " elements\n";
          throw TerminalException{1};
        }
      }
    }
  }
  return elts;
}

/*
  A subgroup is identified by the sorted hashes of its elements. The hash
  must see the signs: std::hash<mpz_class> is mpz_get_ui, the absolute
  value, so that diag(-1, 1) and diag(1, -1) would collide.
 */
template <typename Tint> size_t SignedMatrixHash(MyMatrix<Tint> const &x) {
  std::hash<Tint> hasher;
  size_t seed = 0x1b873540;
  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      Tint const &v = x(i, j);
      size_t h = hasher(v) * 2 + ((v < 0) ? 1 : 0);
      seed ^= h + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
  }
  return seed;
}

template <typename Tint>
std::vector<size_t> SubgroupKey(std::vector<MyMatrix<Tint>> const &elts) {
  std::vector<size_t> key;
  key.reserve(elts.size());
  for (auto const &x : elts) {
    key.push_back(SignedMatrixHash(x));
  }
  std::sort(key.begin(), key.end());
  return key;
}

template <typename Tint>
std::vector<MyMatrix<Tint>>
ConjugateElements(MyMatrix<Tint> const &x, MyMatrix<Tint> const &x_inv,
                  std::vector<MyMatrix<Tint>> const &elts) {
  std::vector<MyMatrix<Tint>> res;
  res.reserve(elts.size());
  for (auto const &e : elts) {
    res.push_back(x * e * x_inv);
  }
  return res;
}

template <typename Tint> MyMatrix<Tint> IntegralInverse(MyMatrix<Tint> const &x) {
  using T = typename overlying_field<Tint>::field_type;
  MyMatrix<T> x_T = UniversalMatrixConversion<T, Tint>(x);
  MyMatrix<T> inv_T = Inverse(x_T);
  if (!IsIntegralMatrix(inv_T)) {
    std::cerr << "FGN: the matrix is not in GL_n(Z)\n";
    throw TerminalException{1};
  }
  return UniversalMatrixConversion<Tint, T>(inv_T);
}

/*
  The orbit of the finite subgroup H (by its elements) under conjugation by
  the group generated by gens (the orbit is assumed finite, so the forward
  closure is the orbit). transversal[i] maps H to the i-th conjugate. If
  target is given, the search stops when the target is met and its
  transversal element is returned in found. Otherwise the Schreier
  generators of the stabilizer of H are computed.
 */
template <typename Tint> struct SubgroupOrbitResult {
  std::vector<MyMatrix<Tint>> transversal;
  std::vector<MyMatrix<Tint>> stabilizer_gens;
  std::optional<MyMatrix<Tint>> found;
};

template <typename Tint>
SubgroupOrbitResult<Tint>
OrbitSubgroupConjugation(std::vector<MyMatrix<Tint>> const &H, int n,
                         std::vector<MyMatrix<Tint>> const &gens,
                         std::optional<std::vector<size_t>> const &target,
                         size_t max_orbit) {
  SubgroupOrbitResult<Tint> res;
  std::vector<MyMatrix<Tint>> gens_inv;
  for (auto const &g : gens) {
    gens_inv.push_back(IntegralInverse(g));
  }
  std::unordered_map<std::vector<size_t>, size_t,
                     std::function<size_t(std::vector<size_t> const &)>>
      map(16, [](std::vector<size_t> const &v) -> size_t {
        size_t seed = v.size();
        for (size_t x : v) {
          seed ^= x + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
      });
  std::vector<size_t> key0 = SubgroupKey(H);
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  res.transversal.push_back(Id);
  map[key0] = 0;
  if (target && *target == key0) {
    res.found = Id;
    return res;
  }
  std::unordered_set<MyMatrix<Tint>> stab_set;
  for (size_t pos = 0; pos < res.transversal.size(); pos++) {
    MyMatrix<Tint> t = res.transversal[pos];
    MyMatrix<Tint> t_inv = IntegralInverse(t);
    std::vector<MyMatrix<Tint>> Hi = ConjugateElements(t, t_inv, H);
    for (size_t u = 0; u < gens.size(); u++) {
      std::vector<MyMatrix<Tint>> Himg =
          ConjugateElements(gens[u], gens_inv[u], Hi);
      std::vector<size_t> key = SubgroupKey(Himg);
      MyMatrix<Tint> st = gens[u] * t;
      auto iter = map.find(key);
      if (iter == map.end()) {
        size_t idx = res.transversal.size();
        map[key] = idx;
        res.transversal.push_back(st);
        if (target && *target == key) {
          res.found = st;
          return res;
        }
        if (res.transversal.size() > max_orbit) {
          std::cerr << "FGN: the conjugation orbit has more than " << max_orbit
                    << " elements\n";
          throw TerminalException{1};
        }
      } else if (!target) {
        // t_j^{-1} s t_i stabilizes H
        MyMatrix<Tint> sch = IntegralInverse(res.transversal[iter->second]) * st;
        if (sch != Id && stab_set.insert(sch).second) {
          res.stabilizer_gens.push_back(sch);
        }
      }
    }
  }
  return res;
}

template <typename Tint>
bool NormalizesFiniteGroup(MyMatrix<Tint> const &x,
                           std::vector<MyMatrix<Tint>> const &gens_G,
                           std::unordered_set<MyMatrix<Tint>> const &set_G) {
  MyMatrix<Tint> x_inv = IntegralInverse(x);
  for (auto const &g : gens_G) {
    if (set_G.count(x * g * x_inv) == 0) {
      return false;
    }
  }
  return true;
}

template <typename Tint>
std::vector<MyMatrix<Tint>> ReduceGenerators(std::vector<MyMatrix<Tint>> const &gens,
                                             int n, std::ostream &os) {
  OnlineHierarchicalMatrixReduction<Tint> red(n, os);
  for (auto const &g : gens) {
    red.insert_generator(g);
  }
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  std::vector<MyMatrix<Tint>> res;
  for (auto const &g : red.get_current_matrix_t()) {
    if (g != Id) {
      res.push_back(g);
    }
  }
  return res;
}

/*
  The T-space F(G) of a finite group, its Bravais group as the pointwise
  stabilizer, the G-perfect forms up to N(B) and generators of N(B).
 */
template <typename T, typename Tint> struct BravaisNormalizerData {
  LinSpaceMatrix<T> LinSpa;
  std::vector<MyMatrix<Tint>> bravais_gens;
  std::vector<MyMatrix<Tint>> normalizer_gens;
  std::vector<MyMatrix<T>> perfect_forms;
};

/*
  The T-space F(G) of a finite group, with its Bravais group as the
  pointwise stabilizer. The positive definite element is the average
  sum_{g in G} g g^T of the identity, which is invariant and positive
  definite: the general search of BuildLinSpaceMatrix for a positive
  element of a T-space hardly converges on groups with large entries.
 */
template <typename T, typename Tint, typename Tgroup>
LinSpaceMatrix<T> InvariantFormSpace(std::vector<MyMatrix<Tint>> const &gens,
                                     int n, std::ostream &os) {
  std::vector<MyMatrix<T>> gens_T =
      UniversalStdVectorMatrixConversion<T, Tint>(gens);
  std::vector<MyMatrix<T>> ListMat = BasisInvariantForm(n, gens_T, os);
  if (ListMat.empty()) {
    std::cerr << "FGN: the group has no invariant form, it is not finite\n";
    throw TerminalException{1};
  }
  FiniteGroupNormalizerBounds bounds;
  std::vector<MyMatrix<Tint>> elts =
      FiniteMatrixGroupElements(gens, n, bounds.max_group);
  MyMatrix<T> SuperMat = ZeroMatrix<T>(n, n);
  for (auto const &g : elts) {
    MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
    SuperMat += g_T * g_T.transpose();
  }
  LinSpaceMatrix<T> LinSpa;
  LinSpa.ListMat = ListMat;
  reset_paperwork(LinSpa);
  LinSpa.SuperMat = SuperMat;
  reset_pt_stab_gens<T, Tint, Tgroup>(LinSpa, os);
  LinSpa.isBravais =
      IsBravaisSpace(LinSpa.n, LinSpa.ListMat, LinSpa.PtStabGens, os);
  return LinSpa;
}

template <typename T>
bool IsSymmetricMatrixInSpace(LinSpaceMatrix<T> const &LinSpa,
                              MyMatrix<T> const &M) {
  // the same vectorization as ListMatAsBigMat
  MyVector<T> V = SymmetricMatrixToVector(M);
  return SolutionMat(LinSpa.ListMatAsBigMat, V).has_value();
}

template <typename T, typename Tint, typename Tgroup>
BravaisNormalizerData<T, Tint>
ComputeBravaisNormalizer(LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
#ifdef TIMINGS_FINITE_GROUP_NORMALIZER
  MicrosecondTime time;
#endif
  using TintGroup = typename Tgroup::Tint;
  BravaisNormalizerData<T, Tint> res;
  res.LinSpa = LinSpa;
  int n = LinSpa.n;
  res.bravais_gens = UniversalStdVectorMatrixConversion<Tint, T>(LinSpa.PtStabGens);
  int dimEXT = LinSpa.ListMat.size();
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  bool keep_generators = false;
  bool reduce_gram_matrix = false;
  DataPerfectTspace<T, Tint, Tgroup> data{
      LinSpa, OnlineHierarchicalMatrixReduction<Tint>(n, os), keep_generators,
      reduce_gram_matrix, std::move(rddo)};
  using Tdata = DataPerfectTspaceFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::optional<std::vector<Tout>> opt_l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  std::vector<Tout> l_tot =
      unfold_opt(opt_l_tot, "EnumerateAndStore_Serial (Bravais normalizer)");
  for (auto const &entry : l_tot) {
    res.perfect_forms.push_back(entry.x.Gram);
  }
  // Opgenorth: N(B) is generated by the stabilizers of the perfect forms
  // and the equivalences of their neighbours to the representatives. They
  // are taken raw, each oriented on its own, before any product is formed.
  std::vector<MyMatrix<Tint>> gens;
  for (auto const &entry : l_tot) {
    for (auto const &g : entry.x.GRP_matr) {
      gens.push_back(g);
    }
    for (auto const &adj : entry.ListAdj) {
      gens.push_back(adj.x.eBigMat);
    }
  }
  for (auto const &g : res.bravais_gens) {
    gens.push_back(g);
  }
  // The convention of the generators is fixed here: g acts by M -> g M
  // g^T, and a generator acting by M -> g^T M g is transposed.
  auto preserves = [&](MyMatrix<Tint> const &g) -> bool {
    MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
    for (auto const &M : LinSpa.ListMat) {
      MyMatrix<T> img = g_T * M * g_T.transpose();
      if (!IsSymmetricMatrixInSpace(LinSpa, img)) {
        return false;
      }
    }
    return true;
  };
  std::vector<MyMatrix<Tint>> gens_oriented;
  for (auto const &g : gens) {
    if (preserves(g)) {
      gens_oriented.push_back(g);
    } else if (preserves(g.transpose())) {
      gens_oriented.push_back(g.transpose());
    } else {
      std::cerr << "FGN: a generator does not preserve the T-space\n";
      throw TerminalException{1};
    }
  }
  res.normalizer_gens = ReduceGenerators(gens_oriented, n, os);
#ifdef TIMINGS_FINITE_GROUP_NORMALIZER
  os << "|FGN: ComputeBravaisNormalizer|=" << time << "\n";
#endif
#ifdef DEBUG_FINITE_GROUP_NORMALIZER
  os << "FGN: dim F=" << dimEXT << " |perfect forms|=" << res.perfect_forms.size()
     << " |gens N(B)|=" << res.normalizer_gens.size() << "\n";
#endif
  return res;
}

/*
  N(G) as the stabilizer of G in the conjugation action of N(B) on the
  subgroups of B: the orbit length [N(B) : N(G)] and the generators.
 */
template <typename Tint>
std::pair<size_t, std::vector<MyMatrix<Tint>>>
NormalizerFromBravaisNormalizer(std::vector<MyMatrix<Tint>> const &elts_G,
                                std::vector<MyMatrix<Tint>> const &gens_G, int n,
                                std::vector<MyMatrix<Tint>> const &gens_NB,
                                FiniteGroupNormalizerBounds const &bounds,
                                std::ostream &os) {
  SubgroupOrbitResult<Tint> orb =
      OrbitSubgroupConjugation(elts_G, n, gens_NB, {}, bounds.max_orbit);
  std::vector<MyMatrix<Tint>> sch = orb.stabilizer_gens;
  for (auto const &g : gens_G) {
    sch.push_back(g);
  }
  return {orb.transversal.size(), ReduceGenerators(sch, n, os)};
}

template <typename Tint> struct FiniteGroupNormalizerResult {
  int n;
  size_t order_G;
  size_t order_B;
  int dim_space;
  size_t n_perfect;
  size_t index_in_NB; // [N(B) : N(G)], the orbit length of G
  std::vector<MyMatrix<Tint>> gens;
};

template <typename T, typename Tint, typename Tgroup>
FiniteGroupNormalizerResult<Tint>
FiniteGroupNormalizer(std::vector<MyMatrix<Tint>> const &gens_G,
                      FiniteGroupNormalizerBounds const &bounds,
                      std::ostream &os) {
  if (gens_G.empty()) {
    std::cerr << "FGN: the group needs at least one generator\n";
    throw TerminalException{1};
  }
  int n = gens_G[0].rows();
  FiniteGroupNormalizerResult<Tint> res;
  res.n = n;
  std::vector<MyMatrix<Tint>> elts_G =
      FiniteMatrixGroupElements(gens_G, n, bounds.max_group);
  res.order_G = elts_G.size();
  LinSpaceMatrix<T> LinSpa = InvariantFormSpace<T, Tint, Tgroup>(gens_G, n, os);
  res.dim_space = LinSpa.ListMat.size();
  BravaisNormalizerData<T, Tint> bnd =
      ComputeBravaisNormalizer<T, Tint, Tgroup>(LinSpa, os);
  res.n_perfect = bnd.perfect_forms.size();
  std::vector<MyMatrix<Tint>> elts_B =
      FiniteMatrixGroupElements(bnd.bravais_gens, n, bounds.max_group);
  res.order_B = elts_B.size();
  if (res.order_B == res.order_G) {
    // G is its own Bravais group
    res.index_in_NB = 1;
    res.gens = bnd.normalizer_gens;
  } else {
    std::pair<size_t, std::vector<MyMatrix<Tint>>> pair =
        NormalizerFromBravaisNormalizer(elts_G, gens_G, n, bnd.normalizer_gens,
                                        bounds, os);
    res.index_in_NB = pair.first;
    res.gens = pair.second;
  }
  std::unordered_set<MyMatrix<Tint>> set_G(elts_G.begin(), elts_G.end());
  for (auto const &x : res.gens) {
    if (!NormalizesFiniteGroup(x, gens_G, set_G)) {
      std::cerr << "FGN: a computed generator does not normalize G, x=\n";
      WriteMatrix(std::cerr, x);
      throw TerminalException{1};
    }
  }
  return res;
}

/*
  The isometries of forms, in the convention x A x^T = B: the equivalence
  code returns some matrix, and the convention is fixed here.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> OrientIsometry(MyMatrix<Tint> const &M,
                                             MyMatrix<T> const &A,
                                             MyMatrix<T> const &B) {
  std::vector<MyMatrix<Tint>> cands{M, M.transpose(), IntegralInverse(M),
                                    IntegralInverse(M).transpose()};
  for (auto const &c : cands) {
    MyMatrix<T> c_T = UniversalMatrixConversion<T, Tint>(c);
    if (c_T * A * c_T.transpose() == B) {
      return c;
    }
  }
  return {};
}

template <typename Tint> struct FiniteGroupConjugacyResult {
  std::optional<MyMatrix<Tint>> conjugator;
  std::string reason;
  // the generators of N(B2) when they were computed, for reuse
  std::vector<MyMatrix<Tint>> normalizer_B2_gens;
};

template <typename T, typename Tint, typename Tgroup>
FiniteGroupConjugacyResult<Tint>
FiniteGroupConjugacy(std::vector<MyMatrix<Tint>> const &gens_G1,
                     std::vector<MyMatrix<Tint>> const &gens_G2,
                     FiniteGroupNormalizerBounds const &bounds,
                     std::ostream &os) {
  FiniteGroupConjugacyResult<Tint> res;
  if (gens_G1.empty() || gens_G2.empty()) {
    std::cerr << "FGN: the groups need at least one generator\n";
    throw TerminalException{1};
  }
  int n = gens_G1[0].rows();
  if (gens_G2[0].rows() != n) {
    res.reason = "different dimensions";
    return res;
  }
  std::vector<MyMatrix<Tint>> elts_G1 =
      FiniteMatrixGroupElements(gens_G1, n, bounds.max_group);
  std::vector<MyMatrix<Tint>> elts_G2 =
      FiniteMatrixGroupElements(gens_G2, n, bounds.max_group);
  if (elts_G1.size() != elts_G2.size()) {
    res.reason = "different orders";
    return res;
  }
  LinSpaceMatrix<T> LinSpa1 = InvariantFormSpace<T, Tint, Tgroup>(gens_G1, n, os);
  LinSpaceMatrix<T> LinSpa2 = InvariantFormSpace<T, Tint, Tgroup>(gens_G2, n, os);
  if (LinSpa1.ListMat.size() != LinSpa2.ListMat.size()) {
    res.reason = "different dimensions of the spaces of invariant forms";
    return res;
  }
  std::vector<MyMatrix<Tint>> gens_B1 =
      UniversalStdVectorMatrixConversion<Tint, T>(LinSpa1.PtStabGens);
  std::vector<MyMatrix<Tint>> gens_B2 =
      UniversalStdVectorMatrixConversion<Tint, T>(LinSpa2.PtStabGens);
  std::vector<MyMatrix<Tint>> elts_B1 =
      FiniteMatrixGroupElements(gens_B1, n, bounds.max_group);
  std::vector<MyMatrix<Tint>> elts_B2 =
      FiniteMatrixGroupElements(gens_B2, n, bounds.max_group);
  if (elts_B1.size() != elts_B2.size()) {
    res.reason = "Bravais groups of different orders";
    return res;
  }
  // Step 1: a matrix mapping B1 to B2
  BravaisNormalizerData<T, Tint> bnd2 =
      ComputeBravaisNormalizer<T, Tint, Tgroup>(LinSpa2, os);
  res.normalizer_B2_gens = bnd2.normalizer_gens;
  MyMatrix<T> P1 = GetOnePerfectForm<T, Tint>(LinSpa1, os).first;
  std::vector<size_t> key_B2 = SubgroupKey(elts_B2);
  std::optional<MyMatrix<Tint>> opt_X1;
  for (auto const &P2 : bnd2.perfect_forms) {
    std::optional<MyMatrix<Tint>> opt_iso =
        ArithmeticEquivalence<T, Tint, Tgroup>(P1, P2, os);
    if (!opt_iso) {
      continue;
    }
    std::optional<MyMatrix<Tint>> opt_X0 = OrientIsometry(*opt_iso, P1, P2);
    if (!opt_X0) {
      std::cerr << "FGN: the isometry of the perfect forms is not oriented\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> X0 = *opt_X0;
    MyMatrix<Tint> X0_inv = IntegralInverse(X0);
    // B1' = X0^{-1} B2 X0, a subgroup of O(P1) as B1 is
    std::vector<MyMatrix<Tint>> elts_B1p = ConjugateElements(X0_inv, X0, elts_B2);
    IntAutomorphismGroup<Tint, Tgroup> aut =
        ArithmeticAutomorphismGroup<T, Tint, Tgroup>(P1, os);
    std::vector<MyMatrix<Tint>> gens_O;
    for (auto const &g : aut.ListGen) {
      std::optional<MyMatrix<Tint>> opt_g = OrientIsometry(g, P1, P1);
      if (!opt_g) {
        std::cerr << "FGN: an automorphism of the perfect form is wrong\n";
        throw TerminalException{1};
      }
      gens_O.push_back(*opt_g);
    }
    SubgroupOrbitResult<Tint> orb = OrbitSubgroupConjugation(
        elts_B1, n, gens_O, SubgroupKey(elts_B1p), bounds.max_orbit);
    if (orb.found) {
      opt_X1 = X0 * *orb.found;
      break;
    }
  }
  if (!opt_X1) {
    res.reason = "the Bravais groups are not Z-equivalent";
    return res;
  }
  MyMatrix<Tint> X1 = *opt_X1;
  MyMatrix<Tint> X1_inv = IntegralInverse(X1);
#ifdef SANITY_CHECK_FINITE_GROUP_NORMALIZER
  if (SubgroupKey(ConjugateElements(X1, X1_inv, elts_B1)) != key_B2) {
    std::cerr << "FGN: X1 does not map B1 to B2\n";
    throw TerminalException{1};
  }
#endif
  // Step 2: X1 G1 X1^{-1} and G2 under N(B2)
  std::vector<MyMatrix<Tint>> elts_G1p = ConjugateElements(X1, X1_inv, elts_G1);
  SubgroupOrbitResult<Tint> orb = OrbitSubgroupConjugation(
      elts_G1p, n, bnd2.normalizer_gens, SubgroupKey(elts_G2), bounds.max_orbit);
  if (!orb.found) {
    res.reason = "same Bravais flock, but not conjugate under N(B2)";
    return res;
  }
  MyMatrix<Tint> Y = *orb.found * X1;
  MyMatrix<Tint> Y_inv = IntegralInverse(Y);
  std::unordered_set<MyMatrix<Tint>> set_G2(elts_G2.begin(), elts_G2.end());
  for (auto const &g : gens_G1) {
    if (set_G2.count(Y * g * Y_inv) == 0) {
      std::cerr << "FGN: the conjugating matrix is wrong\n";
      throw TerminalException{1};
    }
  }
  res.conjugator = Y;
  res.reason = "conjugate";
  return res;
}

/*
  Conjugacy of two elements g, h of finite order of GL_n(Z): X g X^{-1} = h.

  If X g X^{-1} = h then X <g> X^{-1} = <h>, so the cyclic groups are
  Z-equivalent, by some X1. Every X with X <g> X^{-1} = <h> is then Y X1
  with Y in N(<h>), so g and h are conjugate iff h is in the orbit of
  X1 g X1^{-1}, an element of <h>, under the conjugation action of N(<h>).
  The orbit is finite (it lies in <h>). N(<h>) is obtained from the
  normalizer of the Bravais group of <h>, which the group conjugacy has
  already computed.

  The characteristic polynomial of g is prod_d Phi_d^{r_d}, determined by
  the ranks of the g^k - I for the k dividing the order, which are
  compared first.
 */
template <typename T, typename Tint>
std::vector<int> FixedSpaceRanks(MyMatrix<Tint> const &g, int m) {
  int n = g.rows();
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
  MyMatrix<T> pw = IdentityMat<T>(n);
  std::vector<int> res;
  for (int k = 1; k <= m; k++) {
    pw = pw * g_T;
    if (m % k == 0) {
      res.push_back(n - RankMat(MyMatrix<T>(pw - IdentityMat<T>(n))));
    }
  }
  return res;
}

template <typename Tint> int FiniteMatrixOrder(MyMatrix<Tint> const &g, size_t max_order) {
  int n = g.rows();
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  MyMatrix<Tint> pw = g;
  for (size_t k = 1; k <= max_order; k++) {
    if (pw == Id) {
      return static_cast<int>(k);
    }
    pw = pw * g;
  }
  std::cerr << "FGN: the matrix is not of finite order\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Tgroup>
FiniteGroupConjugacyResult<Tint>
FiniteElementConjugacy(MyMatrix<Tint> const &g, MyMatrix<Tint> const &h,
                       FiniteGroupNormalizerBounds const &bounds,
                       std::ostream &os) {
  FiniteGroupConjugacyResult<Tint> res;
  int n = g.rows();
  if (h.rows() != n) {
    res.reason = "different dimensions";
    return res;
  }
  int m = FiniteMatrixOrder(g, bounds.max_group);
  if (FiniteMatrixOrder(h, bounds.max_group) != m) {
    res.reason = "different orders";
    return res;
  }
  if (FixedSpaceRanks<T, Tint>(g, m) != FixedSpaceRanks<T, Tint>(h, m)) {
    res.reason = "different characteristic polynomials";
    return res;
  }
  if (g == h) {
    res.conjugator = IdentityMat<Tint>(n);
    res.reason = "conjugate";
    return res;
  }
  std::vector<MyMatrix<Tint>> gens_g{g};
  std::vector<MyMatrix<Tint>> gens_h{h};
  FiniteGroupConjugacyResult<Tint> rg =
      FiniteGroupConjugacy<T, Tint, Tgroup>(gens_g, gens_h, bounds, os);
  if (!rg.conjugator) {
    res.reason = "the cyclic groups are not Z-equivalent (" + rg.reason + ")";
    return res;
  }
  MyMatrix<Tint> X1 = *rg.conjugator;
  MyMatrix<Tint> g1 = X1 * g * IntegralInverse(X1);
  std::vector<MyMatrix<Tint>> elts_h =
      FiniteMatrixGroupElements(gens_h, n, bounds.max_group);
  std::vector<MyMatrix<Tint>> gens_N =
      NormalizerFromBravaisNormalizer(elts_h, gens_h, n, rg.normalizer_B2_gens,
                                      bounds, os)
          .second;
  // the orbit of g1 under the conjugation by N(<h>), with a transversal
  std::vector<MyMatrix<Tint>> gens_N_inv;
  for (auto const &x : gens_N) {
    gens_N_inv.push_back(IntegralInverse(x));
  }
  std::unordered_map<MyMatrix<Tint>, MyMatrix<Tint>> transversal;
  std::vector<MyMatrix<Tint>> work{g1};
  transversal[g1] = IdentityMat<Tint>(n);
  std::optional<MyMatrix<Tint>> found;
  if (g1 == h) {
    found = IdentityMat<Tint>(n);
  }
  for (size_t pos = 0; pos < work.size() && !found; pos++) {
    MyMatrix<Tint> x = work[pos];
    MyMatrix<Tint> t = transversal[x];
    for (size_t u = 0; u < gens_N.size(); u++) {
      MyMatrix<Tint> y = gens_N[u] * x * gens_N_inv[u];
      if (transversal.count(y) == 0) {
        MyMatrix<Tint> ty = gens_N[u] * t;
        transversal[y] = ty;
        work.push_back(y);
        if (y == h) {
          found = ty;
          break;
        }
      }
    }
  }
  if (!found) {
    res.reason = "the cyclic groups are Z-equivalent, but not the elements";
    return res;
  }
  MyMatrix<Tint> Y = *found * X1;
  if (Y * g != h * Y) {
    std::cerr << "FGN: the conjugating matrix of the elements is wrong\n";
    throw TerminalException{1};
  }
  res.conjugator = Y;
  res.reason = "conjugate";
  return res;
}

// clang-format off
#endif  // SRC_PERFECT_FINITEGROUPNORMALIZER_H_
// clang-format on
