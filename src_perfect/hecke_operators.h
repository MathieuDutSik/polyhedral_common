// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_OPERATORS_H_
#define SRC_PERFECT_HECKE_OPERATORS_H_

// clang-format off
#include "perfect_complex.h"
#include "MAT_SparseMatrixMod.h"
#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HECKE_OPERATORS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_HECKE_OPERATORS
#endif

#ifdef TIMINGS
#define TIMINGS_HECKE_OPERATORS
#endif

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
  A clean failure: a method exhausted its budget (memory, iterations,
  cells). It is caught by the driver, which reports it and stops,
  instead of letting the memory or the time grow without bound.
 */
struct HeckeFailure {
  std::string message;
};

// The budgets of the methods. They bound the memory and the running time
// of every step; exceeding one of them is a clean failure.
struct HeckeBudget {
  size_t max_cosets = 20000000;         // cosets of G / Gamma
  size_t max_hecke_cosets = 1000000;    // cosets of a Hecke double coset
  size_t max_closure = 1000000;         // elements of a finite closure
  size_t max_walk_steps = 1000000;      // steps of a walk through the cones
  size_t max_peel_moves = 20000;        // moves of the peeling solver
  size_t max_chain = 50000;             // cells of a chain in the peeling
  size_t max_unknowns = 60000;          // unknowns of a linear system
  size_t max_rings = 12;                // rings of unknowns
  size_t max_entries_modular = 6000000; // fill-in of the modular elimination
  size_t max_entries_exact = 3000000;   // fill-in of the rational elimination
  size_t max_nnz_modular = 400000;      // nonzeros for trying the modular test
  int max_yall1_iter = 20000;           // iterations of the L1 solver
  int n_reweight = 1;                   // reweighting rounds of the L1 solver
  size_t max_sharbly_symbols = 200000;  // symbols of an Ash-Rudolph reduction
};

//
// The ambient group G of the T-space and its generators
//

template <typename T>
bool is_in_ambient_group(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &g) {
  if (!IsIntegralMatrix(g)) {
    return false;
  }
  T det = DeterminantMat(g);
  if (det != T(1) && det != T(-1)) {
    return false;
  }
  return is_stab_space(g, LinSpa);
}

// The stabilizers of the perfect cones together with the adjacency
// elements generate the group of a connected complex.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
get_ambient_group_generators(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                             std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  std::unordered_set<MyMatrix<Tint>> set_gen;
  std::vector<MyMatrix<Tint>> l_gens;
  auto f_insert = [&](MyMatrix<Tint> const &M) -> void {
    if (M == id) {
      return;
    }
    if (set_gen.insert(M).second) {
      l_gens.push_back(M);
    }
  };
  for (auto &cone : fce.pctdi.l_perfect) {
    for (auto &ePerm : cone.GRP_ext.GeneratorsOfGroup()) {
      f_insert(cone.find_matrix(ePerm, os));
    }
    for (auto &adj : cone.l_sing_adj) {
      f_insert(adj.eMat);
    }
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: get_ambient_group_generators |l_gens|=" << l_gens.size() << "\n";
#endif
  return l_gens;
}

//
// Finite index subgroups Gamma of G
//

/*
  The supported subgroups all contain the principal congruence subgroup
  of level N, so that the coset g Gamma is determined by g modulo N:
  - "Full": Gamma = G (N = 1).
  - "Principal": g = Id mod N.
  - "Gamma0": the last row of g is (0, ..., 0, *) mod N, i.e. Gamma is
    the stabilizer of the line spanned by e_n modulo N for the right
    action on row vectors.
  - "Gamma1": the last row of g is (0, ..., 0, 1) mod N.
  - "SL", "PrincipalSL", "Gamma0SL", "Gamma1SL": the intersection of
    "Full", "Principal", "Gamma0", "Gamma1" with det(g) = 1. For a group
    containing elements of determinant -1, they can reverse the
    orientation of the cells and kill them in the quotient complex, so
    the cohomology of the subgroups of SL_n(Z), which is the one of the
    literature, needs those. (The determinant is the one of the integral
    matrix, so this is for the classic T-space.)
  get_key(g) returns a complete invariant of the coset g Gamma, which
  is what the coset table is keyed by.
 */
/*
  The keys of the cosets g Gamma are vectors of residues modulo N, computed
  with int64_t arithmetic: they are formed millions of times for the large
  levels. The determinants of the elements of G are +-1, so the inverses
  are given by the adjugates.
 */
using GammaKey = std::vector<int64_t>;

struct GammaKeyHash {
  size_t operator()(GammaKey const &key) const {
    size_t h = key.size();
    for (auto &x : key) {
      h = h * 1000003 ^ std::hash<int64_t>()(x);
    }
    return h;
  }
};

// The determinant modulo m of the n x n matrix a (row major), by the
// cofactor expansion (n is small).
inline int64_t determinant_mod(std::vector<int64_t> const &a, int const &n, int64_t const &m) {
  if (n == 1) {
    return ResidueModPrime(a[0], m);
  }
  int64_t det = 0;
  std::vector<int64_t> minor((n - 1) * (n - 1));
  for (int j = 0; j < n; j++) {
    if (a[j] == 0) {
      continue;
    }
    int pos = 0;
    for (int i2 = 1; i2 < n; i2++) {
      for (int j2 = 0; j2 < n; j2++) {
        if (j2 != j) {
          minor[pos] = a[i2 * n + j2];
          pos += 1;
        }
      }
    }
    int64_t term = (a[j] * determinant_mod(minor, n - 1, m)) % m;
    det = (j % 2 == 0) ? (det + term) % m : (det - term) % m;
  }
  return ResidueModPrime(det, m);
}

// The inverse modulo m of a matrix of determinant det_mod (a unit mod m):
// inv(i, j) = det^{-1} (-1)^{i+j} minor(j, i).
inline std::vector<int64_t> inverse_mod(std::vector<int64_t> const &a, int const &n,
                                        int64_t const &m, int64_t const &det_mod) {
  std::vector<int64_t> inv(n * n, 0);
  if (m == 1) {
    return inv;
  }
  int64_t det_inv = InverseModPrime(det_mod, m);
  std::vector<int64_t> minor((n - 1) * (n - 1));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // The minor deleting the row j and the column i.
      int pos = 0;
      for (int i2 = 0; i2 < n; i2++) {
        if (i2 == j) {
          continue;
        }
        for (int j2 = 0; j2 < n; j2++) {
          if (j2 != i) {
            minor[pos] = a[i2 * n + j2];
            pos += 1;
          }
        }
      }
      int64_t val = (n == 1) ? 1 : determinant_mod(minor, n - 1, m);
      if ((i + j) % 2 == 1) {
        val = m - val;
      }
      inv[i * n + j] = (ResidueModPrime(val, m) * det_inv) % m;
    }
  }
  return inv;
}

// What the action of an element s of G on the keys needs: its residues,
// the ones of its inverse and the sign of its determinant.
struct GammaAction {
  std::vector<int64_t> s_mod;
  std::vector<int64_t> s_inv_mod;
  int64_t det_sign;
};

template <typename Tint> struct FiniteIndexSubgroup {
  std::string type;
  Tint N;
  int n;
  bool is_sl() const {
    return type == "SL" || type == "PrincipalSL" || type == "Gamma0SL" || type == "Gamma1SL";
  }
  // The type without the determinant condition.
  std::string base_type() const {
    if (type == "SL") {
      return "Full";
    }
    if (is_sl()) {
      return type.substr(0, type.size() - 2);
    }
    return type;
  }
  bool is_in_gamma(MyMatrix<Tint> const &g) const {
    if (is_sl()) {
      if (DeterminantMat(g) != 1) {
        return false;
      }
    }
    std::string type = base_type();
    if (type == "Full") {
      return true;
    }
    if (type == "Principal") {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Tint val = g(i, j);
          if (i == j) {
            val -= 1;
          }
          if (ResInt(val, N) != 0) {
            return false;
          }
        }
      }
      return true;
    }
    if (type == "Gamma0" || type == "Gamma1") {
      for (int j = 0; j < n - 1; j++) {
        if (ResInt(g(n - 1, j), N) != 0) {
          return false;
        }
      }
      if (type == "Gamma1") {
        Tint val = g(n - 1, n - 1) - 1;
        if (ResInt(val, N) != 0) {
          return false;
        }
      }
      return true;
    }
    std::cerr << "HECKE: Unsupported subgroup type=" << type << "\n";
    throw TerminalException{1};
  }
  int64_t get_n_mod() const { return UniversalScalarConversion<int64_t, Tint>(N); }
  std::vector<int64_t> get_residues(MyMatrix<Tint> const &g) const {
    std::vector<int64_t> r(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        r[i * n + j] = UniversalScalarConversion<int64_t, Tint>(ResInt(g(i, j), N));
      }
    }
    return r;
  }
  // The sign of det(g) = +-1, from its residue when N > 2.
  int64_t get_det_sign(MyMatrix<Tint> const &g, std::vector<int64_t> const &r) const {
    int64_t m = get_n_mod();
    if (m > 2) {
      return (determinant_mod(r, n, m) == 1) ? 1 : -1;
    }
    Tint det = DeterminantMat(g);
    return (det > 0) ? 1 : -1;
  }
  // The canonical form of the row of a Gamma0 key under the units mod N:
  // the first nonzero entry is 1 for N prime, the smallest multiple in the
  // lexicographic order otherwise.
  void canonicalize_gamma0(std::vector<int64_t> &row) const {
    int64_t m = get_n_mod();
    if (is_prime_n()) {
      for (int j = 0; j < n; j++) {
        if (row[j] != 0) {
          int64_t inv = InverseModPrime(row[j], m);
          for (int k = 0; k < n; k++) {
            row[k] = (row[k] * inv) % m;
          }
          return;
        }
      }
      return;
    }
    std::vector<int64_t> best = row;
    for (int64_t u = 2; u < m; u++) {
      if (std::gcd(u, m) != 1) {
        continue;
      }
      std::vector<int64_t> cand(n);
      for (int k = 0; k < n; k++) {
        cand[k] = (u * row[k]) % m;
      }
      if (cand < best) {
        best = cand;
      }
    }
    row = best;
  }
  bool is_prime_n() const {
    int64_t m = get_n_mod();
    if (m < 2) {
      return false;
    }
    for (int64_t d = 2; d * d <= m; d++) {
      if (m % d == 0) {
        return false;
      }
    }
    return true;
  }
  /*
    A complete invariant of the coset g Gamma, as residues modulo N:
    - "Full": empty; "Principal": g mod N;
    - "Gamma0" / "Gamma1": g Gamma = h Gamma iff e_n h^{-1} g = lambda e_n
      mod N, that is iff the rows e_n g^{-1} and e_n h^{-1} agree modulo N
      up to a unit lambda (lambda = 1 for Gamma1);
    - for the SL variants the sign of det(g) is appended.
   */
  GammaKey get_key(MyMatrix<Tint> const &g) const {
    std::vector<int64_t> r = get_residues(g);
    int64_t m = get_n_mod();
    std::string type = base_type();
    GammaKey key;
    if (type == "Principal") {
      key = r;
    } else if (type == "Gamma0" || type == "Gamma1") {
      int64_t det_sign = get_det_sign(g, r);
      int64_t det_mod = ResidueModPrime(det_sign, m);
      std::vector<int64_t> inv = inverse_mod(r, n, m, (m == 1) ? 0 : det_mod);
      key = std::vector<int64_t>(inv.begin() + (n - 1) * n, inv.end());
      if (type == "Gamma0") {
        canonicalize_gamma0(key);
      }
    } else if (type != "Full") {
      std::cerr << "HECKE: Unsupported subgroup type=" << type << "\n";
      throw TerminalException{1};
    }
    if (is_sl()) {
      key.push_back(get_det_sign(g, r));
    }
    return key;
  }
  GammaAction get_action(MyMatrix<Tint> const &s) const {
    GammaAction a;
    int64_t m = get_n_mod();
    a.s_mod = get_residues(s);
    a.det_sign = get_det_sign(s, a.s_mod);
    a.s_inv_mod = inverse_mod(a.s_mod, n, m, (m == 1) ? 0 : ResidueModPrime(a.det_sign, m));
    return a;
  }
  // The key of s g Gamma from the key of g Gamma: s g mod N for
  // "Principal", e_n g^{-1} s^{-1} for "Gamma0" / "Gamma1".
  GammaKey act(GammaKey const &key, GammaAction const &a) const {
    int64_t m = get_n_mod();
    std::string type = base_type();
    GammaKey key_ret;
    if (type == "Principal") {
      key_ret = std::vector<int64_t>(n * n, 0);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          int64_t sum = 0;
          for (int k = 0; k < n; k++) {
            sum = (sum + a.s_mod[i * n + k] * key[k * n + j]) % m;
          }
          key_ret[i * n + j] = sum;
        }
      }
    } else if (type == "Gamma0" || type == "Gamma1") {
      key_ret = std::vector<int64_t>(n, 0);
      for (int j = 0; j < n; j++) {
        int64_t sum = 0;
        for (int k = 0; k < n; k++) {
          sum = (sum + key[k] * a.s_inv_mod[k * n + j]) % m;
        }
        key_ret[j] = sum;
      }
      if (type == "Gamma0") {
        canonicalize_gamma0(key_ret);
      }
    }
    if (is_sl()) {
      key_ret.push_back(a.det_sign * key.back());
    }
    return key_ret;
  }
};

/*
  The cosets G = union g_i Gamma with the map g -> (i, gamma) such that
  g = g_i gamma.
 */
template <typename Tint> struct CosetTableGamma {
  FiniteIndexSubgroup<Tint> gamma;
  std::vector<MyMatrix<Tint>> l_coset;
  std::vector<MyMatrix<Tint>> l_coset_inv;
  std::vector<GammaKey> l_key;
  std::unordered_map<GammaKey, int, GammaKeyHash> map_key;
  size_t size() const { return l_coset.size(); }
  std::optional<int> find_index(MyMatrix<Tint> const &g) const {
    GammaKey key = gamma.get_key(g);
    auto iter = map_key.find(key);
    if (iter == map_key.end()) {
      return {};
    }
    return iter->second;
  }
  // The index of the coset s g_i Gamma, from the key of g_i Gamma.
  int index_of_action(int const &i, GammaAction const &a) const {
    auto iter = map_key.find(gamma.act(l_key[i], a));
    if (iter == map_key.end()) {
      std::cerr << "HECKE: The coset of s g was not found in the coset table\n";
      throw TerminalException{1};
    }
    return iter->second;
  }
  int index_of(MyMatrix<Tint> const &g) const {
    std::optional<int> opt = find_index(g);
    if (!opt) {
      std::cerr << "HECKE: The coset of g was not found in the coset table\n";
      throw TerminalException{1};
    }
    return *opt;
  }
  std::pair<int, MyMatrix<Tint>> decompose(MyMatrix<Tint> const &g) const {
    int i = index_of(g);
    MyMatrix<Tint> gam = l_coset_inv[i] * g;
#ifdef SANITY_CHECK_HECKE_OPERATORS
    if (!gamma.is_in_gamma(gam)) {
      std::cerr << "HECKE: The coset key does not match the membership\n";
      throw TerminalException{1};
    }
#endif
    return {i, gam};
  }
};

template <typename Tint>
CosetTableGamma<Tint>
enumerate_cosets_gamma(FiniteIndexSubgroup<Tint> const &gamma,
                       std::vector<MyMatrix<Tint>> const &l_gens, int n,
                       [[maybe_unused]] std::ostream &os) {
  CosetTableGamma<Tint> ct;
  ct.gamma = gamma;
  // The keys of the new cosets are obtained by the action on the keys,
  // and the inverses as products.
  std::vector<GammaAction> l_act;
  std::vector<MyMatrix<Tint>> l_gens_inv;
  for (auto &s : l_gens) {
    l_act.push_back(gamma.get_action(s));
    l_gens_inv.push_back(Inverse(s));
  }
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  ct.l_key.push_back(gamma.get_key(id));
  ct.map_key[ct.l_key[0]] = 0;
  ct.l_coset.push_back(id);
  ct.l_coset_inv.push_back(id);
  size_t pos = 0;
  HeckeBudget budget;
  while (pos < ct.l_coset.size()) {
    if (ct.l_coset.size() > budget.max_cosets) {
      throw HeckeFailure{"enumerate_cosets_gamma: more than " + std::to_string(budget.max_cosets) + " cosets (is the index finite?)"};
    }
    for (size_t i_gen = 0; i_gen < l_gens.size(); i_gen++) {
      GammaKey key = gamma.act(ct.l_key[pos], l_act[i_gen]);
      if (ct.map_key.contains(key)) {
        continue;
      }
      int idx = ct.l_coset.size();
      ct.map_key[key] = idx;
      ct.l_key.push_back(key);
      // The vectors may be reallocated by the push_back, so we copy.
      MyMatrix<Tint> g = ct.l_coset[pos];
      MyMatrix<Tint> g_inv = ct.l_coset_inv[pos];
      ct.l_coset.push_back(l_gens[i_gen] * g);
      ct.l_coset_inv.push_back(g_inv * l_gens_inv[i_gen]);
    }
    pos += 1;
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: enumerate_cosets_gamma type=" << gamma.type << " N=" << gamma.N
     << " index=" << ct.l_coset.size() << "\n";
#endif
#ifdef SANITY_CHECK_HECKE_OPERATORS
  for (size_t i = 0; i < ct.l_coset.size(); i++) {
    for (size_t j = 0; j < ct.l_coset.size(); j++) {
      MyMatrix<Tint> prod = ct.l_coset_inv[i] * ct.l_coset[j];
      bool test = gamma.is_in_gamma(prod);
      if (test != (i == j)) {
        std::cerr << "HECKE: Inconsistent coset table at i=" << i << " j=" << j << "\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return ct;
}

// The Schreier generators of Gamma from the coset table.
template <typename Tint>
std::vector<MyMatrix<Tint>>
get_schreier_generators(CosetTableGamma<Tint> const &ct,
                        std::vector<MyMatrix<Tint>> const &l_gens, int n,
                        [[maybe_unused]] std::ostream &os) {
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  std::unordered_set<MyMatrix<Tint>> set_gen;
  std::vector<MyMatrix<Tint>> l_gens_gamma;
  for (auto &g : ct.l_coset) {
    for (auto &s : l_gens) {
      MyMatrix<Tint> prod = s * g;
      std::pair<int, MyMatrix<Tint>> pair = ct.decompose(prod);
      if (pair.second != id) {
        if (set_gen.insert(pair.second).second) {
          l_gens_gamma.push_back(pair.second);
        }
      }
    }
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: get_schreier_generators |l_gens_gamma|=" << l_gens_gamma.size() << "\n";
#endif
  return l_gens_gamma;
}

//
// The cosets of a Hecke double coset
//

/*
  H x H = union x_k H for H a group given by a membership test and a
  generating set. decompose(y) for y in H x H returns (k, h) with
  y = x_k h.
 */
template <typename T> struct HeckeCosets {
  MyMatrix<T> x;
  std::vector<MyMatrix<T>> l_coset;
  std::vector<MyMatrix<T>> l_coset_inv;
  std::function<bool(MyMatrix<T> const &)> f_member;
  size_t size() const { return l_coset.size(); }
  std::optional<std::pair<int, MyMatrix<T>>>
  find_decompose(MyMatrix<T> const &y) const {
    size_t n_coset = l_coset.size();
    for (size_t k = 0; k < n_coset; k++) {
      MyMatrix<T> h = l_coset_inv[k] * y;
      if (f_member(h)) {
        std::pair<int, MyMatrix<T>> pair{static_cast<int>(k), h};
        return pair;
      }
    }
    return {};
  }
  std::pair<int, MyMatrix<T>> decompose(MyMatrix<T> const &y) const {
    std::optional<std::pair<int, MyMatrix<T>>> opt = find_decompose(y);
    if (!opt) {
      std::cerr << "HECKE: The matrix y does not belong to the double coset\n";
      throw TerminalException{1};
    }
    return *opt;
  }
};

template <typename T, typename Tint>
HeckeCosets<T>
enumerate_hecke_cosets(MyMatrix<T> const &x,
                       std::vector<MyMatrix<Tint>> const &l_gens,
                       std::function<bool(MyMatrix<T> const &)> f_member,
                       [[maybe_unused]] std::ostream &os,
                       size_t const &n_expected = 0) {
  HeckeCosets<T> hc;
  hc.x = x;
  hc.f_member = f_member;
  auto f_insert = [&](MyMatrix<T> const &y) -> void {
    if (hc.find_decompose(y)) {
      return;
    }
    hc.l_coset_inv.push_back(Inverse(y));
    hc.l_coset.push_back(y);
  };
  f_insert(x);
  size_t pos = 0;
  HeckeBudget budget;
  while (pos < hc.l_coset.size()) {
    if (hc.l_coset.size() > budget.max_hecke_cosets) {
      throw HeckeFailure{"enumerate_hecke_cosets: more than " + std::to_string(budget.max_hecke_cosets) + " cosets"};
    }
    MyMatrix<T> y = hc.l_coset[pos];
    for (auto &s : l_gens) {
      // When the number of cosets is known, the enumeration stops as soon
      // as they are all found.
      if (n_expected > 0 && hc.l_coset.size() == n_expected) {
        break;
      }
      MyMatrix<T> s_T = UniversalMatrixConversion<T, Tint>(s);
      MyMatrix<T> prod = s_T * y;
      f_insert(prod);
    }
    pos += 1;
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: enumerate_hecke_cosets |l_coset|=" << hc.l_coset.size() << "\n";
#endif
  return hc;
}

//
// Chains: elementary operations
//

template <typename T, typename Tint>
std::vector<PerfectFaceEntry<T, Tint>>
chain_right_multiply(std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                     T const &scal, MyMatrix<Tint> const &M) {
  std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
  for (auto &fe : chain) {
    T value = fe.value * scal;
    MyMatrix<Tint> Mprod = fe.M * M;
    PerfectFaceEntry<T, Tint> fe_new{fe.iOrb, value, std::move(Mprod)};
    chain_ret.emplace_back(std::move(fe_new));
  }
  return chain_ret;
}

// The vectors of all the cells of a chain.
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
chain_vector_family(int const &index,
                    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                    FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  int n = fce.pctdi.LinSpa.n;
  std::vector<MyVector<Tint>> l_vect;
  for (auto &fe : chain) {
    MyMatrix<Tint> EXT = fce.levels[index].l_faces[fe.iOrb].EXT * fe.M;
    for (int i_row = 0; i_row < EXT.rows(); i_row++) {
      l_vect.push_back(GetMatrixRow(EXT, i_row));
    }
  }
  return MatrixFromVectorFamilyDim(n, l_vect);
}

/*
  The filter selecting the cells whose vectors belong to the span of the
  vectors of the chain. The span is saturated by the matrices of ListComm
  (e.g. the multiplication by omega for Hermitian T-spaces): the cells at
  infinity of the complex live in subspaces stable under them and the
  plain rational span of the vectors would contain no cell at all. If the
  span is everything then no filtering.
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<std::function<bool(MyMatrix<Tint> const &)>, bool>
get_subspace_filter(int const &index,
                    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                    FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                    [[maybe_unused]] std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> V_int = chain_vector_family(index, chain, fce);
  MyMatrix<T> V = UniversalMatrixConversion<T, Tint>(V_int);
  std::vector<MyMatrix<T>> const &ListComm = fce.pctdi.LinSpa.ListComm;
  int rnk = RankMat(V);
  while (true) {
    std::vector<MyVector<T>> l_vect;
    for (int i_row = 0; i_row < V.rows(); i_row++) {
      l_vect.push_back(GetMatrixRow(V, i_row));
    }
    for (auto &C : ListComm) {
      MyMatrix<T> Vimg = V * C;
      for (int i_row = 0; i_row < Vimg.rows(); i_row++) {
        l_vect.push_back(GetMatrixRow(Vimg, i_row));
      }
    }
    MyMatrix<T> V_new = MatrixFromVectorFamilyDim(n, l_vect);
    int rnk_new = RankMat(V_new);
    if (rnk_new == rnk) {
      break;
    }
    V = V_new;
    rnk = rnk_new;
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: get_subspace_filter |V|=" << V.rows() << " rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk == n) {
    std::function<bool(MyMatrix<Tint> const &)> f_all = [](MyMatrix<Tint> const &) -> bool { return true; };
    return {f_all, false};
  }
  // The rows e of Equa satisfy V e = 0.
  MyMatrix<T> Equa = NullspaceTrMat(V);
  MyMatrix<T> EquaT = Equa.transpose();
  std::function<bool(MyMatrix<Tint> const &)> f_sub = [EquaT](MyMatrix<Tint> const &EXT) -> bool {
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    MyMatrix<T> prod = EXT_T * EquaT;
    return IsZeroMatrix(prod);
  };
  return {f_sub, true};
}

//
// Walking through the top dimensional cones to reach a cell
//

/*
  The generic search of a cell (fce_face_search) goes through the
  computation of a positive definite form with the cell as set of
  shortest vectors, which is not adapted to the cells at infinity. Here
  we instead walk from a starting cone along a segment towards the
  barycenter of the cell: at every step the exit facet of the segment is
  computed and the adjacent cone is entered. The walk ends with a cone
  whose closure contains the barycenter, which forces the cell to be a
  face of that cone. The starting point is chosen randomly in the
  interior of the starting cone so that the segment is generic.

  All computations are done in the frame of the reference cone: the cone
  EXT_i M is mapped back to EXT_i by w -> w M^{-1} on the vectors, and a
  vector w corresponds to the point of coordinates (w ListMat_k w^T)_k
  of the T-space.
 */
template <typename T, typename Tint, typename Tgroup> struct ConeWalker {
private:
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  std::vector<MyMatrix<T>> l_scal;
  std::vector<std::vector<MyVector<T>>> ll_facet_fct;
  std::ostream &os;
  MyVector<T> get_coords(MyMatrix<T> const &W,
                         std::vector<T> const &weights) const {
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyVector<T> coords = ZeroVector<T>(n_mat);
    int n_row = W.rows();
    for (int i_row = 0; i_row < n_row; i_row++) {
      MyVector<T> V = GetMatrixRow(W, i_row);
      for (int k = 0; k < n_mat; k++) {
        coords(k) += weights[i_row] * EvaluationQuadForm<T, T>(ListMat[k], V);
      }
    }
    return coords;
  }
public:
  ConeWalker(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, std::ostream &_os)
      : fce(_fce), os(_os) {
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    for (auto &cone : fce.pctdi.l_perfect) {
      MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, cone.EXT);
      l_scal.push_back(ScalMat);
      std::vector<MyVector<T>> l_fct;
      for (auto &adj : cone.full_adjacencies(os)) {
        MyMatrix<T> ScalFacet = SelectRow(ScalMat, adj.f_ext);
        MyMatrix<T> NSP = NullspaceTrMat(ScalFacet);
        if (NSP.rows() != 1) {
          std::cerr << "HECKE: The facet should define a functional of dimension 1, |NSP|=" << NSP.rows() << "\n";
          throw TerminalException{1};
        }
        MyVector<T> fct = GetMatrixRow(NSP, 0);
        // Orientation: positive on the vectors outside of the facet.
        int n_row = ScalMat.rows();
        T sign_val(0);
        for (int i_row = 0; i_row < n_row; i_row++) {
          if (adj.f_ext[i_row] == 0) {
            MyVector<T> V = GetMatrixRow(ScalMat, i_row);
            T val = ScalarProduct(V, fct);
            if (val != 0) {
              sign_val = val;
              break;
            }
          }
        }
        if (sign_val == 0) {
          std::cerr << "HECKE: Failed to orient the facet functional\n";
          throw TerminalException{1};
        }
        if (sign_val < 0) {
          fct = -fct;
        }
        l_fct.push_back(fct);
      }
      ll_facet_fct.push_back(l_fct);
    }
  }
  /*
    A top dimensional cone (i_perfect, M) whose closure contains the
    barycenter of the rays of the (rational) vectors of P.
   */
  /*
    Whether the closure of the cone (i_perfect, M) contains the barycenter
    of the rays of the vectors of P. If yes, returns the number of facets
    containing it.
   */
  std::optional<int> n_facet_containing(int const &i_perfect, MyMatrix<Tint> const &M,
                                        MyMatrix<T> const &P) const {
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights_p(P.rows(), T(1));
    MyVector<T> cp = get_coords(P * Minv, weights_p);
    int n_tight = 0;
    for (auto &fct : ll_facet_fct[i_perfect]) {
      T b = ScalarProduct(fct, cp);
      if (b < 0) {
        return {};
      }
      if (b == 0) {
        n_tight += 1;
      }
    }
    return n_tight;
  }
  TopPerfectCone<Tint> find_containing_cone(MyMatrix<Tint> const &EXT_target) const {
    MyMatrix<T> P = UniversalMatrixConversion<T, Tint>(EXT_target);
    return find_containing_cone_T(P);
  }
  TopPerfectCone<Tint> find_containing_cone_T(MyMatrix<T> const &P) const {
    int n = fce.pctdi.LinSpa.n;
    std::vector<T> weights_p(P.rows(), T(1));
    MyMatrix<T> Q0 = UniversalMatrixConversion<T, Tint>(fce.pctdi.l_perfect[0].EXT);
    size_t n_attempt = 0;
    uint64_t seed = 12345;
    auto next_random = [&]() -> uint64_t {
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      return seed >> 33;
    };
    while (true) {
      n_attempt += 1;
      if (n_attempt > 20) {
        throw HeckeFailure{"find_containing_cone: the walk cycled in 20 attempts"};
      }
      std::vector<T> weights_q(Q0.rows());
      for (int i_row = 0; i_row < Q0.rows(); i_row++) {
        int w_val = 1 + static_cast<int>(next_random() % 1000);
        weights_q[i_row] = T(w_val);
      }
      int i_cone = 0;
      MyMatrix<Tint> M = IdentityMat<Tint>(n);
      std::unordered_set<MyMatrix<Tint>> set_visited;
      bool cycle = false;
      HeckeBudget budget;
      while (true) {
        if (set_visited.size() > budget.max_walk_steps) {
          throw HeckeFailure{"find_containing_cone: walk longer than " + std::to_string(budget.max_walk_steps) + " steps"};
        }
        MyMatrix<Tint> EXT_cone = fce.pctdi.l_perfect[i_cone].EXT * M;
        MyMatrix<Tint> key = tot_set(EXT_cone);
        if (!set_visited.insert(key).second) {
          cycle = true;
          break;
        }
        MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
        MyMatrix<T> Minv = Inverse(M_T);
        MyVector<T> cq = get_coords(Q0 * Minv, weights_q);
        MyVector<T> cp = get_coords(P * Minv, weights_p);
        std::optional<T> best_t;
        int best_f = -1;
        std::vector<MyVector<T>> const &l_fct = ll_facet_fct[i_cone];
        int n_facet = l_fct.size();
        for (int i_f = 0; i_f < n_facet; i_f++) {
          T a = ScalarProduct(l_fct[i_f], cq);
          T b = ScalarProduct(l_fct[i_f], cp);
          if (b < 0) {
            T t_f = a / (a - b);
            if (!best_t || t_f < *best_t) {
              best_t = t_f;
              best_f = i_f;
            }
          }
        }
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: walk i_cone=" << i_cone << " best_f=" << best_f;
        if (best_t) {
          os << " best_t=" << *best_t;
        }
        os << " cp=" << StringVectorGAP(cp) << " cq=" << StringVectorGAP(cq) << "\n";
#endif
        if (best_f == -1) {
          return {i_cone, M};
        }
        sing_adj<Tint> const &adj = fce.pctdi.l_perfect[i_cone].full_adjacencies(os)[best_f];
        M = adj.eMat * M;
        i_cone = adj.jCone;
      }
      if (cycle) {
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: find_containing_cone, cycle detected, new attempt\n";
#endif
      }
    }
  }
  /*
    The vertex of the complex on the ray of the functional
    Q -> sum_w Q[w] for the rational vectors w of W. The vectors of the
    vertex are the vectors of the containing cone whose functional is
    proportional to it (they are not necessarily proportional to the
    rows of W over Z, e.g. for Hermitian T-spaces).
   */
  FceFaceSearch<Tint> search_vertex(MyMatrix<T> const &W) const {
    TopPerfectCone<Tint> tpc = find_containing_cone_T(W);
    int i_cone = tpc.i_perfect;
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights(W.rows(), T(1));
    MyVector<T> p = get_coords(W * Minv, weights);
    int n_mat = p.size();
    int l_nz = -1;
    for (int k = 0; k < n_mat; k++) {
      if (p(k) != 0) {
        l_nz = k;
        break;
      }
    }
    if (l_nz == -1) {
      std::cerr << "HECKE: The functional of the vertex should be nonzero\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &ScalMat = l_scal[i_cone];
    int n_row = ScalMat.rows();
    Face f_ext(n_row);
    for (int i_row = 0; i_row < n_row; i_row++) {
      bool is_prop = true;
      for (int k = 0; k < n_mat; k++) {
        if (ScalMat(i_row, k) * p(l_nz) != p(k) * ScalMat(i_row, l_nz)) {
          is_prop = false;
          break;
        }
      }
      if (is_prop) {
        f_ext[i_row] = 1;
      }
    }
    if (f_ext.count() == 0) {
      std::cerr << "HECKE: No vector of the containing cone is on the ray of the vertex\n";
      throw TerminalException{1};
    }
    triple<Tint> t{static_cast<size_t>(i_cone), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    int vlevel = fce.levels.size() - 1;
    int n_orb = fce.levels[vlevel].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[vlevel].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {vlevel, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the vertex among the orbits\n";
    throw TerminalException{1};
  }
  /*
    The smallest face of the complex whose relative interior contains the
    barycenter of the rays of the rational vectors of W: the walk gives a
    cone containing it, the tight facets of the cone cut out the face.
   */
  FceFaceSearch<Tint> search_face_containing(MyMatrix<T> const &W) const {
    TopPerfectCone<Tint> tpc = find_containing_cone_T(W);
    int i_cone = tpc.i_perfect;
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights(W.rows(), T(1));
    MyVector<T> p = get_coords(W * Minv, weights);
    std::vector<sing_adj<Tint>> const &l_adj = fce.pctdi.l_perfect[i_cone].full_adjacencies(os);
    std::vector<MyVector<T>> const &l_fct = ll_facet_fct[i_cone];
    int n_row = fce.pctdi.l_perfect[i_cone].EXT.rows();
    Face f_ext(n_row);
    f_ext.set();
    for (size_t i_f = 0; i_f < l_fct.size(); i_f++) {
      T val = ScalarProduct(l_fct[i_f], p);
      if (val == 0) {
        f_ext &= l_adj[i_f].f_ext;
      }
    }
    triple<Tint> t{static_cast<size_t>(i_cone), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    MyMatrix<Tint> EXT_face = SelectRow(fce.pctdi.l_perfect[i_cone].EXT, f_ext);
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, EXT_face);
    int index = n_mat - RankMat(ScalMat);
    int n_orb = fce.levels[index].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[index].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {index, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the face among the orbits of level " << index << "\n";
    throw TerminalException{1};
  }
  /*
    The position of a cell in the complex, by walking to a cone
    containing it.
   */
  FceFaceSearch<Tint> search_cell(MyMatrix<Tint> const &EXT_cell) const {
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: search_cell EXT_cell=" << StringMatrixGAP(EXT_cell) << "\n";
#endif
    TopPerfectCone<Tint> tpc = find_containing_cone(EXT_cell);
    MyMatrix<Tint> EXT_cone = fce.pctdi.l_perfect[tpc.i_perfect].EXT * tpc.M;
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: search_cell i_perfect=" << tpc.i_perfect << " M=" << StringMatrixGAP(tpc.M) << " EXT_cone=" << StringMatrixGAP(EXT_cone) << "\n";
#endif
    ContainerMatrix<Tint> cont(EXT_cone);
    Face f_ext(EXT_cone.rows());
    for (int i_row = 0; i_row < EXT_cell.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(EXT_cell, i_row);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      if (!opt) {
        std::cerr << "HECKE: The vector of the cell is not a vector of the containing cone\n";
        throw TerminalException{1};
      }
      f_ext[*opt] = 1;
    }
    triple<Tint> t{static_cast<size_t>(tpc.i_perfect), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, EXT_cell);
    int index = n_mat - RankMat(ScalMat);
    int n_orb = fce.levels[index].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[index].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {index, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the cell among the orbits of level " << index << "\n";
    throw TerminalException{1};
  }
};

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

//
// The dual construction: induction from the perfect cones
//

/*
  The induction from the vertices goes through the cells at infinity,
  whose stabilizers are infinite and whose stars are infinite; this is
  what makes the coverings explode for GL_4(Z) and GL_3(O_K). The dual
  construction runs the induction in the other direction, on the well
  rounded cells only: the images of the perfect cones (the vertices of
  the well rounded retract W) are found by the walk to the image of the
  barycenter; a facet (an edge of W, between two cones) is mapped to a
  path of adjacent cones; a cell of codimension k is mapped to a chain
  of cells of codimension k whose "coboundary" (the signed sum of the
  cells of codimension k-1 containing it) equals the image of its own
  coboundary. The chain complex is the transposed one, the cellular chain
  complex of W, and everything in it is finite. The harmonic
  representatives of the homology are the same as for the original
  complex, so step B applies unchanged.
 */

// The matrix A of the action Q -> M Q M^T on the T-space: the image of
// ListMat[k] is sum_l A(k,l) ListMat[l].
template <typename T>
MyMatrix<T> tspace_action_matrix(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &M) {
  int n_mat = LinSpa.ListMat.size();
  MyMatrix<T> A(n_mat, n_mat);
  for (int k = 0; k < n_mat; k++) {
    MyMatrix<T> Qimg = M * LinSpa.ListMat[k] * M.transpose();
    MyVector<T> V = SymmetricMatrixToVector(Qimg);
    std::optional<MyVector<T>> opt = SolutionMat(LinSpa.ListMatAsBigMat, V);
    if (!opt) {
      std::cerr << "HECKE: The matrix does not preserve the T-space\n";
      throw TerminalException{1};
    }
    for (int l = 0; l < n_mat; l++) {
      A(k, l) = (*opt)(l);
    }
  }
  return A;
}

// The sign of the determinant of the action Q -> M Q M^T on the T-space.
template <typename T>
int tspace_orientation_sign(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &M) {
  MyMatrix<T> A = tspace_action_matrix(LinSpa, M);
  T det = DeterminantMat(A);
  if (det > 0) {
    return 1;
  }
  if (det < 0) {
    return -1;
  }
  std::cerr << "HECKE: The action on the T-space is degenerate\n";
  throw TerminalException{1};
}

/*
  The cofaces of the representative cells, with the sign of the cell in
  the boundary of the coface, cached per orbit. Only for well rounded
  cells (finite stabilizers).
 */
template <typename T, typename Tint, typename Tgroup> struct UpperFaceCache {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  std::ostream &os;
  // cofaces[index][iOrb] : (coface (jOrb, P) at level index-1, sign)
  std::map<std::pair<int, int>, std::vector<std::pair<PerfectFace<Tint>, int>>> cofaces;
  // stars[index][iOrb] : the perfect cones containing the representative
  std::map<std::pair<int, int>, std::vector<PerfectFace<Tint>>> stars;
  UpperFaceCache(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, std::ostream &_os)
      : fce(_fce), os(_os) {}
  std::vector<std::pair<PerfectFace<Tint>, int>> const &get_cofaces(int const &index, int const &iOrb) {
    std::pair<int, int> key{index, iOrb};
    auto iter = cofaces.find(key);
    if (iter != cofaces.end()) {
      return iter->second;
    }
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    MyMatrix<Tint> key_cell = tot_set(face.EXT);
    std::pair<std::vector<MyMatrix<Tint>>, std::vector<PerfectFace<Tint>>> pair =
        get_all_upper_faces(fce, index, iOrb, os);
    std::vector<std::pair<PerfectFace<Tint>, int>> l_ret;
    for (auto &pf : pair.second) {
      // The sign of F_iOrb in d(F_jOrb P)
      int sign_tot = 0;
      int n_match = 0;
      for (auto &e : fce.boundaries[index - 1].ll_bound[pf.iOrb].l_bound) {
        if (e.iOrb == iOrb) {
          MyMatrix<Tint> Mcell = e.M * pf.M;
          MyMatrix<Tint> EXT_prod = face.EXT * Mcell;
          MyMatrix<Tint> EXT2 = tot_set(EXT_prod);
          if (EXT2 == key_cell) {
            MyMatrix<Tint> t = Inverse(Mcell);
            int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
            sign_tot += e.sign * sign;
            n_match += 1;
          }
        }
      }
      if (n_match != 1) {
        std::cerr << "HECKE: The cell should appear exactly once in the boundary of its coface, n_match=" << n_match << "\n";
        throw TerminalException{1};
      }
      l_ret.emplace_back(pf, sign_tot);
    }
    auto ret = cofaces.emplace(key, std::move(l_ret));
    return ret.first->second;
  }
  std::vector<PerfectFace<Tint>> const &get_star(int const &index, int const &iOrb) {
    std::pair<int, int> key{index, iOrb};
    auto iter = stars.find(key);
    if (iter != stars.end()) {
      return iter->second;
    }
    std::vector<PerfectFace<Tint>> l_star;
    if (index == 0) {
      int n = fce.pctdi.LinSpa.n;
      l_star.push_back({iOrb, IdentityMat<Tint>(n)});
    } else {
      std::unordered_set<MyMatrix<Tint>> set_key;
      for (auto &pair : get_cofaces(index, iOrb)) {
        std::vector<PerfectFace<Tint>> const &l_up = get_star(index - 1, pair.first.iOrb);
        // A copy: l_up is a reference into the map that the recursive
        // calls may not invalidate (std::map), but we keep it simple.
        std::vector<PerfectFace<Tint>> l_up_copy = l_up;
        for (auto &pf : l_up_copy) {
          MyMatrix<Tint> M = pf.M * pair.first.M;
          MyMatrix<Tint> EXT_cone = fce.levels[0].l_faces[pf.iOrb].EXT * M;
          MyMatrix<Tint> key_cone = tot_set(EXT_cone);
          if (set_key.insert(key_cone).second) {
            l_star.push_back({pf.iOrb, M});
          }
        }
      }
    }
    auto ret = stars.emplace(key, std::move(l_star));
    return ret.first->second;
  }
};

// The coboundary of a chain: the signed sum of the cofaces.
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
chain_coboundary(int const &index, std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                 UpperFaceCache<T, Tint, Tgroup> &ufc, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  ChainBuilder<T, Tint, Tgroup> cb(index - 1, fce, os);
  for (auto &fe : chain) {
    for (auto &pair : ufc.get_cofaces(index, fe.iOrb)) {
      T value = fe.value * T(pair.second);
      MyMatrix<Tint> M = pair.first.M * fe.M;
      cb.f_insert(value, pair.first.iOrb, M);
    }
  }
  return cb.get_faces();
}

/*
  The L1 solution of x A = b in floating point, sparsified by a few
  rounds of reweighting (weights 1/(|x_i| + eps)).
 */
inline MyVector<double> l1_solution_reweighted(MySparseMatrix<double> const &SpMat_d,
                                               MyVector<double> const &Bvect_d, int const &maxit,
                                               int const &n_round, [[maybe_unused]] std::ostream &os) {
  RecSparse<double> rec = AMP_linear_operators(SpMat_d);
  RecOptSparse<double> opt;
  opt.delta = 0;
  opt.mu = 0;
  opt.nu = 0;
  opt.eps = 2.2204e-16;
  opt.rho = 0;
  opt.tol = 1e-10;
  opt.gamma = 0;
  opt.nonneg = false;
  opt.nonorth = true;
  opt.print = false;
  opt.UseWeight = false;
  opt.basis = false;
  opt.stepfreq = 1;
  opt.maxit = maxit;
  opt.xs = -1;
  [[maybe_unused]] auto count_above = [](MyVector<double> const &v, double thr) -> int {
    int n_nz = 0;
    for (int i = 0; i < v.size(); i++) {
      if (std::abs(v(i)) > thr) {
        n_nz += 1;
      }
    }
    return n_nz;
  };
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  OutSolver<double> out = AMP_yall1(rec, Bvect_d, opt);
  MyVector<double> x = out.x;
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: yall1 round 0 iter=" << out.iter << " " << out.exit
     << " support(1e-4)=" << count_above(x, 1e-4) << "\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: yall1 round 0|=" << time << "\n";
#endif
  int n = x.size();
  for (int i_round = 0; i_round < n_round; i_round++) {
    double max_abs = 0;
    for (int i = 0; i < n; i++) {
      max_abs = std::max(max_abs, std::abs(x(i)));
    }
    double eps = 1e-3 * max_abs + 1e-12;
    for (int i = 0; i < n; i++) {
      rec.weights(i) = 1.0 / (std::abs(x(i)) + eps);
    }
    out = AMP_yall1(rec, Bvect_d, opt);
    x = out.x;
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: yall1 round " << (i_round + 1) << " iter=" << out.iter << " " << out.exit
       << " support(1e-4)=" << count_above(x, 1e-4) << "\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: yall1 round " << (i_round + 1) << "|=" << time << "\n";
#endif
  }
  return x;
}

/*
  Solution of x A = b for the sparse matrix given by its triplets.
  When is_graph is set (the facet level), the rows with two entries are
  the edges of a graph on the columns and b is a 0-chain: the solution
  is built as a flow, by shortest paths between the nodes with nonzero
  residual, each path being solved exactly. Otherwise the consistency is
  tested modulo a prime when the system is small enough, the L1 solver in
  floating point gives a sparse approximate solution which is solved
  exactly on its support and the columns it touches, and the exact
  rational elimination is the last resort for the small systems. Every
  step is bounded by the budget; the status says whether the system was
  solved, is inconsistent, or could not be decided within the budget.
 */
template <typename T>
SparseSolveResult<T>
solve_chain_system(std::vector<Eigen::Triplet<T>> const &tripletList, int const &nbRow,
                   int const &nbCol, MyVector<T> const &Bvect, bool const &is_graph,
                   HeckeBudget const &budget, std::ostream &os) {
  MyVector<T> empty;
  [[maybe_unused]] auto count_nz = [](MyVector<T> const &v) -> int {
    int n_nz = 0;
    for (int i = 0; i < v.size(); i++) {
      if (v(i) != 0) {
        n_nz += 1;
      }
    }
    return n_nz;
  };
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: solve_chain_system nbRow=" << nbRow << " nbCol=" << nbCol
     << " nnz=" << tripletList.size() << " is_graph=" << is_graph << "\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time_solve;
#endif
  if (is_graph) {
    std::vector<std::vector<std::pair<int, int>>> adj(nbCol);
    std::vector<std::vector<int>> row_cols(nbRow);
    for (auto &t : tripletList) {
      row_cols[t.row()].push_back(t.col());
    }
    for (int i1 = 0; i1 < nbRow; i1++) {
      std::vector<int> const &lv = row_cols[i1];
      if (lv.size() == 2 && lv[0] != lv[1]) {
        adj[lv[0]].push_back({lv[1], i1});
        adj[lv[1]].push_back({lv[0], i1});
      }
    }
    MyVector<T> resid = Bvect;
    MyVector<T> sol = ZeroVector<T>(nbRow);
    size_t n_path = 0;
    while (true) {
      if (n_path > static_cast<size_t>(nbCol)) {
        // Each path zeroes one node, so this cannot happen.
        return {SparseSolveStatus::Aborted, empty};
      }
      int v_start = -1;
      for (int v = 0; v < nbCol; v++) {
        if (resid(v) != 0) {
          if (v_start == -1 || T_abs(resid(v)) > T_abs(resid(v_start))) {
            v_start = v;
          }
        }
      }
      if (v_start == -1) {
        break;
      }
      std::vector<int> prev_edge(nbCol, -1), prev_vert(nbCol, -1);
      std::vector<uint8_t> visited(nbCol, 0);
      std::vector<int> queue{v_start};
      visited[v_start] = 1;
      size_t pos = 0;
      int v_end = -1;
      while (pos < queue.size() && v_end == -1) {
        int v = queue[pos];
        pos += 1;
        for (auto &e : adj[v]) {
          if (!visited[e.first]) {
            visited[e.first] = 1;
            prev_edge[e.first] = e.second;
            prev_vert[e.first] = v;
            queue.push_back(e.first);
            if (resid(e.first) != 0) {
              v_end = e.first;
              break;
            }
          }
        }
      }
      if (v_end == -1) {
        return {SparseSolveStatus::Inconsistent, empty};
      }
      std::vector<int> l_edge;
      std::vector<int> l_node;
      int v = v_end;
      while (v != v_start) {
        l_edge.push_back(prev_edge[v]);
        v = prev_vert[v];
        l_node.push_back(v);
      }
      std::unordered_map<int, int> map_edge, map_node;
      for (size_t u = 0; u < l_edge.size(); u++) {
        map_edge[l_edge[u]] = u;
        map_node[l_node[u]] = u;
      }
      MyMatrix<T> Mpath = ZeroMatrix<T>(l_edge.size(), l_node.size());
      for (auto &t : tripletList) {
        auto iter = map_edge.find(t.row());
        if (iter != map_edge.end()) {
          auto iter2 = map_node.find(t.col());
          if (iter2 != map_node.end()) {
            Mpath(iter->second, iter2->second) += t.value();
          }
        }
      }
      MyVector<T> target = ZeroVector<T>(l_node.size());
      target(map_node[v_start]) = resid(v_start);
      std::optional<MyVector<T>> opt_path = SolutionMat(Mpath, target);
      if (!opt_path) {
        return {SparseSolveStatus::Aborted, empty};
      }
      for (size_t u = 0; u < l_edge.size(); u++) {
        sol(l_edge[u]) += (*opt_path)(u);
      }
      for (auto &t : tripletList) {
        auto iter = map_edge.find(t.row());
        if (iter != map_edge.end()) {
          resid(t.col()) -= (*opt_path)(iter->second) * t.value();
        }
      }
      n_path += 1;
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: solve_chain_system, flow with " << n_path << " paths, n_nz=" << count_nz(sol) << "\n";
#endif
    return {SparseSolveStatus::Solved, sol};
  }
  bool small = (tripletList.size() <= budget.max_nnz_modular);
  bool known_consistent = false;
  if (small) {
    SparseSolveStatus status = SparseSystemConsistent_Mod_Budget(SpMat, Bvect, budget.max_entries_modular, os);
    if (status == SparseSolveStatus::Inconsistent) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: solve_chain_system nbRow=" << nbRow << " nbCol=" << nbCol << " inconsistent (mod p)\n";
#endif
      return {SparseSolveStatus::Inconsistent, empty};
    }
    known_consistent = (status == SparseSolveStatus::Solved);
  }
  // The L1 solution and the exact solution on its support
  MySparseMatrix<double> SpMat_d = UniversalSparseMatrixConversion<double, T>(SpMat);
  MyVector<double> Bvect_d = UniversalVectorConversion<double, T>(Bvect);
  MyVector<double> x_out = l1_solution_reweighted(SpMat_d, Bvect_d, budget.max_yall1_iter, budget.n_reweight, os);
  std::vector<std::vector<std::pair<int, T>>> row_entries(nbRow);
  for (auto &t : tripletList) {
    row_entries[t.row()].emplace_back(t.col(), t.value());
  }
  std::vector<double> l_thr{1e-4, 1e-6, 1e-8};
  for (auto &thr : l_thr) {
    std::vector<int> l_rows;
    for (int i_row = 0; i_row < nbRow; i_row++) {
      if (std::abs(x_out(i_row)) > thr) {
        l_rows.push_back(i_row);
      }
    }
    std::vector<int> col_index(nbCol, -1);
    std::vector<int> l_cols;
    for (auto &i_row : l_rows) {
      for (auto &pair : row_entries[i_row]) {
        if (col_index[pair.first] == -1) {
          col_index[pair.first] = l_cols.size();
          l_cols.push_back(pair.first);
        }
      }
    }
    bool is_ok = true;
    for (int j = 0; j < nbCol; j++) {
      if (col_index[j] == -1 && Bvect(j) != 0) {
        is_ok = false;
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: solve_chain_system thr=" << thr << " support=" << l_rows.size()
       << " touched columns=" << l_cols.size() << " is_ok=" << is_ok << "\n";
#endif
    if (!is_ok) {
      continue;
    }
    MyVector<T> Bred(l_cols.size());
    for (size_t v = 0; v < l_cols.size(); v++) {
      Bred(v) = Bvect(l_cols[v]);
    }
    std::vector<Eigen::Triplet<T>> trip_red;
    for (size_t u = 0; u < l_rows.size(); u++) {
      for (auto &pair : row_entries[l_rows[u]]) {
        trip_red.emplace_back(u, col_index[pair.first], pair.second);
      }
    }
    MySparseMatrix<T> SpRed(l_rows.size(), l_cols.size());
    SpRed.setFromTriplets(trip_red.begin(), trip_red.end());
#ifdef TIMINGS_HECKE_OPERATORS
    MicrosecondTime time_exact;
#endif
    SparseSolveResult<T> res_red = SparseSolutionMat_Exact_Budget(SpRed, Bred, budget.max_entries_exact, os);
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: exact solve on support " << l_rows.size() << " x " << l_cols.size()
       << " status=" << static_cast<int>(res_red.status) << "|=" << time_exact << "\n";
#endif
    if (res_red.status == SparseSolveStatus::Solved) {
      MyVector<T> sol = ZeroVector<T>(nbRow);
      for (size_t u = 0; u < l_rows.size(); u++) {
        sol(l_rows[u]) = res_red.sol(u);
      }
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: solve_chain_system nbRow=" << nbRow << " nbCol=" << nbCol
         << " L1 support=" << l_rows.size() << " thr=" << thr << " n_nz=" << count_nz(sol) << "\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
      os << "|HECKE: solve_chain_system total|=" << time_solve << "\n";
#endif
      return {SparseSolveStatus::Solved, sol};
    }
    if (static_cast<int>(l_rows.size()) == nbRow) {
      break;
    }
  }
  if (known_consistent) {
    // Consistent but the L1 support does not carry a solution: the
    // rational elimination, within its budget.
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: solve_chain_system nbRow=" << nbRow << " nbCol=" << nbCol << " exact elimination fallback\n";
#endif
    return SparseSolutionMat_Exact_Budget(SpMat, Bvect, budget.max_entries_exact, os);
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: solve_chain_system nbRow=" << nbRow << " nbCol=" << nbCol << " undecided within the budget\n";
#endif
  return {SparseSolveStatus::Aborted, empty};
}

/*
  The generic potential used by the peeling: a positive definite form of
  the T-space, the Gram matrix of the first perfect form slightly
  perturbed (so that the vectors of a cell have distinct heights).
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> get_peeling_potential(FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  MyMatrix<T> Q0 = fce.pctdi.l_perfect[0].gram;
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  uint64_t seed = 987654321;
  auto next_random = [&]() -> int {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((seed >> 33) % 11) - 5;
  };
  for (int i_try = 0; i_try < 20; i_try++) {
    MyMatrix<T> Q = Q0;
    T scale = T(1) / T(1000 * (i_try + 1));
    for (auto &M : ListMat) {
      int r = next_random();
      Q += (scale * T(r)) * M;
    }
    if (IsPositiveDefinite(Q, os)) {
      return Q;
    }
  }
  return Q0;
}

/*
  The peeling solver: the chain c (level k-1) is reduced by repeatedly
  taking its highest cell G (for the potential, the largest height of a
  vector of G) and subtracting the multiple of the coboundary of a level
  k face E of G that cancels G; the face is the one that gives the
  lexicographically smallest (maximal height, number of cells at that
  height, number of cells). Each accepted move decreases this objective,
  so the chain is pushed towards the low cells until it vanishes; the
  accumulated faces form the filling. The move fails cleanly when no
  face improves the objective or when a budget is exceeded.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
peel_partial_contracting_homotopy(int const &k, std::vector<PerfectFaceEntry<T, Tint>> const &rhs,
                          UpperFaceCache<T, Tint, Tgroup> &ufc, MyMatrix<T> const &Q0,
                          HeckeBudget const &budget, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
  auto cell_height = [&](int const &index, int const &iOrb, MyMatrix<Tint> const &M) -> T {
    MyMatrix<Tint> EXT = fce.levels[index].l_faces[iOrb].EXT * M;
    T h(0);
    for (int i_row = 0; i_row < EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(EXT, i_row);
      T val = EvaluationQuadForm<T, Tint>(Q0, V);
      if (i_row == 0 || val > h) {
        h = val;
      }
    }
    return h;
  };
  struct Objective {
    T max_height;
    size_t n_top;
    size_t n_cell;
    bool operator<(Objective const &o) const {
      if (max_height != o.max_height) {
        return max_height < o.max_height;
      }
      if (n_top != o.n_top) {
        return n_top < o.n_top;
      }
      return n_cell < o.n_cell;
    }
  };
  // The coefficient of the cell (iOrb, M) in a chain.
  auto coefficient_of = [&](std::vector<PerfectFaceEntry<T, Tint>> const &c, int const &index,
                            int const &iOrb, MyMatrix<Tint> const &M) -> T {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    MyMatrix<Tint> EXT_M = face.EXT * M;
    MyMatrix<Tint> key = tot_set(EXT_M);
    for (auto &fe : c) {
      if (fe.iOrb == iOrb) {
        MyMatrix<Tint> EXT_fe = face.EXT * fe.M;
        MyMatrix<Tint> key_fe = tot_set(EXT_fe);
        if (key_fe == key) {
          MyMatrix<Tint> t = fe.M * Inverse(M);
          int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
          return fe.value * T(sign);
        }
      }
    }
    return T(0);
  };
  std::vector<PerfectFaceEntry<T, Tint>> c = rhs;
  std::vector<PerfectFaceEntry<T, Tint>> filling;
  // The heights of the cells of c, the multiset of heights and the map
  // from the cells (as sets of vectors) to their position, refreshed
  // after each accepted move so that the candidates are evaluated
  // incrementally.
  std::vector<T> h_c;
  std::multiset<T> heights;
  std::unordered_map<MyMatrix<Tint>, size_t> map_c;
  auto refresh = [&]() -> void {
    h_c.clear();
    heights.clear();
    map_c.clear();
    for (size_t i = 0; i < c.size(); i++) {
      T h = cell_height(k - 1, c[i].iOrb, c[i].M);
      h_c.push_back(h);
      heights.insert(h);
      MyMatrix<Tint> EXT_c = fce.levels[k - 1].l_faces[c[i].iOrb].EXT * c[i].M;
      map_c[tot_set(EXT_c)] = i;
    }
  };
  auto get_top = [&]() -> int {
    int i_top = -1;
    for (size_t i = 0; i < c.size(); i++) {
      if (i_top == -1 || h_c[i] > h_c[i_top]) {
        i_top = i;
      }
    }
    return i_top;
  };
  auto objective_of_heights = [&](std::multiset<T> const &hs) -> Objective {
    if (hs.empty()) {
      return Objective{T(0), 0, 0};
    }
    T top = *hs.rbegin();
    return Objective{top, hs.count(top), hs.size()};
  };
  refresh();
  Objective obj = objective_of_heights(heights);
  size_t n_move = 0;
  while (!c.empty()) {
    if (n_move >= budget.max_peel_moves) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, budget of moves exhausted, |c|=" << c.size() << "\n";
#endif
      return {};
    }
#ifdef DEBUG_HECKE_OPERATORS
    if (n_move % 100 == 0 && n_move > 0) {
      os << "HECKE: peel, move " << n_move << " |c|=" << c.size() << " max_height=" << obj.max_height << "\n";
    }
#endif
    int i_top = get_top();
    PerfectFaceEntry<T, Tint> const G = c[i_top];
    std::optional<Objective> best_obj;
    std::vector<PerfectFaceEntry<T, Tint>> best_cob;
    PerfectFaceEntry<T, Tint> best_E{0, T(0), G.M};
    for (auto &e : fce.boundaries[k - 1].ll_bound[G.iOrb].l_bound) {
      if (!fce.levels[k].l_faces[e.iOrb].is_well_rounded) {
        continue;
      }
      MyMatrix<Tint> Mface = e.M * G.M;
      std::vector<PerfectFaceEntry<T, Tint>> chain_E{{e.iOrb, T(1), Mface}};
      std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(k, chain_E, ufc, os);
      T s_G = coefficient_of(cob, k - 1, G.iOrb, G.M);
      if (s_G == 0) {
        continue;
      }
      T eps = G.value / s_G;
      // The heights after the move: G disappears, the touched cells of c
      // may vanish, the other cofaces of E appear.
      std::multiset<T> hs = heights;
      hs.erase(hs.find(h_c[i_top]));
      for (auto &fe : cob) {
        FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[k - 1].l_faces[fe.iOrb];
        MyMatrix<Tint> EXT_fe = face.EXT * fe.M;
        MyMatrix<Tint> key = tot_set(EXT_fe);
        auto iter = map_c.find(key);
        if (iter != map_c.end()) {
          size_t idx = iter->second;
          if (static_cast<int>(idx) == i_top) {
            continue;
          }
          MyMatrix<Tint> t = c[idx].M * Inverse(fe.M);
          int sign = get_face_orientation(face.EXT, ListMat, face.or_info, t);
          T new_coef = c[idx].value - eps * fe.value * T(sign);
          if (new_coef == 0) {
            hs.erase(hs.find(h_c[idx]));
          }
        } else {
          T h = cell_height(k - 1, fe.iOrb, fe.M);
          hs.insert(h);
        }
      }
      Objective obj_new = objective_of_heights(hs);
      if (!best_obj || obj_new < *best_obj) {
        best_obj = obj_new;
        best_cob = cob;
        best_E = PerfectFaceEntry<T, Tint>{e.iOrb, eps, Mface};
      }
    }
    if (!best_obj || !(*best_obj < obj)) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, no improving move at move " << n_move << " |c|=" << c.size()
         << " max_height=" << obj.max_height << " n_top=" << obj.n_top << "\n";
#endif
      return {};
    }
    ChainBuilder<T, Tint, Tgroup> cb(k - 1, fce, os);
    for (auto &fe : c) {
      cb.f_insert(fe.value, fe.iOrb, fe.M);
    }
    for (auto &fe : best_cob) {
      cb.f_insert(-best_E.value * fe.value, fe.iOrb, fe.M);
    }
    c = cb.get_faces();
    if (c.size() > budget.max_chain) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: peel, chain larger than the budget\n";
#endif
      return {};
    }
    filling.push_back(best_E);
    refresh();
    obj = objective_of_heights(heights);
    n_move += 1;
  }
  // The filling, with the duplicates merged
  ChainBuilder<T, Tint, Tgroup> cb(k, fce, os);
  for (auto &fe : filling) {
    cb.f_insert(fe.value, fe.iOrb, fe.M);
  }
  std::vector<PerfectFaceEntry<T, Tint>> sol = cb.get_faces();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: peel, solved in " << n_move << " moves, |sol|=" << sol.size() << "\n";
#endif
  return sol;
}

/*
  The dual partial contracting homotopy: given a chain b of level k-1 (which must
  be a coboundary), find a chain c of level k with coboundary(c) = b.
  The methods are tried in turn: the peeling, then the linear systems on
  the rings of unknowns around b (the faces of its cells first, then the
  cells sharing a coface with a current unknown, and so on) with all the
  cofaces of the unknowns as equations. Each method has its budget and
  the failure is reported cleanly (empty optional).
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
dual_partial_contracting_homotopy(int const &k, std::vector<PerfectFaceEntry<T, Tint>> const &rhs,
                          UpperFaceCache<T, Tint, Tgroup> &ufc, MyMatrix<T> const &Q0,
                          HeckeBudget const &budget, std::ostream &os) {
  FullComplexEnumeration<T, Tint, Tgroup> const &fce = ufc.fce;
  if (k > 1) {
    std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_peel =
        peel_partial_contracting_homotopy(k, rhs, ufc, Q0, budget, os);
    if (opt_peel) {
      return opt_peel;
    }
  }
  ComplexBuilder<T, Tint, Tgroup> cb_k(k, fce, os);
  auto insert_faces_of = [&](int const &iOrb_up, MyMatrix<Tint> const &M_up) -> void {
    for (auto &e : fce.boundaries[k - 1].ll_bound[iOrb_up].l_bound) {
      if (fce.levels[k].l_faces[e.iOrb].is_well_rounded) {
        MyMatrix<Tint> M = e.M * M_up;
        cb_k.f_insert(e.iOrb, M);
      }
    }
  };
  for (auto &fe : rhs) {
    insert_faces_of(fe.iOrb, fe.M);
  }
  size_t start = 0;
  for (size_t ring = 0; ring < budget.max_rings; ring++) {
    std::vector<PerfectFace<Tint>> faces_k = cb_k.get_faces();
    if (faces_k.size() > budget.max_unknowns) {
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: dual_partial_contracting_homotopy k=" << k << " ring=" << ring
         << " unknowns=" << faces_k.size() << " beyond the budget\n";
#endif
      return {};
    }
    ComplexBuilder<T, Tint, Tgroup> cb_km1(k - 1, fce, os);
    using T2 = Eigen::Triplet<T>;
    std::vector<T2> tripletList;
    int nbRow = faces_k.size();
    for (int u = 0; u < nbRow; u++) {
      PerfectFace<Tint> const &F = faces_k[u];
      for (auto &pair : ufc.get_cofaces(k, F.iOrb)) {
        MyMatrix<Tint> M = pair.first.M * F.M;
        cb_km1.f_insert(pair.first.iOrb, M);
        std::pair<size_t, int> p = cb_km1.get_position(pair.first.iOrb, M);
        T value = T(pair.second * p.second);
        tripletList.emplace_back(u, p.first, value);
      }
    }
    int nbCol = cb_km1.get_faces().size();
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: dual_partial_contracting_homotopy k=" << k << " ring=" << ring
       << " unknowns=" << nbRow << " equations=" << nbCol << "\n";
#endif
    MyVector<T> Bvect = ZeroVector<T>(nbCol);
    bool rhs_present = true;
    for (auto &fe : rhs) {
      std::optional<std::pair<size_t, int>> opt_p = cb_km1.find_position(fe.iOrb, fe.M);
      if (!opt_p) {
        rhs_present = false;
        break;
      }
      Bvect(opt_p->first) += T(opt_p->second) * fe.value;
    }
    bool is_graph = (k == 1);
    if (rhs_present) {
      SparseSolveResult<T> res = solve_chain_system(tripletList, nbRow, nbCol, Bvect, is_graph, budget, os);
      if (res.status == SparseSolveStatus::Solved) {
        std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
        for (int u = 0; u < nbRow; u++) {
          T val = res.sol(u);
          if (val != 0) {
            PerfectFace<Tint> const &face = faces_k[u];
            chain_ret.push_back({face.iOrb, val, face.M});
          }
        }
        return chain_ret;
      }
    }
    // The next ring: the cells sharing a coface with a current unknown.
    size_t len = faces_k.size();
    for (size_t u = start; u < len; u++) {
      PerfectFace<Tint> const &F = faces_k[u];
      for (auto &pair : ufc.get_cofaces(k, F.iOrb)) {
        MyMatrix<Tint> M = pair.first.M * F.M;
        insert_faces_of(pair.first.iOrb, M);
      }
      if (cb_k.get_faces().size() > budget.max_unknowns) {
        break;
      }
    }
    start = len;
    if (cb_k.get_faces().size() == len) {
      break;
    }
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: dual_partial_contracting_homotopy k=" << k << " no filling found within the budgets\n";
#endif
  return {};
}

/*
  The image of a perfect cone C under y: a cone whose closure contains
  the image of the barycenter of C, found by the walk. When the point
  lies on a lower face the choice is not canonical; the equivariance is
  restored afterwards by the averaging over the coset stabilizer and the
  propagation along its orbit, as for the other levels. The sign is the
  orientation sign of y and of the matrix of the cone.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
get_cone_image(int const &iOrb, MyMatrix<T> const &y,
               FullComplexEnumeration<T, Tint, Tgroup> const &fce,
               ConeWalker<T, Tint, Tgroup> const &walker, int const &eps_y,
               [[maybe_unused]] std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  MyMatrix<Tint> const &EXT = fce.levels[0].l_faces[iOrb].EXT;
  MyMatrix<T> W = UniversalMatrixConversion<T, Tint>(EXT) * y;
  TopPerfectCone<Tint> tpc = walker.find_containing_cone_T(W);
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
  int eps_M = tspace_orientation_sign(LinSpa, M_T);
  PerfectFaceEntry<T, Tint> fe{tpc.i_perfect, T(eps_y * eps_M), tpc.M};
  return {fe};
}

// The right hand side coboundary(F) mapped by y_j, from the images of
// the level index-1.
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>> get_hecke_coboundary_image(
    int const &index, int const &iOrb, int const &j,
    FullComplexEnumeration<T, Tint, Tgroup> const &fce, HeckeCosets<T> const &cosets,
    std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>> const &images_prev,
    UpperFaceCache<T, Tint, Tgroup> &ufc, std::ostream &os) {
  ChainBuilder<T, Tint, Tgroup> cb(index - 1, fce, os);
  MyMatrix<T> const &y = cosets.l_coset[j];
  for (auto &pair : ufc.get_cofaces(index, iOrb)) {
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(pair.first.M);
    MyMatrix<T> z = M_T * y;
    std::pair<int, MyMatrix<T>> dec = cosets.decompose(z);
    MyMatrix<Tint> M_conj = UniversalMatrixConversion<Tint, T>(dec.second);
    T e_sign = T(pair.second);
    for (auto &fe : images_prev[pair.first.iOrb][dec.first]) {
      T value = fe.value * e_sign;
      MyMatrix<Tint> M = fe.M * M_conj;
      cb.f_insert(value, fe.iOrb, M);
    }
  }
  return cb.get_faces();
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map_dual(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                             MyMatrix<T> const &x, std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  PerfectComplexOptions const &pco = fce.pctdi.pco;
  if (!pco.compute_boundary || !pco.compute_partial_contracting_homotopy) {
    std::cerr << "HECKE: The dual method requires ComputeBoundary = T and ComputePartialContractingHomotopy = T\n";
    throw TerminalException{1};
  }
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
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images(n_levels);
  UpperFaceCache<T, Tint, Tgroup> ufc(fce, os);
  ConeWalker<T, Tint, Tgroup> walker(fce, os);
  HeckeBudget budget;
  MyMatrix<T> Q0 = get_peeling_potential(fce, os);
  auto is_wr = [&](int const &index, int const &iOrb) -> bool {
    return fce.levels[index].l_faces[iOrb].is_well_rounded;
  };
  for (int index = 0; index < n_levels; index++) {
    int n_orb = fce.levels[index].l_faces.size();
    images[index] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      if (!is_wr(index, iOrb)) {
        continue;
      }
      FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
      std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
          get_stab_coset_action(face.l_gens, cosets);
      std::vector<int> status(n_coset, 0);
      for (size_t j = 0; j < n_coset; j++) {
        if (status[j] == 1) {
          continue;
        }
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: dual index=" << index << " iOrb=" << iOrb << "/" << n_orb
           << " j=" << j << "/" << n_coset << "\n";
#endif
        std::vector<PerfectFaceEntry<T, Tint>> rhs;
        std::vector<PerfectFaceEntry<T, Tint>> sol;
        if (index == 0) {
          int eps_y = tspace_orientation_sign(LinSpa, cosets.l_coset[j]);
          sol = get_cone_image(iOrb, cosets.l_coset[j], fce, walker, eps_y, os);
        } else {
          rhs = get_hecke_coboundary_image(index, iOrb, j, fce, cosets, images[index - 1], ufc, os);
          if (!rhs.empty()) {
            std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_sol =
                dual_partial_contracting_homotopy(index, rhs, ufc, Q0, budget, os);
            if (!opt_sol) {
              throw HeckeFailure{"no filling found for the cell of level " + std::to_string(index) +
                                 " orbit " + std::to_string(iOrb) + " coset " + std::to_string(j) +
                                 " (right hand side of " + std::to_string(rhs.size()) + " cells)"};
            }
            sol = *opt_sol;
          }
        }
        CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
        std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
            get_coset_stabilizer_schreier(co, face.l_gens, act, n);
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (index > 0 && !is_chain_coset_stabilizer_invariant(index, iOrb, index - 1, rhs, l_schreier, fce, os)) {
          std::cerr << "HECKE: dual, the right hand side should be invariant under the coset stabilizer\n";
          throw TerminalException{1};
        }
#endif
        std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
            average_chain_coset_stabilizer(index, iOrb, sol, l_schreier, fce, os);
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: dual index=" << index << " iOrb=" << iOrb << " j=" << j
           << " |rhs|=" << rhs.size() << " |sol|=" << sol.size()
           << " |sol_avg|=" << sol_avg.size() << "\n";
#endif
#ifdef SANITY_CHECK_HECKE_OPERATORS
        if (index == 0) {
          // The coset stabilizer fixes the image of the barycenter, so the
          // averaged image is made of cones all containing it.
          MyMatrix<T> P = UniversalMatrixConversion<T, Tint>(face.EXT) * cosets.l_coset[j];
          for (auto &fe : sol_avg) {
            std::optional<int> opt = walker.n_facet_containing(fe.iOrb, fe.M, P);
            if (!opt) {
              std::cerr << "HECKE: dual, a cone of the averaged image of iOrb=" << iOrb
                        << " j=" << j << " does not contain the image of the barycenter\n";
              throw TerminalException{1};
            }
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: dual index=0 iOrb=" << iOrb << " j=" << j
               << " n_facet_containing=" << *opt << "\n";
#endif
          }
#ifdef DEBUG_HECKE_OPERATORS
          // The dimension of the space of forms fixed by the coset stabilizer.
          // If it is 1 then the image of the barycenter is the only fixed point
          // and the averaged image cannot be smaller.
          int n_mat = LinSpa.ListMat.size();
          MyMatrix<T> Equa = ZeroMatrix<T>(n_mat, n_mat * l_schreier.size());
          for (size_t i_s = 0; i_s < l_schreier.size(); i_s++) {
            MyMatrix<T> h = UniversalMatrixConversion<T, Tint>(l_schreier[i_s].second);
            MyMatrix<T> A = tspace_action_matrix(LinSpa, h);
            for (int k = 0; k < n_mat; k++) {
              for (int l = 0; l < n_mat; l++) {
                T val = A(k, l);
                if (k == l) {
                  val -= 1;
                }
                Equa(k, i_s * n_mat + l) = val;
              }
            }
          }
          os << "HECKE: dual index=0 iOrb=" << iOrb << " j=" << j
             << " dim_fixed_forms=" << NullspaceMat(Equa).rows() << "\n";
#endif
        }
        {
          if (index > 0) {
            std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(index, sol_avg, ufc, os);
            if (!is_equal_chain(cob, rhs, index - 1, fce, os)) {
              std::cerr << "HECKE: dual, the averaged chain is not a solution\n";
              throw TerminalException{1};
            }
          }
          if (!is_chain_coset_stabilizer_invariant(index, iOrb, index, sol_avg, l_schreier, fce, os)) {
            std::cerr << "HECKE: dual, the averaged chain should be invariant\n";
            throw TerminalException{1};
          }
        }
#endif
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
      for (size_t j = 0; j < n_coset && index > 0; j++) {
        std::vector<PerfectFaceEntry<T, Tint>> rhs =
            get_hecke_coboundary_image(index, iOrb, j, fce, cosets, images[index - 1], ufc, os);
        std::vector<PerfectFaceEntry<T, Tint>> cob = chain_coboundary(index, images[index][iOrb][j], ufc, os);
        if (!is_equal_chain(cob, rhs, index - 1, fce, os)) {
          std::cerr << "HECKE: dual, the chain map property fails at index=" << index
                    << " iOrb=" << iOrb << " j=" << j << "\n";
          throw TerminalException{1};
        }
      }
#endif
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: dual level " << index << " n_orb=" << n_orb << " done\n";
#endif
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: dual level " << index << "|=" << time << "\n";
#endif
  }
  return {x, std::move(cosets), std::move(images), {}, "Dual", true};
}

//
// The sharbly method: Ash-Rudolph reduction at the symbol level
//

/*
  The sharbly method computes the images at the symbol level: the
  well rounded cells whose cones have dimension n, which carry the top
  degree cohomology H^nu, with nu the virtual cohomological dimension.
  The image of a symbol cell F under a coset y is the modular symbol of
  the rows of EXT * y, rewritten as a sum of unimodular symbols by the
  continued fraction algorithm of Ash and Rudolph: for a symbol
  [v_1, ..., v_n] with |det| > 1 and a lattice point w = sum_i q_i v_i,
  0 <= q_i < 1, not all q_i integral,
    [v_1, ..., v_n] = sum_{q_i <> 0} [v_1, .., v_{i-1}, w, v_{i+1}, .., v_n]
  and each term has determinant |q_i| |det| < |det|. The identity is a
  relation of the Steinberg module, so the output is homologous to the
  input in the well rounded quotient complex: the operator induced on the
  homology of the symbol level is the Hecke operator. The level below
  is obtained by filling (see compute_sharbly_level_below) when the
  coverings are available (ComputePartialContractingHomotopy = T). The
  images of the other levels are left out (is_level_computed); the
  homology level sanity checks (the image of a cycle is a cycle and
  projects onto the harmonic space modulo boundaries) verify the result.

  The reducing points w are the lattice points of the box
  {sum_i q_i v_i : |q_i| < 1} minimizing max |q_i| (the largest new
  determinant), and the reduction is averaged over them. The box only
  depends on the symbol, so the reduction Phi is a G-equivariant function
  of the symbol: Phi(W g) = Phi(W) g. This makes the images of all the
  symbols consistent, which the level below relies on; the averaging over
  the coset stabilizer is kept as a safety net and the invariance is
  checked.

  A unimodular symbol is a cell of the unimodular frame orbit iOrb_0 (the
  orbit of the standard frame; GL_n(Z) is transitive on the frames). The
  orientation character of a frame cell on its stabilizer is the sign of
  the permutation of the rays, which is also the sign of the symbol
  relations, so inserting the symbol W as the cell (iOrb_0, E_0^{-1} W)
  is sign-consistent across all insertions.
 */

// The symbol level: the level of the cones of dimension n, which is
// dim T - n since level k holds the cones of codimension k.
template <typename T, typename Tint, typename Tgroup>
int get_symbol_level(FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  int dimT = fce.pctdi.LinSpa.ListMat.size();
  int n = fce.pctdi.LinSpa.n;
  return dimT - n;
}

// One vector per pair {v, -v} of the rows, each with its first non-zero
// entry positive. The EXT of a face holds both v and -v.
template <typename Tint>
MyMatrix<Tint> half_vector_set(MyMatrix<Tint> const &EXT) {
  std::set<MyVector<Tint>> set;
  for (int i_row = 0; i_row < EXT.rows(); i_row++) {
    MyVector<Tint> V = GetMatrixRow(EXT, i_row);
    set.insert(CanonicalizeVectorToInvertible(V));
  }
  MyMatrix<Tint> M(set.size(), EXT.cols());
  int pos = 0;
  for (auto &V : set) {
    AssignMatrixRow(M, pos, V);
    pos += 1;
  }
  return M;
}

template <typename T, typename Tint, typename Tgroup>
struct SharblyReducer {
private:
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  int symbol_level;
  int iOrb0;
  MyMatrix<Tint> E0_inv;
  HeckeBudget budget;
  std::ostream &os;
  // x - floor(x), in [0, 1).
  static T frac(T const &x) {
    Tint fl = UniversalFloorScalarInteger<Tint, T>(x);
    return x - T(fl);
  }
  /*
    The reducing points of the symbol of the rows of W: the nonzero lattice
    points w = sum_i q_i v_i with all |q_i| < 1 for which max_i |q_i| is the
    smallest, then the number of nonzero q_i, then sum_i q_i^2. They are returned as the
    vectors q, one for each pair {w, -w}. The lattice points of this box
    only depend on the symbol (the box is invariant under the permutations
    and the sign changes of the rows), so the reduction, which averages
    over all of them, is G-equivariant: Phi(W g) = Phi(W) g.
   */
  std::vector<MyVector<T>> get_reducing_points(MyMatrix<T> const &W_inv) const {
    int n = W_inv.rows();
    // The coefficients modulo 1 of the lattice points form the group
    // Z^n / L(W) of order |det W|, generated by those of the e_j, which are
    // the rows of W^{-1}.
    std::vector<MyVector<T>> l_gen;
    for (int j = 0; j < n; j++) {
      MyVector<T> g(n);
      for (int i = 0; i < n; i++) {
        g(i) = frac(W_inv(j, i));
      }
      l_gen.push_back(g);
    }
    std::set<MyVector<T>> set_class;
    std::vector<MyVector<T>> l_class;
    MyVector<T> zero = ZeroVector<T>(n);
    set_class.insert(zero);
    l_class.push_back(zero);
    size_t pos = 0;
    while (pos < l_class.size()) {
      MyVector<T> f = l_class[pos];
      for (auto &g : l_gen) {
        MyVector<T> h(n);
        for (int i = 0; i < n; i++) {
          h(i) = frac(f(i) + g(i));
        }
        if (set_class.insert(h).second) {
          l_class.push_back(h);
        }
      }
      pos += 1;
    }
    // In a class f, the smallest |q_i| is min(f_i, 1 - f_i). The classes
    // are ordered by max |q_i|, then the number of nonzero q_i, then
    // sum q_i^2 = w (W^T W)^{-1} w^T, all invariant under G.
    T half = T(1) / T(2);
    std::optional<std::tuple<T, int, T>> best;
    std::vector<MyVector<T>> l_best;
    for (auto &f : l_class) {
      if (IsZeroVector(f)) {
        continue;
      }
      T score(0);
      T sum_sqr(0);
      int n_nonzero = 0;
      for (int i = 0; i < n; i++) {
        if (f(i) != 0) {
          n_nonzero += 1;
          T a = (f(i) <= half) ? f(i) : T(1) - f(i);
          sum_sqr += a * a;
          if (a > score) {
            score = a;
          }
        }
      }
      std::tuple<T, int, T> key{score, n_nonzero, sum_sqr};
      if (!best || key < *best) {
        best = key;
        l_best.clear();
      }
      if (key == *best) {
        l_best.push_back(f);
      }
    }
    // The vectors q of the classes: f_i = 1/2 gives the two choices +-1/2.
    std::set<MyVector<T>> set_q;
    for (auto &f : l_best) {
      std::vector<MyVector<T>> l_q{ZeroVector<T>(n)};
      for (int i = 0; i < n; i++) {
        std::vector<MyVector<T>> l_q_new;
        for (auto &q : l_q) {
          auto f_push = [&](T const &val) -> void {
            MyVector<T> q2 = q;
            q2(i) = val;
            l_q_new.push_back(q2);
          };
          if (f(i) == half) {
            f_push(half);
            f_push(-half);
          } else if (f(i) < half) {
            f_push(f(i));
          } else {
            f_push(f(i) - T(1));
          }
        }
        l_q = std::move(l_q_new);
      }
      for (auto &q : l_q) {
        // w and -w give the same symbols: one of them is kept.
        MyVector<T> q_can = q;
        for (int i = 0; i < n; i++) {
          if (q_can(i) != 0) {
            if (q_can(i) < 0) {
              q_can = -q_can;
            }
            break;
          }
        }
        set_q.insert(q_can);
      }
    }
    return std::vector<MyVector<T>>(set_q.begin(), set_q.end());
  }
public:
  SharblyReducer(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, int const &_symbol_level,
                 int const &_iOrb0, MyMatrix<Tint> const &_E0_inv, std::ostream &_os)
      : fce(_fce), symbol_level(_symbol_level), iOrb0(_iOrb0), E0_inv(_E0_inv), os(_os) {
  }
  // The reducing points of the symbol of the rows of W (|det W| > 1), as
  // primitive integral vectors.
  std::vector<MyVector<Tint>> reducing_points(MyMatrix<Tint> const &W) const {
    MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
    MyMatrix<T> W_inv = Inverse(W_T);
    std::vector<MyVector<Tint>> l_w;
    for (auto &q : get_reducing_points(W_inv)) {
      MyVector<T> w_T = W_T.transpose() * q;
      l_w.push_back(UniversalVectorConversion<Tint, T>(RemoveFractionVector(w_T)));
    }
    return l_w;
  }
  // Make the rows primitive integral vectors: the symbol is invariant
  // under the scaling of its rows.
  MyMatrix<Tint> primitive_rows(MyMatrix<T> const &W_T) const {
    int n_row = W_T.rows();
    int n_col = W_T.cols();
    MyMatrix<Tint> W(n_row, n_col);
    for (int i_row = 0; i_row < n_row; i_row++) {
      MyVector<T> V = GetMatrixRow(W_T, i_row);
      MyVector<T> Vred = RemoveFractionVector(V);
      MyVector<Tint> Vint = UniversalVectorConversion<Tint, T>(Vred);
      AssignMatrixRow(W, i_row, Vint);
    }
    return W;
  }
  // The Ash-Rudolph reduction of the symbol of the rows of W, written
  // as a chain of cells of the symbol level. If l_sample is not null, the
  // symbols met and the 1-sharblies [w, v_1, ..., v_n] of the steps are
  // appended to it: they form a filling of the difference between the
  // symbol and its reduction in the sharbly complex.
  std::vector<PerfectFaceEntry<T, Tint>> reduce(MyMatrix<Tint> const &W_input,
                                                std::vector<MyMatrix<Tint>> *l_sample = nullptr) const {
    ChainBuilder<T, Tint, Tgroup> cb(symbol_level, fce, os);
    std::vector<std::pair<T, MyMatrix<Tint>>> l_work;
    l_work.push_back({T(1), W_input});
    size_t n_symbol = 0;
    while (!l_work.empty()) {
      n_symbol += 1;
      if (n_symbol > budget.max_sharbly_symbols) {
        throw HeckeFailure{"sharbly: more than " + std::to_string(budget.max_sharbly_symbols) +
                           " symbols in one reduction"};
      }
      std::pair<T, MyMatrix<Tint>> pair = l_work.back();
      l_work.pop_back();
      T const &coef = pair.first;
      MyMatrix<Tint> const &W = pair.second;
      Tint det = DeterminantMat(W);
      if (det == 0) {
        // A degenerate symbol vanishes in the quotient complex.
        continue;
      }
      if (T_abs(det) == 1) {
        MyMatrix<Tint> M = E0_inv * W;
        cb.f_insert(coef, iOrb0, M);
        continue;
      }
      MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
      MyMatrix<T> W_inv = Inverse(W_T);
      std::vector<MyVector<T>> l_q = get_reducing_points(W_inv);
      if (l_q.empty()) {
        throw HeckeFailure{"sharbly: no reducing point for a symbol of determinant > 1"};
      }
      int n_q = l_q.size();
      T coef_q = coef / T(n_q);
      int n = W.rows();
      if (l_sample) {
        l_sample->push_back(W);
      }
      for (auto &q : l_q) {
        MyVector<T> w_T = W_T.transpose() * q;
        MyVector<Tint> w = UniversalVectorConversion<Tint, T>(RemoveFractionVector(w_T));
        if (l_sample) {
          MyMatrix<Tint> Wext(n + 1, n);
          AssignMatrixRow(Wext, 0, w);
          for (int i = 0; i < n; i++) {
            AssignMatrixRow(Wext, i + 1, GetMatrixRow(W, i));
          }
          l_sample->push_back(std::move(Wext));
        }
        for (int i = 0; i < n; i++) {
          if (q(i) != 0) {
            MyMatrix<Tint> W_new = W;
            AssignMatrixRow(W_new, i, w);
            l_work.push_back({coef_q, std::move(W_new)});
          }
        }
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: sharbly reduce n_symbol=" << n_symbol << " |chain|=" << cb.chain.size() << "\n";
#endif
    return cb.get_faces();
  }
};

// The orbit of the unimodular frames at the symbol level. GL_n(Z) is
// transitive on the frames so there is exactly one orbit; the standard
// frame is a face of the perfect cone of A_n, so the orbit is a cell.
template <typename T, typename Tint, typename Tgroup>
int get_unimodular_frame_orbit(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                               int const &symbol_level,
                               [[maybe_unused]] std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  std::vector<int> l_found;
  int n_orb = fce.levels[symbol_level].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[symbol_level].l_faces[iOrb];
    // The EXT of a face holds both v and -v: the frame is the half set.
    MyMatrix<Tint> half = half_vector_set(face.EXT);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: sharbly symbol_level=" << symbol_level << " iOrb=" << iOrb
       << " |half|=" << half.rows() << " wr=" << face.is_well_rounded << "\n";
#endif
    if (face.is_well_rounded && half.rows() == n) {
      Tint det = DeterminantMat(half);
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: sharbly iOrb=" << iOrb << " det=" << det << "\n";
#endif
      if (T_abs(det) == 1) {
        l_found.push_back(iOrb);
      }
    }
  }
  if (l_found.size() != 1) {
    throw HeckeFailure{"sharbly: found " + std::to_string(l_found.size()) +
                       " orbits of unimodular frames at the symbol level instead of 1"};
  }
  return l_found[0];
}

/*
  The level below the symbol level (degree nu - 1 of the cohomology, the
  1-sharblies). The well rounded quotient complex (the cells modulo the
  cells which are not well rounded) is the compactly supported cochain
  complex of the well rounded retract W. Since Gamma is a duality group
  with dualizing module the Steinberg module, it is exact except at the
  symbol level, where its homology is St. For a cell F of the level below
  and a coset y, the chain R = Psi(d F, y) of the symbol level represents
  [d F] y = 0 in St, so R is a boundary in the quotient complex and the
  image of F is a filling X of it: d X = R modulo the cells which are not
  well rounded. Then d Psi = Psi d on the level below, which is all that
  is needed for the operator on its homology: the extension to the lower
  levels exists by exactness and does not change the induced map.

  The filling is searched by the covering method among the well rounded
  cells, which for the classic T-space are the cells whose vectors span
  R^n. The solution is averaged over the coset stabilizer as for the
  other methods.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
chain_boundary_well_rounded(int const &index, std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                            FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary(index, chain, fce, os);
  std::vector<PerfectFaceEntry<T, Tint>> bnd_wr;
  for (auto &fe : bnd) {
    if (fce.levels[index + 1].l_faces[fe.iOrb].is_well_rounded) {
      bnd_wr.push_back(fe);
    }
  }
  return bnd_wr;
}

/*
  The perfect cones containing points of the simplex of the rank one forms
  of the rows of P: the points weighted towards a subset of its vertices
  (the vectors of the subset doubled, so their rank one forms have weight
  4), all the subsets for the small simplices, the singletons otherwise.
 */
template <typename T, typename Tint, typename Tgroup>
void locate_simplex_cones(MyMatrix<Tint> const &P, ConeWalker<T, Tint, Tgroup> const &walker,
                          std::vector<TopPerfectCone<Tint>> &l_top) {
  int max_row_all_subsets = 6;
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  int n_row = P_T.rows();
  auto f_locate = [&](std::vector<int> const &l_row) -> void {
    MyMatrix<T> Q_T = P_T;
    for (auto &i_row : l_row) {
      for (int i_col = 0; i_col < Q_T.cols(); i_col++) {
        Q_T(i_row, i_col) *= 2;
      }
    }
    l_top.push_back(walker.find_containing_cone_T(Q_T));
  };
  f_locate({});
  if (n_row <= max_row_all_subsets) {
    for (int mask = 1; mask < (1 << n_row) - 1; mask++) {
      std::vector<int> l_row;
      for (int i_row = 0; i_row < n_row; i_row++) {
        if ((mask >> i_row) & 1) {
          l_row.push_back(i_row);
        }
      }
      f_locate(l_row);
    }
  } else {
    for (int i_row = 0; i_row < n_row; i_row++) {
      f_locate({i_row});
    }
  }
}

/*
  The seed of the covering for the filling of the image of F under y: the
  perfect cones containing the points of a filling in the sharbly complex,
  namely the cone F y with its facets and the symbols and the 1-sharblies
  met by the Ash-Rudolph reduction of the facets. The covering made of the
  cones containing the cells of the right hand side alone is too far from
  the filling.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<TopPerfectCone<Tint>>
get_sharbly_covering_seed(int const &index, int const &iOrb, MyMatrix<T> const &y,
                          FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                          ConeWalker<T, Tint, Tgroup> const &walker,
                          SharblyReducer<T, Tint, Tgroup> const &reducer) {
  std::vector<MyMatrix<Tint>> l_sample;
  MyMatrix<Tint> const &EXT_F = fce.levels[index].l_faces[iOrb].EXT;
  l_sample.push_back(reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(half_vector_set(EXT_F)) * y));
  for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
    FacePerfectComplex<T, Tint, Tgroup> const &face_G = fce.levels[index + 1].l_faces[ebnd.iOrb];
    if (!face_G.is_well_rounded) {
      continue;
    }
    MyMatrix<Tint> EXT_G = half_vector_set(face_G.EXT) * ebnd.M;
    MyMatrix<Tint> W = reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(EXT_G) * y);
    l_sample.push_back(W);
    (void)reducer.reduce(W, &l_sample);
  }
  std::vector<TopPerfectCone<Tint>> l_top;
  for (auto &P : l_sample) {
    locate_simplex_cones(P, walker, l_top);
  }
  return l_top;
}

/*
  The local filling of the image of a cell F of the level below the symbol
  level under a coset y (the 1-sharbly reduction).

  The 1-sharblies D = [u_0, ..., u_n] (n+1 vectors spanning R^n) are
  mapped to the chains Phi(dD) of the symbol level, with Phi the reduction
  of the symbols and dD = sum_k (-1)^k [u_0, .., u_k^, .., u_n]. Two
  identities are used:
  - The right hand side R of F y is sum_S c_S Phi(dD_S) for the
    1-sharblies D_S made of n+1 of the vectors of F y; the coefficients
    c_S are found by a small exact system, which also settles the
    orientations.
  - For any vector w, D = sum_k D[u_k <- w] + d[w, u_0, ..., u_n], so
    Phi(dD) = sum_k Phi(dD[u_k <- w]): a 1-sharbly can be split.
  A filling of R is then the sum of fillings of the pieces. A piece whose
  n-subsets have |det| <= 1 (its faces are reduced symbols, so Phi leaves
  them in place) is filled on a covering seeded by the cones met by its
  simplex. Otherwise it is split along its worst face sigma with the
  reducing points of sigma used by Phi: this follows the reduction of
  sigma by Phi one degree up, as the lifts of Gunnells do. The splitting
  is not equivariant, which does not matter since the result is averaged
  afterwards. Its termination is not proven (as for Gunnells), so the
  depth and the number of pieces are bounded, a piece at the maximal
  depth being filled directly.
 */
template <typename T, typename Tint, typename Tgroup>
struct SharblyLocalFiller {
private:
  int index;
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  ConeWalker<T, Tint, Tgroup> const &walker;
  SharblyReducer<T, Tint, Tgroup> const &reducer;
  std::function<bool(MyMatrix<Tint> const &)> f_filter;
  CoveringOptions options;
  std::ostream &os;
  // The pieces are filled from the depth max_depth, and split further up to
  // max_depth_split when their covering is beyond the budget.
  int max_depth = 1;
  int max_depth_split = 3;
  size_t max_piece = 2000;
  size_t n_piece = 0;
  std::unordered_map<MyMatrix<Tint>, std::vector<PerfectFaceEntry<T, Tint>>> cache_phi;
  // The rows of D except the row k.
  static MyMatrix<Tint> remove_row(MyMatrix<Tint> const &D, int const &k) {
    MyMatrix<Tint> W(D.rows() - 1, D.cols());
    int pos = 0;
    for (int i = 0; i < D.rows(); i++) {
      if (i != k) {
        AssignMatrixRow(W, pos, GetMatrixRow(D, i));
        pos += 1;
      }
    }
    return W;
  }
  std::vector<PerfectFaceEntry<T, Tint>> const &phi(MyMatrix<Tint> const &W) {
    auto iter = cache_phi.find(W);
    if (iter != cache_phi.end()) {
      return iter->second;
    }
    std::vector<PerfectFaceEntry<T, Tint>> chain = reducer.reduce(W);
    return cache_phi.emplace(W, std::move(chain)).first->second;
  }
  // The largest |det| of the n-subsets of the rows of D, and their sum.
  static std::pair<Tint, Tint> det_size(MyMatrix<Tint> const &D) {
    Tint max_val(0), sum_val(0);
    for (int k = 0; k < D.rows(); k++) {
      Tint det = T_abs(DeterminantMat(remove_row(D, k)));
      sum_val += det;
      if (det > max_val) {
        max_val = det;
      }
    }
    return {max_val, sum_val};
  }
  std::vector<PerfectFaceEntry<T, Tint>> fill_piece(MyMatrix<Tint> const &D, int const &depth) {
    n_piece += 1;
    if (n_piece > max_piece) {
      throw HeckeFailure{"sharbly local filling: more than " + std::to_string(max_piece) + " pieces"};
    }
    std::vector<PerfectFaceEntry<T, Tint>> rhs = phi_boundary(D);
    if (rhs.empty()) {
      return {};
    }
    MyMatrix<T> D_T = UniversalMatrixConversion<T, Tint>(D);
    TopPerfectCone<Tint> tpc = walker.find_containing_cone_T(D_T);
    std::pair<Tint, Tint> size_D = det_size(D);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local piece depth=" << depth << " |D|=" << D.rows() << " max_det=" << size_D.first
       << " sum_det=" << size_D.second << " |rhs|=" << rhs.size() << "\n";
#endif
    if (size_D.first <= 1 || depth >= max_depth) {
      std::vector<TopPerfectCone<Tint>> l_top_seed{tpc};
      locate_simplex_cones(D, walker, l_top_seed);
      try {
        return partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options, os);
      } catch (CoveringFailure const &e) {
        // Too large for the budget of a piece: it is split further.
        if (size_D.first <= 1 || depth >= max_depth_split) {
          throw HeckeFailure{"sharbly local filling of a piece: " + e.message};
        }
      }
    }
    // The split along the worst face sigma = D minus u_j: with the reducing
    // points w of sigma used by Phi (averaged as Phi does),
    //   D = (1/|W|) sum_w sum_k D[u_k <- w]   modulo boundaries.
    // The pieces D[u_j <- w] have the faces sigma and its children in the
    // reduction of Phi, so their Phi(d .) sum to zero over w: they are
    // dropped (this is checked). The other pieces hold the children of
    // sigma, of smaller determinants, and are split further.
    int n_row = D.rows();
    int j_worst = -1;
    Tint det_worst(0);
    for (int k = 0; k < n_row; k++) {
      Tint det = T_abs(DeterminantMat(remove_row(D, k)));
      if (det > det_worst) {
        det_worst = det;
        j_worst = k;
      }
    }
    std::vector<MyVector<Tint>> l_w = reducer.reducing_points(remove_row(D, j_worst));
    if (l_w.empty()) {
      throw HeckeFailure{"sharbly local filling: no reducing point for a face"};
    }
    T coef_w = T(1) / T(static_cast<int>(l_w.size()));
    std::vector<MyMatrix<Tint>> l_piece_j;
    {
      ChainBuilder<T, Tint, Tgroup> cb_check(index + 1, fce, os);
      for (auto &w : l_w) {
        MyMatrix<Tint> Dj = D;
        AssignMatrixRow(Dj, j_worst, w);
        l_piece_j.push_back(Dj);
        for (auto &fe : phi_boundary(Dj)) {
          cb_check.f_insert(fe.value, fe.iOrb, fe.M);
        }
      }
      if (!cb_check.is_zero_chain()) {
        // Not expected; the pieces are then kept.
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: local piece, the pieces through the worst face do not cancel\n";
#endif
      } else {
        l_piece_j.clear();
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local piece split along the face " << j_worst << " of det " << det_worst
       << " with " << l_w.size() << " reducing points\n";
#endif
    ChainBuilder<T, Tint, Tgroup> cb(index, fce, os);
    auto f_add = [&](MyMatrix<Tint> const &Dk) -> void {
      if (RankMat(Dk) < Dk.cols()) {
        return;
      }
      for (auto &fe : fill_piece(Dk, depth + 1)) {
        cb.f_insert(coef_w * fe.value, fe.iOrb, fe.M);
      }
    };
    for (auto &w : l_w) {
      for (int k = 0; k < n_row; k++) {
        if (k == j_worst) {
          continue;
        }
        MyMatrix<Tint> Dk = D;
        AssignMatrixRow(Dk, k, w);
        f_add(Dk);
      }
    }
    for (auto &Dj : l_piece_j) {
      f_add(Dj);
    }
    return cb.get_faces();
  }
public:
  SharblyLocalFiller(int const &_index, FullComplexEnumeration<T, Tint, Tgroup> const &_fce,
                     ConeWalker<T, Tint, Tgroup> const &_walker,
                     SharblyReducer<T, Tint, Tgroup> const &_reducer,
                     std::function<bool(MyMatrix<Tint> const &)> const &_f_filter,
                     CoveringOptions const &_options, std::ostream &_os)
      : index(_index), fce(_fce), walker(_walker), reducer(_reducer), f_filter(_f_filter),
        options(_options), os(_os) {
  }
  size_t get_n_piece() const { return n_piece; }
  // Phi(dD) for the (n+1) x n matrix D.
  std::vector<PerfectFaceEntry<T, Tint>> phi_boundary(MyMatrix<Tint> const &D) {
    ChainBuilder<T, Tint, Tgroup> cb(index + 1, fce, os);
    for (int k = 0; k < D.rows(); k++) {
      T sign = (k % 2 == 0) ? T(1) : T(-1);
      for (auto &fe : phi(remove_row(D, k))) {
        cb.f_insert(sign * fe.value, fe.iOrb, fe.M);
      }
    }
    return cb.get_faces();
  }
  // A filling of rhs, the image of the boundary of the cell whose vectors
  // (one per pair {v, -v}) mapped by y are the rows of V. Nothing if rhs
  // is not a combination of the Phi(dD_S).
  std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
  fill(MyMatrix<Tint> const &V, std::vector<PerfectFaceEntry<T, Tint>> const &rhs) {
    n_piece = 0;
    int n = V.cols();
    int m = V.rows();
    // The (n+1)-subsets of the rows spanning R^n.
    std::vector<MyMatrix<Tint>> l_D;
    std::vector<int> subset(n + 1);
    std::function<void(int, int)> f_enum = [&](int start, int depth) -> void {
      if (depth == n + 1) {
        MyMatrix<Tint> D(n + 1, n);
        for (int k = 0; k <= n; k++) {
          AssignMatrixRow(D, k, GetMatrixRow(V, subset[k]));
        }
        if (RankMat(D) == n) {
          l_D.push_back(D);
        }
        return;
      }
      for (int i = start; i < m; i++) {
        subset[depth] = i;
        f_enum(i + 1, depth + 1);
      }
    };
    f_enum(0, 0);
    // The system sum_S c_S Phi(dD_S) = rhs, on the cells met.
    ComplexBuilder<T, Tint, Tgroup> cb2(index + 1, fce, os);
    std::vector<std::vector<PerfectFaceEntry<T, Tint>>> l_bnd;
    for (auto &D : l_D) {
      std::vector<PerfectFaceEntry<T, Tint>> bnd = phi_boundary(D);
      for (auto &fe : bnd) {
        cb2.f_insert(fe.iOrb, fe.M);
      }
      l_bnd.emplace_back(std::move(bnd));
    }
    for (auto &fe : rhs) {
      cb2.f_insert(fe.iOrb, fe.M);
    }
    int n_D = l_D.size();
    int n_cell = cb2.dimension();
    MyMatrix<T> A = ZeroMatrix<T>(n_D, n_cell);
    for (int i_D = 0; i_D < n_D; i_D++) {
      for (auto &fe : l_bnd[i_D]) {
        std::pair<size_t, int> p = cb2.get_position(fe.iOrb, fe.M);
        A(i_D, p.first) += fe.value * T(p.second);
      }
    }
    MyVector<T> b = ZeroVector<T>(n_cell);
    for (auto &fe : rhs) {
      std::pair<size_t, int> p = cb2.get_position(fe.iOrb, fe.M);
      b(p.first) += fe.value * T(p.second);
    }
    std::optional<MyVector<T>> opt = SolutionMat(A, b);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local fill |V|=" << m << " n_D=" << n_D << " n_cell=" << n_cell
       << " solvable=" << (opt ? 1 : 0) << "\n";
#endif
    if (!opt) {
      return {};
    }
    ChainBuilder<T, Tint, Tgroup> cb(index, fce, os);
    for (int i_D = 0; i_D < n_D; i_D++) {
      T c = (*opt)(i_D);
      if (c == 0) {
        continue;
      }
      for (auto &fe : fill_piece(l_D[i_D], 0)) {
        cb.f_insert(c * fe.value, fe.iOrb, fe.M);
      }
    }
    return cb.get_faces();
  }
};

template <typename T, typename Tint, typename Tgroup>
void compute_sharbly_level_below(
    int const &index, FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeCosets<T> const &cosets, ConeWalker<T, Tint, Tgroup> const &walker,
    SharblyReducer<T, Tint, Tgroup> const &reducer,
    std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> &images,
    std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  size_t n_coset = cosets.size();
  std::function<bool(MyMatrix<Tint> const &)> f_filter =
      [n](MyMatrix<Tint> const &EXT) -> bool {
    return RankMat(EXT) == n;
  };
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  HeckeBudget budget;
  int n_orb = fce.levels[index].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    if (!face.is_well_rounded) {
      continue;
    }
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
        get_stab_coset_action(face.l_gens, cosets);
    std::vector<int> status(n_coset, 0);
    for (size_t j = 0; j < n_coset; j++) {
      if (status[j] == 1) {
        continue;
      }
      std::vector<PerfectFaceEntry<T, Tint>> rhs =
          get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
      std::vector<PerfectFaceEntry<T, Tint>> sol;
      if (!rhs.empty()) {
#ifdef SANITY_CHECK_HECKE_OPERATORS
        for (auto &fe : rhs) {
          MyMatrix<Tint> EXT = fce.levels[index + 1].l_faces[fe.iOrb].EXT * fe.M;
          if (!f_filter(EXT) || !fce.levels[index + 1].l_faces[fe.iOrb].is_well_rounded) {
            std::cerr << "HECKE: sharbly, the right hand side should be made of well rounded cells\n";
            throw TerminalException{1};
          }
        }
#endif
        CoveringOptions options;
        options.guided_growth = true;
        options.clean_failure = true;
        options.max_entries_exact = budget.max_entries_exact;
        options.l1_first = true;
        options.max_l1_iter = budget.max_yall1_iter;
        // Three attempts: the filling of the whole right hand side on a
        // small covering, accepted if it is not bloated; the local filling
        // (the 1-sharbly reduction) with the same small coverings for the
        // pieces; the filling of the whole right hand side on a covering of
        // any size.
        CoveringOptions options_small = options;
        options_small.max_top = 100;
        options_small.max_unknowns = 50000;
        options_small.max_entries_exact = 500000;
        std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_sol;
        {
          std::vector<TopPerfectCone<Tint>> l_top_seed =
              get_sharbly_covering_seed(index, iOrb, cosets.l_coset[j], fce, walker, reducer);
          try {
            std::vector<PerfectFaceEntry<T, Tint>> sol_small =
                partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options_small, os);
            if (sol_small.size() <= 20 * rhs.size() + 50) {
              opt_sol = sol_small;
            }
          } catch (CoveringFailure const &e) {
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: sharbly small global filling failed: " << e.message << "\n";
#endif
          }
#ifdef TIMINGS_HECKE_OPERATORS
          os << "|HECKE: sharbly small global filling iOrb=" << iOrb << " j=" << j
             << " success=" << (opt_sol ? 1 : 0) << "|=" << time << "\n";
#endif
        }
        if (!opt_sol) {
          SharblyLocalFiller<T, Tint, Tgroup> filler(index, fce, walker, reducer, f_filter, options_small, os);
          MyMatrix<Tint> V = reducer.primitive_rows(
              UniversalMatrixConversion<T, Tint>(half_vector_set(face.EXT)) * cosets.l_coset[j]);
          try {
            opt_sol = filler.fill(V, rhs);
          } catch (HeckeFailure const &e) {
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: sharbly local filling failed: " << e.message << "\n";
#endif
          }
#ifdef TIMINGS_HECKE_OPERATORS
          os << "|HECKE: sharbly local filling iOrb=" << iOrb << " j=" << j
             << " success=" << (opt_sol ? 1 : 0) << " n_piece=" << filler.get_n_piece() << "|=" << time << "\n";
#endif
        }
        if (opt_sol) {
          sol = *opt_sol;
        } else {
          std::vector<TopPerfectCone<Tint>> l_top_seed =
              get_sharbly_covering_seed(index, iOrb, cosets.l_coset[j], fce, walker, reducer);
          try {
            sol = partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options, os);
          } catch (CoveringFailure const &e) {
            throw HeckeFailure{"sharbly, filling of the cell " + std::to_string(iOrb) + " of level " +
                               std::to_string(index) + " for the coset " + std::to_string(j) + ": " + e.message};
          }
        }
      }
#ifdef TIMINGS_HECKE_OPERATORS
      os << "|HECKE: sharbly level below iOrb=" << iOrb << "/" << n_orb << " j=" << j << "/" << n_coset
         << " |rhs|=" << rhs.size() << " |sol|=" << sol.size() << "|=" << time << "\n";
#endif
      CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
      std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
          get_coset_stabilizer_schreier(co, face.l_gens, act, n);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      if (!is_chain_coset_stabilizer_invariant(index, iOrb, index + 1, rhs, l_schreier, fce, os)) {
        std::cerr << "HECKE: sharbly, the right hand side should be invariant under the coset stabilizer\n";
        throw TerminalException{1};
      }
#endif
      std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
          average_chain_coset_stabilizer(index, iOrb, sol, l_schreier, fce, os);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      {
        std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary_well_rounded(index, sol_avg, fce, os);
        if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
          std::cerr << "HECKE: sharbly, the averaged chain is not a filling\n";
          throw TerminalException{1};
        }
        if (!is_chain_coset_stabilizer_invariant(index, iOrb, index, sol_avg, l_schreier, fce, os)) {
          std::cerr << "HECKE: sharbly, the averaged filling should be invariant\n";
          throw TerminalException{1};
        }
      }
#endif
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
    // The chain map property on all the cosets, after the propagation.
    for (size_t j = 0; j < n_coset; j++) {
      std::vector<PerfectFaceEntry<T, Tint>> rhs =
          get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
      std::vector<PerfectFaceEntry<T, Tint>> bnd =
          chain_boundary_well_rounded(index, images[index][iOrb][j], fce, os);
      if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
        std::cerr << "HECKE: sharbly, the chain map property fails at index=" << index
                  << " iOrb=" << iOrb << " j=" << j << "\n";
        throw TerminalException{1};
      }
    }
#endif
  }
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map_sharbly(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                MyMatrix<T> const &x, std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  PerfectComplexOptions const &pco = fce.pctdi.pco;
  if (!pco.compute_boundary) {
    std::cerr << "HECKE: The sharbly method requires the boundary, please set ComputeBoundary = T\n";
    throw TerminalException{1};
  }
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
  int symbol_level = get_symbol_level(fce);
  if (symbol_level < 0 || symbol_level >= n_levels) {
    std::cerr << "HECKE: sharbly, the symbol level " << symbol_level
              << " is not in the complex (n_levels=" << n_levels << ")\n";
    throw TerminalException{1};
  }
  int iOrb0 = get_unimodular_frame_orbit(fce, symbol_level, os);
  MyMatrix<Tint> E0 = half_vector_set(fce.levels[symbol_level].l_faces[iOrb0].EXT);
  MyMatrix<Tint> E0_inv = Inverse(E0);
  ConeWalker<T, Tint, Tgroup> walker(fce, os);
  SharblyReducer<T, Tint, Tgroup> reducer(fce, symbol_level, iOrb0, E0_inv, os);
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images(n_levels);
  std::vector<uint8_t> l_level_computed(n_levels, 0);
  l_level_computed[symbol_level] = 1;
  for (int index = 0; index < n_levels; index++) {
    int n_orb = fce.levels[index].l_faces.size();
    images[index] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
  }
  int n_orb = fce.levels[symbol_level].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[symbol_level].l_faces[iOrb];
    if (!face.is_well_rounded) {
      continue;
    }
    MyMatrix<Tint> half = half_vector_set(face.EXT);
    if (half.rows() != n) {
      throw HeckeFailure{"sharbly: the cell " + std::to_string(iOrb) +
                         " of the symbol level has " + std::to_string(half.rows()) +
                         " vectors instead of " + std::to_string(n) +
                         ": not a symbol, use another method"};
    }
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(half);
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
        get_stab_coset_action(face.l_gens, cosets);
    std::vector<int> status(n_coset, 0);
    for (size_t j = 0; j < n_coset; j++) {
      if (status[j] == 1) {
        continue;
      }
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: sharbly iOrb=" << iOrb << "/" << n_orb << " j=" << j << "/" << n_coset << "\n";
#endif
      MyMatrix<T> W_T = EXT_T * cosets.l_coset[j];
      MyMatrix<Tint> W = reducer.primitive_rows(W_T);
      std::vector<PerfectFaceEntry<T, Tint>> sol = reducer.reduce(W);
      CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
      std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
          get_coset_stabilizer_schreier(co, face.l_gens, act, n);
      std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
          average_chain_coset_stabilizer(symbol_level, iOrb, sol, l_schreier, fce, os);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      if (!is_chain_coset_stabilizer_invariant(symbol_level, iOrb, symbol_level, sol_avg, l_schreier, fce, os)) {
        std::cerr << "HECKE: sharbly, the averaged chain should be invariant\n";
        throw TerminalException{1};
      }
#endif
      size_t n_orbit = co.l_orbit.size();
      for (size_t pos = 0; pos < n_orbit; pos++) {
        int j2 = co.l_orbit[pos];
        int sign = get_face_orientation(face.EXT, LinSpa.ListMat, face.or_info, co.transversal[pos]);
        MyMatrix<Tint> Minv = Inverse(co.transversal_conj[pos]);
        images[symbol_level][iOrb][j2] = chain_right_multiply(sol_avg, T(sign), Minv);
        status[j2] = 1;
      }
    }
  }
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: sharbly images|=" << time << "\n";
#endif
  // The level below needs the coverings by the perfect cones.
  if (symbol_level > 0 && pco.compute_partial_contracting_homotopy) {
    compute_sharbly_level_below(symbol_level - 1, fce, cosets, walker, reducer, images, os);
    l_level_computed[symbol_level - 1] = 1;
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: sharbly images level below|=" << time << "\n";
#endif
  }
  HeckeChainMap<T, Tint> hcm;
  hcm.x = x;
  hcm.cosets = std::move(cosets);
  hcm.images = std::move(images);
  hcm.l_level_computed = std::move(l_level_computed);
  hcm.method = "Sharbly";
  hcm.transposed = false;
  return hcm;
}

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
        size_t n_rep = gos.l_rep.size();
        for (size_t i_rep = 0; i_rep < n_rep; i_rep++) {
          if (gos.l_orientable[i_rep]) {
            int rep = gos.l_rep[i_rep];
            gc.basis_index[index][iOrb][rep] = gc.basis[index].size();
            gc.basis[index].push_back({iOrb, rep});
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
  The Hecke operator on the basis cell (iOrb, i) of level index:
    T(F g_i) = sum_k Phi(F, y_j) g'   with g_i x_k = y_j g'
  expressed in the basis of the level.
 */
template <typename T, typename Tint>
std::map<int, T> get_hecke_image_cell_sparse(
    int const &index, int const &i_cell, GammaComplex<T, Tint> const &gc,
    HeckeChainMap<T, Tint> const &hcm, HeckeCosets<T> const &cosets_gamma) {
  std::map<int, T> V;
  int iOrb = gc.basis[index][i_cell].first;
  int i_coset = gc.basis[index][i_cell].second;
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(gc.ct.l_coset[i_coset]);
  size_t n_hecke = cosets_gamma.size();
  for (size_t k = 0; k < n_hecke; k++) {
    MyMatrix<T> z = g_T * cosets_gamma.l_coset[k];
    std::pair<int, MyMatrix<T>> pair = hcm.cosets.decompose(z);
    MyMatrix<Tint> g_conj = UniversalMatrixConversion<Tint, T>(pair.second);
    for (auto &fe : hcm.images[index][iOrb][pair.first]) {
      MyMatrix<Tint> M = fe.M * g_conj;
      std::optional<std::pair<int, int>> opt = gc.cell_position(index, fe.iOrb, M);
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
    [[maybe_unused]] FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  int n_cell = gc.basis[index].size();
  MyVector<T> V = ZeroVector<T>(n_cell);
  for (auto &kv : get_hecke_image_cell_sparse(index, i_cell, gc, hcm, cosets_gamma)) {
    V(kv.first) += kv.second;
  }
  return V;
}

/*
  The Hecke operator on the homology of a level modulo a prime p, by
  sparse linear algebra: no matrix of the size of the complex is formed,
  which is what limits the rational computation. With B the space of the
  boundaries and Z the one of the cycles, H = Z / B is computed in the
  coordinates of the cells outside a set S whose cycle condition rows are
  a basis of the cycle conditions (see below); a basis of H is given by
  the free cells of the echelon form of the boundaries there, lifted to
  cycles, and the image of a basis vector is restricted to the complement
  of S and reduced modulo the boundaries, its normal form giving its
  coordinates.
  For a chain map commuting with the boundary (Sharbly, Vertex) the cycles
  are the v with v D_index = 0 and the boundaries the rows of D_{index-1};
  for the transposed one (Dual), the v with D_{index-1} v = 0 and the
  columns of D_index. The entries of the Hecke matrix are residues mod p.
 */
template <typename T, typename Tint>
HeckeHomologyLevel<T> get_hecke_level_mod(int const &index, GammaComplex<T, Tint> const &gc,
                                 HeckeChainMap<T, Tint> const &hcm,
                                 HeckeCosets<T> const &cosets_gamma, int64_t const &p,
                                 [[maybe_unused]] std::ostream &os) {
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
  std::vector<int> l_label_next = f_label(l_cycle_row, n_next);
  SparseEchelonMod echelon_cyc(p);
  std::vector<uint8_t> in_S(n_cell, 0);
  for (int c = 0; c < n_cell; c++) {
    if (!l_cycle_row[c].empty() && echelon_cyc.insert(f_relabel(l_cycle_row[c], l_label_next))) {
      in_S[c] = 1;
    }
  }
  // The same rows with markers at n_next + c: the reduction of the cycle
  // condition row of a cell gives its expression in the rows of S.
  SparseEchelonMod echelon_cyc_mark(p);
  for (int c = 0; c < n_cell; c++) {
    if (in_S[c] == 1) {
      SparseVectorMod v = f_relabel(l_cycle_row[c], l_label_next);
      v.push_back({n_next + c, 1});
      echelon_cyc_mark.insert(v);
    }
  }
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
  std::vector<int> l_label = f_label(l_bnd_res, n_cell);
  std::sort(l_bnd_res.begin(), l_bnd_res.end(),
            [](SparseVectorMod const &a, SparseVectorMod const &b) -> bool { return a.size() < b.size(); });
  SparseRrefMod echelon_bnd(p);
  for (auto &v : l_bnd_res) {
    echelon_bnd.insert(f_relabel(v, l_label));
  }
  // The free cells of S^c give the basis of H; their positions.
  std::vector<int> l_free;
  std::unordered_map<int, int> map_free;
  for (int c = 0; c < n_cell; c++) {
    if (in_S[c] == 0 && !echelon_bnd.is_pivot(l_label[c])) {
      map_free[l_label[c]] = l_free.size();
      l_free.push_back(c);
    }
  }
  int dim_hom = l_free.size();
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: mod p index=" << index << " n_cell=" << n_cell << " |S|=" << echelon_cyc.rank()
     << " rank(B)=" << echelon_bnd.rank() << " dim_hom=" << dim_hom << "\n";
#endif
  MyMatrix<T> hecke_matrix = ZeroMatrix<T>(dim_hom, dim_hom);
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
    // Its image, restricted to S^c and reduced modulo the boundaries.
    std::map<int, int64_t> acc;
    for (auto &e : v) {
      for (auto &kv : get_hecke_image_cell_sparse(index, e.first, gc, hcm, cosets_gamma)) {
        if (in_S[kv.first] == 1) {
          continue;
        }
        std::optional<int64_t> opt = residue_mod_p(kv.second, p);
        if (!opt) {
          // The homology is known, but not the Hecke matrix.
          return {index, n_cell, MyMatrix<T>(dim_hom, 0), MyMatrix<T>(0, 0),
                  "a denominator of an image is divisible by p=" + std::to_string(p)};
        }
        int64_t &val = acc[l_label[kv.first]];
        val = ResidueModPrime(val + e.second * (*opt), p);
      }
    }
    SparseVectorMod w;
    for (auto &kv : acc) {
      if (kv.second != 0) {
        w.push_back(kv);
      }
    }
    for (auto &e : echelon_bnd.reduce(w)) {
      auto iter = map_free.find(e.first);
      if (iter == map_free.end()) {
        throw HeckeFailure{"Hecke matrix mod p: a normal form is not supported on the free cells, index=" +
                           std::to_string(index)};
      }
      hecke_matrix(i, iter->second) = UniversalScalarConversion<T, int64_t>(e.second);
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
    int64_t const &prime, std::ostream &os) {
  [[maybe_unused]] bool const &dual = hcm.transposed;
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> l_gens_G = get_ambient_group_generators(fce, os);
  CosetTableGamma<Tint> ct = enumerate_cosets_gamma(gamma, l_gens_G, n, os);
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
#ifdef SANITY_CHECK_HECKE_OPERATORS
  {
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
    if (hcm.is_level_computed(index)) {
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
    if (!hcm.is_level_computed(index)) {
      continue;
    }
    int n_cell = gc.basis[index].size();
    if (prime > 0) {
      result.l_level.emplace_back(get_hecke_level_mod(index, gc, hcm, cosets_gamma, prime, os));
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
    MyMatrix<T> hecke_matrix = ZeroMatrix<T>(dim_hom, dim_hom);
    if (dim_hom > 0) {
      // The images of the basis cells, computed once.
      std::vector<MyVector<T>> l_img(n_cell);
      for (int i_cell = 0; i_cell < n_cell; i_cell++) {
        l_img[i_cell] = get_hecke_image_cell(index, i_cell, gc, hcm, cosets_gamma, fce);
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
#endif  // SRC_PERFECT_HECKE_OPERATORS_H_
// clang-format on
