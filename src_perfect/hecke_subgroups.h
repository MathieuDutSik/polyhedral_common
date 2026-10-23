// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_SUBGROUPS_H_
#define SRC_PERFECT_HECKE_SUBGROUPS_H_

// clang-format off
#include "hecke_common.h"
#include "MAT_SparseMatrixMod.h"
#include "MAT_MatrixModInt64.h"
#include "ImagQuadRing.h"
// clang-format on

/*
  The finite index subgroups Gamma of G (Hecke operators, step B): the
  congruence subgroups Principal, Gamma0, Gamma1 and their intersections
  with SL, of GL_n(Z) or of GL_m(O_K) for a Hermitian T-space, given by
  the keys of the cosets g Gamma; the coset table and the Schreier
  generators.
 */

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

// What the action of an element s of G on the keys needs: its residues,
// the ones of its inverse and the sign of its determinant.
struct GammaAction {
  std::vector<int64_t> s_mod;
  std::vector<int64_t> s_inv_mod;
  int64_t det_sign;
  // For a Hermitian T-space, the determinant over O_K (a unit).
  std::pair<int64_t, int64_t> det_k;
};

template <typename Tint> struct FiniteIndexSubgroup {
  std::string type;
  Tint N;
  int n;
  /*
    For a Hermitian T-space over O_K = Z[omega], omega^2 = s omega - p,
    the matrices are the realifications of the m x m matrices over O_K
    (n = 2m, the entry (i, j) over O_K being (g(i, j), g(i, m + j))), and
    the subgroups are those of GL_m(O_K): the congruences are on the
    entries over O_K modulo N O_K and the SL condition is det_K = 1.
   */
  int m_field = 0;
  int64_t eSum_field = 0;
  int64_t eProd_field = 0;
  // The units of O_K / N O_K, for the canonical form of the Gamma0 keys.
  std::vector<std::pair<int64_t, int64_t>> l_unit_mod;
  bool is_hermitian() const { return m_field > 0; }
  void set_field(int const &m, int64_t const &eSum, int64_t const &eProd) {
    m_field = m;
    eSum_field = eSum;
    eProd_field = eProd;
    int64_t mod = get_n_mod();
    l_unit_mod.clear();
    for (int64_t a = 0; a < mod; a++) {
      for (int64_t b = 0; b < mod; b++) {
        int64_t nrm = ResidueModPrime(a * a + eSum * a * b + eProd * b * b, mod);
        if (std::gcd(nrm, mod) == 1 && !(a == 1 && b == 0)) {
          l_unit_mod.push_back({a, b});
        }
      }
    }
  }
  template <typename Tval>
  std::pair<Tval, Tval> mul_k(std::pair<Tval, Tval> const &x, std::pair<Tval, Tval> const &y) const {
    ImagQuadRing<Tval> ring(UniversalScalarConversion<Tval, int64_t>(eSum_field),
                            UniversalScalarConversion<Tval, int64_t>(eProd_field));
    return ring.mul(x, y);
  }

  // The determinant over O_K of the realification g, by the expansion
  // along the first row.
  std::pair<Tint, Tint> determinant_k(MyMatrix<Tint> const &g) const {
    int m = m_field;
    std::function<std::pair<Tint, Tint>(std::vector<int> const &, std::vector<int> const &)> f_det =
        [&](std::vector<int> const &l_row, std::vector<int> const &l_col) -> std::pair<Tint, Tint> {
      if (l_row.empty()) {
        return {Tint(1), Tint(0)};
      }
      std::pair<Tint, Tint> det{Tint(0), Tint(0)};
      int i = l_row[0];
      std::vector<int> l_row_sub(l_row.begin() + 1, l_row.end());
      for (size_t pos = 0; pos < l_col.size(); pos++) {
        int j = l_col[pos];
        std::vector<int> l_col_sub;
        for (size_t pos2 = 0; pos2 < l_col.size(); pos2++) {
          if (pos2 != pos) {
            l_col_sub.push_back(l_col[pos2]);
          }
        }
        std::pair<Tint, Tint> entry{g(i, j), g(i, m + j)};
        std::pair<Tint, Tint> term = mul_k(entry, f_det(l_row_sub, l_col_sub));
        if (pos % 2 == 0) {
          det = {det.first + term.first, det.second + term.second};
        } else {
          det = {det.first - term.first, det.second - term.second};
        }
      }
      return det;
    };
    std::vector<int> l_idx(m);
    for (int i = 0; i < m; i++) {
      l_idx[i] = i;
    }
    return f_det(l_idx, l_idx);
  }
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
      if (is_hermitian()) {
        if (determinant_k(g) != std::pair<Tint, Tint>{Tint(1), Tint(0)}) {
          return false;
        }
      } else if (DeterminantMat(g) != 1) {
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
    if ((type == "Gamma0" || type == "Gamma1") && is_hermitian()) {
      // The last row over O_K is the row m - 1 of the realification.
      int m = m_field;
      for (int j = 0; j < m - 1; j++) {
        if (ResInt(g(m - 1, j), N) != 0 || ResInt(g(m - 1, m + j), N) != 0) {
          return false;
        }
      }
      if (type == "Gamma1") {
        Tint val = g(m - 1, m - 1) - 1;
        if (ResInt(val, N) != 0 || ResInt(g(m - 1, 2 * m - 1), N) != 0) {
          return false;
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
      return (DeterminantMatModInt64(r, n, m, is_prime_n()) == 1) ? 1 : -1;
    }
    Tint det = DeterminantMat(g);
    return (det > 0) ? 1 : -1;
  }
  // The canonical form of the row of a Gamma0 key under the units mod N:
  // the first nonzero entry is 1 for N prime, the smallest multiple in the
  // lexicographic order otherwise.
  void canonicalize_gamma0(std::vector<int64_t> &row) const {
    int64_t m = get_n_mod();
    if (is_hermitian()) {
      // The row is an element of (O_K / N)^m, in the realified layout.
      // The smallest of the u row in the lexicographic order, u a unit
      // of O_K / N (the product written out, this is the inner loop of
      // the coset keys).
      int mf = m_field;
      std::vector<int64_t> best = row;
      std::vector<int64_t> cand(2 * mf);
      for (auto &u : l_unit_mod) {
        for (int j = 0; j < mf; j++) {
          int64_t a = row[j], b = row[mf + j];
          int64_t ub_b = (u.second * b) % m;
          cand[j] = ResidueModPrime(u.first * a - eProd_field * ub_b, m);
          cand[mf + j] = ResidueModPrime(u.first * b + u.second * a + eSum_field * ub_b, m);
        }
        if (cand < best) {
          best = cand;
        }
      }
      row = best;
      return;
    }
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
  // Whether N is prime, computed once (it is needed for every key).
  mutable int8_t is_prime_cache = -1;
  bool is_prime_n() const {
    if (is_prime_cache == -1) {
      int64_t m = get_n_mod();
      bool is_prime = (m >= 2);
      for (int64_t d = 2; d * d <= m; d++) {
        if (m % d == 0) {
          is_prime = false;
          break;
        }
      }
      is_prime_cache = is_prime ? 1 : 0;
    }
    return is_prime_cache == 1;
  }
  // The inverse modulo N of the residues r of an element of G.
  std::vector<int64_t> inverse_residues(std::vector<int64_t> const &r) const {
    std::optional<std::vector<int64_t>> opt = InverseMatModInt64(r, n, get_n_mod(), is_prime_n());
    if (!opt) {
      std::cerr << "HECKE: an element of G is not invertible modulo N\n";
      throw TerminalException{1};
    }
    return *opt;
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
    std::string type = base_type();
    GammaKey key;
    if (type == "Principal") {
      key = r;
    } else if (type == "Gamma0" || type == "Gamma1") {
      std::vector<int64_t> inv = inverse_residues(r);
      // The last row of g^{-1} (over O_K: the row m - 1 of the realification).
      int i_last = is_hermitian() ? m_field - 1 : n - 1;
      key = std::vector<int64_t>(inv.begin() + i_last * n, inv.begin() + (i_last + 1) * n);
      if (type == "Gamma0") {
        canonicalize_gamma0(key);
      }
    } else if (type != "Full") {
      std::cerr << "HECKE: Unsupported subgroup type=" << type << "\n";
      throw TerminalException{1};
    }
    if (is_sl()) {
      if (is_hermitian()) {
        std::pair<Tint, Tint> det = determinant_k(g);
        key.push_back(UniversalScalarConversion<int64_t, Tint>(det.first));
        key.push_back(UniversalScalarConversion<int64_t, Tint>(det.second));
      } else {
        key.push_back(get_det_sign(g, r));
      }
    }
    return key;
  }
  GammaAction get_action(MyMatrix<Tint> const &s) const {
    GammaAction a;
    a.s_mod = get_residues(s);
    a.det_sign = get_det_sign(s, a.s_mod);
    a.s_inv_mod = inverse_residues(a.s_mod);
    if (is_hermitian()) {
      std::pair<Tint, Tint> det = determinant_k(s);
      a.det_k = {UniversalScalarConversion<int64_t, Tint>(det.first),
                 UniversalScalarConversion<int64_t, Tint>(det.second)};
    }
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
      if (is_hermitian()) {
        size_t len = key.size();
        std::pair<int64_t, int64_t> det = mul_k(a.det_k, std::pair<int64_t, int64_t>{key[len - 2], key[len - 1]});
        key_ret.push_back(det.first);
        key_ret.push_back(det.second);
      } else {
        key_ret.push_back(a.det_sign * key.back());
      }
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
  // The keys against the membership, linearly in the index: each coset
  // representative has its own key and g_i^{-1} g_i is in Gamma.
  for (size_t i = 0; i < ct.l_coset.size(); i++) {
    if (gamma.get_key(ct.l_coset[i]) != ct.l_key[i] ||
        !gamma.is_in_gamma(ct.l_coset_inv[i] * ct.l_coset[i])) {
      std::cerr << "HECKE: Inconsistent coset table at i=" << i << "\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE
  // All the pairs: quadratic in the index.
  std::cerr << "HECKE: SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE, checking all the pairs of the coset table\n";
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
  std::cerr << "HECKE: SANITY_CHECK_HECKE_OPERATORS_EXTENSIVE, done\n";
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

// clang-format off
#endif  // SRC_PERFECT_HECKE_SUBGROUPS_H_
// clang-format on
