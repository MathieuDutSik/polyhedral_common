// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PROJECTIVE_SYSTEMS_PROJECTIVE_COMPLEX_H_
#define SRC_PROJECTIVE_SYSTEMS_PROJECTIVE_COMPLEX_H_

// clang-format off
#include "hecke_operators.h"
#include "ImagQuadRing.h"
#include "PolytopeEquiStab.h"
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PROJECTIVE_COMPLEX
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PROJECTIVE_COMPLEX
#endif

#ifdef TIMINGS
#define TIMINGS_PROJECTIVE_COMPLEX
#endif

/*
  The cell complex of a projective system (port of ProjectiveSystem.g of
  MyPolyhedral, made independent of the group).

  The group G is the group of the matrices of GL_m(R) preserving a form h
  on R^m, for R = Z or R = O_K = Z[omega] (omega^2 = s omega - p, K
  imaginary quadratic), acting on the row vectors on the right:
    G = { g : g h g^* = h }.
  The cells are the configurations of vectors of R^m taken up to the units
  of R, given by a description file (see README.md): the representative
  configurations of each level, the top cells first. The cells of a level
  are the subsets of the vectors of the top cells that are equivalent
  under G to one of the representatives of the level, as in
  SP4Z_GetCellComplex; the faces of a cell are the cells of the next level
  contained in it. The orientation of a cell is combinatorial: the sign
  character of its stabilizer, determined from the signs of its boundary.

  The vectors and matrices over O_K are realified as in src_perfect (see
  FiniteIndexSubgroup): the vector (x_1, ..., x_m), x_j = a_j + b_j omega,
  is (a_1, ..., a_m, b_1, ..., b_m), and a matrix of GL_m(O_K) is the
  matrix of the induced Z-linear map of Z^{2m}. The form h is then
  replaced by the real form Tr_{K/Q} h(u, v), which a matrix commuting with
  the multiplication by omega preserves iff it preserves h.

  The levels are indexed as in src_perfect: the top cells at the index 0
  and the faces of a cell of the index k at the index k + 1.
 */

namespace psys {

using Tint = mpz_class;
using T = mpq_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

// The subsets of k elements of [0, n).
inline std::vector<std::vector<int>> combinations(int const &n, int const &k) {
  std::vector<std::vector<int>> l_comb;
  std::vector<int> sel(k);
  std::function<void(int, int)> f_rec = [&](int start, int depth) -> void {
    if (depth == k) {
      l_comb.push_back(sel);
      return;
    }
    for (int i = start; i < n; i++) {
      sel[depth] = i;
      f_rec(i + 1, depth + 1);
    }
  };
  f_rec(0, 0);
  return l_comb;
}

inline MyMatrix<Tint> select_rows(MyMatrix<Tint> const &EXT, std::vector<int> const &l_row) {
  MyMatrix<Tint> M(l_row.size(), EXT.cols());
  for (size_t i = 0; i < l_row.size(); i++) {
    for (int j = 0; j < EXT.cols(); j++) {
      M(i, j) = EXT(l_row[i], j);
    }
  }
  return M;
}

/*
  The ring R and the form: the dimension n = d m of the realification
  (d = 1 for Z, 2 for O_K), the multiplication by omega W, the units (as
  the matrices of the scalar multiplications) and the real form B.
 */
struct ProjectiveRing {
  int d = 1;
  int m = 0;
  long s = 0;
  long p = 0;
  MyMatrix<Tint> W;
  std::vector<MyMatrix<Tint>> l_unit;
  std::vector<MyMatrix<Tint>> l_ring_basis;
  MyMatrix<Tint> B;
  int n() const { return d * m; }
  // Whether g is in G: integral, R-linear and preserving B.
  template <typename Tf> bool in_group(MyMatrix<Tf> const &g) const {
    if (!IsIntegralMatrix(g)) {
      return false;
    }
    MyMatrix<Tf> B_f = UniversalMatrixConversion<Tf, Tint>(B);
    if (d == 2) {
      MyMatrix<Tf> W_f = UniversalMatrixConversion<Tf, Tint>(W);
      if (g * W_f != W_f * g) {
        return false;
      }
    }
    return g * B_f * g.transpose() == B_f;
  }
  // The element a + b omega of O_K for the coordinate j of v.
  std::pair<T, T> entry(MyVector<Tint> const &v, int const &j) const {
    return {T(v(j)), T(v(m + j))};
  }
  /*
    The primitive vector of the line of v: v divided by the gcd of its
    entries (over R).
   */
  MyVector<Tint> primitive(MyVector<T> const &v_in) const {
    MyVector<T> v = RemoveFractionVector(v_in);
    if (d == 1) {
      return UniversalVectorConversion<Tint, T>(v);
    }
    ImagQuadRing<T> ring{T(s), T(p)};
    std::vector<std::pair<T, T>> l_x;
    for (int j = 0; j < m; j++) {
      std::pair<T, T> x{v(j), v(m + j)};
      if (!ring.is_zero(x)) {
        l_x.push_back(x);
      }
    }
    std::optional<std::pair<T, T>> opt = ring.gcd(l_x);
    if (!opt) {
      std::cerr << "PSYS: primitive, the gcd of the entries is not principal\n";
      throw TerminalException{1};
    }
    MyVector<Tint> w(n());
    for (int j = 0; j < m; j++) {
      std::pair<T, T> y = ring.div({v(j), v(m + j)}, *opt);
      w(j) = UniversalScalarConversion<Tint, T>(y.first);
      w(m + j) = UniversalScalarConversion<Tint, T>(y.second);
    }
    return w;
  }
};

/*
  The description of a projective system: the ring and the form, the
  representative configurations of each level (top cells first), which
  cells are at infinity and how the Hecke operators locate a vector.
 */
struct ProjectiveDescription {
  ProjectiveRing ring;
  std::vector<std::vector<MyMatrix<Tint>>> l_conf;
  // "rank": the cells whose vectors do not span R^m are at infinity;
  // "single": only the cells of one vector are at infinity.
  std::string at_infinity;
  // "voronoi": the walk across the facets of the cones of rank one forms
  // of the top cells; "potential": the walk decreasing min_x |h(x, v)|.
  std::string walk;
};

/*
  The tokens of a description file, the comments (from # to the end of
  the line) removed.
 */
inline std::vector<std::string> read_tokens(std::string const &file) {
  std::ifstream is(file);
  if (!is.good()) {
    std::cerr << "PSYS: the description file " << file << " cannot be read\n";
    throw TerminalException{1};
  }
  std::vector<std::string> l_tok;
  std::string line;
  while (std::getline(is, line)) {
    size_t pos = line.find('#');
    if (pos != std::string::npos) {
      line = line.substr(0, pos);
    }
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) {
      l_tok.push_back(tok);
    }
  }
  return l_tok;
}

struct TokenReader {
  std::vector<std::string> l_tok;
  size_t pos = 0;
  std::string next() {
    if (pos >= l_tok.size()) {
      std::cerr << "PSYS: the description file ends too early\n";
      throw TerminalException{1};
    }
    return l_tok[pos++];
  }
  void expect(std::string const &key) {
    std::string tok = next();
    if (tok != key) {
      std::cerr << "PSYS: expected " << key << " in the description file, found " << tok << "\n";
      throw TerminalException{1};
    }
  }
  int next_int() { return ParseScalar<int>(next()); }
  MyMatrix<Tint> next_matrix() {
    int n_row = next_int();
    int n_col = next_int();
    MyMatrix<Tint> M(n_row, n_col);
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j < n_col; j++) {
        M(i, j) = ParseScalar<Tint>(next());
      }
    }
    return M;
  }
};

/*
  The real form B: for R = Z the form h itself, for R = O_K the trace form
  Tr_{K/Q} h(e, f) on the real basis e_j, omega e_j of O_K^m.
 */
inline MyMatrix<Tint> get_real_form(int const &m, long const &s, long const &p, MyMatrix<Tint> const &h_a,
                                    MyMatrix<Tint> const &h_b) {
  ImagQuadRing<T> ring{T(s), T(p)};
  std::vector<std::pair<T, T>> l_basis{{T(1), T(0)}, {T(0), T(1)}};
  MyMatrix<Tint> B(2 * m, 2 * m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      std::pair<T, T> h_ij{T(h_a(i, j)), T(h_b(i, j))};
      for (int ki = 0; ki < 2; ki++) {
        for (int kj = 0; kj < 2; kj++) {
          std::pair<T, T> val = ring.mul(ring.mul(l_basis[ki], h_ij), ring.conj(l_basis[kj]));
          B(ki * m + i, kj * m + j) = UniversalScalarConversion<Tint, T>(ring.trace(val));
        }
      }
    }
  }
  return B;
}

inline ProjectiveDescription read_description(std::string const &file) {
  TokenReader tr{read_tokens(file)};
  ProjectiveDescription desc;
  ProjectiveRing &ring = desc.ring;
  tr.expect("ring");
  std::string ring_name = tr.next();
  if (ring_name == "Z") {
    ring.d = 1;
  } else if (ring_name == "QuadraticImag") {
    ring.d = 2;
    ring.s = tr.next_int();
    ring.p = tr.next_int();
  } else {
    std::cerr << "PSYS: the ring should be Z or QuadraticImag, not " << ring_name << "\n";
    throw TerminalException{1};
  }
  tr.expect("dimension");
  ring.m = tr.next_int();
  int n = ring.n();
  tr.expect("form");
  if (ring.d == 1) {
    ring.B = tr.next_matrix();
  } else {
    MyMatrix<Tint> h_a = tr.next_matrix();
    MyMatrix<Tint> h_b = tr.next_matrix();
    ring.B = get_real_form(ring.m, ring.s, ring.p, h_a, h_b);
  }
  if (ring.B.rows() != n || ring.B.cols() != n) {
    std::cerr << "PSYS: the form should be of dimension " << n << "\n";
    throw TerminalException{1};
  }
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  ring.l_ring_basis.push_back(id);
  if (ring.d == 1) {
    ring.l_unit = {id, MyMatrix<Tint>(-id)};
  } else {
    // (a, b) -> omega (a + b omega) = -p b + (a + s b) omega.
    ring.W = ZeroMatrix<Tint>(n, n);
    for (int j = 0; j < ring.m; j++) {
      ring.W(j, ring.m + j) = 1;
      ring.W(ring.m + j, j) = -ring.p;
      ring.W(ring.m + j, ring.m + j) = ring.s;
    }
    ring.l_ring_basis.push_back(ring.W);
    ImagQuadRing<T> iq{T(ring.s), T(ring.p)};
    for (auto &u : iq.units()) {
      MyMatrix<Tint> U = UniversalScalarConversion<Tint, T>(u.first) * id +
                         UniversalScalarConversion<Tint, T>(u.second) * ring.W;
      ring.l_unit.push_back(U);
    }
  }
  tr.expect("at_infinity");
  desc.at_infinity = tr.next();
  tr.expect("walk");
  desc.walk = tr.next();
  tr.expect("levels");
  int n_levels = tr.next_int();
  for (int index = 0; index < n_levels; index++) {
    tr.expect("level");
    int n_conf = tr.next_int();
    std::vector<MyMatrix<Tint>> l_conf;
    for (int i = 0; i < n_conf; i++) {
      MyMatrix<Tint> EXT = tr.next_matrix();
      if (EXT.cols() != n) {
        std::cerr << "PSYS: the configurations should have " << n << " columns\n";
        throw TerminalException{1};
      }
      l_conf.push_back(EXT);
    }
    desc.l_conf.push_back(l_conf);
  }
  return desc;
}

/*
  The configurations up to the units.
 */

inline bool lex_less(MyVector<Tint> const &v1, MyVector<Tint> const &v2) {
  for (int i = 0; i < v1.size(); i++) {
    if (v1(i) != v2(i)) {
      return v1(i) < v2(i);
    }
  }
  return false;
}

// The canonical representative of v up to the units: the smallest of the
// u v for the lexicographic order.
inline MyVector<Tint> unit_normalize(ProjectiveRing const &ring, MyVector<Tint> const &v) {
  MyVector<Tint> best = v;
  for (auto &U : ring.l_unit) {
    MyVector<Tint> w = U.transpose() * v;
    if (lex_less(w, best)) {
      best = w;
    }
  }
  return best;
}

// The key of a configuration taken up to the units of its vectors and to
// their order: the sorted list of the normalized vectors.
using ConfigKey = std::vector<int64_t>;

inline ConfigKey config_key(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  std::vector<std::vector<int64_t>> l_row;
  for (int i = 0; i < n_row; i++) {
    MyVector<Tint> v = unit_normalize(ring, GetMatrixRow(EXT, i));
    std::vector<int64_t> w(n_col);
    for (int j = 0; j < n_col; j++) {
      w[j] = v(j).get_si();
    }
    l_row.emplace_back(std::move(w));
  }
  std::sort(l_row.begin(), l_row.end());
  ConfigKey key;
  for (auto &w : l_row) {
    key.insert(key.end(), w.begin(), w.end());
  }
  return key;
}

// The position of a unit multiple of v among the rows of EXT.
inline std::optional<int> position_unit(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT,
                                        MyVector<Tint> const &v) {
  MyVector<Tint> v_norm = unit_normalize(ring, v);
  for (int i = 0; i < EXT.rows(); i++) {
    if (unit_normalize(ring, GetMatrixRow(EXT, i)) == v_norm) {
      return i;
    }
  }
  return {};
}

// The permutation of the rows of EXT (up to the units) induced by M, if M
// maps the configuration to itself.
inline std::optional<std::vector<Tidx>> permutation_unit(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT,
                                                         MyMatrix<Tint> const &M) {
  MyMatrix<Tint> EXTimg = EXT * M;
  std::vector<Tidx> l_pos;
  for (int i = 0; i < EXT.rows(); i++) {
    std::optional<int> opt = position_unit(ring, EXT, GetMatrixRow(EXTimg, i));
    if (!opt) {
      return {};
    }
    l_pos.push_back(*opt);
  }
  return l_pos;
}

// The unit multiples of the rows: EXT u_0, EXT u_1, ... (the row k n_ext + i
// is the row i times the unit k).
inline MyMatrix<Tint> totality(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  int n_row = EXT.rows();
  int n_unit = ring.l_unit.size();
  MyMatrix<Tint> EXTtot(n_unit * n_row, EXT.cols());
  for (int k = 0; k < n_unit; k++) {
    MyMatrix<Tint> EXTu = EXT * ring.l_unit[k];
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j < EXT.cols(); j++) {
        EXTtot(k * n_row + i, j) = EXTu(i, j);
      }
    }
  }
  return EXTtot;
}

// The rank over R of the span of the rows.
inline int ring_rank(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  return RankMat(totality(ring, EXT)) / ring.d;
}

// The rows v rho (rho in the basis 1, omega of R) of the rows l_row of EXT.
inline MyMatrix<Tint> real_expansion(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT,
                                     std::vector<int> const &l_row) {
  int d = ring.d;
  MyMatrix<Tint> M(d * l_row.size(), EXT.cols());
  for (size_t k = 0; k < l_row.size(); k++) {
    MyVector<Tint> v = GetMatrixRow(EXT, l_row[k]);
    for (int a = 0; a < d; a++) {
      MyVector<Tint> w = ring.l_ring_basis[a].transpose() * v;
      for (int j = 0; j < EXT.cols(); j++) {
        M(d * k + a, j) = w(j);
      }
    }
  }
  return M;
}

// The positions of rows of EXT forming a basis over R of their span
// (greedy).
inline std::vector<int> ring_basis_rows(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  std::vector<int> l_idx;
  int rnk = 0;
  for (int i = 0; i < EXT.rows(); i++) {
    std::vector<int> l_try = l_idx;
    l_try.push_back(i);
    if (RankMat(real_expansion(ring, EXT, l_try)) > rnk) {
      l_idx = l_try;
      rnk += ring.d;
    }
  }
  return l_idx;
}

// The positions of rows of EXT forming a basis of its span over Q.
inline std::vector<int> basis_rows(MyMatrix<Tint> const &EXT) {
  std::vector<int> l_idx;
  int rnk = 0;
  for (int i = 0; i < EXT.rows(); i++) {
    std::vector<int> l_try = l_idx;
    l_try.push_back(i);
    if (RankMat(select_rows(EXT, l_try)) > rnk) {
      l_idx = l_try;
      rnk += 1;
    }
  }
  return l_idx;
}

// Whether the rows of M (independent) span a saturated sublattice of Z^n.
inline bool is_saturated(MyMatrix<Tint> const &M) {
  MyMatrix<Tint> Msat = IntegralSpaceSaturation(M);
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
  for (int i = 0; i < Msat.rows(); i++) {
    MyVector<T> v = UniversalVectorConversion<T, Tint>(GetMatrixRow(Msat, i));
    std::optional<MyVector<T>> opt = SolutionMat(M_T, v);
    if (!opt || !IsIntegralVector(*opt)) {
      return false;
    }
  }
  return true;
}

/*
  The elements of G preserving a configuration spanning R^m, as the list of
  all of them (the group is finite) and a generating set.
 */
struct FiniteStabilizer {
  std::vector<MyMatrix<Tint>> l_elt;
  std::vector<MyMatrix<Tint>> l_gens;
};

inline FiniteStabilizer projective_stabilizer(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT,
                                              std::ostream &os) {
  MyMatrix<Tint> EXTtot = totality(ring, EXT);
  MyMatrix<T> EXTtot_T = UniversalMatrixConversion<T, Tint>(EXTtot);
  Tgroup GRP = LinPolytope_Automorphism<T, Tgroup>(EXTtot_T, os);
  FiniteStabilizer fs;
  for (auto &elt : GRP) {
    MyMatrix<T> M_T = FindTransformation(EXTtot_T, EXTtot_T, elt);
    if (ring.in_group(M_T)) {
      fs.l_elt.push_back(UniversalMatrixConversion<Tint, T>(M_T));
    }
  }
  // A generating set: the elements not in the group generated by the
  // previous ones (the groups are small).
  std::set<MyMatrix<Tint>> set_generated;
  set_generated.insert(IdentityMat<Tint>(ring.n()));
  for (auto &M : fs.l_elt) {
    if (set_generated.count(M) > 0) {
      continue;
    }
    fs.l_gens.push_back(M);
    std::vector<MyMatrix<Tint>> l_new(set_generated.begin(), set_generated.end());
    size_t pos = 0;
    while (pos < l_new.size()) {
      for (auto &g : fs.l_gens) {
        MyMatrix<Tint> h = l_new[pos] * g;
        if (set_generated.insert(h).second) {
          l_new.push_back(h);
        }
      }
      pos += 1;
    }
  }
  return fs;
}

/*
  The maps sending a basis over R of the span of eFace1 (rows of eFace1) to
  rows of the unit multiples of eFace2 with the same Gram matrix for B,
  by backtracking (the rows of eFace2 are chosen one at a time with the
  pairings of the previous ones). f_leaf receives the images S of the real
  expansion of the basis.
 */
inline void enumerate_isometric_images(ProjectiveRing const &ring, MyMatrix<Tint> const &eFace1,
                                       MyMatrix<Tint> const &eFace2,
                                       std::function<void(std::vector<int> const &, MyMatrix<Tint> const &)> const &f_leaf) {
  int d = ring.d;
  std::vector<int> l_basis = ring_basis_rows(ring, eFace1);
  int r = l_basis.size();
  MyMatrix<Tint> Bexp = real_expansion(ring, eFace1, l_basis);
  MyMatrix<Tint> GB = Bexp * ring.B * Bexp.transpose();
  MyMatrix<Tint> eFace2tot = totality(ring, eFace2);
  int nb2tot = eFace2tot.rows();
  std::vector<int> l_all(nb2tot);
  for (int c = 0; c < nb2tot; c++) {
    l_all[c] = c;
  }
  MyMatrix<Tint> Eexp = real_expansion(ring, eFace2tot, l_all);
  MyMatrix<Tint> G2 = Eexp * ring.B * Eexp.transpose();
  std::vector<int> sel(r);
  auto f_block_equal = [&](int k, int l, int c, int c2) -> bool {
    for (int a = 0; a < d; a++) {
      for (int b = 0; b < d; b++) {
        if (G2(d * c + a, d * c2 + b) != GB(d * k + a, d * l + b)) {
          return false;
        }
      }
    }
    return true;
  };
  std::function<void(int)> f_rec = [&](int k) -> void {
    if (k == r) {
      MyMatrix<Tint> S(d * r, ring.n());
      for (int l = 0; l < r; l++) {
        for (int a = 0; a < d; a++) {
          for (int j = 0; j < ring.n(); j++) {
            S(d * l + a, j) = Eexp(d * sel[l] + a, j);
          }
        }
      }
      f_leaf(l_basis, S);
      return;
    }
    for (int c = 0; c < nb2tot; c++) {
      if (!f_block_equal(k, k, c, c)) {
        continue;
      }
      bool is_ok = true;
      for (int l = 0; l < k; l++) {
        if (sel[l] % eFace2.rows() == c % eFace2.rows() || !f_block_equal(k, l, c, sel[l])) {
          is_ok = false;
          break;
        }
      }
      if (is_ok) {
        sel[k] = c;
        f_rec(k + 1);
      }
    }
  };
  f_rec(0);
}

/*
  The elements M of G with eFace1 * M contained in eFace2 (up to the
  units), for eFace1 spanning R^m, one per image (ListEmbedding): with B a
  basis of eFace1 and S its image, M = B^{-1} S is in G iff the Gram
  matrices of B and S for the form are equal, then M must be integral and
  map eFace1 into eFace2.
 */
inline std::vector<MyMatrix<Tint>> list_embedding(ProjectiveRing const &ring, MyMatrix<Tint> const &eFace1,
                                                  MyMatrix<Tint> const &eFace2) {
  std::vector<MyMatrix<Tint>> l_sol;
  std::set<ConfigKey> set_img;
  MyMatrix<T> Bexp_inv;
  bool has_inv = false;
  auto f_leaf = [&](std::vector<int> const &l_basis, MyMatrix<Tint> const &S) -> void {
    if (!has_inv) {
      Bexp_inv = Inverse(UniversalMatrixConversion<T, Tint>(real_expansion(ring, eFace1, l_basis)));
      has_inv = true;
    }
    MyMatrix<T> M_T = Bexp_inv * UniversalMatrixConversion<T, Tint>(S);
    if (!IsIntegralMatrix(M_T)) {
      return;
    }
    MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
    MyMatrix<Tint> EXTimg = eFace1 * M;
    for (int i = 0; i < EXTimg.rows(); i++) {
      if (!position_unit(ring, eFace2, GetMatrixRow(EXTimg, i))) {
        return;
      }
    }
    if (set_img.insert(config_key(ring, EXTimg)).second) {
      l_sol.push_back(M);
    }
  };
  enumerate_isometric_images(ring, eFace1, eFace2, f_leaf);
  return l_sol;
}

// The matrices M with eFace1 contained in eFace2 * M, one per eFace2 * M
// (ListEmbeddingRev).
inline std::vector<MyMatrix<Tint>> list_embedding_rev(ProjectiveRing const &ring, MyMatrix<Tint> const &eFace1,
                                                      MyMatrix<Tint> const &eFace2) {
  std::vector<MyMatrix<Tint>> l_ret;
  std::set<ConfigKey> set_key;
  for (auto &M : list_embedding(ring, eFace1, eFace2)) {
    MyMatrix<Tint> M_inv = Inverse(M);
    if (set_key.insert(config_key(ring, eFace2 * M_inv)).second) {
      l_ret.push_back(M_inv);
    }
  }
  return l_ret;
}

/*
  The subsets of the rows of eFace2 (a top cell) that are images of eFace1
  by an element of G (SP4Z_GetAllEmbeddings), as sorted lists of positions.
  For eFace1 spanning R^m, by list_embedding. Otherwise a basis over R of
  the span of eFace1 is mapped to rows of eFace2 (up to the units) with the
  same Gram matrix; the linear map must map eFace1 into eFace2 and the
  saturated lattice of the span of eFace1 onto a saturated lattice. Such
  an isometry of saturated sublattices is assumed to extend to G (true for
  the primitive isotropic sublattices of Sp_4(Z), and for the sublattices
  of rank 2 spanned by isotropic vectors of SU(2,1; O_K) of class number 1,
  whose classes are given by the pairings).
 */
inline std::vector<std::vector<int>> get_all_embeddings(ProjectiveRing const &ring, MyMatrix<Tint> const &eFace1,
                                                        MyMatrix<Tint> const &eFace2) {
  std::set<std::vector<int>> set_img;
  auto f_image_set = [&](MyMatrix<Tint> const &EXTimg) -> std::optional<std::vector<int>> {
    std::vector<int> l_pos;
    for (int i = 0; i < EXTimg.rows(); i++) {
      std::optional<int> opt = position_unit(ring, eFace2, GetMatrixRow(EXTimg, i));
      if (!opt) {
        return {};
      }
      l_pos.push_back(*opt);
    }
    std::sort(l_pos.begin(), l_pos.end());
    return l_pos;
  };
  if (ring_rank(ring, eFace1) == ring.m) {
    for (auto &M : list_embedding(ring, eFace1, eFace2)) {
      set_img.insert(*f_image_set(eFace1 * M));
    }
    return std::vector<std::vector<int>>(set_img.begin(), set_img.end());
  }
  std::vector<int> l_basis = ring_basis_rows(ring, eFace1);
  MyMatrix<T> Bexp = UniversalMatrixConversion<T, Tint>(real_expansion(ring, eFace1, l_basis));
  MyMatrix<Tint> EXT1tot = totality(ring, eFace1);
  MyMatrix<Tint> EXTsat = IntegralSpaceSaturation(select_rows(EXT1tot, basis_rows(EXT1tot)));
  // The coordinates of the rows of eFace1 and of the saturated basis in the
  // real expansion of the basis.
  auto f_coord = [&](MyMatrix<Tint> const &M) -> MyMatrix<T> {
    MyMatrix<T> C(M.rows(), Bexp.rows());
    for (int i = 0; i < M.rows(); i++) {
      std::optional<MyVector<T>> opt = SolutionMat(Bexp, UniversalVectorConversion<T, Tint>(GetMatrixRow(M, i)));
      MyVector<T> c = unfold_opt(opt, "PSYS: the rows should be in the span of the basis");
      for (int j = 0; j < Bexp.rows(); j++) {
        C(i, j) = c(j);
      }
    }
    return C;
  };
  MyMatrix<T> coord = f_coord(eFace1);
  MyMatrix<T> coord_sat = f_coord(EXTsat);
  auto f_leaf = [&]([[maybe_unused]] std::vector<int> const &l_b, MyMatrix<Tint> const &S) -> void {
    MyMatrix<T> S_T = UniversalMatrixConversion<T, Tint>(S);
    MyMatrix<T> EXTimg_T = coord * S_T;
    MyMatrix<T> SATimg_T = coord_sat * S_T;
    if (!IsIntegralMatrix(EXTimg_T) || !IsIntegralMatrix(SATimg_T)) {
      return;
    }
    if (!is_saturated(UniversalMatrixConversion<Tint, T>(SATimg_T))) {
      return;
    }
    std::optional<std::vector<int>> opt = f_image_set(UniversalMatrixConversion<Tint, T>(EXTimg_T));
    if (opt) {
      set_img.insert(*opt);
    }
  };
  enumerate_isometric_images(ring, eFace1, eFace2, f_leaf);
  return std::vector<std::vector<int>>(set_img.begin(), set_img.end());
}

/*
  The neighbors of the representative top cells across their ridges
  (SP4Z_GetFlippingInfo): for the top cell iPerf, the ridge (of the index
  1) given by eSet, the neighbor is the top cell jPerf mapped by eMatB.
 */
struct FlipInfo {
  int iRidge;
  std::vector<int> eSet;
  int jPerf;
  MyMatrix<Tint> eMatB;
};

inline std::vector<std::vector<FlipInfo>> get_flipping_info(ProjectiveRing const &ring,
                                                            std::vector<std::vector<MyMatrix<Tint>>> const &l_conf,
                                                            [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrix<Tint>> const &l_perf = l_conf[0];
  std::vector<MyMatrix<Tint>> const &l_ridge = l_conf[1];
  std::vector<std::vector<FlipInfo>> ll_adj;
  for (size_t iPerf = 0; iPerf < l_perf.size(); iPerf++) {
    MyMatrix<Tint> const &ePerf = l_perf[iPerf];
    ConfigKey key_perf = config_key(ring, ePerf);
    std::vector<FlipInfo> l_adj;
    for (size_t iRidge = 0; iRidge < l_ridge.size(); iRidge++) {
      for (auto &eSet : get_all_embeddings(ring, l_ridge[iRidge], ePerf)) {
        MyMatrix<Tint> eFace1 = select_rows(ePerf, eSet);
        std::vector<std::pair<int, MyMatrix<Tint>>> l_cont;
        for (size_t jPerf = 0; jPerf < l_perf.size(); jPerf++) {
          for (auto &eMatB : list_embedding_rev(ring, eFace1, l_perf[jPerf])) {
            if (config_key(ring, l_perf[jPerf] * eMatB) != key_perf) {
              l_cont.push_back({jPerf, eMatB});
            }
          }
        }
        if (l_cont.size() != 1) {
          std::cerr << "PSYS: get_flipping_info, the ridge should have one neighbor, |l_cont|=" << l_cont.size()
                    << "\n";
          throw TerminalException{1};
        }
        l_adj.push_back({static_cast<int>(iRidge), eSet, l_cont[0].first, l_cont[0].second});
      }
    }
#ifdef DEBUG_PROJECTIVE_COMPLEX
    os << "PSYS: iPerf=" << iPerf << " |l_adj|=" << l_adj.size() << "\n";
#endif
    ll_adj.emplace_back(std::move(l_adj));
  }
  return ll_adj;
}

/*
  A top cell of the complex: the representative top cell iPerf mapped by
  g, with its stabilizer in G, the permutation group on the rows of the
  unit multiples of EXT and its neighbors.
 */
struct TopCell {
  int iPerf;
  MyMatrix<Tint> g;
  MyMatrix<Tint> EXT;
  std::vector<MyMatrix<Tint>> l_stab_gens;
  Tgroup perm_grp;
  struct Adj {
    int iCell;
    MyMatrix<Tint> eEquiv;
    std::vector<int> eSet;
  };
  std::vector<Adj> l_adj;
  // Rows of EXTtot forming a basis, and the inverse of that basis.
  std::vector<int> l_basis;
  MyMatrix<T> basis_inv;
  MyMatrix<Tint> EXTtot;
  // For the walk "potential": a vector of the cell with B(X, X) < 0, the
  // sum of unit multiples of its vectors.
  MyVector<Tint> center;
  // The matrix of a permutation of the rows of EXTtot in perm_grp: the
  // basis rows go to their images.
  MyMatrix<Tint> matrix_of(Telt const &elt) const {
    std::vector<int> l_img;
    for (auto &i : l_basis) {
      l_img.push_back(elt.at(i));
    }
    MyMatrix<T> M_T = basis_inv * UniversalMatrixConversion<T, Tint>(select_rows(EXTtot, l_img));
    return UniversalMatrixConversion<Tint, T>(M_T);
  }
};

// The permutation of the rows of EXTtot induced by M.
inline Telt permutation_totality(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &M) {
  MyMatrix<Tint> EXTtot = totality(ring, EXT);
  MyMatrix<Tint> EXTimg = EXTtot * M;
  int n = EXTtot.rows();
  std::map<MyVector<Tint>, int, std::function<bool(MyVector<Tint> const &, MyVector<Tint> const &)>> map_row(
      [](MyVector<Tint> const &a, MyVector<Tint> const &b) -> bool { return lex_less(a, b); });
  for (int k = 0; k < n; k++) {
    map_row[GetMatrixRow(EXTtot, k)] = k;
  }
  std::vector<Tidx> l_pos(n);
  for (int i = 0; i < n; i++) {
    auto iter = map_row.find(GetMatrixRow(EXTimg, i));
    if (iter == map_row.end()) {
      std::cerr << "PSYS: permutation_totality, M does not preserve the configuration\n";
      throw TerminalException{1};
    }
    l_pos[i] = iter->second;
  }
  return Telt(l_pos);
}

/*
  A vector X = sum_k u_k x_k of the cell (x_k its vectors, u_k units) with
  B(X, X) minimal, for the walk "potential": for a form of signature
  (2, 1) it must be negative, a point of the symmetric space (the negative
  lines) attached to the cell.
 */
inline MyVector<Tint> get_cell_center(ProjectiveRing const &ring, MyMatrix<Tint> const &EXT) {
  int n_ext = EXT.rows();
  int n_unit = ring.l_unit.size();
  std::vector<int> choice(n_ext, 0);
  MyVector<Tint> best;
  Tint val_best(0);
  bool has_best = false;
  while (true) {
    MyVector<Tint> X = ZeroVector<Tint>(ring.n());
    for (int k = 0; k < n_ext; k++) {
      X += ring.l_unit[choice[k]].transpose() * GetMatrixRow(EXT, k);
    }
    Tint val = X.dot(ring.B * X);
    if (!has_best || val < val_best) {
      best = X;
      val_best = val;
      has_best = true;
    }
    // The next choice, the first unit fixed to 1.
    int k = 1;
    while (k < n_ext && choice[k] == n_unit - 1) {
      choice[k] = 0;
      k += 1;
    }
    if (k >= n_ext) {
      break;
    }
    choice[k] += 1;
  }
  if (val_best >= 0) {
    std::cerr << "PSYS: get_cell_center, no vector of negative norm in the cell\n";
    throw TerminalException{1};
  }
  return best;
}

/*
  The orbits of top cells under G (SP4Z_GetTopCells with the group G
  itself), by flipping across the ridges, and the generators of G found on
  the way (the stabilizers and the equivalences).
 */
struct TopCellResult {
  std::vector<TopCell> l_top;
  std::vector<MyMatrix<Tint>> l_gens_G;
};

inline TopCellResult get_top_cells(ProjectiveRing const &ring, std::vector<std::vector<MyMatrix<Tint>>> const &l_conf,
                                   std::vector<std::vector<FlipInfo>> const &ll_adj, std::ostream &os) {
  std::vector<MyMatrix<Tint>> const &l_perf = l_conf[0];
  std::vector<FiniteStabilizer> l_stab;
  for (auto &ePerf : l_perf) {
    l_stab.push_back(projective_stabilizer(ring, ePerf, os));
  }
  TopCellResult res;
  std::set<MyMatrix<Tint>> set_gens;
  MyMatrix<Tint> id = IdentityMat<Tint>(ring.n());
  auto f_insert_gen = [&](MyMatrix<Tint> const &M) -> void {
    if (M != id && set_gens.insert(M).second) {
      res.l_gens_G.push_back(M);
    }
  };
  struct InsertResult {
    int iCell;
    MyMatrix<Tint> eEquiv;
  };
  // Two top cells of the same type are equivalent: c1.g^{-1} g.
  auto f_insert = [&](int const &iPerf, MyMatrix<Tint> const &g) -> InsertResult {
    for (size_t iCell = 0; iCell < res.l_top.size(); iCell++) {
      if (res.l_top[iCell].iPerf == iPerf) {
        MyMatrix<Tint> test = Inverse(res.l_top[iCell].g) * g;
        f_insert_gen(test);
        return {static_cast<int>(iCell), test};
      }
    }
    TopCell cell;
    cell.iPerf = iPerf;
    cell.g = g;
    cell.EXT = l_perf[iPerf] * g;
    MyMatrix<Tint> g_inv = Inverse(g);
    std::vector<Telt> l_perm;
    for (auto &s : l_stab[iPerf].l_gens) {
      MyMatrix<Tint> s_conj = g_inv * s * g;
      cell.l_stab_gens.push_back(s_conj);
      f_insert_gen(s_conj);
      l_perm.push_back(permutation_totality(ring, cell.EXT, s_conj));
    }
    cell.EXTtot = totality(ring, cell.EXT);
    Telt id_perm(cell.EXTtot.rows());
    cell.perm_grp = Tgroup(l_perm, id_perm);
    cell.l_basis = basis_rows(cell.EXTtot);
    cell.basis_inv = Inverse(UniversalMatrixConversion<T, Tint>(select_rows(cell.EXTtot, cell.l_basis)));
    res.l_top.emplace_back(std::move(cell));
    return {static_cast<int>(res.l_top.size()) - 1, id};
  };
  f_insert(0, id);
  size_t pos = 0;
  while (pos < res.l_top.size()) {
    int iPerf = res.l_top[pos].iPerf;
    MyMatrix<Tint> g = res.l_top[pos].g;
    std::vector<TopCell::Adj> l_adj;
    for (auto &adj : ll_adj[iPerf]) {
      InsertResult ins = f_insert(adj.jPerf, adj.eMatB * g);
      l_adj.push_back({ins.iCell, ins.eEquiv, adj.eSet});
    }
    res.l_top[pos].l_adj = l_adj;
    pos += 1;
  }
#ifdef DEBUG_PROJECTIVE_COMPLEX
  os << "PSYS: |l_top|=" << res.l_top.size() << " |l_gens_G|=" << res.l_gens_G.size() << "\n";
#endif
  return res;
}

// A cell given as the subset eSet of the rows of the top cell iTopCell,
// mapped by eMat: its vectors are l_top[iTopCell].EXT{eSet} * eMat.
struct Triple {
  int iTopCell;
  std::vector<int> eSet;
  MyMatrix<Tint> eMat;
};

// The face of EXTtot of the unit multiples of the rows eSet.
inline Face face_totality(TopCell const &cell, std::vector<int> const &eSet) {
  int n_ext = cell.EXT.rows();
  int n_unit = cell.EXTtot.rows() / n_ext;
  Face f(cell.EXTtot.rows());
  for (auto &i : eSet) {
    for (int k = 0; k < n_unit; k++) {
      f[k * n_ext + i] = 1;
    }
  }
  return f;
}

/*
  The equivalence of two triples (SP4Z_TestEquivalenceTriple): an element
  M of G with EXT1 * M = EXT2 (up to the units) obtained from the
  stabilizer of their top cell.
 */
inline std::optional<MyMatrix<Tint>> test_equivalence_triple(std::vector<TopCell> const &l_top, Triple const &t1,
                                                             Triple const &t2) {
  if (t1.iTopCell != t2.iTopCell || t1.eSet.size() != t2.eSet.size()) {
    return {};
  }
  TopCell const &cell = l_top[t1.iTopCell];
  Face f1 = face_totality(cell, t1.eSet);
  Face f2 = face_totality(cell, t2.eSet);
  std::optional<Telt> opt = cell.perm_grp.RepresentativeAction_OnSets(f1, f2);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> M = Inverse(t1.eMat) * cell.matrix_of(*opt) * t2.eMat;
  return M;
}

/*
  An orbit of cells: its triples (the ways the cell appears in the top
  cells, up to the equivalence), the generators of its stabilizer, its
  vectors, its boundary (the faces, with the elements mapping the
  representative of their orbit to them and the signs) and its
  orientation character: the stabilizer acts on the vectors (up to the
  units) by permutations with signs; the group perm_or on the vectors and
  two extra points, generated by the permutations extended by the
  transposition of the extra points for the sign -1, gives the sign of any
  element of the stabilizer.
 */
struct CellOrbit {
  std::vector<Triple> l_triple;
  std::vector<MyMatrix<Tint>> l_stab_gens;
  MyMatrix<Tint> EXT;
  std::vector<BoundEntry<Tint>> l_bound;
  std::vector<std::vector<int>> l_set_m1;
  Tgroup perm_or;
  Tgroup perm_grp;
  // The order of the stabilizer in G, 0 when it is infinite (more than
  // max_finite_order elements).
  int64_t stab_order = 0;
};

/*
  The order of the group generated by the matrices l_gens, by the closure
  of the products, or 0 when it has more than max_order elements (the
  stabilizer of a cell at infinity is infinite). The kernel of the action
  on the vectors up to the units can be larger than the scalar units, for
  a cell not spanning R^m, so the order is not the one of perm_grp times
  the number of units.
 */
inline int64_t matrix_group_order(std::vector<MyMatrix<Tint>> const &l_gens, int const &n, size_t const &max_order) {
  std::set<MyMatrix<Tint>> set_elt;
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  set_elt.insert(id);
  std::vector<MyMatrix<Tint>> l_elt{id};
  size_t pos = 0;
  while (pos < l_elt.size()) {
    for (auto &g : l_gens) {
      MyMatrix<Tint> h = l_elt[pos] * g;
      if (set_elt.insert(h).second) {
        l_elt.push_back(h);
        if (l_elt.size() > max_order) {
          return 0;
        }
      }
    }
    pos += 1;
  }
  return l_elt.size();
}

// The sign of an element g of the stabilizer of the cell.
inline int orbit_orientation(ProjectiveRing const &ring, CellOrbit const &orbit, MyMatrix<Tint> const &g) {
  std::optional<std::vector<Tidx>> opt = permutation_unit(ring, orbit.EXT, g);
  if (!opt) {
    std::cerr << "PSYS: orbit_orientation, the element does not stabilize the cell\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> l_pos = *opt;
  int n = orbit.EXT.rows();
  l_pos.push_back(n);
  l_pos.push_back(n + 1);
  Telt elt(l_pos);
  if (orbit.perm_or.isin(elt)) {
    return 1;
  }
  return -1;
}

// The sign changes of the faces (UntangleAllSigns): each codimension two
// face occurs in two faces, which gives a relation between their signs.
struct SignOccurrence {
  int sign;
  int idx;
};

inline std::vector<int> untangle_all_signs(std::vector<std::vector<SignOccurrence>> const &l_case, int const &n_var) {
  std::vector<int> v_sign(n_var, 0);
  std::vector<int> assigned(n_var, 0);
  v_sign[0] = 1;
  assigned[0] = 1;
  std::vector<int> used(l_case.size(), 0);
  while (true) {
    int n_unused = 0;
    for (size_t i = 0; i < l_case.size(); i++) {
      std::vector<SignOccurrence> const &c = l_case[i];
      if (c.size() != 2) {
        std::cerr << "PSYS: untangle_all_signs, a face of codimension two occurs " << c.size() << " times\n";
        throw TerminalException{1};
      }
      int sum = assigned[c[0].idx] + assigned[c[1].idx];
      if (used[i] == 0 && sum > 0) {
        n_unused += 1;
        if (sum == 2) {
          if (c[0].sign * v_sign[c[0].idx] + c[1].sign * v_sign[c[1].idx] != 0) {
            std::cerr << "PSYS: untangle_all_signs, inconsistent signs\n";
            throw TerminalException{1};
          }
        } else {
          int ref = (assigned[c[0].idx] == 1) ? c[0].idx : c[1].idx;
          int oth = (assigned[c[0].idx] == 1) ? c[1].idx : c[0].idx;
          v_sign[oth] = -c[0].sign * c[1].sign * v_sign[ref];
          assigned[oth] = 1;
        }
        used[i] = 1;
      }
    }
    if (n_unused == 0) {
      break;
    }
  }
  for (auto &s : v_sign) {
    if (s != 1 && s != -1) {
      std::cerr << "PSYS: untangle_all_signs, a sign is not determined\n";
      throw TerminalException{1};
    }
  }
  return v_sign;
}

/*
  The cell complex (SP4Z_GetCellComplex): the orbits of cells of each
  level with their stabilizers, boundaries and orientations, and for each
  top cell its faces of each level (the triples with the orbit and the
  element mapping the representative to them).
 */
struct FaceIndex {
  std::vector<int> eSet;
  int iOrbit;
  MyMatrix<Tint> eEquiv;
};

struct ProjectiveCellComplex {
  ProjectiveDescription desc;
  std::vector<TopCell> l_top;
  std::vector<MyMatrix<Tint>> l_gens_G;
  // levels[index][iOrbit]
  std::vector<std::vector<CellOrbit>> levels;
  // face_index[index][iTopCell]
  std::vector<std::vector<std::vector<FaceIndex>>> face_index;
  ProjectiveRing const &ring() const { return desc.ring; }
};

inline ProjectiveCellComplex compute_projective_cell_complex(ProjectiveDescription const &desc, std::ostream &os) {
#ifdef TIMINGS_PROJECTIVE_COMPLEX
  MicrosecondTime time;
#endif
  ProjectiveCellComplex cx;
  cx.desc = desc;
  ProjectiveRing const &ring = cx.desc.ring;
  std::vector<std::vector<MyMatrix<Tint>>> const &l_conf = cx.desc.l_conf;
  int n = ring.n();
  int n_levels = l_conf.size();
  std::vector<std::vector<FlipInfo>> ll_adj = get_flipping_info(ring, l_conf, os);
  TopCellResult tcr = get_top_cells(ring, l_conf, ll_adj, os);
  cx.l_top = std::move(tcr.l_top);
  cx.l_gens_G = std::move(tcr.l_gens_G);
#ifdef TIMINGS_PROJECTIVE_COMPLEX
  os << "|PSYS: top cells|=" << time << "\n";
#endif
  if (cx.desc.walk == "potential") {
    for (auto &cell : cx.l_top) {
      cell.center = get_cell_center(ring, cell.EXT);
    }
  }
  std::vector<TopCell> const &l_top = cx.l_top;
  int n_top = l_top.size();
  // The saturation of a triple: the triples of the same cell in the top
  // cells containing it, and the stabilizer.
  auto f_saturation = [&](Triple const &t0) -> CellOrbit {
    CellOrbit orbit;
    orbit.EXT = select_rows(l_top[t0.iTopCell].EXT, t0.eSet) * t0.eMat;
    std::set<MyMatrix<Tint>> set_gens;
    auto f_insert_gen = [&](MyMatrix<Tint> const &M) -> void {
#ifdef SANITY_CHECK_PROJECTIVE_COMPLEX
      if (!permutation_unit(ring, orbit.EXT, M)) {
        std::cerr << "PSYS: a generator of the stabilizer does not stabilize the cell\n";
        throw TerminalException{1};
      }
#endif
      if (set_gens.insert(M).second) {
        orbit.l_stab_gens.push_back(M);
      }
    };
    TopCell const &cell0 = l_top[t0.iTopCell];
    Face f = face_totality(cell0, t0.eSet);
    MyMatrix<Tint> eMat_inv = Inverse(t0.eMat);
    for (auto &elt : cell0.perm_grp.Stabilizer_OnSets(f).GeneratorsOfGroup()) {
      f_insert_gen(eMat_inv * cell0.matrix_of(elt) * t0.eMat);
    }
    auto f_insert_triple = [&](Triple const &t) -> void {
      for (auto &t_old : orbit.l_triple) {
        std::optional<MyMatrix<Tint>> opt = test_equivalence_triple(l_top, t, t_old);
        if (opt) {
          f_insert_gen(*opt);
          return;
        }
      }
      orbit.l_triple.push_back(t);
    };
    f_insert_triple(t0);
    size_t pos = 0;
    while (pos < orbit.l_triple.size()) {
      Triple t = orbit.l_triple[pos];
      for (auto &adj : l_top[t.iTopCell].l_adj) {
        if (!std::includes(adj.eSet.begin(), adj.eSet.end(), t.eSet.begin(), t.eSet.end())) {
          continue;
        }
        MyMatrix<Tint> fMat = adj.eEquiv * t.eMat;
        MyMatrix<Tint> EXTimg = l_top[adj.iCell].EXT * fMat;
        std::vector<int> fSet;
        for (int i = 0; i < orbit.EXT.rows(); i++) {
          std::optional<int> opt = position_unit(ring, EXTimg, GetMatrixRow(orbit.EXT, i));
          if (!opt) {
            std::cerr << "PSYS: the cell is not in the adjacent top cell\n";
            throw TerminalException{1};
          }
          fSet.push_back(*opt);
        }
        std::sort(fSet.begin(), fSet.end());
        f_insert_triple({adj.iCell, fSet, fMat});
      }
      pos += 1;
    }
    return orbit;
  };
  cx.levels = std::vector<std::vector<CellOrbit>>(n_levels);
  cx.face_index = std::vector<std::vector<std::vector<FaceIndex>>>(n_levels);
  for (int index = 0; index < n_levels; index++) {
    std::vector<CellOrbit> &l_orbit = cx.levels[index];
    auto f_insert_orbit = [&](Triple const &t) -> std::pair<int, MyMatrix<Tint>> {
      for (size_t iOrbit = 0; iOrbit < l_orbit.size(); iOrbit++) {
        for (auto &t_orb : l_orbit[iOrbit].l_triple) {
          std::optional<MyMatrix<Tint>> opt = test_equivalence_triple(l_top, t_orb, t);
          if (opt) {
            return {static_cast<int>(iOrbit), *opt};
          }
        }
      }
      l_orbit.push_back(f_saturation(t));
      return {static_cast<int>(l_orbit.size()) - 1, IdentityMat<Tint>(n)};
    };
    cx.face_index[index] = std::vector<std::vector<FaceIndex>>(n_top);
    for (int iTopCell = 0; iTopCell < n_top; iTopCell++) {
      MyMatrix<Tint> const &ePerf = l_conf[0][l_top[iTopCell].iPerf];
      for (auto &eFace : l_conf[index]) {
        for (auto &eSet : get_all_embeddings(ring, eFace, ePerf)) {
          Triple t{iTopCell, eSet, IdentityMat<Tint>(n)};
          std::pair<int, MyMatrix<Tint>> pair = f_insert_orbit(t);
          cx.face_index[index][iTopCell].push_back({eSet, pair.first, pair.second});
        }
      }
    }
#ifdef DEBUG_PROJECTIVE_COMPLEX
    os << "PSYS: index=" << index << " |l_orbit|=" << l_orbit.size() << "\n";
#endif
  }
#ifdef TIMINGS_PROJECTIVE_COMPLEX
  os << "|PSYS: orbits of cells|=" << time << "\n";
#endif
  /*
    The boundaries, from the single vectors (the last index) upwards, with
    the orientation characters.
   */
  for (int index = n_levels - 1; index >= 0; index--) {
    for (auto &orbit : cx.levels[index]) {
      Triple const &t0 = orbit.l_triple[0];
      std::vector<std::vector<Tidx>> l_perm_gens;
      for (auto &g : orbit.l_stab_gens) {
        l_perm_gens.push_back(*permutation_unit(ring, orbit.EXT, g));
      }
      int n_ext = orbit.EXT.rows();
      std::vector<int> l_sign_gens(orbit.l_stab_gens.size(), 1);
      if (index < n_levels - 1) {
        // The faces: the cells of the next index contained in the cell, in
        // its top cell.
        std::vector<int> l_iface;
        std::vector<MyMatrix<Tint>> l_elt;
        for (auto &fi : cx.face_index[index + 1][t0.iTopCell]) {
          if (std::includes(t0.eSet.begin(), t0.eSet.end(), fi.eSet.begin(), fi.eSet.end())) {
            l_iface.push_back(fi.iOrbit);
            l_elt.push_back(fi.eEquiv);
            std::vector<int> set_m1;
            for (auto &x : fi.eSet) {
              set_m1.push_back(std::lower_bound(t0.eSet.begin(), t0.eSet.end(), x) - t0.eSet.begin());
            }
            orbit.l_set_m1.push_back(set_m1);
          }
        }
        int n_face = l_iface.size();
        std::vector<int> l_sign;
        if (index == n_levels - 2) {
          if (n_face != 2) {
            std::cerr << "PSYS: a cell of two vectors should have two faces\n";
            throw TerminalException{1};
          }
          l_sign = {1, -1};
        } else {
          // The faces of codimension two, each occurring in two faces.
          std::vector<std::vector<int>> l_set_m2;
          std::vector<MyMatrix<Tint>> l_elt_m2;
          std::vector<std::vector<SignOccurrence>> l_case;
          for (int i_face = 0; i_face < n_face; i_face++) {
            CellOrbit const &orbit_m1 = cx.levels[index + 1][l_iface[i_face]];
            for (auto &eb : orbit_m1.l_bound) {
              MyMatrix<Tint> add_elt = eb.M * l_elt[i_face];
              CellOrbit const &orbit_m2 = cx.levels[index + 2][eb.iOrb];
              MyMatrix<Tint> EXTimg = orbit_m2.EXT * add_elt;
              std::vector<int> set_img;
              for (int i = 0; i < EXTimg.rows(); i++) {
                set_img.push_back(*position_unit(ring, orbit.EXT, GetMatrixRow(EXTimg, i)));
              }
              std::sort(set_img.begin(), set_img.end());
              auto iter = std::find(l_set_m2.begin(), l_set_m2.end(), set_img);
              if (iter == l_set_m2.end()) {
                l_set_m2.push_back(set_img);
                l_elt_m2.push_back(add_elt);
                l_case.push_back({{eb.sign, i_face}});
              } else {
                size_t pos = iter - l_set_m2.begin();
                MyMatrix<Tint> quot = l_elt_m2[pos] * Inverse(add_elt);
                int mul_sign = orbit_orientation(ring, orbit_m2, quot);
                l_case[pos].push_back({eb.sign * mul_sign, i_face});
              }
            }
          }
          l_sign = untangle_all_signs(l_case, n_face);
        }
        for (int i_face = 0; i_face < n_face; i_face++) {
          orbit.l_bound.push_back({l_iface[i_face], l_sign[i_face], l_elt[i_face]});
        }
        // The signs of the generators of the stabilizer.
        for (size_t i_gen = 0; i_gen < orbit.l_stab_gens.size(); i_gen++) {
          MyMatrix<Tint> const &g = orbit.l_stab_gens[i_gen];
          if (index == n_levels - 2) {
            // The sign of the permutation of the two vectors.
            l_sign_gens[i_gen] = (l_perm_gens[i_gen][0] == 0) ? 1 : -1;
            continue;
          }
          std::vector<int> const &set1 = orbit.l_set_m1[0];
          MyMatrix<Tint> EXT1img = select_rows(orbit.EXT, set1) * g;
          std::vector<int> set1_img;
          for (int i = 0; i < EXT1img.rows(); i++) {
            set1_img.push_back(*position_unit(ring, orbit.EXT, GetMatrixRow(EXT1img, i)));
          }
          std::sort(set1_img.begin(), set1_img.end());
          auto iter = std::find(orbit.l_set_m1.begin(), orbit.l_set_m1.end(), set1_img);
          if (iter == orbit.l_set_m1.end()) {
            std::cerr << "PSYS: the image of a face is not a face\n";
            throw TerminalException{1};
          }
          size_t pos = iter - orbit.l_set_m1.begin();
          MyMatrix<Tint> quot = orbit.l_bound[0].M * g * Inverse(orbit.l_bound[pos].M);
          CellOrbit const &orbit_m1 = cx.levels[index + 1][orbit.l_bound[0].iOrb];
          int sign = orbit_orientation(ring, orbit_m1, quot);
          l_sign_gens[i_gen] = orbit.l_bound[0].sign * orbit.l_bound[pos].sign * sign;
        }
      }
      // The orientation group: the permutations extended by the sign.
      std::vector<Telt> l_ext_perm, l_perm;
      for (size_t i_gen = 0; i_gen < orbit.l_stab_gens.size(); i_gen++) {
        std::vector<Tidx> v = l_perm_gens[i_gen];
        l_perm.push_back(Telt(v));
        if (l_sign_gens[i_gen] == 1) {
          v.push_back(n_ext);
          v.push_back(n_ext + 1);
        } else {
          v.push_back(n_ext + 1);
          v.push_back(n_ext);
        }
        l_ext_perm.push_back(Telt(v));
      }
      orbit.perm_or = Tgroup(l_ext_perm, Telt(n_ext + 2));
      orbit.perm_grp = Tgroup(l_perm, Telt(n_ext));
      // The scalar units are in the stabilizer, generated by the elements
      // found from the top cells only if the cell spans R^m.
      std::vector<MyMatrix<Tint>> l_gens_unit = orbit.l_stab_gens;
      l_gens_unit.insert(l_gens_unit.end(), ring.l_unit.begin(), ring.l_unit.end());
      orbit.stab_order = matrix_group_order(l_gens_unit, n, 20000);
    }
  }
#ifdef TIMINGS_PROJECTIVE_COMPLEX
  os << "|PSYS: boundaries|=" << time << "\n";
#endif
  return cx;
}

/*
  The complex as step B of src_perfect sees it (see PerfectComplexView in
  hecke_homology.h).
 */
struct ProjectiveComplexView {
  ProjectiveCellComplex const &cx;
  int dim() const { return cx.ring().n(); }
  int n_levels() const { return cx.levels.size(); }
  int n_orbit(int const &index) const { return cx.levels[index].size(); }
  std::vector<MyMatrix<Tint>> const &stabilizer_generators(int const &index, int const &iOrb) const {
    return cx.levels[index][iOrb].l_stab_gens;
  }
  int orientation(int const &index, int const &iOrb, MyMatrix<Tint> const &g) const {
    return orbit_orientation(cx.ring(), cx.levels[index][iOrb], g);
  }
  bool is_well_rounded(int const &index, int const &iOrb) const {
    MyMatrix<Tint> const &EXT = cx.levels[index][iOrb].EXT;
    if (cx.desc.at_infinity == "single") {
      return EXT.rows() > 1;
    }
    return ring_rank(cx.ring(), EXT) == cx.ring().m;
  }
  int64_t stabilizer_order(int const &index, int const &iOrb) const {
    int64_t order = cx.levels[index][iOrb].stab_order;
    if (order == 0) {
      std::cerr << "PSYS: the stabilizer of a well rounded cell should be finite\n";
      throw TerminalException{1};
    }
    return order;
  }
  std::vector<BoundEntry<Tint>> const &boundary(int const &index, int const &iOrb) const {
    return cx.levels[index][iOrb].l_bound;
  }
  std::vector<MyMatrix<Tint>> ambient_generators([[maybe_unused]] std::ostream &os) const { return cx.l_gens_G; }
  bool in_ambient_group(MyMatrix<T> const &g) const { return cx.ring().in_group(g); }
  bool has_cells_at_infinity() const { return true; }
};

}  // namespace psys

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_PROJECTIVE_COMPLEX_H_
// clang-format on
