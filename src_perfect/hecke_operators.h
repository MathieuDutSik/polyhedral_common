// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_OPERATORS_H_
#define SRC_PERFECT_HECKE_OPERATORS_H_

// clang-format off
#include "perfect_complex.h"
#include <functional>
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
  obtained by the contracting homotopy. The solution is searched among
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
  get_key(g) returns a complete invariant of the coset g Gamma, which
  is what the coset table is keyed by.
 */
template <typename Tint> struct FiniteIndexSubgroup {
  std::string type;
  Tint N;
  int n;
  bool is_in_gamma(MyMatrix<Tint> const &g) const {
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
  // A complete invariant of the coset g Gamma.
  MyMatrix<Tint> get_key(MyMatrix<Tint> const &g) const {
    if (type == "Full") {
      return ZeroMatrix<Tint>(1, 1);
    }
    if (type == "Principal") {
      MyMatrix<Tint> key(n, n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          key(i, j) = ResInt(g(i, j), N);
        }
      }
      return key;
    }
    if (type == "Gamma0" || type == "Gamma1") {
      // g Gamma = h Gamma iff e_n h^{-1} g = lambda e_n mod N, that is
      // iff the rows e_n g^{-1} and e_n h^{-1} agree modulo N up to a
      // unit lambda (lambda = 1 for Gamma1).
      MyMatrix<Tint> gInv = Inverse(g);
      MyMatrix<Tint> row(1, n);
      for (int j = 0; j < n; j++) {
        row(0, j) = ResInt(gInv(n - 1, j), N);
      }
      if (type == "Gamma1") {
        return row;
      }
      // For Gamma0 we take the lexicographically smallest multiple by a
      // unit modulo N.
      MyMatrix<Tint> best = row;
      bool has_best = false;
      for (Tint u = 1; u < N; u++) {
        Tint g_val = GcdPair(u, N);
        if (g_val != 1) {
          continue;
        }
        MyMatrix<Tint> cand(1, n);
        for (int j = 0; j < n; j++) {
          Tint val = u * row(0, j);
          cand(0, j) = ResInt(val, N);
        }
        auto is_smaller = [&]() -> bool {
          for (int j = 0; j < n; j++) {
            if (cand(0, j) < best(0, j)) {
              return true;
            }
            if (cand(0, j) > best(0, j)) {
              return false;
            }
          }
          return false;
        };
        if (!has_best || is_smaller()) {
          best = cand;
          has_best = true;
        }
      }
      return best;
    }
    std::cerr << "HECKE: Unsupported subgroup type=" << type << "\n";
    throw TerminalException{1};
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
  std::unordered_map<MyMatrix<Tint>, int> map_key;
  size_t size() const { return l_coset.size(); }
  std::optional<int> find_index(MyMatrix<Tint> const &g) const {
    MyMatrix<Tint> key = gamma.get_key(g);
    auto iter = map_key.find(key);
    if (iter == map_key.end()) {
      return {};
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
  auto f_insert = [&](MyMatrix<Tint> const &g) -> void {
    MyMatrix<Tint> key = gamma.get_key(g);
    if (ct.map_key.contains(key)) {
      return;
    }
    int idx = ct.l_coset.size();
    ct.map_key[key] = idx;
    ct.l_coset_inv.push_back(Inverse(g));
    ct.l_coset.push_back(g);
  };
  f_insert(IdentityMat<Tint>(n));
  size_t pos = 0;
  while (pos < ct.l_coset.size()) {
    // The vector may be reallocated by f_insert, so we copy.
    MyMatrix<Tint> g = ct.l_coset[pos];
    for (auto &s : l_gens) {
      MyMatrix<Tint> prod = s * g;
      f_insert(prod);
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
                       [[maybe_unused]] std::ostream &os) {
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
  while (pos < hc.l_coset.size()) {
    MyMatrix<T> y = hc.l_coset[pos];
    for (auto &s : l_gens) {
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
std::function<bool(MyMatrix<Tint> const &)>
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
    return [](MyMatrix<Tint> const &) -> bool { return true; };
  }
  // The rows e of Equa satisfy V e = 0.
  MyMatrix<T> Equa = NullspaceTrMat(V);
  MyMatrix<T> EquaT = Equa.transpose();
  return [EquaT](MyMatrix<Tint> const &EXT) -> bool {
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    MyMatrix<T> prod = EXT_T * EquaT;
    return IsZeroMatrix(prod);
  };
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
        std::cerr << "HECKE: The walk towards the cell failed to terminate\n";
        throw TerminalException{1};
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
      while (true) {
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

template <typename T, typename Tint> struct HeckeChainMap {
  MyMatrix<T> x;
  HeckeCosets<T> cosets;
  // images[index][iOrb][j] is the chain of level index representing the
  // image of the cell F_{iOrb} under the coset y_j.
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images;
};

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
      return;
    }
    map_key[key] = l_elt.size();
    int sign = get_face_orientation(EXT, ListMat, face.or_info, h);
    l_elt.emplace_back(h, h_conj);
    l_sign.push_back(sign);
  };
  f_insert(id, id);
  size_t pos = 0;
  while (pos < l_elt.size()) {
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
  T inv_size = T(1) / T(n_elt);
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
  if (!pco.compute_contracting_homotopy) {
    std::cerr << "HECKE: The Hecke operators require ComputeContractingHomotopy = T\n";
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
          std::function<bool(MyMatrix<Tint> const &)> f_filter =
              get_subspace_filter(index + 1, rhs, fce, os);
          sol = contracting_homotopy(index + 1, rhs, fce, f_filter, os);
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
  return {x, std::move(cosets), std::move(images)};
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
    for (size_t i = 0; i < n_coset; i++) {
      MyMatrix<Tint> prod = s * ct.l_coset[i];
      perm[i_gen][i] = ct.index_of(prod);
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
  std::vector<MyMatrix<T>> boundary;
  // The position and sign of the cell F_{iOrb} M of level index in the
  // basis, if it is present.
  std::optional<std::pair<int, int>>
  cell_position(int const &index, int const &iOrb, MyMatrix<Tint> const &M) const {
    int i0 = ct.index_of(M);
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
                      bool const &only_well_rounded, std::ostream &os) {
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
    int n_row = gc.basis[index].size();
    int n_col = gc.basis[index + 1].size();
    MyMatrix<T> D = ZeroMatrix<T>(n_row, n_col);
    for (int i_row = 0; i_row < n_row; i_row++) {
      int iOrb = gc.basis[index][i_row].first;
      int i_coset = gc.basis[index][i_row].second;
      MyMatrix<Tint> const &g = ct.l_coset[i_coset];
      for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
        MyMatrix<Tint> M = ebnd.M * g;
        std::optional<std::pair<int, int>> opt = gc.cell_position(index + 1, ebnd.iOrb, M);
        if (opt) {
          D(i_row, opt->first) += T(ebnd.sign * opt->second);
        }
      }
    }
    gc.boundary.push_back(D);
  }
#ifdef SANITY_CHECK_HECKE_OPERATORS
  for (int index = 0; index < n_levels - 2; index++) {
    MyMatrix<T> prod = gc.boundary[index] * gc.boundary[index + 1];
    if (!IsZeroMatrix(prod)) {
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
};

template <typename T> struct HeckeHomologyResult {
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
template <typename T, typename Tint, typename Tgroup>
MyVector<T> get_hecke_image_cell(
    int const &index, int const &i_cell, GammaComplex<T, Tint> const &gc,
    HeckeChainMap<T, Tint> const &hcm, HeckeCosets<T> const &cosets_gamma,
    [[maybe_unused]] FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  int n_cell = gc.basis[index].size();
  MyVector<T> V = ZeroVector<T>(n_cell);
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
        V(opt->first) += fe.value * T(opt->second);
      }
    }
  }
  return V;
}

template <typename T, typename Tint, typename Tgroup>
HeckeHomologyResult<T> compute_hecke_homology(
    FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeChainMap<T, Tint> const &hcm,
    FiniteIndexSubgroup<Tint> const &gamma, bool const &only_well_rounded,
    std::ostream &os) {
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
  std::vector<MyMatrix<Tint>> l_gens_gamma = get_schreier_generators(ct, l_gens_G, n, os);
  std::function<bool(MyMatrix<T> const &)> f_member_gamma =
      [&LinSpa, &gamma](MyMatrix<T> const &g) -> bool {
    if (!is_in_ambient_group(LinSpa, g)) {
      return false;
    }
    MyMatrix<Tint> g_int = UniversalMatrixConversion<Tint, T>(g);
    return gamma.is_in_gamma(g_int);
  };
  HeckeCosets<T> cosets_gamma =
      enumerate_hecke_cosets<T, Tint>(hcm.x, l_gens_gamma, f_member_gamma, os);
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: enumerate_hecke_cosets (Gamma)|=" << time << "\n";
#endif
  GammaComplex<T, Tint> gc = compute_gamma_complex(fce, ct, only_well_rounded, os);
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: compute_gamma_complex|=" << time << "\n";
#endif
  HeckeHomologyResult<T> result;
  result.subgroup_type = gamma.type;
  result.only_well_rounded = only_well_rounded;
  result.n_coset_gamma = ct.size();
  result.n_hecke_coset = cosets_gamma.size();
  int n_levels = fce.levels.size();
  for (int index = 0; index < n_levels; index++) {
    int n_cell = gc.basis[index].size();
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
        if (index < n_levels - 1) {
          MyVector<T> wD = gc.boundary[index].transpose() * w;
          if (!IsZeroVector(wD)) {
            std::cerr << "HECKE: The image of a cycle should be a cycle, index=" << index << "\n";
            throw TerminalException{1};
          }
        }
        MyVector<T> diff = w - harmonic.transpose() * proj;
        if (!IsZeroVector(diff)) {
          if (index == 0) {
            std::cerr << "HECKE: The difference should be zero at index 0\n";
            throw TerminalException{1};
          }
          std::optional<MyVector<T>> opt = SolutionMat(gc.boundary[index - 1], diff);
          if (!opt) {
            std::cerr << "HECKE: The difference should be a boundary, index=" << index << "\n";
            throw TerminalException{1};
          }
        }
#endif
      }
    }
    HeckeHomologyLevel<T> level{index, n_cell, std::move(harmonic), std::move(hecke_matrix)};
    result.l_level.emplace_back(std::move(level));
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: homology level " << index << "|=" << time << "\n";
#endif
  }
  return result;
}

template <typename T>
void WriteHeckeHomologyResultGAP(std::ostream &os_out, HeckeHomologyResult<T> const &result) {
  os_out << "rec(subgroup_type:=\"" << result.subgroup_type << "\"";
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
    os_out << ")";
  }
  os_out << "])";
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_OPERATORS_H_
// clang-format on
