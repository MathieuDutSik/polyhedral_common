// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_POLYTOPE_ENUMERATION_H_
#define SRC_ERDAHL_ERDAHL_POLYTOPE_ENUMERATION_H_

// clang-format off
#include "erdahl_enumeration.h"
#include <algorithm>
#include <optional>
#include <random>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_POLYTOPE_ENUMERATION
#endif

#ifdef TIMINGS
#define TIMINGS_ERDAHL_POLYTOPE_ENUMERATION
#endif

/*
  Enumeration of the perfect Delaunay polytopes relative to a space W by
  restricting to the polytopes (degeneracy rank 0), as in the paper of
  Dutour, Erdahl and Rybnikov on perfect Delaunay polytopes in low
  dimensions.

  The perfect Delaunay polyhedra of rank 1 are the extreme rays of
  Erdahl(n) cap W of Delaunay type, and two of them are adjacent if they
  contain a common sub-polyhedron of rank 2. For a perfect polytope D, the
  sub-polyhedra of rank 2 are the facets of its cone of evaluations
  (erdahl_polytope.h) and the flip with D3 = Z^n gives the adjacent perfect
  polyhedra. Only those that are polytopes are kept: the recursion through
  the degenerate polyhedra such as {0,1} x Z^{n-1} is avoided, at the price
  of the guarantee of completeness (the graph restricted to the polytopes
  may not be connected).

  More generally the perfect polyhedra D = P + L with dim L <= max_dim_L
  are kept. For dim L > 0 the sub-polyhedra of rank 2 are obtained by the
  recursive method (erdahl_sub_delaunay), which only goes through
  polyhedra of degeneracy rank < dim L. For max_dim_L = 1 this connects
  for example the perfect polytopes through the P x Z with P a perfect
  polytope of dimension n-1.

  A starting perfect polytope is obtained from any Delaunay polytope by
  moving its function f along a direction g of Space_W(D) until new zeros
  appear, which lowers the rank, keeping only the moves that give a
  polytope.
 */

// Z^n as a Delaunay polyhedron: the zero function.
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint> erdahl_whole_lattice(int n, std::ostream &os) {
  return erdahl_polyhedron_from_function<T, Tint>(ZeroMatrix<T>(n + 1, n + 1),
                                                  os);
}

// The simplex {0, e_1, ..., e_n}: sum_i x_i (x_i - 1) + s (s - 1) with
// s = sum_i x_i is nonnegative on Z^n and vanishes exactly on it.
template <typename T> MyMatrix<T> erdahl_simplex_function(int n) {
  MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      F(i + 1, j + 1) = (i == j) ? T(2) : T(1);
    }
    F(0, i + 1) = T(-1);
    F(i + 1, 0) = T(-1);
  }
  return F;
}

// The cube {0,1}^n: sum_i x_i (x_i - 1), of center (1/2, ..., 1/2).
template <typename T> MyMatrix<T> erdahl_cube_function(int n) {
  MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
  for (int i = 0; i < n; i++) {
    F(i + 1, i + 1) = 1;
    F(0, i + 1) = T(-1) / T(2);
    F(i + 1, 0) = T(-1) / T(2);
  }
  return F;
}

/*
  One move from the polytope D of rank r > 1 along the direction g of
  Space_W(D): the function f + lambda g with lambda maximal such that it
  stays nonnegative on Z^n. Its zero set contains D strictly and has a
  smaller rank. Returns it if it is a polytope.
 */
template <typename T, typename Tint>
std::optional<DelaunayPolyhedron<T, Tint>>
erdahl_move_along(ErdahlFunctionSpace<T> const &W,
                  DelaunayPolyhedron<T, Tint> const &D,
                  DelaunayPolyhedron<T, Tint> const &Zn, MyMatrix<T> const &g,
                  std::vector<MyVector<Tint>> const &l_near,
                  std::ostream &os) {
  // f + lambda g = f - lambda h with h = -g, lambda = min f(v) / h(v).
  MyMatrix<T> h = -g;
  std::vector<MyVector<Tint>> l_cand;
  for (auto &e : l_near) {
    if (EvaluationQuadForm<T, Tint>(h, e) > 0) {
      l_cand.push_back(e);
    }
  }
  std::optional<MyMatrix<T>> opt_G =
      erdahl_minimal_ratio_function<T, Tint>(D.F, h, Zn, l_cand, os);
  if (!opt_G) {
#ifdef DEBUG_ERDAHL_POLYTOPE_ENUMERATION
    os << "ERDAHL: move, no candidate with g < 0 among |l_near|="
       << l_near.size() << "\n";
#endif
    return {};
  }
  auto res = erdahl_zero_set<T, Tint>(*opt_G, os);
  if (res.is_err()) {
    std::cerr << "ERDAHL: the moved function should be nonnegative\n";
    throw TerminalException{1};
  }
  ErdahlLatticeSet<Tint> const &zs = res.get_ok();
#ifdef DEBUG_ERDAHL_POLYTOPE_ENUMERATION
  os << "ERDAHL: move, zero set |EXT|=" << zs.EXT.rows()
     << " dim L=" << zs.L.rows()
     << " full_dim=" << erdahl_is_full_dimensional(zs) << "\n";
#endif
  if (zs.L.rows() > 0 || !erdahl_is_full_dimensional(zs)) {
    return {};
  }
  ErdahlLatticeSet<Tint> can = erdahl_canonical_lattice_set(zs);
  MyMatrix<T> F = erdahl_canonical_function<T, Tint>(W, can, nullptr, os);
  return DelaunayPolyhedron<T, Tint>{can.EXT, can.L, F};
}

/*
  The polytope D of center C moved to the center c (both in (1/2) Z^n but
  not in Z^n) by x -> x A + b with A unimodular mapping 2C to 2c modulo 2.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_move_to_center(DelaunayPolyhedron<T, Tint> const &D,
                      MyVector<T> const &C, MyVector<T> const &c) {
  int n = c.size();
  auto get_parity = [&](MyVector<T> const &v) -> MyVector<Tint> {
    MyVector<Tint> u(n);
    for (int i = 0; i < n; i++) {
      T two_v = 2 * v(i);
      if (!IsInteger(two_v)) {
        std::cerr << "ERDAHL: the center should be half-integral\n";
        throw TerminalException{1};
      }
      u(i) = ResInt(UniversalScalarConversion<Tint, T>(two_v), Tint(2));
    }
    if (IsZeroVector(u)) {
      std::cerr << "ERDAHL: the center should not be a lattice point\n";
      throw TerminalException{1};
    }
    return u;
  };
  MyMatrix<Tint> A =
      Inverse(erdahl_unimodular_with_first_row(get_parity(C))) *
      erdahl_unimodular_with_first_row(get_parity(c));
  MyMatrix<T> A_T = UniversalMatrixConversion<T, Tint>(A);
  MyVector<T> img = A_T.transpose() * C;
  MyMatrix<Tint> g = IdentityMat<Tint>(n + 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g(1 + i, 1 + j) = A(i, j);
    }
    g(0, 1 + i) = UniversalScalarConversion<Tint, T>(c(i) - img(i));
  }
  return erdahl_apply_transformation(D, g);
}

/*
  The Erdahl-Rybnikov function on Z^n, n >= 6, whose zero set is the
  asymmetric perfect polytope P_ER(n) (the Gosset polytope 2_21 for n = 6,
  the 35-tope for n = 7). It is x^t Q x + l^t x with
    Q_ii = n^2 - 5n + 2, Q_ij = n^2 - 7n + 12 for 1 <= i < j < n,
    Q_nn = n^2 - 3n - 2, Q_in = n^2 - 5n + 4,
    l_i = -Q_ii, l_n = -(Q_nn - 4),
  see V. P. Grishukhin, Infinite series of extreme Delaunay polytopes,
  European J. Combin. 27 (2006), equation (6).
 */
template <typename T> MyMatrix<T> erdahl_er_function(int n) {
  T n_T(n);
  T q_ii = n_T * n_T - 5 * n_T + 2;
  T q_ij = n_T * n_T - 7 * n_T + 12;
  T q_nn = n_T * n_T - 3 * n_T - 2;
  T q_in = n_T * n_T - 5 * n_T + 4;
  MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - 1; j++) {
      F(1 + i, 1 + j) = (i == j) ? q_ii : q_ij;
    }
    F(1 + i, n) = q_in;
    F(n, 1 + i) = q_in;
    F(0, 1 + i) = -q_ii / 2;
    F(1 + i, 0) = -q_ii / 2;
  }
  F(n, n) = q_nn;
  F(0, n) = -(q_nn - 4) / 2;
  F(n, 0) = -(q_nn - 4) / 2;
  return F;
}

/*
  The symmetrization of Lemma 15.3.7 of Deza and Laurent, Geometry of Cuts
  and Metrics, in the form of section 3 of the paper of Grishukhin above.
  For the function f(x) = q(x - c_P) - r^2 of a Delaunay polytope P of Z^m
  and a point c_F of (1/2) Z^m, the function on Z^{m+1}
    f(x, t) = q(x - c_F - (t - 1/2) a) + gamma t (t - 1) - r^2
  with a = 2 (c_F - c_P) vanishes on the layers P x {0} and
  (2 c_F - P) x {1}, is symmetric around (c_F, 1/2), and is positive on the
  other layers for gamma > r^2 / 2.
 */
template <typename T>
MyMatrix<T> erdahl_symmetrized_function(MyMatrix<T> const &F,
                                        MyVector<T> const &c_F) {
  int m = F.rows() - 1;
  MyMatrix<T> Q(m, m);
  MyVector<T> b(m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Q(i, j) = F(1 + i, 1 + j);
    }
    b(i) = F(1 + i, 0);
  }
  MyMatrix<T> Qinv = Inverse(Q);
  MyVector<T> c_P = -(Qinv * b);
  T r2 = b.dot(Qinv * b) - F(0, 0);
  MyVector<T> a = 2 * (c_F - c_P);
  T gamma = r2 + 1;
  // y = x - c_F - (t - 1/2) a as a linear function of (1, x, t).
  MyMatrix<T> M = ZeroMatrix<T>(m, m + 2);
  for (int i = 0; i < m; i++) {
    M(i, 0) = -c_F(i) + a(i) / 2;
    M(i, 1 + i) = 1;
    M(i, m + 1) = -a(i);
  }
  MyMatrix<T> Fs = M.transpose() * Q * M;
  Fs(m + 1, m + 1) += gamma;
  Fs(0, m + 1) -= gamma / 2;
  Fs(m + 1, 0) -= gamma / 2;
  Fs(0, 0) -= r2;
  return Fs;
}

/*
  A centrally symmetric perfect polytope of Z^{m+1} with a section equal
  to the perfect polytope P of function F on Z^m. The two layer polytope
  of erdahl_symmetrized_function has perfection rank 2 relative to the
  functions of center (c_F, 1/2): its space is spanned by its function and
  t (t - 1). Decreasing gamma until new zeros appear gives a perfect
  polytope if they are not on a line. The points c_F tried are the
  midpoints of the pairs of vertices of P, i.e. the centers of the
  centrally symmetric faces and more, by increasing q(c_F - c_P).
 */
template <typename T, typename Tint>
std::optional<std::pair<DelaunayPolyhedron<T, Tint>, MyVector<T>>>
erdahl_symmetrized_perfect_polytope(MyMatrix<T> const &F, std::ostream &os) {
  int m = F.rows() - 1;
  DelaunayPolyhedron<T, Tint> P =
      erdahl_polyhedron_from_function<T, Tint>(F, os);
  int n_vert = P.EXT.rows();
  MyMatrix<T> Q(m, m);
  MyVector<T> b(m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Q(i, j) = F(1 + i, 1 + j);
    }
    b(i) = F(1 + i, 0);
  }
  MyVector<T> c_P = -(Inverse(Q) * b);
  std::vector<std::pair<T, MyVector<T>>> l_cand;
  std::set<MyVector<T>> set_cF;
  for (int i = 0; i < n_vert; i++) {
    for (int j = i; j < n_vert; j++) {
      MyVector<T> c_F(m);
      for (int k = 0; k < m; k++) {
        c_F(k) = UniversalScalarConversion<T, Tint>(P.EXT(i, 1 + k) +
                                                    P.EXT(j, 1 + k)) /
                 2;
      }
      if (set_cF.insert(c_F).second) {
        MyVector<T> d = c_F - c_P;
        l_cand.push_back({d.dot(Q * d), c_F});
      }
    }
  }
  std::sort(l_cand.begin(), l_cand.end(),
            [](auto const &x, auto const &y) { return x.first < y.first; });
  int n = m + 1;
  DelaunayPolyhedron<T, Tint> Zn = erdahl_whole_lattice<T, Tint>(n, os);
  // -t (t - 1): the decrease of gamma.
  MyMatrix<T> g = ZeroMatrix<T>(n + 1, n + 1);
  g(n, n) = -1;
  g(0, n) = T(1) / T(2);
  g(n, 0) = T(1) / T(2);
  for (auto &cand : l_cand) {
    MyVector<T> const &c_F = cand.second;
    MyVector<T> C(n);
    for (int k = 0; k < m; k++) {
      C(k) = c_F(k);
    }
    C(m) = T(1) / T(2);
    ErdahlFunctionSpace<T> W_C = erdahl_centered_space(C);
    DelaunayPolyhedron<T, Tint> D = erdahl_polyhedron_from_function<T, Tint>(
        erdahl_symmetrized_function(F, c_F), os);
    // The points of the layers t = -1 and t = 2 above the vertices.
    std::vector<MyVector<Tint>> l_near;
    for (int i = 0; i < D.EXT.rows(); i++) {
      MyVector<Tint> e = GetMatrixRow(D.EXT, i);
      MyVector<Tint> e_up = e;
      e_up(n) += 2;
      l_near.push_back(e_up);
      MyVector<Tint> e_down = e;
      e_down(n) -= 2;
      l_near.push_back(e_down);
    }
    std::optional<DelaunayPolyhedron<T, Tint>> opt =
        erdahl_move_along(W_C, D, Zn, g, l_near, os);
#ifdef DEBUG_ERDAHL_POLYTOPE_ENUMERATION
    os << "ERDAHL: symmetrization, q(c_F - c_P)=" << cand.first
       << " |EXT|=" << D.EXT.rows() << " polytope=" << opt.has_value();
    if (opt) {
      os << " |EXT(move)|=" << opt->EXT.rows()
         << " rank=" << erdahl_perfection_rank(W_C, *opt);
    }
    os << "\n";
#endif
    if (opt && erdahl_perfection_rank(W_C, *opt) == 1) {
      return std::make_pair(*opt, C);
    }
  }
  return {};
}

/*
  A perfect Delaunay polytope of W to start from, if there is one:
  * n = 1: the segment {0, 1}.
  * Full space, n >= 6: the polytope P_ER(n). There are none for
    2 <= n <= 5 (Erdahl).
  * Space of functions of center c, n >= 7: the symmetrization of
    P_ER(n-1), moved to the center c (the Gosset polytope 3_21 for n = 7,
    the 72-tope of the 35-tope for n = 8). A centrally symmetric perfect
    polytope is perfect in the full space, so there are none for
    2 <= n <= 6 (2_21 is the only perfect polytope of dimension 6).
 */
template <typename T, typename Tint>
std::optional<DelaunayPolyhedron<T, Tint>>
erdahl_initial_polytope(ErdahlFunctionSpace<T> const &W, std::ostream &os) {
  int n = W.n;
  auto get_polytope =
      [&](DelaunayPolyhedron<T, Tint> const &D) -> DelaunayPolyhedron<T, Tint> {
    if (!erdahl_is_in_space(W, D.F)) {
      std::cerr << "ERDAHL: the initial polytope is not in the space\n";
      throw TerminalException{1};
    }
    MyMatrix<T> F = erdahl_canonical_function<T, Tint>(
        W, erdahl_lattice_set(D), nullptr, os);
    return DelaunayPolyhedron<T, Tint>{D.EXT, D.L, F};
  };
  if (!W.center) {
    erdahl_check_supported_space(W);
    if (n == 1) {
      return get_polytope(erdahl_polyhedron_from_function<T, Tint>(
          erdahl_simplex_function<T>(n), os));
    }
    if (n <= 5) {
      return {};
    }
    return get_polytope(
        erdahl_polyhedron_from_function<T, Tint>(erdahl_er_function<T>(n), os));
  }
  MyVector<T> const &c = *W.center;
  if (n == 1) {
    MyVector<T> C(1);
    C(0) = T(1) / T(2);
    DelaunayPolyhedron<T, Tint> D = erdahl_polyhedron_from_function<T, Tint>(
        erdahl_cube_function<T>(1), os);
    return get_polytope(erdahl_move_to_center(D, C, c));
  }
  if (n <= 6) {
    return {};
  }
  auto opt = erdahl_symmetrized_perfect_polytope<T, Tint>(
      erdahl_er_function<T>(n - 1), os);
  if (!opt) {
    std::cerr << "ERDAHL: the symmetrization of P_ER(n-1) failed\n";
    throw TerminalException{1};
  }
  return get_polytope(erdahl_move_to_center(opt->first, opt->second, c));
}

/*
  A perfect polytope relative to W containing the polytope D, obtained by
  successive moves. The directions tried are the basis of Space_W(D) and
  their opposites, then random integral combinations of them.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_perfect_polytope_from(ErdahlFunctionSpace<T> const &W,
                             DelaunayPolyhedron<T, Tint> D, std::ostream &os) {
  int n = W.n;
  DelaunayPolyhedron<T, Tint> Zn = erdahl_whole_lattice<T, Tint>(n, os);
  std::mt19937 rng(1);
  while (true) {
    if (D.L.rows() > 0) {
      std::cerr << "ERDAHL: the starting point should be a polytope\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<T>> basis =
        erdahl_vanishing_space(W, erdahl_lattice_set(D));
    int r = basis.size();
#ifdef DEBUG_ERDAHL_POLYTOPE_ENUMERATION
    os << "ERDAHL: perfection, |EXT|=" << D.EXT.rows() << " rank=" << r << "\n";
#endif
    if (r == 1) {
      return D;
    }
    // The points near D where the directions are evaluated first.
    std::vector<MyVector<Tint>> l_near;
    for (int i = 0; i < D.EXT.rows(); i++) {
      MyVector<Tint> e = GetMatrixRow(D.EXT, i);
      for (int k = 0; k < n; k++) {
        for (int sign = -1; sign <= 1; sign += 2) {
          MyVector<Tint> f = e;
          f(1 + k) += sign;
          l_near.push_back(f);
        }
      }
    }
    std::optional<std::vector<std::vector<int64_t>>> pool =
        erdahl_triangle_pool(erdahl_lattice_set(D));
    if (pool) {
      for (auto &pt : *pool) {
        MyVector<Tint> e(n + 1);
        for (int a = 0; a <= n; a++) {
          e(a) = UniversalScalarConversion<Tint, int64_t>(pt[a]);
        }
        l_near.push_back(e);
      }
    }
    std::vector<MyMatrix<T>> l_dir;
    for (auto &g : basis) {
      l_dir.push_back(g);
      l_dir.push_back(-g);
    }
    std::uniform_int_distribution<int> dist(-3, 3);
    for (int i_rand = 0; i_rand < 20 * r; i_rand++) {
      MyMatrix<T> g = ZeroMatrix<T>(n + 1, n + 1);
      for (auto &b : basis) {
        g += T(dist(rng)) * b;
      }
      l_dir.push_back(g);
    }
    bool progress = false;
    for (auto &g : l_dir) {
      // A direction proportional to the function of D does not move it.
      MyMatrix<T> Mpair(2, (n + 1) * (n + 2) / 2);
      AssignMatrixRow(Mpair, 0, SymmetricMatrixToVector(D.F));
      AssignMatrixRow(Mpair, 1, SymmetricMatrixToVector(g));
      if (RankMat(Mpair) < 2) {
        continue;
      }
      std::optional<DelaunayPolyhedron<T, Tint>> opt =
          erdahl_move_along(W, D, Zn, g, l_near, os);
      if (opt && erdahl_perfection_rank(W, *opt) < r) {
        D = *opt;
        progress = true;
        break;
      }
    }
    if (!progress) {
      std::cerr << "ERDAHL: no move from the polytope gives a polytope of "
                   "smaller rank\n";
      throw TerminalException{1};
    }
  }
}

template <typename T, typename Tint> struct ErdahlPerfectPolytopes {
  std::vector<DelaunayPolyhedron<T, Tint>> l_perfect;
  // The flips leading to a perfect polyhedron with dim L > max_dim_L,
  // which are discarded.
  size_t n_degenerate;
  size_t n_flip;
};

/*
  The orbits of perfect Delaunay polyhedra relative to W with
  dim L <= max_dim_L that are connected to the perfect polytope Dstart
  through such polyhedra.
 */
template <typename T, typename Tint, typename Tgroup>
ErdahlPerfectPolytopes<T, Tint>
erdahl_enumerate_perfect_polytopes(ErdahlFunctionSpace<T> const &W,
                                   DelaunayPolyhedron<T, Tint> const &Dstart,
                                   int max_dim_L,
                                   std::string const &FileDualDesc,
                                   std::ostream &os) {
  int n = W.n;
  if (erdahl_perfection_rank(W, Dstart) != 1 || Dstart.L.rows() > 0) {
    std::cerr << "ERDAHL: the starting point should be a perfect polytope\n";
    throw TerminalException{1};
  }
  DelaunayPolyhedron<T, Tint> Zn = erdahl_whole_lattice<T, Tint>(n, os);
  // The sub-polyhedra computed by erdahl_sub_delaunay, shared between the
  // perfect polyhedra.
  ErdahlBank<T, Tint> bank{{}, FileDualDesc};
  struct Entry {
    DelaunayPolyhedron<T, Tint> D;
    bool done;
    ErdahlChainConfig<T, Tint> cfg;
  };
  std::vector<Entry> l_entry;
  std::unordered_map<size_t, std::vector<size_t>> map_hash;
  std::vector<ErdahlSuperData<T, Tint>> l_sd;
  auto insert = [&](DelaunayPolyhedron<T, Tint> Dnew) -> void {
    ErdahlChainConfig<T, Tint> cfg = erdahl_chain_config(Dnew, l_sd);
    size_t hash = erdahl_invariant_hash(cfg, os);
    for (auto &idx : map_hash[hash]) {
      if (erdahl_equivalence<T, Tint, Tgroup>(W, l_entry[idx].D,
                                               l_entry[idx].cfg, Dnew, cfg,
                                               l_sd, os)) {
        return;
      }
    }
    erdahl_ensure_function(W, Dnew, os);
    map_hash[hash].push_back(l_entry.size());
    l_entry.push_back({Dnew, false, cfg});
#ifdef DEBUG_ERDAHL_POLYTOPE_ENUMERATION
    os << "ERDAHL: new perfect polyhedron, |EXT|=" << Dnew.EXT.rows()
       << " dim L=" << Dnew.L.rows() << " n_orbit=" << l_entry.size() << "\n";
#endif
  };
  insert(Dstart);
  size_t n_degenerate = 0;
  size_t n_flip = 0;
  while (true) {
    int i_sel = -1;
    for (size_t i = 0; i < l_entry.size(); i++) {
      if (!l_entry[i].done) {
        i_sel = i;
        break;
      }
    }
    if (i_sel == -1) {
      break;
    }
    l_entry[i_sel].done = true;
    DelaunayPolyhedron<T, Tint> D = l_entry[i_sel].D;
#ifdef TIMINGS_ERDAHL_POLYTOPE_ENUMERATION
    MicrosecondTime time;
#endif
    std::vector<DelaunayPolyhedron<T, Tint>> l_sub =
        erdahl_sub_delaunay<T, Tint, Tgroup>(W, D, bank, os);
#ifdef TIMINGS_ERDAHL_POLYTOPE_ENUMERATION
    size_t n_deg_loc = 0;
#endif
    for (auto &S : l_sub) {
      DelaunayPolyhedron<T, Tint> Dadj = erdahl_flip(W, S, D, Zn, false, os);
      n_flip++;
      if (Dadj.L.rows() > max_dim_L) {
        n_degenerate++;
#ifdef TIMINGS_ERDAHL_POLYTOPE_ENUMERATION
        n_deg_loc++;
#endif
        continue;
      }
      insert(Dadj);
    }
#ifdef TIMINGS_ERDAHL_POLYTOPE_ENUMERATION
    os << "|ERDAHL: perfect polyhedron |EXT|=" << D.EXT.rows()
       << " dim L=" << D.L.rows() << " |sub|=" << l_sub.size() << " n_degenerate=" << n_deg_loc
       << " n_orbit=" << l_entry.size() << "|=" << time << "\n";
#endif
  }
  std::vector<DelaunayPolyhedron<T, Tint>> l_perfect;
  for (auto &entry : l_entry) {
    l_perfect.push_back(entry.D);
  }
  return {l_perfect, n_degenerate, n_flip};
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_POLYTOPE_ENUMERATION_H_
// clang-format on
