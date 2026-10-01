// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_ENUMERATION_H_
#define SRC_ERDAHL_ERDAHL_ENUMERATION_H_

// clang-format off
#include "erdahl_flip.h"
#include "erdahl_group.h"
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_ENUMERATION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ERDAHL_ENUMERATION
#endif

#ifdef TIMINGS
#define TIMINGS_ERDAHL_ENUMERATION
#endif

/*
  The recursive adjacency decomposition method of the paper (Algorithm
  "Enumeration_inequivalent", the GAP function InfDel_SubDelaunayPolytopes)
  relative to a space W.

  Input: a Delaunay polyhedron D of perfection rank r in W.
  Output: the orbits under Aut_W(D) of the Delaunay polyhedra D' subset D
  of rank r+1 (their functions being in W).

  * If L(D) = 0, they are the facets of the cone of evaluations
    (erdahl_polytope.h). The facets whose zero sets are not full
    dimensional are not Delaunay polyhedra and are left out.
  * Otherwise, start from D_init = D cap { phi in {0,1} } for phi an
    integral affine function nonconstant on L(D) with phi(phi - 1) in W,
    which is the image of the D_init of the paper by an element of Aff(D)
    and so has rank r+1 and degeneracy rank d-1. Then, as long as there is
    an untreated orbit D' of degeneracy rank < d:
    - compute the orbits of the sub-polyhedra D1 of D' of rank r+2 under
      Aut_W(D') (recursive call),
    - split them into orbits under Stab_W(D', D),
    - for each of them, the flip gives the other D'' with D1 subset D''
      subset D of rank r+1, inserted if new up to Aut_W(D).
    By Theorem "ErdahlConnectivityResult" this gives all the orbits.

  The results of the recursive calls are kept in a bank, indexed by the
  polyhedron up to the transformations preserving W.

  The perfect Delaunay polyhedra relative to W are the rank 1
  sub-polyhedra of Z^n = Z(0), which has rank 0.
 */

template <typename T, typename Tint> struct ErdahlBankEntry {
  DelaunayPolyhedron<T, Tint> D;
  // erdahl_invariant_hash of D without supers.
  size_t hash;
  std::vector<DelaunayPolyhedron<T, Tint>> l_sub;
};

template <typename T, typename Tint> struct ErdahlBank {
  std::vector<ErdahlBankEntry<T, Tint>> l_entry;
  // The heuristics file of the dual descriptions ("unset" for the default).
  std::string FileDualDesc;
};

template <typename T, typename Tint>
std::pair<int, int> erdahl_size_key(DelaunayPolyhedron<T, Tint> const &D) {
  return {static_cast<int>(D.L.rows()), static_cast<int>(D.EXT.rows())};
}

// The function phi (phi - 1) for phi(x) = a + w.x.
template <typename T, typename Tint>
MyMatrix<T> erdahl_slab_function(T const &a, MyVector<Tint> const &w) {
  MyVector<T> w_T = UniversalVectorConversion<T, Tint>(w);
  MyMatrix<T> Q = w_T * w_T.transpose();
  MyVector<T> lin = ((2 * a - 1) / 2) * w_T;
  T cst = a * a - a;
  return erdahl_assemble_function(cst, lin, Q);
}

/*
  D_init = D cap { phi in {0,1} } with phi(x) = a + w.x, w.L(D) = Z and
  phi(phi - 1) in W. For the centered spaces this means phi(c) = 1/2.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_initial_sub_polyhedron(ErdahlFunctionSpace<T> const &W,
                              DelaunayPolyhedron<T, Tint> const &D,
                              std::ostream &os) {
  int n = erdahl_dimension(D);
  int d = D.L.rows();
  int r = erdahl_perfection_rank(W, D);
  MyVector<Tint> x0(n);
  for (int i = 0; i < n; i++) {
    x0(i) = D.EXT(0, i + 1);
  }
  auto try_vector = [&](MyVector<Tint> const &w)
      -> std::optional<DelaunayPolyhedron<T, Tint>> {
    Tint gcd(0);
    for (int i = 0; i < d; i++) {
      Tint scal(0);
      for (int j = 0; j < n; j++) {
        scal += w(j) * D.L(i, j);
      }
      gcd = GcdPair(gcd, scal);
    }
    if (T_abs(gcd) != 1) {
      return {};
    }
    T a;
    if (W.center) {
      MyVector<T> w_T = UniversalVectorConversion<T, Tint>(w);
      a = T(1) / T(2) - w_T.dot(*W.center);
      if (!IsInteger(a)) {
        return {};
      }
    } else {
      Tint scal(0);
      for (int j = 0; j < n; j++) {
        scal += w(j) * x0(j);
      }
      a = UniversalScalarConversion<T, Tint>(-scal);
    }
    MyMatrix<T> Fslab = erdahl_slab_function<T, Tint>(a, w);
    if (!erdahl_is_in_space(W, Fslab)) {
      return {};
    }
    MyMatrix<T> F = D.F + Fslab;
    DelaunayPolyhedron<T, Tint> Dinit =
        erdahl_polyhedron_from_function<T, Tint>(F, os);
    if (Dinit.L.rows() != d - 1 || erdahl_perfection_rank(W, Dinit) != r + 1) {
      std::cerr << "ERDAHL: the initial sub-polyhedron does not have the "
                   "expected degeneracy rank or perfection rank\n";
      throw TerminalException{1};
    }
    Dinit.F = erdahl_canonical_function<T, Tint>(
        W, erdahl_lattice_set(Dinit), nullptr, os);
    return Dinit;
  };
  for (int bound = 1; bound < 10; bound++) {
    // The vectors of Z^n of sup norm exactly bound.
    std::vector<int> v(n, -bound);
    while (true) {
      int max_abs = 0;
      for (int i = 0; i < n; i++) {
        max_abs = std::max(max_abs, std::abs(v[i]));
      }
      if (max_abs == bound) {
        MyVector<Tint> w(n);
        for (int i = 0; i < n; i++) {
          w(i) = v[i];
        }
        std::optional<DelaunayPolyhedron<T, Tint>> opt = try_vector(w);
        if (opt) {
          return *opt;
        }
      }
      int pos = 0;
      while (pos < n && v[pos] == bound) {
        v[pos] = -bound;
        pos++;
      }
      if (pos == n) {
        break;
      }
      v[pos]++;
    }
  }
  std::cerr << "ERDAHL: failed to find an initial sub-polyhedron\n";
  throw TerminalException{1};
}

/*
  The sub-polyhedra of a polytope D as Delaunay polyhedra, orbit
  representatives under Aut_W(D).
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<DelaunayPolyhedron<T, Tint>>
erdahl_sub_delaunay_polytope(ErdahlFunctionSpace<T> const &W,
                             DelaunayPolyhedron<T, Tint> const &D,
                             std::string const &FileDualDesc,
                             std::ostream &os) {
  ErdahlSubPolyhedraPolytope<Tint> sub =
      erdahl_sub_polyhedra_polytope<T, Tint, Tgroup>(W, D, FileDualDesc, os);
  std::vector<DelaunayPolyhedron<T, Tint>> l_sub;
  // The pool and the vertices of D, shared by the canonical functions of
  // all its sub-polyhedra.
  ErdahlCanonicalHint<Tint> hint = erdahl_canonical_hint(erdahl_lattice_set(D));
  for (auto &eFace : sub.l_face) {
    l_sub.push_back(erdahl_sub_polyhedron_of_face(W, D, eFace, &hint, os));
  }
  return l_sub;
}

/*
  The orbits under Aut_W(Dp) of the list are split into orbits under
  Stab_W(Dp, D). Since Aff(Dp) is contained in Stab(Dp, D), the finite
  quotients G = Aut(P(Dp)) and H = G_1 are what matters: the orbit of S
  splits into the H-orbits of the g S for g running over the cosets. The
  representatives of both the left and the right cosets are used, the
  duplicates being removed by the triple equivalence.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<DelaunayPolyhedron<T, Tint>>
erdahl_orbit_splitting(ErdahlFunctionSpace<T> const &W,
                       DelaunayPolyhedron<T, Tint> const &D,
                       DelaunayPolyhedron<T, Tint> const &Dp,
                       std::vector<DelaunayPolyhedron<T, Tint>> const &l_orbit,
                       std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  ErdahlFiniteGroup<T, Tint, Tgroup> G =
      erdahl_finite_group<T, Tint, Tgroup>(Dp, {}, os);
  ErdahlFiniteGroup<T, Tint, Tgroup> H =
      erdahl_finite_group<T, Tint, Tgroup>(Dp, {D}, os);
#ifdef DEBUG_ERDAHL_ENUMERATION
  os << "ERDAHL: orbit splitting |G|=" << G.grp.size()
     << " |H|=" << H.grp.size() << " |l_orbit|=" << l_orbit.size() << "\n";
#endif
#ifdef TIMINGS_ERDAHL_ENUMERATION
  MicrosecondTime time;
#endif
  if (G.grp.size() == H.grp.size()) {
    return l_orbit;
  }
  std::set<Telt> set_rep;
  for (auto &eCos : G.grp.left_cosets(H.grp)) {
    set_rep.insert(eCos);
  }
  for (auto &eCos : G.grp.right_cosets(H.grp)) {
    set_rep.insert(eCos);
  }
  std::vector<MyMatrix<Tint>> l_g;
  for (auto &elt : set_rep) {
    MyMatrix<Tint> g = erdahl_lift_permutation<T, Tint, Tgroup>(G.cfg, elt);
    l_g.push_back(erdahl_correct_for_space(W, g, Dp));
  }
  std::vector<DelaunayPolyhedron<T, Tint>> supers{Dp, D};
  std::vector<DelaunayPolyhedron<T, Tint>> l_ret;
#ifdef TIMINGS_ERDAHL_ENUMERATION
  // The quality of the hash invariant: equivalence tests run after a hash
  // match, and those that fail.
  size_t n_test = 0, n_fail = 0;
#endif
  // The data of the super polyhedra and the configuration of each candidate
  // are computed once, for the hash and for the equivalence tests.
  std::vector<ErdahlSuperData<T, Tint>> l_sd = erdahl_list_super_data(supers);
  for (auto &S : l_orbit) {
    std::vector<DelaunayPolyhedron<T, Tint>> l_part;
    std::vector<ErdahlChainConfig<T, Tint>> l_part_cfg;
    std::unordered_map<size_t, std::vector<size_t>> map_hash;
    for (auto &g : l_g) {
      DelaunayPolyhedron<T, Tint> Simg = erdahl_apply_transformation(S, g);
      ErdahlChainConfig<T, Tint> cfg = erdahl_chain_config(Simg, l_sd);
      size_t hash = erdahl_invariant_hash(cfg, os);
      bool is_new = true;
      for (auto &idx : map_hash[hash]) {
        bool test = erdahl_equivalence<T, Tint, Tgroup>(W, l_part[idx],
                                                         l_part_cfg[idx], Simg,
                                                         cfg, l_sd, os)
                        .has_value();
#ifdef TIMINGS_ERDAHL_ENUMERATION
        n_test++;
        if (!test) {
          n_fail++;
        }
#endif
        if (test) {
          is_new = false;
          break;
        }
      }
      if (is_new) {
        map_hash[hash].push_back(l_part.size());
        l_part.push_back(Simg);
        l_part_cfg.push_back(cfg);
      }
    }
    for (auto &Snew : l_part) {
      l_ret.push_back(Snew);
    }
  }
#ifdef DEBUG_ERDAHL_ENUMERATION
  os << "ERDAHL: orbit splitting |l_orbit|=" << l_orbit.size()
     << " -> |l_ret|=" << l_ret.size() << "\n";
#endif
#ifdef TIMINGS_ERDAHL_ENUMERATION
  os << "|ERDAHL: orbit splitting |G|/|H|=" << G.grp.size() / H.grp.size()
     << " |l_orbit|=" << l_orbit.size() << " |l_ret|=" << l_ret.size()
     << " |coset reps|=" << l_g.size() << " n_test=" << n_test
     << " n_fail=" << n_fail << "|=" << time << "\n";
#endif
  return l_ret;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<DelaunayPolyhedron<T, Tint>>
erdahl_sub_delaunay(ErdahlFunctionSpace<T> const &W,
                    DelaunayPolyhedron<T, Tint> const &D,
                    ErdahlBank<T, Tint> &bank, std::ostream &os) {
#ifdef SANITY_CHECK_ERDAHL_ENUMERATION
  erdahl_check_polyhedron(W, D, os);
#endif
#ifdef TIMINGS_ERDAHL_ENUMERATION
  MicrosecondTime time_total;
  MicrosecondTime time;
#endif
  // The bank.
  size_t hash_D = erdahl_invariant_hash<T, Tint>(D, {}, os);
  for (auto &entry : bank.l_entry) {
    if (entry.hash == hash_D) {
      std::optional<MyMatrix<Tint>> opt =
          erdahl_equivalence<T, Tint, Tgroup>(W, entry.D, D, {}, os);
#ifdef TIMINGS_ERDAHL_ENUMERATION
      if (!opt) {
        os << "|ERDAHL: bank, hash match but not equivalent|=0\n";
      }
#endif
      if (opt) {
        std::vector<DelaunayPolyhedron<T, Tint>> l_sub;
        for (auto &S : entry.l_sub) {
          l_sub.push_back(erdahl_apply_transformation(S, *opt));
        }
#ifdef TIMINGS_ERDAHL_ENUMERATION
        os << "|ERDAHL: bank hit, |bank|=" << bank.l_entry.size()
           << "|=" << time << "\n";
#endif
        return l_sub;
      }
    }
  }
#ifdef TIMINGS_ERDAHL_ENUMERATION
  os << "|ERDAHL: bank miss, |bank|=" << bank.l_entry.size() << "|=" << time
     << "\n";
#endif
  std::vector<DelaunayPolyhedron<T, Tint>> l_sub;
  int d = D.L.rows();
  if (d == 0) {
    l_sub = erdahl_sub_delaunay_polytope<T, Tint, Tgroup>(W, D,
                                                          bank.FileDualDesc, os);
  } else {
    erdahl_check_supported_space(W);
    struct Entry {
      DelaunayPolyhedron<T, Tint> Dp;
      bool done;
      ErdahlChainConfig<T, Tint> cfg;
    };
    std::vector<Entry> l_entry;
    // The entries by erdahl_invariant_hash relative to D.
    std::unordered_map<size_t, std::vector<size_t>> map_hash;
    std::vector<DelaunayPolyhedron<T, Tint>> supers{D};
    std::vector<ErdahlSuperData<T, Tint>> l_sd = erdahl_list_super_data(supers);
    auto insert = [&](DelaunayPolyhedron<T, Tint> Dnew) -> void {
#ifdef SANITY_CHECK_ERDAHL_ENUMERATION
      if (!erdahl_is_subset(Dnew, D)) {
        std::cerr << "ERDAHL: the inserted polyhedron is not contained in D\n";
        throw TerminalException{1};
      }
#endif
#ifdef TIMINGS_ERDAHL_ENUMERATION
      MicrosecondTime time_insert;
      size_t n_test = 0;
#endif
      ErdahlChainConfig<T, Tint> cfg = erdahl_chain_config(Dnew, l_sd);
      size_t hash = erdahl_invariant_hash(cfg, os);
      for (auto &idx : map_hash[hash]) {
#ifdef TIMINGS_ERDAHL_ENUMERATION
        n_test++;
#endif
        if (erdahl_equivalence<T, Tint, Tgroup>(W, l_entry[idx].Dp,
                                                 l_entry[idx].cfg, Dnew, cfg,
                                                 l_sd, os)) {
#ifdef TIMINGS_ERDAHL_ENUMERATION
          os << "|ERDAHL: insert, old, n_test=" << n_test
             << "|=" << time_insert << "\n";
#endif
          return;
        }
      }
#ifdef TIMINGS_ERDAHL_ENUMERATION
      os << "|ERDAHL: insert, new, n_test=" << n_test << "|=" << time_insert
         << "\n";
#endif
      erdahl_ensure_function(W, Dnew, os);
      map_hash[hash].push_back(l_entry.size());
      l_entry.push_back({Dnew, false, cfg});
#ifdef DEBUG_ERDAHL_ENUMERATION
      os << "ERDAHL: n=" << erdahl_dimension(D) << " d=" << d
         << " new orbit, |EXT|=" << Dnew.EXT.rows()
         << " d'=" << Dnew.L.rows() << " n_orbit=" << l_entry.size() << "\n";
#endif
    };
    insert(erdahl_initial_sub_polyhedron(W, D, os));
    while (true) {
      int i_sel = -1;
      for (size_t i = 0; i < l_entry.size(); i++) {
        Entry const &entry = l_entry[i];
        if (!entry.done && entry.Dp.L.rows() < d) {
          if (i_sel == -1) {
            i_sel = i;
          } else {
            std::pair<int, int> inv_sel = erdahl_size_key(l_entry[i_sel].Dp);
            if (erdahl_size_key(entry.Dp) < inv_sel) {
              i_sel = i;
            }
          }
        }
      }
      if (i_sel == -1) {
        break;
      }
      l_entry[i_sel].done = true;
      DelaunayPolyhedron<T, Tint> Dp = l_entry[i_sel].Dp;
      std::vector<DelaunayPolyhedron<T, Tint>> l_sub_p =
          erdahl_sub_delaunay<T, Tint, Tgroup>(W, Dp, bank, os);
      std::vector<DelaunayPolyhedron<T, Tint>> l_split =
          erdahl_orbit_splitting<T, Tint, Tgroup>(W, D, Dp, l_sub_p, os);
      for (auto &D1 : l_split) {
        // The function of Dpp is only computed if it is a new orbit.
        DelaunayPolyhedron<T, Tint> Dpp = erdahl_flip(W, D1, Dp, D, false, os);
        insert(Dpp);
      }
    }
    for (auto &entry : l_entry) {
      l_sub.push_back(entry.Dp);
    }
  }
#ifdef TIMINGS_ERDAHL_ENUMERATION
  os << "|ERDAHL: sub_delaunay n=" << erdahl_dimension(D) << " d=" << d
     << " |EXT|=" << D.EXT.rows() << " n_orbit=" << l_sub.size()
     << "|=" << time_total << "\n";
#endif
  bank.l_entry.push_back({D, hash_D, l_sub});
  return l_sub;
}

// The perfect Delaunay polyhedra relative to W, up to the transformations
// preserving W.
template <typename T, typename Tint, typename Tgroup>
std::vector<DelaunayPolyhedron<T, Tint>>
erdahl_enumerate_perfect(ErdahlFunctionSpace<T> const &W,
                         std::string const &FileDualDesc, std::ostream &os) {
  int n = W.n;
  MyMatrix<T> Fzero = ZeroMatrix<T>(n + 1, n + 1);
  DelaunayPolyhedron<T, Tint> D0 =
      erdahl_polyhedron_from_function<T, Tint>(Fzero, os);
  if (erdahl_perfection_rank(W, D0) != 0) {
    std::cerr << "ERDAHL: a nonzero function of W vanishes on Z^n\n";
    throw TerminalException{1};
  }
  ErdahlBank<T, Tint> bank{{}, FileDualDesc};
  return erdahl_sub_delaunay<T, Tint, Tgroup>(W, D0, bank, os);
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_ENUMERATION_H_
// clang-format on
