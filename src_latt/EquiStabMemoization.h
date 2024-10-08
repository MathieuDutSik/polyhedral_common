// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_EQUISTABMEMOIZATION_H_
#define SRC_INDEFINITE_MODELS_EQUISTABMEMOIZATION_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <utility>
#include <vector>
// clang-format on

// This is a simple process for memoizing equivalence and stabilizers
// in order to get better speed.
// We keep track of
// ---Isomorphisms
// ---Non isomorphisms
// ---Stabilizers
// Strategies used:
// ---Connected isomorphism:
//    ---If x is isomorphic to y
//    ---If y is isomorphic to z
//    ---Then an isomorphsism between x and z can be deduced
// ---Connected non-isomorphism:
//    ---If x is isomorphic to y
//    ---If y is not-isomorphic to z
//    ---Then x and z are not isomorphic.
// ---Connected stabilizers
//    ---If x is isomorphic to y
//    ---If the stabilizer of y is known
//    ---Then the stabilizer of x can be computed.
// The above examples are miniature examples. In the code, connected components
// are computed using known isomorphisms.
//
// Assumptions:
// ---The stabilizer is represented by a finite number of equivalence.
// ---The action is (equiv, x1) = x2.
// ---The composition rule (equi1, (equi2, x)) = (equi1 * equi2, x)
// ---In terms of template function
//    ---The composition is represented by the "operator*" function.
//    ---The inverse is represented by the "Inverse" function
//    ---The identity function is represented by the "IdentityObject" template
//    function
//
// Consequences and examples
// ---if (equi, x1) = x2
//    then any matrix g2 satisfying (g2, x2) = x2 give rise
//    to a matrix g1 = equi^{-1} g2 equi
// ---The arithmetic group action is our basic one:
//    (P, A) = P A P^T

#ifdef DEBUG
#define DEBUG_EQUI_STAB_MEMOIZATION
#endif

template <typename T, typename Tint>
void IDENTITY_OBJECT(stc<MyMatrix<T>> const &x, MyMatrix<Tint> &equiv) {
  int n = x.val.rows();
  equiv = IdentityMat<Tint>(n);
}

template <typename Tequiv, typename Tdata>
Tequiv IdentityObject(Tdata const &x) {
  stc<Tdata> x_in{x};
  Tequiv equiv;
  IDENTITY_OBJECT(x_in, equiv);
  return equiv;
}

template <typename Tdata, typename Tequiv> struct DatabaseResultEquiStab {
private:
  // If an entry (x2, eq) is in list_iso[x1], that is (eq, x1) = x2
  std::unordered_map<Tdata, std::vector<std::pair<Tdata, Tequiv>>> list_iso;
  // If (x,y) is in list_non_iso then x and y are not isomorphic
  std::unordered_set<std::pair<Tdata, Tdata>> list_non_iso;
  // If eq is in list_stab[x] then this means that (eq, x) = x
  std::unordered_map<Tdata, std::vector<Tequiv>> list_stab;

public:
  void insert_equi(Tdata const &x1, Tdata const &x2,
                   std::optional<Tequiv> const &result) {
    if (result) {
      Tequiv const &res = *result;
      Tequiv res_inv = Inverse(res);
      //
      std::pair<Tdata, Tequiv> pair2{x2, res};
      list_iso[x1].push_back(pair2);
      //
      std::pair<Tdata, Tequiv> pair1{x1, res_inv};
      list_iso[x2].push_back(pair1);
    } else {
      std::pair<Tdata, Tdata> pair{x1, x2};
      list_non_iso.insert(pair);
    }
  }
  void insert_stab(Tdata const &x, std::vector<Tequiv> const &ListGens) {
    list_stab[x] = ListGens;
  }
  template <typename Fterminate>
  std::vector<std::pair<Tdata, Tequiv>>
  connected_component(Tdata const &x, Fterminate f_terminate) const {
    std::vector<std::pair<Tdata, Tequiv>> l_vertices;
    std::set<Tdata> set_vert;
    auto f_insert = [&](std::pair<Tdata, Tequiv> const &pair) -> void {
      if (set_vert.count(pair.first) == 1) {
        return;
      }
      set_vert.insert(pair.first);
      l_vertices.push_back(pair);
    };
    Tequiv id = IdentityObject<Tequiv, Tdata>(x);
    std::pair<Tdata, Tequiv> start_pair{x, id};
    if (f_terminate(start_pair)) {
      return {};
    }
    f_insert(start_pair);
    if (list_iso.find(x) == list_iso.end()) {
      // It is an isolated vertices, so no need to get into the loop.
      return l_vertices;
    }
    //
    size_t pos_start = 0;
    while (true) {
      size_t pos_end = l_vertices.size();
#ifdef DEBUG_EQUI_STAB_MEMOIZATION
      std::cerr << "ESM: connected_components, pos_start=" << pos_start
                << " pos_end=" << pos_end << "\n";
#endif
      for (size_t pos = pos_start; pos < pos_end; pos++) {
        std::pair<Tdata, Tequiv> pair = l_vertices[pos];
        // pair = (y, eqA) which gets us (eqA, x) = y
        // epair = (z, eqB) in list_iso[y] means that (eqB, y) = z
        // Therefore z = (eqB, (eqA, x)) = (eqB * eqA, x)
        std::vector<std::pair<Tdata, Tequiv>> const &l_pair =
            list_iso.at(pair.first);
        for (auto &epair : l_pair) {
          Tequiv new_equiv = epair.second * pair.second;
          std::pair<Tdata, Tequiv> new_pair{epair.first, new_equiv};
          if (f_terminate(new_pair)) {
            return {};
          }
          f_insert(new_pair);
        }
      }
      pos_start = pos_end;
      if (pos_end == l_vertices.size()) {
        break;
      }
    }
    return l_vertices;
  }
  bool has_noniso_conn_edge(
      std::vector<std::pair<Tdata, Tequiv>> const &eList1,
      std::vector<std::pair<Tdata, Tequiv>> const &eList2) const {
    auto has_noniso_edge = [&](Tdata const &x, Tdata const &y) -> bool {
      std::pair<Tdata, Tdata> pair1{x, y};
      if (list_non_iso.count(pair1) == 1) {
        return true;
      }
      std::pair<Tdata, Tdata> pair2{y, x};
      if (list_non_iso.count(pair2) == 1) {
        return true;
      }
      return false;
    };
    for (auto &pair1 : eList1) {
      for (auto &pair2 : eList2) {
        if (has_noniso_edge(pair1.first, pair2.first)) {
          return true;
        }
      }
    }
    return false;
  }
  std::optional<std::vector<Tequiv>> attempt_stabilizer(Tdata const &x) const {
    if (list_stab.size() == 0) {
      // No stabilizer computed so no point in computing anything.
      return {};
    }
    std::optional<std::vector<Tequiv>> opt;
    auto f_terminate = [&](std::pair<Tdata, Tequiv> const &pair) -> bool {
      auto iter = list_stab.find(pair.first);
      if (iter == list_stab.end()) {
        return false;
      }
      // We have pair = (y, eqA)  with (eqA, x) = y
      // And we have a set of generator eq s.t. (eqB, y) = y
      // therefore we have (eqA^{-1} eqB eqA, x) = x
      Tequiv const &eq = pair.second;
      Tequiv eq_inv = Inverse(eq);
      std::vector<Tequiv> New_list_gens;
      std::vector<Tequiv> const &list_gens = iter->second;
      for (auto &eGen : list_gens) {
        Tequiv NewGen = eq_inv * eGen * eq;
        New_list_gens.push_back(NewGen);
      }
      opt = New_list_gens;
      return true;
    };
    (void)connected_component<decltype(f_terminate)>(x, f_terminate);
    return opt;
  }
  std::optional<std::optional<Tequiv>> attempt_equiv(Tdata const &x,
                                                     Tdata const &y) const {
    std::optional<std::optional<Tequiv>> opt;
    auto f_terminate_x = [&](std::pair<Tdata, Tequiv> const &pair) -> bool {
      if (pair.first == y) {
        std::optional<Tequiv> new_opt = pair.second;
        opt = new_opt;
        return true;
      }
      return false;
    };
    std::vector<std::pair<Tdata, Tequiv>> l_pair_x =
        connected_component<decltype(f_terminate_x)>(x, f_terminate_x);
    if (opt) {
#ifdef DEBUG_EQUI_STAB_MEMOIZATION
      if (l_pair_x.size() != 0) {
        std::cerr << "l_pair_x should be of zero length\n";
        throw TerminalException{1};
      }
#endif
      return opt;
    }
    auto f_terminate_y =
        [&]([[maybe_unused]] std::pair<Tdata, Tequiv> const &pair) -> bool {
      return false;
    };
    std::vector<std::pair<Tdata, Tequiv>> l_pair_y =
        connected_component<decltype(f_terminate_y)>(y, f_terminate_y);
    bool test = has_noniso_conn_edge(l_pair_x, l_pair_y);
    if (test) {
      std::optional<Tequiv> ret_opt = {};
      return ret_opt;
    } else {
      return {};
    }
  }
};

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_EQUISTABMEMOIZATION_H_
// clang-format on
