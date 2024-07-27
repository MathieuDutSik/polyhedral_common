// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_EQUISTABMEMOIZATION_H_
#define SRC_INDEFINITE_MODELS_EQUISTABMEMOIZATION_H_


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
// Of course, this is combined with connected components and other iterative algorithms
//
// Assumptions:
// ---The stabilizer is represented by a finite number of equivalence.
// ---The action is (equiv, x1) = x2.
// ---In terms of template function
//    ---The composition is represented by the "operator*" function.
//    ---The inverse is represented by the "Inverse" function
//    ---The identity function is represented by the "IdentityObject" template function

#ifdef DEBUG
#define DEBUG_EQUI_STAB_MEMOIZATION
#endif


template<typename Tdata, template Tequiv>
Tequiv IdentityObject(Tdata const& x) {




// Used for the memoization process for isomorphism
template<typename Tdata, typename Tequiv>
struct ResultEquivalence {
  Tdata x1;
  Tdata x2;
  std::optional<Tequiv> result;
};

template<typename Tdata, typename Tequiv>
struct ResultStabilizer {
  Tdata x;
  std::vector<Tequiv> ListGens;
};

template<typename Tdata, typename Tequiv>
struct DatabaseResultEquiStab {
private:
  std::unordered_map<Tdata, std::vector<std::pair<Tdata, Tequiv>>> list_iso;
  std::unordered_set<std::pair<Tdata,Tdata>> list_non_iso;
  std::unordered_map<Tdata, std::vector<Tequiv>> list_stab;

public:
  void insert_equi(ResultEquivalence<Tdata,Tequiv> const& result_equiv) {
    if (result_equiv.result) {
      Tequiv const& res = *result_equiv.result;
      std::pair<Tdata, Tequiv> pair1{result_equiv.x1, res};
      list_iso[result_equiv.x2].push_back(pair1);
      //
      Tequiv res_inv = Inverse(res);
      std::pair<Tdata, Tequiv> pair2{result_equiv.x2, res_inv};
      list_iso[result_equiv.x1].push_back(pair2);
    } else {
      std::pair<Tdata,Tdata> pair{result_equiv.x1, result_equiv.x2};
      list_non_iso.insert(pair);
    }
  }
  void insert_stab(ResultStabilizer<Tdata, Tequiv> const& result_stab) {
    list_stab[result_stab.x] = result_stab.ListGens;
  }
  template<typename Tterminate>
  std::vector<std::pair<Tdata,Tequiv>> connected_component(MyMatrix<T> const& x, Fterminate f_terminate) const {
    std::vector<std::pair<Tdata,Tequiv>> l_vertices;
    std::set<Tdata> set_vert;
    auto f_insert=[&](std::pair<Tdata,Tequiv> const& pair) -> void {
      if (set_vert.count(pair.first) == 1) {
        return;
      }
      set_vert.insert(pair.first);
      l_vertices.push_back(pair);
    };
    Tequiv id = IdentityObject<Tequiv,T>(x);
    std::pair<Tdata,Tequiv> pair{x, id};
    f_insert(pair);
    //
    size_t pos_start = 0;
    while(true) {
      size_t pos_end = l_vertices.size();
      for (size_t pos=pos_start; pos<pos_end; pos++) {
        std::pair<Tdata,Tequiv> pair = l_vertices[pos];
        for (auto & epair : list_iso[pair.first]) {
          Tequiv new_equiv = epair.second * pair.second;
          std::pair<Tdata,Tequiv> new_pair{epair.first, new_equiv};
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
  bool has_noniso_conn_edge(std::vector<std::pair<Tdata,Tequiv>> const& eList1, std::vector<std::pair<Tdata,Tequiv>> const& eList2) const {
    auto has_noniso_edge=[&](Tdata const&x, Tdata const& y) -> bool {
      std::pair<Tdata,Tdata> pair1{x,y};
      if (list_non_iso.count(pair1) == 1) {
        return true;
      }
      std::pair<Tdata,Tdata> pair2{y,x};
      if (list_non_iso.count(pair2) == 1) {
        return true;
      }
      return false;
    };

    for (auto &pair1 : eList1) {
      for (auto & pair2 : eList2) {
        if (has_noniso_edge(pair1.first, pair2.first)) {
          return true;
        }
      }
    }
    return false;
  }
  std::optional<std::vector<Tequiv>> attempt_stabilizer(Tdata const& x) const {
    if (list_stab.size() == 0) {
      return {};
    }
    std::optional<std::vector<Tequiv>> opt;
    auto f_terminate=[&](std::pair<Tdata,Tequiv> const& pair) -> bool {
      auto iter = list_stab.find(pair.first);
      if (iter == list_stab.end()) {
        return false;
      }
      Tequiv const& eq = pair.second
      Tequiv eq_inv = Inverse(eq);
      std::vector<Tequiv> l_gens;
      for (auto& eGen : *iter) {
        Tequiv NewGen = eq * eGen * eq_inv;
        l_gens.push_back(NewGen);
      }
      opt = l_gens;
      return true;
    };
    (void)connected_component<decltype(f_terminate)>(x);
    return opt;
  }
  std::optional<std::optional<Tequiv>> attempt_equiv(Tdata const& x, Tdata const& y) const {
    std::optional<std::optional<Tequiv>> opt;
    auto f_terminate_x=[&](std::pair<Tdata,Tequiv> const& pair) -> bool {
      if (pair.first == y) {
        std::optional<Tequiv> new_opt = pair.second;
        opt = new_opt;
        return true;
      }
      return false;
    };
    std::vector<std::pair<Tdata,Tequiv>> l_pair_x = connected_component<decltype(f_terminate_x)>(x);
    if (opt) {
#ifdef DEBUG_EQUI_STAB_MEMOIZATION
      if (l_pair_x.size() != 0) {
        std::cerr << "l_pair_x should be of zero length\n";
        throw TerminalException{1};
      }
#endif
      return opt;
    }
    auto f_terminate_y=[&](std::pair<Tdata,Tequiv> const& pair) -> bool {
      return false;
    };
    std::vector<std::pair<Tdata,Tequiv>> l_pair_y = connected_component<decltype(f_terminate_y)>(y);
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
