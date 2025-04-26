// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_
#define SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_


template<typename T>
struct ComplexityMeasure {
  T ell1;
  T ellinfinity;
};

template<typename T>
ComplexityMeasure<T> get_complexity_measure(MyMatrix<T> const& M) {
  int n = M.rows();
  T ell1(0);
  T ellinfinity(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = M(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
      if (abs_val > ellinfinity) {
        ellinfinity = abs_val;
      }
    }
  }
  return {ell1, ellinfinity};
}


template<typename T, typename Fcomplexity>
std::vector<MyMatrix<T>> ExhaustiveReductionComplexity(std::vector<MyMatrix<T>> const& ListM, Fcomplexity f_complexity) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  using Tcomb = std::pair<MyMatrix<T>, T>;
  auto f_comp=[](Tcomb const& a, Tcomb const& b) -> bool {
    if (a.second < b.second) {
      return true;
    }
    if (a.second > b.second) {
      return false;
    }
    return a.first < b.first;
  };
  auto get_pair=[&](MyMatrix<T> const& eM) -> Tcomb {
    T comp = f_complexity(eM);
    return {eM, comp};
  };
  std::set<Tcomb, decltype(f_comp)> set(f_comp);
  for (auto & eM: ListM) {
    set.insert(get_pair(eM));
  }
  auto f_generate_candidate=[&](Tcomb const& a, Tcomb const& b) -> std::vector<Tcomb> {
    MyMatrix<T> a_inv = Inverse(a.first);
    MyMatrix<T> b_inv = Inverse(b.first);
    MyMatrix<T> prod1 = a.first * b.first;
    Tcomb pair1 = get_pair(prod1);
    MyMatrix<T> prod2 = a_inv * b.first;
    Tcomb pair2 = get_pair(prod1);
    MyMatrix<T> prod3 = a.first * b_inv;
    Tcomb pair3 = get_pair(prod1);
    MyMatrix<T> prod4 = a_inv * b_inv;
    Tcomb pair4 = get_pair(prod1);
    return {pair1, pair2, pair3, pair4};
  };
  auto f_get_best_candidate=[&](Tcomb const& a, Tcomb const& b) -> Tcomb {
    std::vector<Tcomb> l_comb = f_generate_candidate(a, b);
    Tcomb best_comp = l_comb[0];
    for (size_t i=1; i<l_comb.size(); i++) {
      if (l_comb[i].second < best_comp.second) {
        best_comp = l_comb[i];
      }
    }
    return best_comp;
  };
  auto f_reduce=[&](Tcomb const& a, Tcomb const& b) -> std::pair<size_t, std::vector<Tcomb>> {
    Tcomb a_work = a;
    Tcomb b_work = b;
    size_t n_change = 0;
    while(true) {
      Tcomb cand = f_get_best_candidate(a_work, b_work);
      T a_norm = a_work.second;
      T b_norm = b_work.second;
      T cand_norm = cand.second;
      bool do_something = true;
      if (cand_norm < a_norm && cand_norm < b_norm) {
        if (a_norm < b_norm) {
          b_work = cand;
        } else {
          a_work = cand;
        }
      } else {
        if (cand_norm < b_norm) {
          b_work = cand;
        } else {
          if (cand_norm < a_norm) {
            a_work = cand;
          } else {
            do_something = false;
          }
        }
      }
      if (!do_something) {
        std::vector<Tcomb> npair{a_work, b_work};
        return {n_change, npair};
      } else {
        n_change += 1;
      }
    }
  };
  auto look_for_simplification=[&]() -> void {
    // Iterating over the elements and looking for simplifications.
    //
    // The iterators are unstable, so everytime we change the state,
    // the iterators are invalidated and need to be recomputed.
    size_t u = 0;
    size_t v = 1;
    while(true) {
      auto iter1 = set.begin();
      std::advance(iter1, u);
      Tcomb a = *iter1;
      //
      auto iter2 = set.begin();
      std::advance(iter2, v);
      Tcomb b = *iter2;
      //
      std::pair<size_t, std::vector<Tcomb>> pair = f_reduce(a, b);
      if (pair.first > 0) {
        std::vector<size_t> att;
        std::vector<Tcomb> new_elt;
        size_t min_distance = miss_val;
        bool a_attained = false;
        bool b_attained = false;
        for (auto & ent: pair.second) {
          auto iter = set.find(ent);
          if (iter == set.end()) {
            new_elt.push_back(ent);
          } else {
            size_t distance = std::distance(set.begin(), iter);
            if (distance < min_distance) {
              min_distance = distance;
            }
            if (distance == u) {
              a_attained = true;
            }
            if (distance == v) {
              b_attained = true;
            }
            att.push_back(distance);
          }
        }
        if (!a_attained) {
          auto iter = set.find(a);
          if (iter == set.end()) {
            std::cerr << "SIMP: a shoulf be attained (to be removed)\n";
            throw TerminalException{1};
          }
          set.delete(iter);
        }
        if (!b_attained) {
          auto iter = set.find(b);
          if (iter == set.end()) {
            std::cerr << "SIMP: b shoulf be attained (to be removed)\n";
            throw TerminalException{1};
          }
          set.delete(iter);
        }
        if (a_attained && b_attained) {
          std::cerr << "SIMP: if n_changes > 0 then either a or b is removed\n";
          throw TerminalException{1};
        }
        size_t min_distance_bis = miss_val;
        for (auto & elt: new_elt) {
          set.insert(elt);
          auto iter = set.find(elt);
          size_t distance = std::distance(set.first, iter);
          if (distance < min_distance_bis) {
            min_distance_bis = distance;
          }
        }
        if (min_distance_bis == miss_val) {
          std::cerr << "SIMP: min_distance_bis should be new\n";
          throw TerminalException{1};
        }
        if (!a_attained) {
          u = min_distance_bis;
          v = min_distance_bis + 1;
        }
        if (!b_attained) {
          v = min_distance_bis;
        }
      } else {
        if (v < set.size() - 1) {
          v += 1;
        } else {
          u += 1;
          v = u + 1;
          if (u == set.size() - 1) {
            break;
          }
        }
      }
    }
  };
  look_for_simplification();
  std::vector<MyMatrix<T>> new_list_gens;
  for (auto & elt: set) {
    new_list_gens.push_back(elt.first);
  }
  return new_list_gens;
}







// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
