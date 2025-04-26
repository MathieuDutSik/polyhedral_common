// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_
#define SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_

#include "MAT_Matrix.h"


#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_SIMPLIFICATION
#endif


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
std::vector<MyMatrix<T>> ExhaustiveReductionComplexityKernel(std::vector<MyMatrix<T>> const& ListM, Fcomplexity f_complexity, [[maybe_unused]] std::ostream& os) {
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
    Tcomb pair2 = get_pair(prod2);
    MyMatrix<T> prod3 = a.first * b_inv;
    Tcomb pair3 = get_pair(prod3);
    MyMatrix<T> prod4 = a_inv * b_inv;
    Tcomb pair4 = get_pair(prod4);
    return {pair1, pair2, pair3, pair4};
  };
  auto f_get_best_candidate=[&](Tcomb const& a, Tcomb const& b) -> Tcomb {
    std::vector<Tcomb> l_comb = f_generate_candidate(a, b);
    Tcomb best_comp = l_comb[0];
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    size_t i_choice = 0;
#endif
    for (size_t i=1; i<l_comb.size(); i++) {
      if (l_comb[i].second < best_comp.second) {
        best_comp = l_comb[i];
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        i_choice = i;
#endif
      }
    }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: f_get_best_candidate, choosing i_choice=" << i_choice << "\n";
#endif
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
  auto erase_entry=[&](Tcomb const& val) -> void {
    auto iter = set.find(val);
    if (iter == set.end()) {
      std::cerr << "SIMP: val shoulf be present in order to be removed\n";
      throw TerminalException{1};
    }
    set.erase(iter);
  };
  auto get_pos=[&](Tcomb const& val) -> size_t {
    auto iter = set.find(val);
    if (iter == set.end()) {
      std::cerr << "SIMP: val shoulf be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t distance = std::distance(set.begin(), iter);
    return distance;
  };
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  T total_complexity(0);
  for (auto & pair: set) {
    total_complexity += pair.second;
  }
  os << "SIMP: total_complexity=" << total_complexity << "\n";
#endif
  auto look_for_simplification=[&]() -> void {
    // Iterating over the elements and looking for simplifications.
    //
    // The iterators are unstable, so everytime we change the state,
    // the iterators are invalidated and need to be recomputed.
    //
    // As the iteration is being run, the number of elements may decrease.
    // All of this forces a single loop and to handle all the special
    // cases separately.
    size_t u = 0;
    size_t v = 1;
    auto increment_uv=[&]() -> bool {
      if (v < set.size() - 1) {
        v += 1;
      } else {
        u += 1;
        v = u + 1;
        if (u == set.size() - 1) {
          return false;
        }
      }
      return true;
    };
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
    size_t pos = 0;
    os << "SIMP: starting with the following matrices\n";
    for (auto & pair: set) {
      os << "SIMP: pos=" << pos << " comp=" << pair.second << " eM=\n";
      WriteMatrix(os, pair.first);
      pos += 1;
    }
#endif
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    size_t n_operation = 0;
#endif
    while(true) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP: Addressing u=" << u << " v=" << v << " n_operation=" << n_operation << " |set|=" << set.size() << "\n";
#endif
      auto iter1 = set.begin();
      std::advance(iter1, u);
      Tcomb a = *iter1;
      //
      auto iter2 = set.begin();
      std::advance(iter2, v);
      Tcomb b = *iter2;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP:   Complexities, a.second=" << a.second << " b.second=" << b.second << "\n";
#endif
      //
      std::pair<size_t, std::vector<Tcomb>> pair = f_reduce(a, b);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP:   n_changes=" << pair.first << "\n";
#endif
      if (pair.first > 0) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        n_operation += 1;
#endif
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
        os << "SIMP:   Working with a=\n";
        WriteMatrix(os, a.first);
        os << "SIMP:   and b=\n";
        WriteMatrix(os, b.first);
#endif
        std::vector<size_t> att;
        std::vector<Tcomb> new_elt;
        size_t min_distance = miss_val;
        bool a_attained = false;
        bool b_attained = false;
        for (auto & ent: pair.second) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  ent.comp=" << ent.second << " elt.eM=\n";
          WriteMatrix(os, ent.first);
#endif
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
            if (distance != u && distance != v) {
              att.push_back(distance);
            }
          }
        }
        if (!a_attained) {
          erase_entry(a);
        }
        if (!b_attained) {
          erase_entry(b);
        }
        if (a_attained && b_attained) {
          std::cerr << "SIMP: if n_changes > 0 then either a or b is removed\n";
          throw TerminalException{1};
        }
        size_t min_distance_bis = miss_val;
        for (auto & elt: new_elt) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  elt.comp=" << elt.second << " elt.eM=\n";
          WriteMatrix(os, elt.first);
#endif
          set.insert(elt);
          size_t distance = get_pos(elt);
          if (distance < min_distance_bis) {
            min_distance_bis = distance;
          }
        }
        if (min_distance_bis == miss_val) {
          // This scenario occurs if the new found generators are already present
          // Two scenarios
          if (!a_attained) {
            // a was removed,
            if (u == set.size()) { // This can happen if a and b are removed and there is nothing after.
              break;
            }
            if (u == set.size() - 1) {
              // This can happen if a is removed but either b or something is after and that is it.
              // So, no more operation possible
              break;
            }
            // Left u to the same value as if u remains the same, we are in the next element.
            v = u + 1; // This is a valid position, continuing
          } else {
            // a still exists
            if (v == set.size()) {
              // v is not invalid, going to the next if possible.
              if (u < set.size() - 2) { // Enough space to make something work
                u += 1;
                v = u + 1;
              } else {
                break;
              }
            } else {
              // No change in value of v, but since b is dropped, we go to the next one.
            }
          }
        } else {
          // We have a new generator, adjusting accordingly.
          if (!a_attained) {
            u = min_distance_bis;
            v = min_distance_bis + 1;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
            os << "SIMP:   setup A\n";
#endif
          } else {
            size_t pos_a = get_pos(a);
            if (pos_a < min_distance_bis) {
              v = min_distance_bis; // Setting up to where we are and incrementing.
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
              os << "SIMP:   setup B(1)\n";
#endif
              bool test = increment_uv();
              if (!test) {
                break;
              }
            } else {
              u = min_distance_bis;
              v = min_distance_bis + 1;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
              os << "SIMP:   setup B(2)\n";
#endif
            }
          }
        }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        os << "SIMP:   a/b_attained=" << a_attained << "/" << b_attained << " min_distance_bis=" << min_distance_bis << "\n";
#endif
      } else {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        os << "SIMP:   no change, incrementing u / v\n";
#endif
        bool test = increment_uv();
        if (!test) {
          break;
        }
      }
    }
  };
  look_for_simplification();
  std::vector<MyMatrix<T>> new_list_gens;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  T return_complexity(0);
#endif
  for (auto & elt: set) {
    new_list_gens.push_back(elt.first);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    return_complexity += elt.second;
#endif
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  os << "SIMP: total_complexity=" << total_complexity << " return_complexity=" << return_complexity << "\n";
#endif
  return new_list_gens;
}

template<typename T, typename Fcomplexity>
std::vector<MyMatrix<T>> ExhaustiveReductionComplexity(std::vector<MyMatrix<T>> const& ListM, Fcomplexity f_complexity, [[maybe_unused]] std::ostream& os) {
  std::unordered_set<MyMatrix<T>> SetMred;
  for (auto & eM : ListM) {
    MyMatrix<T> eM_inv = Inverse(eM);
    if (eM < eM_inv) {
      SetMred.insert(eM);
    } else {
      SetMred.insert(eM_inv);
    }
  }
  std::vector<MyMatrix<T>> ListMred;
  for (auto & eM: SetMred) {
    ListMred.push_back(eM);
  }
  return ExhaustiveReductionComplexityKernel(ListM, f_complexity, os);
}






// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
