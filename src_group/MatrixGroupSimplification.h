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

// LocalSimpProduct is a local encapsulation of the product operation
// Supposed to be used only here
template<typename T>
MyMatrix<T> LocalSimpProduct(MyMatrix<T> const& M1, MyMatrix<T> const& M2) {
  return M1 * M2;
}



template<typename Tnorm, typename Ttype, typename Fcomplexity>
std::vector<Ttype> ExhaustiveReductionComplexityKernel(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, [[maybe_unused]] std::ostream& os) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  using Tcomb = std::pair<Ttype, Tnorm>;
  auto f_comp=[](Tcomb const& a, Tcomb const& b) -> bool {
    if (a.second < b.second) {
      return true;
    }
    if (a.second > b.second) {
      return false;
    }
    return a.first < b.first;
  };
  auto get_pair=[&](Ttype const& eM) -> Tcomb {
    Tnorm comp = f_complexity(eM);
    return {eM, comp};
  };
  std::map<Tcomb, size_t, decltype(f_comp)> map(f_comp);
  size_t nonce = 0;
  for (auto & eM: ListM) {
    map[get_pair(eM)] = nonce;
    nonce += 1;
  }
  std::unordered_set<std::pair<size_t,size_t>> set_treated;
  auto f_generate_candidate=[&](Tcomb const& a, Tcomb const& b) -> std::vector<Tcomb> {
    Ttype a_inv = Inverse(a.first);
    Ttype b_inv = Inverse(b.first);
    Ttype prod1 = LocalSimpProduct(a.first, b.first);
    Tcomb pair1 = get_pair(prod1);
    //
    Ttype prod2 = LocalSimpProduct(a_inv, b.first);
    Tcomb pair2 = get_pair(prod2);
    //
    Ttype prod3 = LocalSimpProduct(a.first, b_inv);
    Tcomb pair3 = get_pair(prod3);
    //
    Ttype prod4 = LocalSimpProduct(a_inv, b_inv);
    Tcomb pair4 = get_pair(prod4);
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
      Tnorm a_norm = a_work.second;
      Tnorm b_norm = b_work.second;
      Tnorm cand_norm = cand.second;
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
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val shoulf be present in order to be removed\n";
      throw TerminalException{1};
    }
    map.erase(iter);
  };
  auto get_pos=[&](Tcomb const& val) -> size_t {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val shoulf be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t distance = std::distance(map.begin(), iter);
    return distance;
  };
  Tnorm total_complexity(0);
  for (auto & kv: map) {
    total_complexity += kv.first.second;
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
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
      if (v < map.size() - 1) {
        v += 1;
      } else {
        u += 1;
        v = u + 1;
        if (u == map.size() - 1) {
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
    size_t n_already0 = 0;
    size_t n_already1 = 0;
#endif
    while(true) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP: Addressing u=" << u << " v=" << v << " n_operation=" << n_operation << " |set|=" << map.size() << "\n";
#endif
      auto iter1 = map.begin();
      std::advance(iter1, u);
      Tcomb a = iter1->first;
      size_t nonce_a = iter1->second;
      //
      auto iter2 = map.begin();
      std::advance(iter2, v);
      Tcomb b = iter2->first;
      size_t nonce_b = iter2->second;
      std::pair<size_t, size_t> nonce_pair{nonce_a, nonce_b};
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP:   Complexities, a.second=" << a.second << " b.second=" << b.second << "\n";
#endif
      //
      bool already_treated = false;
      std::pair<size_t, std::vector<Tcomb>> pair{0,{}};
      if (set_treated.find(nonce_pair) != set_treated.end()) {
        already_treated = true;
      } else {
        pair = f_reduce(a, b);
      }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      if (already_treated) {
        n_already1 += 1;
      }
      if (!already_treated) {
        n_already0 += 1;
      }
      os << "SIMP:   n_changes=" << pair.first << " already_treated=" << already_treated << " n_already0=" << n_already0 << " n_already1=" << n_already1 << "\n";
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
          auto iter = map.find(ent);
          if (iter == map.end()) {
            new_elt.push_back(ent);
          } else {
            size_t distance = std::distance(map.begin(), iter);
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
          map[elt] = nonce;
          nonce += 1;
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
            if (u == map.size()) { // This can happen if a and b are removed and there is nothing after.
              break;
            }
            if (u == map.size() - 1) {
              // This can happen if a is removed but either b or something is after and that is it.
              // So, no more operation possible
              break;
            }
            // Left u to the same value as if u remains the same, we are in the next element.
            v = u + 1; // This is a valid position, continuing
          } else {
            // a still exists
            if (v == map.size()) {
              // v is not invalid, going to the next if possible.
              if (u < map.size() - 2) { // Enough space to make something work
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
        if (!already_treated) {
          set_treated.insert(nonce_pair);
        }
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


  while(true) {
    look_for_simplification();
    Tnorm new_complexity(0);
    for (auto & kv: map) {
      new_complexity += kv.first.second;
    }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: total_complexity=" << total_complexity << " new_complexity=" << new_complexity << "\n";
#endif
    if (total_complexity == new_complexity) {
      break;
    }
    total_complexity = new_complexity;
  }
  std::vector<Ttype> new_list_gens;
  for (auto & kv: map) {
    new_list_gens.push_back(kv.first.first);
  }
  return new_list_gens;
}

template<typename Tnorm, typename Ttype, typename Fcomplexity>
std::vector<Ttype> ExhaustiveReductionComplexity(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, [[maybe_unused]] std::ostream& os) {
  std::unordered_set<Ttype> SetMred;
  for (auto & eM : ListM) {
    Ttype eM_inv = Inverse(eM);
    if (eM < eM_inv) {
      SetMred.insert(eM);
    } else {
      SetMred.insert(eM_inv);
    }
  }
  std::vector<Ttype> ListMred;
  for (auto & eM: SetMred) {
    ListMred.push_back(eM);
  }
  return ExhaustiveReductionComplexityKernel<Tnorm,Ttype,Fcomplexity>(ListM, f_complexity, os);
}






// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
