// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_
#define SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_

// clang-format off
#include "MatrixGroupSimplification.h"
#include <stdexcept>
// clang-format on

// Online/incremental version of ExhaustiveReductionComplexityKernelInner_V2
// Allows generators to be inserted one by one rather than providing a complete list
template<typename Ttype, typename Tnorm, typename Fcomplexity, typename Fproduct, typename Fcheck>
class OnlineExhaustiveReductionComplexityKernel {
private:
  std::map<TcombPair<Ttype,Tnorm>, BlockInterval> map;
  std::vector<TcombPair<Ttype,Tnorm>> vect;
  Fcomplexity f_complexity;
  Fproduct f_product;
  Fcheck f_check;
  std::ostream& os;

  // Internal helper functions from the original algorithm
  void delete_entry(TcombPair<Ttype,Tnorm> const& val) {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val should be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t pos = std::distance(map.begin(), iter);
    map.erase(iter);
    vect.erase(vect.begin() + pos);

    for (auto & kv : map) {
      BlockInterval & blk_int = kv.second;
      blk_int.remove_entry_and_shift(pos, os);
    }
  }

  void insert_entry(TcombPair<Ttype,Tnorm> const& val) {
    BlockInterval blk_int;
    auto result = map.insert({val, blk_int});
    if (!result.second) {
      std::cerr << "SIMP: Entry is already present. Unexpected\n";
      throw TerminalException{1};
    }
    auto iter = result.first;
    size_t pos = std::distance(map.begin(), iter);
    vect.insert(vect.begin() + pos, val);

    size_t n_entry = vect.size();
    size_t idx = 0;
    for (auto & kv: map) {
      BlockInterval & blk_int = kv.second;
      if (idx < pos) {
        blk_int.insert_entry_and_shift(pos, os);
      } else {
        if (idx == pos) {
          blk_int.insert_interval(pos + 1, n_entry);
        } else {
          blk_int.noinsert_and_shift(pos, os);
        }
      }
      idx += 1;
    }
  }

  struct FoundImprov {
    std::vector<TcombPair<Ttype,Tnorm>> list_delete;
    std::vector<TcombPair<Ttype,Tnorm>> list_insert;
  };

  std::optional<FoundImprov> search_for_improvement() {
    for (auto & kv: map) {
      TcombPair<Ttype,Tnorm> const& x1 = kv.first;
      BlockInterval & blk_int = kv.second;
      while(true) {
        std::optional<size_t> opt = blk_int.get_first();
        if (opt) {
          size_t idx2 = *opt;
          TcombPair<Ttype,Tnorm> const& x2 = vect[idx2];
          GenResult<Ttype,Tnorm> gen = f_reduce<Ttype,Tnorm,Fcomplexity,Fproduct>(x1, x2, f_complexity, f_product);
          if (gen.do_something) {
            bool x1_attained = false;
            bool x2_attained = false;
            std::vector<TcombPair<Ttype,Tnorm>> list_delete;
            std::vector<TcombPair<Ttype,Tnorm>> list_insert;
            for (auto & ent : gen.l_ent) {
              bool is_x1 = ent == x1;
              bool is_x2 = ent == x2;
              if (is_x1) {
                x1_attained = true;
              }
              if (is_x2) {
                x2_attained = true;
              }
              if (!is_x1 && !is_x2) {
                if (map.find(ent) == map.end()) {
                  list_insert.push_back(ent);
                }
              }
            }
            if (!x1_attained) {
              list_delete.push_back(x1);
            }
            if (!x2_attained) {
              list_delete.push_back(x2);
            }
            return FoundImprov{list_delete, list_insert};
          }
        } else {
          break;
        }
      }
    }
    return {};
  }

  bool run_reduction_pass() {
    while(true) {
      std::optional<FoundImprov> opt = search_for_improvement();
      if (opt) {
        FoundImprov found_improv = *opt;
        // We need to do the insertions before the deletes
        // since if there is a failure in the reduction, we
        // want to avoid losing a generator when we switch to
        // another numerics.
        for (auto & val : found_improv.list_insert) {
          if (!f_check(val.pair.first)) {
            return false;
          }
          insert_entry(val);
        }
        for (auto & val : found_improv.list_delete) {
          delete_entry(val);
        }
      } else {
        break;
      }
    }
    return true;
  }

public:
  OnlineExhaustiveReductionComplexityKernel(Fcomplexity _f_complexity, Fproduct _f_product, Fcheck _f_check, std::ostream& _os)
    : f_complexity(_f_complexity), f_product(_f_product), f_check(_f_check), os(_os) {}

  // Insert a new generator and run reduction
  bool insert_generator(TcombPair<Ttype,Tnorm> const& generator) {
    if (!f_check(generator.pair.first)) {
      return false;
    }

    if (map.count(generator) == 1) {
      return true; // Already present
    }

    insert_entry(generator);
    return run_reduction_pass();
  }

  // Get current simplified set
  std::vector<TcombPair<Ttype,Tnorm>> get_current_set() const {
    return vect;
  }

  // Get total complexity
  Tnorm get_total_complexity() const {
    Tnorm total_complexity(0);
    for (auto & kv: map) {
      total_complexity += kv.first.norm;
    }
    return total_complexity;
  }

  // Get number of generators
  size_t size() const {
    return vect.size();
  }

  // Check if empty
  bool empty() const {
    return vect.empty();
  }

  // Clear all generators
  void clear() {
    map.clear();
    vect.clear();
  }

  // Get current generators as vector (alias for compatibility)
  std::vector<TcombPair<Ttype,Tnorm>> get_generators() const {
    return get_current_set();
  }
};

// clang-format off
#endif  // SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_
// clang-format on
