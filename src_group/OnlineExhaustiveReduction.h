// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_
#define SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_

// clang-format off
#include "MatrixGroupSimplification.h"
#include <stdexcept>
#include <functional>
#include <memory>
#include <limits>
// clang-format on

#ifdef DEBUG
#define DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
#endif

template<typename U>
using PairMatrix = std::pair<MyMatrix<U>,MyMatrix<U>>;


template <typename T2, typename T1>
PairMatrix<T2> UniversalPairMatrixConversion(PairMatrix<T1> const &pair) {
  MyMatrix<T2> m1 = UniversalMatrixConversion<T2,T1>(pair.first);
  MyMatrix<T2> m2 = UniversalMatrixConversion<T2,T1>(pair.second);
  return {std::move(m1), std::move(m2)};
}

template <typename T2, typename T1>
std::optional<PairMatrix<T2>> UniversalPairMatrixConversionCheck(PairMatrix<T1> const &pair) {
  std::optional<MyMatrix<T2>> opt1 = UniversalMatrixConversionCheck<T2,T1>(pair.first);
  if (!opt1) {
    return {};
  }
  std::optional<MyMatrix<T2>> opt2 = UniversalMatrixConversionCheck<T2,T1>(pair.second);
  if (!opt2) {
    return {};
  }
  PairMatrix<T2> p{*opt1, *opt2};
  return p;
}

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
          if (!f_check(val.pair.first) || !f_check(val.pair.second)) {
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
  bool insert_generator(std::pair<Ttype, Ttype> const& pair) {
    if (!f_check(pair.first) || !f_check(pair.second)) {
      return false;
    }
    TcombPair<Ttype,Tnorm> comb = generate_comb_pair<Ttype,Tnorm,Fcomplexity>(pair, f_complexity);
    if (map.count(comb) == 1) {
      return true; // Already present
    }

    insert_entry(comb);
    return run_reduction_pass();
  }

  // Get current simplified set
  std::vector<TcombPair<Ttype,Tnorm>> get_current_set() {
    std::vector<TcombPair<Ttype,Tnorm>> ret_vect = vect;
    vect.clear();
    map.clear();
    return ret_vect;
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
    return vect;
  }
};


// Hierarchical online reduction system that automatically switches between
// int16_t -> int32_t -> int64_t -> T numerics when overflow occurs
template<typename T>
class OnlineHierarchicalMatrixReduction {
private:
  int n_;  // Matrix dimension
  std::ostream& os;

  // Type aliases for matrix-matrix pairs
  using Ttype = std::pair<MyMatrix<T>, MyMatrix<T>>;
  using Tnorm = T;

  template<typename Tin>
  using Tinput = TcombPair<MyMatrix<Tin>,Tin>;

  // Upper bounds for each numeric type
  int16_t max_val_16;
  int32_t max_val_32;
  int64_t max_val_64;

  // Type aliases for function types
  using f_complexity_16_t = std::function<int16_t(MyMatrix<int16_t> const&)>;
  using f_product_16_t = std::function<MyMatrix<int16_t>(MyMatrix<int16_t> const&, MyMatrix<int16_t> const&)>;
  using f_check_16_t = std::function<bool(MyMatrix<int16_t> const&)>;

  using f_complexity_32_t = std::function<int32_t(MyMatrix<int32_t> const&)>;
  using f_product_32_t = std::function<MyMatrix<int32_t>(MyMatrix<int32_t> const&, MyMatrix<int32_t> const&)>;
  using f_check_32_t = std::function<bool(MyMatrix<int32_t> const&)>;

  using f_complexity_64_t = std::function<int64_t(MyMatrix<int64_t> const&)>;
  using f_product_64_t = std::function<MyMatrix<int64_t>(MyMatrix<int64_t> const&, MyMatrix<int64_t> const&)>;
  using f_check_64_t = std::function<bool(MyMatrix<int64_t> const&)>;

  using f_complexity_T_t = std::function<Tnorm(MyMatrix<T> const&)>;
  using f_product_T_t = std::function<MyMatrix<T>(MyMatrix<T> const&, MyMatrix<T> const&)>;
  using f_check_T_t = std::function<bool(MyMatrix<T> const&)>;

  // Online kernels for each numeric type
  std::unique_ptr<OnlineExhaustiveReductionComplexityKernel<
    MyMatrix<int16_t>, int16_t, f_complexity_16_t, f_product_16_t, f_check_16_t
  >> kernel_16;

  std::unique_ptr<OnlineExhaustiveReductionComplexityKernel<
    MyMatrix<int32_t>, int32_t, f_complexity_32_t, f_product_32_t, f_check_32_t
  >> kernel_32;

  std::unique_ptr<OnlineExhaustiveReductionComplexityKernel<
    MyMatrix<int64_t>, int64_t, f_complexity_64_t, f_product_64_t, f_check_64_t
  >> kernel_64;

  std::unique_ptr<OnlineExhaustiveReductionComplexityKernel<
    MyMatrix<T>, Tnorm, f_complexity_T_t, f_product_T_t, f_check_T_t
  >> kernel_T;

  // Current active kernel (0=int16_t, 1=int32_t, 2=int64_t, 3=T)
  int current_level;

  template<typename Tin>
  bool insert_level_32(std::vector<Tinput<Tin>> const& l_gen) {
    for (auto const& comb_pair : l_gen) {
      PairMatrix<int32_t> pair_32 = UniversalPairMatrixConversion<int32_t,Tin>(comb_pair.pair);
      if (!kernel_32->insert_generator(pair_32)) {
        kernel_32->clear();
        return false;
      }
    }
    return true;
  }

  template<typename Tin>
  bool insert_level_64(std::vector<Tinput<Tin>> const& l_gen) {
    for (auto const& comb_pair : l_gen) {
      PairMatrix<int64_t> pair_64 = UniversalPairMatrixConversion<int64_t,Tin>(comb_pair.pair);
      if (!kernel_64->insert_generator(pair_64)) {
        kernel_64->clear();
        return false;
      }
    }
    return true;
  }

  template<typename Tin>
  void insert_level_T(std::vector<Tinput<Tin>> const& l_gen) {
    for (auto const& comb_pair : l_gen) {
      Ttype pair_T = UniversalPairMatrixConversion<T,Tin>(comb_pair.pair);
      if (!kernel_T->insert_generator(pair_T)) {
        std::cerr << "ONL: This should never happen in infinite precision\n";
        throw TerminalException{1};
      }
    }
  }

  // Migrate from current level to next level
  void migrate_to_next_level() {

#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
    os << "SIMP: Migrating from level " << current_level << " to level " << (current_level + 1) << "\n";
#endif

    if (current_level == 0) {
      // Migrate from int16_t to int32_t or int64_t or T
      std::vector<Tinput<int16_t>> l_gen = kernel_16->get_current_set();
      bool result1 = insert_level_32<int16_t>(l_gen);
      if (result1) {
        current_level = 1;
        return;
      }
      bool result2 = insert_level_64<int16_t>(l_gen);
      if (result2) {
        current_level = 2;
        return;
      }
      insert_level_T<int16_t>(l_gen);
      current_level = 3;
    } else if (current_level == 1) {
      // Migrate from int32_t to int64_t or T
      std::vector<Tinput<int32_t>> l_gen = kernel_32->get_current_set();
      bool result1 = insert_level_64<int32_t>(l_gen);
      if (result1) {
        current_level = 2;
        return;
      }
      insert_level_T<int32_t>(l_gen);
      current_level = 3;
    } else if (current_level == 2) {
      // Migrate from int64_t to T
      std::vector<Tinput<int64_t>> l_gen = kernel_64->get_current_set();
      insert_level_T<int64_t>(l_gen);
      current_level = 3;
    }
  }

public:
  OnlineHierarchicalMatrixReduction(int n, std::ostream& _os) : n_(n), os(_os), current_level(0) {
    // Calculate upper bounds for each numeric type
    double max_poss_16 = static_cast<double>(std::numeric_limits<int16_t>::max());
    double max_poss_32 = static_cast<double>(std::numeric_limits<int32_t>::max());
    double max_poss_64 = static_cast<double>(std::numeric_limits<int64_t>::max());

    max_val_16 = static_cast<int16_t>(sqrt(max_poss_16 / (10 * n)));
    max_val_32 = static_cast<int32_t>(sqrt(max_poss_32 / (10 * n)));
    max_val_64 = static_cast<int64_t>(sqrt(max_poss_64 / (10 * n)));

#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
    os << "SIMP: Upper bounds: int16_t=" << max_val_16
       << ", int32_t=" << max_val_32
       << ", int64_t=" << max_val_64 << "\n";
#endif

    // Create function objects
    f_complexity_16_t f_complexity_16 = [](MyMatrix<int16_t> const& M) -> int16_t {
      return get_ell1_complexity_measure(M);
    };
    f_product_16_t f_product_16 = [](MyMatrix<int16_t> const& A, MyMatrix<int16_t> const& B) -> MyMatrix<int16_t> {
      return A * B;
    };
    f_check_16_t f_check_16 = [this](MyMatrix<int16_t> const& M) -> bool {
      return check_matrix_coefficients(M, max_val_16);
    };

    f_complexity_32_t f_complexity_32 = [](MyMatrix<int32_t> const& M) -> int32_t {
      return get_ell1_complexity_measure(M);
    };
    f_product_32_t f_product_32 = [](MyMatrix<int32_t> const& A, MyMatrix<int32_t> const& B) -> MyMatrix<int32_t> {
      return A * B;
    };
    f_check_32_t f_check_32 = [this](MyMatrix<int32_t> const& M) -> bool {
      return check_matrix_coefficients(M, max_val_32);
    };

    f_complexity_64_t f_complexity_64 = [](MyMatrix<int64_t> const& M) -> int64_t {
      return get_ell1_complexity_measure(M);
    };
    f_product_64_t f_product_64 = [](MyMatrix<int64_t> const& A, MyMatrix<int64_t> const& B) -> MyMatrix<int64_t> {
      return A * B;
    };
    f_check_64_t f_check_64 = [this](MyMatrix<int64_t> const& M) -> bool {
      return check_matrix_coefficients(M, max_val_64);
    };

    f_complexity_T_t f_complexity_T = [](MyMatrix<T> const& M) -> Tnorm {
      return get_ell1_complexity_measure(M);
    };
    f_product_T_t f_product_T = [](MyMatrix<T> const& A, MyMatrix<T> const& B) -> MyMatrix<T> {
      return A * B;
    };
    f_check_T_t f_check_T = []([[maybe_unused]] MyMatrix<T> const& M) -> bool {
      return true; // No bounds checking for arbitrary precision
    };

    // Initialize all kernels
    kernel_16 = std::make_unique<typename decltype(kernel_16)::element_type>(
      f_complexity_16, f_product_16, f_check_16, os);
    kernel_32 = std::make_unique<typename decltype(kernel_32)::element_type>(
      f_complexity_32, f_product_32, f_check_32, os);
    kernel_64 = std::make_unique<typename decltype(kernel_64)::element_type>(
      f_complexity_64, f_product_64, f_check_64, os);
    kernel_T = std::make_unique<typename decltype(kernel_T)::element_type>(
      f_complexity_T, f_product_T, f_check_T, os);
  }

  // Insert a single matrix (compute inverse in T first)
  void insert_generator(MyMatrix<T> const& generator) {
    // Compute inverse in T arithmetic first
    MyMatrix<T> generator_inv = Inverse(generator);
    return insert_generator_pair({generator, generator_inv});
  }

  bool simple_insert_generator_pair(Ttype const& generator_pair) {
    if (current_level == 0) {
      // Try int16_t
      std::optional<PairMatrix<int16_t>> opt = UniversalPairMatrixConversionCheck<int16_t,T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_16->insert_generator(*opt);
    } else if (current_level == 1) {
      // Try int32_t
      std::optional<PairMatrix<int32_t>> opt = UniversalPairMatrixConversionCheck<int32_t,T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_32->insert_generator(*opt);
    } else if (current_level == 2) {
      // Try int64_t
      std::optional<PairMatrix<int64_t>> opt = UniversalPairMatrixConversionCheck<int64_t,T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_64->insert_generator(*opt);
    } else {
      // Use T. infinite precision, cannot fail.
      return kernel_T->insert_generator(generator_pair);
    }
  }


  // Insert a matrix pair
  void insert_generator_pair(Ttype const& generator_pair) {
    while (true) {
      bool test = simple_insert_generator_pair(generator_pair);
      if (test) {
        return;
      }
      migrate_to_next_level();
    }
  }

  // Extract the final reduced generators (always in type T)
  std::vector<MyMatrix<T>> get_reduced_generators() {
    std::vector<MyMatrix<T>> result;

    if (current_level == 0) {
      auto current_set = kernel_16->get_current_set();
      for (auto const& comb_pair : current_set) {
        MyMatrix<T> mat_T = UniversalMatrixConversion<T,int16_t>(comb_pair.pair.first);
        result.push_back(mat_T);
      }
    } else if (current_level == 1) {
      auto current_set = kernel_32->get_current_set();
      for (auto const& comb_pair : current_set) {
        MyMatrix<T> mat_T = UniversalMatrixConversion<T,int32_t>(comb_pair.pair.first);
        result.push_back(mat_T);
      }
    } else if (current_level == 2) {
      auto current_set = kernel_64->get_current_set();
      for (auto const& comb_pair : current_set) {
        MyMatrix<T> mat_T = UniversalMatrixConversion<T,int64_t>(comb_pair.pair.first);
        result.push_back(mat_T);
      }
    } else if (current_level == 3) {
      auto current_set = kernel_T->get_current_set();
      for (auto const& comb_pair : current_set) {
        result.push_back(comb_pair.pair.first);
      }
    }

    return result;
  }

  // Get current level information
  int get_current_level() const {
    return current_level;
  }
  size_t get_current_size() const {
    if (current_level == 0) return kernel_16->size();
    if (current_level == 1) return kernel_32->size();
    if (current_level == 2) return kernel_64->size();
    if (current_level == 3) return kernel_T->size();
    std::cerr << "ONL: We do not have the right level for OnlineExhaustiveReduction\n";
    throw TerminalException{1};
  }
};

// clang-format off
#endif  // SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_
// clang-format on
