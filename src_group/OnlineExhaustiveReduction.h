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
#define DEBUG_ONLINE_SIMPLIFICATION_SHIFT_NUMERICS
#endif

#ifdef DISABLE_DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
#undef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
#endif

#ifdef DISABLE_DEBUG_ONLINE_SIMPLIFICATION_SHIFT_NUMERICS
#undef DEBUG_ONLINE_SIMPLIFICATION_SHIFT_NUMERICS
#endif

template <typename U> using PairMatrix = std::pair<MyMatrix<U>, MyMatrix<U>>;

template <typename T2, typename T1>
PairMatrix<T2> UniversalPairMatrixConversion(PairMatrix<T1> const &pair) {
  MyMatrix<T2> m1 = UniversalMatrixConversion<T2, T1>(pair.first);
  MyMatrix<T2> m2 = UniversalMatrixConversion<T2, T1>(pair.second);
  return {std::move(m1), std::move(m2)};
}

template <typename T2, typename T1>
std::optional<PairMatrix<T2>>
UniversalPairMatrixConversionCheck(PairMatrix<T1> const &pair) {
  std::optional<MyMatrix<T2>> opt1 =
      UniversalMatrixConversionCheck<T2, T1>(pair.first);
  if (!opt1) {
    return {};
  }
  std::optional<MyMatrix<T2>> opt2 =
      UniversalMatrixConversionCheck<T2, T1>(pair.second);
  if (!opt2) {
    return {};
  }
  PairMatrix<T2> p{*opt1, *opt2};
  return p;
}

// Online/incremental version of ExhaustiveReductionComplexityKernelInner_V2
// Allows generators to be inserted one by one rather than providing a complete
// list
template <typename Ttype, typename Tnorm, typename Fcomplexity,
          typename Fproduct, typename Fcheck>
class OnlineExhaustiveReductionComplexityKernel {
private:
  std::map<TcombPair<Ttype, Tnorm>, BlockInterval> map;
  std::vector<TcombPair<Ttype, Tnorm>> vect;
  Fcomplexity f_complexity;
  Fproduct f_product;
  Fcheck f_check;
  std::ostream &os;

  // Internal helper functions from the original algorithm
  void delete_entry(TcombPair<Ttype, Tnorm> const &val) {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val should be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t pos = std::distance(map.begin(), iter);
    map.erase(iter);
    vect.erase(vect.begin() + pos);

    for (auto &kv : map) {
      BlockInterval &blk_int = kv.second;
      blk_int.remove_entry_and_shift(pos, os);
    }
  }

  size_t insert_entry_noop(TcombPair<Ttype, Tnorm> const &val) {
    BlockInterval blk_int;
    auto result = map.insert({val, blk_int});
    if (!result.second) {
      std::cerr << "ONL: Entry is already present. Unexpected\n";
      throw TerminalException{1};
    }
    auto iter = result.first;
    size_t pos = std::distance(map.begin(), iter);
    vect.insert(vect.begin() + pos, val);
    return pos;
  }

  void insert_entry(TcombPair<Ttype, Tnorm> const &val) {
    size_t pos = insert_entry_noop(val);
    // Setting the untreated entries
    size_t n_entry = vect.size();
    size_t idx = 0;
    for (auto &kv : map) {
      BlockInterval &blk_int = kv.second;
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
    std::vector<TcombPair<Ttype, Tnorm>> list_delete;
    std::vector<TcombPair<Ttype, Tnorm>> list_insert;
  };

  std::optional<FoundImprov> search_for_improvement() {
    for (auto &kv : map) {
      TcombPair<Ttype, Tnorm> const &x1 = kv.first;
      BlockInterval &blk_int = kv.second;
      while (true) {
        std::optional<size_t> opt = blk_int.get_first();
        if (opt) {
          size_t idx2 = *opt;
          TcombPair<Ttype, Tnorm> const &x2 = vect[idx2];
          GenResult<Ttype, Tnorm> gen =
              f_reduce<Ttype, Tnorm, Fcomplexity, Fproduct>(
                  x1, x2, f_complexity, f_product);
          if (gen.do_something) {
            bool x1_attained = false;
            bool x2_attained = false;
            std::vector<TcombPair<Ttype, Tnorm>> list_delete;
            std::vector<TcombPair<Ttype, Tnorm>> list_insert;
            for (auto &ent : gen.l_ent) {
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
    while (true) {
      std::optional<FoundImprov> opt = search_for_improvement();
      if (opt) {
        FoundImprov found_improv = *opt;
        // We need to do the insertions before the deletes
        // since if there is a failure in the reduction, we
        // want to avoid losing a generator when we switch to
        // another numerics.
        for (auto &val : found_improv.list_insert) {
          if (!f_check(val.pair.first) || !f_check(val.pair.second)) {
            return false;
          }
          insert_entry(val);
        }
        for (auto &val : found_improv.list_delete) {
          delete_entry(val);
        }
      } else {
        break;
      }
    }
    return true;
  }

public:
  OnlineExhaustiveReductionComplexityKernel(Fcomplexity _f_complexity,
                                            Fproduct _f_product,
                                            Fcheck _f_check, std::ostream &_os)
      : f_complexity(_f_complexity), f_product(_f_product), f_check(_f_check),
        os(_os) {}

  // Insert a new generator and run reduction
  bool insert_generator(std::pair<Ttype, Ttype> const &pair) {
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
    os << "gen1=" << compute_complexity_matrix(pair.first)
       << " gen2=" << compute_complexity_matrix(pair.second) << "\n";
#endif
    if (!f_check(pair.first) || !f_check(pair.second)) {
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
      os << "f_check fails at the beginning\n";
#endif
      return false;
    }
    TcombPair<Ttype, Tnorm> comb =
        generate_comb_pair<Ttype, Tnorm, Fcomplexity>(pair, f_complexity);
    if (map.count(comb) == 1) {
      return true; // Already present
    }

    insert_entry(comb);
    return run_reduction_pass();
  }

  bool insert_generator_noop(std::pair<Ttype, Ttype> const &pair) {
    if (!f_check(pair.first) || !f_check(pair.second)) {
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
      os << "f_check fails at the beginning\n";
#endif
      return false;
    }
    TcombPair<Ttype, Tnorm> comb =
        generate_comb_pair<Ttype, Tnorm, Fcomplexity>(pair, f_complexity);
    (void)insert_entry_noop(comb);
    return true;
  }

  std::vector<TcombPair<Ttype, Tnorm>> get_final_set() {
    std::vector<TcombPair<Ttype, Tnorm>> ret_vect = vect;
    vect.clear();
    map.clear();
    return ret_vect;
  }

  std::vector<TcombPair<Ttype, Tnorm>> get_current_set() const { return vect; }

  const std::vector<TcombPair<Ttype, Tnorm>> &get_current_set_ref() const {
    return vect;
  }

  // Get total complexity
  Tnorm get_total_complexity() const {
    Tnorm total_complexity(0);
    for (auto &kv : map) {
      total_complexity += kv.first.norm;
    }
    return total_complexity;
  }

  // Clear all generators
  void clear() {
    map.clear();
    vect.clear();
  }

  // Get number of generators
  size_t size() const { return vect.size(); }
};

template <typename Tfinite>
class OnlineExhaustiveReductionComplexityMatrixFinite {
private:
  OnlineExhaustiveReductionComplexityKernel<
      MyMatrix<Tfinite>, Tfinite,
      std::function<Tfinite(MyMatrix<Tfinite> const &)>,
      std::function<MyMatrix<Tfinite>(MyMatrix<Tfinite> const &,
                                      MyMatrix<Tfinite> const &)>,
      std::function<bool(MyMatrix<Tfinite> const &)>>
      inner;

public:
  OnlineExhaustiveReductionComplexityMatrixFinite(Tfinite max_val,
                                                  std::ostream &_os)
      : inner(
            [](MyMatrix<Tfinite> const &M) -> Tfinite {
              return get_ell1_complexity_measure(M);
            },
            [](MyMatrix<Tfinite> const &x, MyMatrix<Tfinite> const &y)
                -> MyMatrix<Tfinite> { return x * y; },
            [max_val](MyMatrix<Tfinite> const &x) -> bool {
              return check_matrix_coefficients(x, max_val);
            },
            _os) {}
  bool insert_generator(PairMatrix<Tfinite> const &M) {
    return inner.insert_generator(M);
  }
  template <typename Tinput>
  bool insert_generators_tinput(
      std::vector<TcombPair<MyMatrix<Tinput>, Tinput>> const &l_gen) {
    for (auto const &comb_pair : l_gen) {
      PairMatrix<Tfinite> pair =
          UniversalPairMatrixConversion<Tfinite, Tinput>(comb_pair.pair);
      if (!inner.insert_generator_noop(pair)) {
        clear();
        return false;
      }
    }
    return true;
  }
  std::vector<TcombPair<MyMatrix<Tfinite>, Tfinite>> get_final_set() {
    return inner.get_final_set();
  }
  template <typename T> std::vector<MyMatrix<T>> get_current_matrix_t() const {
    std::vector<MyMatrix<T>> ret_vect;
    const std::vector<TcombPair<MyMatrix<Tfinite>, Tfinite>> &vect =
        inner.get_current_set_ref();
    for (auto &pair : vect) {
      MyMatrix<T> M = UniversalMatrixConversion<T, Tfinite>(pair.pair.first);
      ret_vect.push_back(M);
    }
    return ret_vect;
  }
  void print_invariants(std::ostream &os_out) const {
    size_t seed = 10;
    const std::vector<TcombPair<MyMatrix<Tfinite>, Tfinite>> &vect =
        inner.get_current_set_ref();
    std::vector<MyMatrix<Tfinite>> ret_vect;
    os_out << "l_hash = [";
    for (auto &pair : vect) {
      MyMatrix<Tfinite> const &M = pair.pair.first;
      size_t hash = matrix_type_independent_hash(M, seed);
      os_out << " " << hash;
      ret_vect.push_back(M);
    }
    os_out << " ] invs=" << compute_complexity_listmat(ret_vect) << "\n";
  }
  size_t type_independent_hash(size_t const &seed) const {
    size_t ret_hash = seed;
    const std::vector<TcombPair<MyMatrix<Tfinite>, Tfinite>> &vect =
        inner.get_current_set_ref();
    for (auto &pair : vect) {
      MyMatrix<Tfinite> const &M = pair.pair.first;
      ret_hash = matrix_type_independent_hash(M, ret_hash);
    }
    return ret_hash;
  }
  std::vector<TcombPair<MyMatrix<Tfinite>, Tfinite>> get_current_set() const {
    return inner.get_current_set();
  }
  void clear() { inner.clear(); }
  size_t size() { return inner.size(); }
};

template <typename T> class OnlineExhaustiveReductionComplexityMatrixInfinite {
private:
  OnlineExhaustiveReductionComplexityKernel<
      MyMatrix<T>, T, std::function<T(MyMatrix<T> const &)>,
      std::function<MyMatrix<T>(MyMatrix<T> const &, MyMatrix<T> const &)>,
      std::function<bool(MyMatrix<T> const &)>>
      inner;

public:
  OnlineExhaustiveReductionComplexityMatrixInfinite(std::ostream &_os)
      : inner(
            [](MyMatrix<T> const &M) -> T {
              return get_ell1_complexity_measure(M);
            },
            [](MyMatrix<T> const &x, MyMatrix<T> const &y) -> MyMatrix<T> {
              return x * y;
            },
            []([[maybe_unused]] MyMatrix<T> const &x) -> bool { return true; },
            _os) {}
  bool insert_generator(PairMatrix<T> const &pair) {
    return inner.insert_generator(pair);
  }
  template <typename Tinput>
  bool insert_generators_tinput(
      std::vector<TcombPair<MyMatrix<Tinput>, Tinput>> const &l_gen) {
    for (auto const &comb_pair : l_gen) {
      PairMatrix<T> pair =
          UniversalPairMatrixConversion<T, Tinput>(comb_pair.pair);
      inner.insert_generator_noop(pair);
    }
    return true;
  }
  std::vector<TcombPair<MyMatrix<T>, T>> get_final_set() {
    return inner.get_final_set();
  }
  std::vector<MyMatrix<T>> get_current_matrix_t() const {
    std::vector<MyMatrix<T>> ret_vect;
    const std::vector<TcombPair<MyMatrix<T>, T>> &vect =
        inner.get_current_set_ref();
    for (auto &pair : vect) {
      ret_vect.push_back(pair.pair.first);
    }
    return ret_vect;
  }
  void print_invariants(std::ostream &os_out) const {
    size_t seed = 10;
    const std::vector<TcombPair<MyMatrix<T>, T>> &vect =
        inner.get_current_set_ref();
    std::vector<MyMatrix<T>> ret_vect;
    os_out << "l_hash = [";
    for (auto &pair : vect) {
      MyMatrix<T> const &M = pair.pair.first;
      size_t hash = matrix_type_independent_hash(M, seed);
      os_out << " " << hash;
      ret_vect.push_back(M);
    }
    os_out << " ] invs=" << compute_complexity_listmat(ret_vect) << "\n";
  }
  size_t type_independent_hash(size_t const &seed) const {
    size_t ret_hash = seed;
    const std::vector<TcombPair<MyMatrix<T>, T>> &vect =
        inner.get_current_set_ref();
    for (auto &pair : vect) {
      MyMatrix<T> const &M = pair.pair.first;
      ret_hash = matrix_type_independent_hash(M, ret_hash);
    }
    return ret_hash;
  }
  std::vector<TcombPair<MyMatrix<T>, T>> get_current_set() const {
    return inner.get_current_set();
  }
  void clear() { inner.clear(); }
  size_t size() const { return inner.size(); }
};

// Hierarchical online reduction system that automatically switches between
// int16_t -> int32_t -> int64_t -> T numerics when overflow occurs
template <typename T> class OnlineHierarchicalMatrixReduction {
private:
  // Type aliases for matrix-matrix pairs
  using Ttype = std::pair<MyMatrix<T>, MyMatrix<T>>;
  using Tnorm = T;

  template <typename Tin> using Tinput = TcombPair<MyMatrix<Tin>, Tin>;

  int n; // Matrix dimension
  std::ostream &os;
  // Current active kernel (0=int16_t, 1=int32_t, 2=int64_t, 3=T)
  int current_level;

  // Online kernels for each numeric type
  std::unique_ptr<OnlineExhaustiveReductionComplexityMatrixFinite<int16_t>>
      kernel_16;
  std::unique_ptr<OnlineExhaustiveReductionComplexityMatrixFinite<int32_t>>
      kernel_32;
  std::unique_ptr<OnlineExhaustiveReductionComplexityMatrixFinite<int64_t>>
      kernel_64;
  OnlineExhaustiveReductionComplexityMatrixInfinite<T> kernel_T;

  // Migrate from current level to next level
  void migrate_to_next_level() {
#ifdef DEBUG_ONLINE_SIMPLIFICATION_SHIFT_NUMERICS
    os << "SIMP: Migrating from level " << current_level << " to level "
       << (current_level + 1) << "\n";
#endif

    if (current_level == 0) {
      // Migrate from int16_t to int32_t or int64_t or T
      std::vector<Tinput<int16_t>> l_gen = kernel_16->get_current_set();
      bool result1 = kernel_32->insert_generators_tinput<int16_t>(l_gen);
      if (result1) {
        current_level = 1;
        return;
      }
      bool result2 = kernel_64->insert_generators_tinput<int16_t>(l_gen);
      if (result2) {
        current_level = 2;
        return;
      }
      kernel_T.template insert_generators_tinput<int16_t>(l_gen);
      current_level = 3;
    } else if (current_level == 1) {
      // Migrate from int32_t to int64_t or T
      std::vector<Tinput<int32_t>> l_gen = kernel_32->get_current_set();
      bool result1 = kernel_64->insert_generators_tinput<int32_t>(l_gen);
      if (result1) {
        current_level = 2;
        return;
      }
      kernel_T.template insert_generators_tinput<int32_t>(l_gen);
      current_level = 3;
    } else if (current_level == 2) {
      // Migrate from int64_t to T
      std::vector<Tinput<int64_t>> l_gen = kernel_64->get_current_set();
      kernel_T.template insert_generators_tinput<int64_t>(l_gen);
      current_level = 3;
    }
  }

public:
  OnlineHierarchicalMatrixReduction(int _n, std::ostream &_os)
      : n(_n), os(_os), current_level(0),
        kernel_T(OnlineExhaustiveReductionComplexityMatrixInfinite<T>(_os)) {
    // Calculate upper bounds for each numeric type
    double max_poss_16 =
        static_cast<double>(std::numeric_limits<int16_t>::max());
    double max_poss_32 =
        static_cast<double>(std::numeric_limits<int32_t>::max());
    double max_poss_64 =
        static_cast<double>(std::numeric_limits<int64_t>::max());

    int16_t max_val_16 = static_cast<int16_t>(sqrt(max_poss_16 / (10 * n)));
    int32_t max_val_32 = static_cast<int32_t>(sqrt(max_poss_32 / (10 * n)));
    int64_t max_val_64 = static_cast<int64_t>(sqrt(max_poss_64 / (10 * n)));

#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
    os << "SIMP: Upper bounds: int16_t=" << max_val_16
       << ", int32_t=" << max_val_32 << ", int64_t=" << max_val_64 << "\n";
#endif
    // Initialize all kernels
    kernel_16 = std::make_unique<
        OnlineExhaustiveReductionComplexityMatrixFinite<int16_t>>(max_val_16,
                                                                  os);
    kernel_32 = std::make_unique<
        OnlineExhaustiveReductionComplexityMatrixFinite<int32_t>>(max_val_32,
                                                                  os);
    kernel_64 = std::make_unique<
        OnlineExhaustiveReductionComplexityMatrixFinite<int64_t>>(max_val_64,
                                                                  os);
  }

  // Insert a single matrix (compute inverse in T first)
  void insert_generator(MyMatrix<T> const &generator) {
    // Compute inverse in T arithmetic first
    MyMatrix<T> generator_inv = Inverse(generator);
    return insert_generator_pair({generator, generator_inv});
  }

  bool simple_insert_generator_pair(Ttype const &generator_pair) {
    if (current_level == 0) {
      // Try int16_t
      std::optional<PairMatrix<int16_t>> opt =
          UniversalPairMatrixConversionCheck<int16_t, T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_16->insert_generator(*opt);
    } else if (current_level == 1) {
      // Try int32_t
      std::optional<PairMatrix<int32_t>> opt =
          UniversalPairMatrixConversionCheck<int32_t, T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_32->insert_generator(*opt);
    } else if (current_level == 2) {
      // Try int64_t
      std::optional<PairMatrix<int64_t>> opt =
          UniversalPairMatrixConversionCheck<int64_t, T>(generator_pair);
      if (!opt) {
        return false;
      }
      return kernel_64->insert_generator(*opt);
    } else {
      // Use T. infinite precision, cannot fail.
      return kernel_T.insert_generator(generator_pair);
    }
  }
  size_t type_independent_hash(size_t const &seed) const {
    if (current_level == 0) {
      return kernel_16->type_independent_hash(seed);
    } else if (current_level == 1) {
      return kernel_32->type_independent_hash(seed);
    } else if (current_level == 2) {
      return kernel_64->type_independent_hash(seed);
    } else {
      return kernel_T.type_independent_hash(seed);
    }
  }
  // Insert a matrix pair
  void insert_generator_pair(Ttype const &generator_pair) {
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
    size_t seed = 60;
#endif
    while (true) {
      bool test = simple_insert_generator_pair(generator_pair);
      if (test) {
        return;
      }
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
      size_t hash1 = type_independent_hash(seed);
      std::vector<MyMatrix<T>> l_mat1 = get_current_matrix_t();
      os << "ONL: Before, hash1=" << hash1 << " |l_mat1|=" << l_mat1.size()
         << "\n";
      print_invariants(os);
      //      os << "ONL: l_mat1(10)=\n";
      //      WriteMatrix(os, l_mat1[10]);
      //      os << "ONL: l_mat1(11)=\n";
      //      WriteMatrix(os, l_mat1[11]);
#endif
      migrate_to_next_level();
#ifdef DEBUG_ONLINE_EXHAUSTIVE_REDUCTION
      size_t hash2 = type_independent_hash(seed);
      std::vector<MyMatrix<T>> l_mat2 = get_current_matrix_t();
      os << "ONL: After, hash2=" << hash2 << " |l_mat2|=" << l_mat2.size()
         << "\n";
      print_invariants(os);
      //      os << "ONL: l_mat2(10)=\n";
      //      WriteMatrix(os, l_mat2[10]);
      //      os << "ONL: l_mat2(11)=\n";
      //      WriteMatrix(os, l_mat2[11]);
      if (hash1 == hash2) {
        os << "hash1 EQUALS hash2\n";
      } else {
        os << "hash1 does NOT equal hash2\n";
      }
#endif
    }
  }

  // Extract the final reduced generators (always in type T)
  std::vector<MyMatrix<T>> get_current_matrix_t() const {
    if (current_level == 0) {
      return kernel_16->get_current_matrix_t<T>();
    } else if (current_level == 1) {
      return kernel_32->get_current_matrix_t<T>();
    } else if (current_level == 2) {
      return kernel_64->get_current_matrix_t<T>();
    } else {
      return kernel_T.get_current_matrix_t();
    }
  }

  void print_invariants(std::ostream &os_out) const {
    if (current_level == 0) {
      return kernel_16->print_invariants(os_out);
    } else if (current_level == 1) {
      return kernel_32->print_invariants(os_out);
    } else if (current_level == 2) {
      return kernel_64->print_invariants(os_out);
    } else {
      return kernel_T.print_invariants(os_out);
    }
  }

  // Get current level information
  int get_current_level() const { return current_level; }

  size_t size() const {
    if (current_level == 0)
      return kernel_16->size();
    if (current_level == 1)
      return kernel_32->size();
    if (current_level == 2)
      return kernel_64->size();
    return kernel_T.size();
  }
};

// clang-format off
#endif  // SRC_GROUP_ONLINEEXHAUSTIVEREDUCTION_H_
// clang-format on
