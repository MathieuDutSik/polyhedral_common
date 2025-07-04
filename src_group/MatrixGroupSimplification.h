// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_
#define SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_

#include "MAT_Matrix.h"
#include "PermutationElt.h"
#include "Indefinite_LLL.h"

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_SIMPLIFICATION
#define DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
#endif

#ifdef TIMINGS
#define TIMINGS_MATRIX_GROUP_SIMPLIFICATION
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP_SIMPLIFICATION
#endif


template<typename T>
std::string compute_complexity_matrix(MyMatrix<T> const& mat) {
  int n = mat.rows();
  T ell1(0);
  T ellinfinity(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = mat(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
      if (abs_val > ellinfinity) {
        ellinfinity = abs_val;
      }
    }
  }
  return "(ell1=" + std::to_string(ell1) + ", ellinf=" + std::to_string(ellinfinity) + ")";
}

template<typename T>
std::string compute_complexity_listmat(std::vector<MyMatrix<T>> const& list_mat) {
  if (list_mat.size() == 0) {
    return "zero generators";
  }
  int n = list_mat[0].rows();
  size_t n_mat = list_mat.size();
  T ell1_global(0);
  T ellinfinite_global(0);
  for (auto & e_mat: list_mat) {
    T ell1(0);
    T ellinfinity(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        T val = e_mat(i,j);
        T abs_val = T_abs(val);
        ell1 += abs_val;
        if (abs_val > ellinfinity) {
          ellinfinity = abs_val;
        }
      }
    }
    ell1_global += ell1;
    if (ellinfinity > ellinfinite_global) {
      ellinfinite_global = ellinfinity;
    }
  }
  return "(n_gen=" + std::to_string(n_mat) + ", ell1_global=" + std::to_string(ell1_global) + ", ellinfinity=" + std::to_string(ellinfinite_global) + ")";
}

template<bool always_equal>
std::string compute_complexity_listseq(std::vector<permutalib::SequenceType<always_equal>> const& list_seq) {
  size_t n_seq = list_seq.size();
  size_t ell1_global = 0;
  size_t ellinfinite_global = 0;
  for (auto & seq: list_seq) {
    std::vector<int64_t> const& ListIdx = seq.getVect();
    size_t len = ListIdx.size();
    ell1_global += len;
    if (ellinfinite_global < len) {
      ellinfinite_global = len;
    }
  }
  return "(n_seq=" + std::to_string(n_seq) + ", ell1_global=" + std::to_string(ell1_global) + ", ellinfinity=" + std::to_string(ellinfinite_global) + ")";
}

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

template<typename T>
T get_ell1_complexity_measure(MyMatrix<T> const& M) {
  int n = M.rows();
  T ell1(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = M(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
    }
  }
  return ell1;
}

void print_vector_val(std::vector<size_t> const& V, std::ostream& os) {
  os << "[";
  for (size_t u=0; u<V.size(); u++) {
    if (u > 0) {
      os << ",";
    }
    os << V[u];
  }
  os << "]\n";
}


struct Interval {
  size_t start;
  size_t end;
};

struct BlockInterval {
  std::vector<Interval> intervals;
  BlockInterval() {
  }
  void insert_interval(size_t start, size_t end) {
    if (start < end) {
      Interval interval{start, end};
      intervals.push_back(interval);
    }
  }
  std::optional<size_t> get_first() {
    if (intervals.size() == 0) {
      return {};
    }
    Interval& interval = intervals[0];
    size_t val = interval.start;
    interval.start += 1;
    if (interval.start == interval.end) {
      intervals.erase(intervals.begin());
    }
    return val;
  }
  void print(std::ostream& os) const {
    size_t n_intervals = intervals.size();
    os << "SIMP: symbolic: ";
    for (size_t i_int=0; i_int<n_intervals; i_int++) {
      if (i_int > 0) {
        os << ", ";
      }
      size_t start = intervals[i_int].start;
      size_t end = intervals[i_int].end;
      os << "[" << start << "," << end << ")";
    }
    os << "\n";
  }
  size_t n_intervals() const {
    return intervals.size();
  }
  bool contains(size_t x) const {
    size_t low = 0;
    size_t high = intervals.size();

    while (low < high) {
      size_t mid = low + (high - low) / 2;
      const Interval& iv = intervals[mid];

      if (x < iv.start) {
        high = mid;
      } else if (x >= iv.end) {
        low = mid + 1;
      } else {
        // x is within [start, end)
        return true;
      }
    }
    return false;
  }
  // Used only for debugging purposes
  bool contains_direct(size_t x) const {
    size_t n_intervals = intervals.size();
    for (size_t i_int=0; i_int<n_intervals; i_int++) {
      if (intervals[i_int].start <= x && x < intervals[i_int].end) {
        return true;
      }
    }
    return false;
  }
  // Testing if the intervals are correct
  void test_correctness() const {
    size_t n_intervals = intervals.size();
    // Checking that the intervals are increasing
    for (size_t i=1; i<n_intervals; i++) {
      if (intervals[i-1].end > intervals[i].start) {
        std::cerr << "The intervals are overlapping\n";
        throw TerminalException{1};
      }
      if (intervals[i-1].end == intervals[i].start) {
        std::cerr << "The intervals should be merged\n";
        throw TerminalException{1};
      }
    }
  }

  std::vector<size_t> list_indices() const {
    std::vector<size_t> indices;
    for (auto & interval : intervals) {
      for (size_t u=interval.start; u<interval.end; u++) {
        indices.push_back(u);
      }
    }
    return indices;
  }


  // Remove an entry from the BlockInterval,
  // and shift the positions. This corresponds to a removal
  // from the underlying map.
  void remove_entry_and_shift_inner(size_t x) {
    size_t n_intervals = intervals.size();
    if (n_intervals == 0) {
      return;
    }
    if (x >= intervals[n_intervals - 1].end) {
      return;
    }
    size_t low = 0;
    size_t high = n_intervals;
    //#define DEBUG_REMOVE_ENTRY_AND_SHIFT
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
    std::cerr << "SIMP: ------------------- START remove_entry_and_shift_inner -------------------\n";
    std::cerr << "SIMP: low=" << low << " high=" << high << "\n";
#endif

    while (low < high) {
      size_t mid = low + (high - low) / 2;
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
      std::cerr << "SIMP: low=" << low << " high=" << high << " mid=" << mid << "\n";
#endif
      Interval& iv = intervals[mid];

      if (x < iv.start) {
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
        std::cerr << "SIMP: case 1\n";
#endif
        high = mid;
      } else if (x >= iv.end) {
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
        std::cerr << "SIMP: case 2\n";
#endif
        low = mid + 1;
      } else {
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
        std::cerr << "SIMP: case 3\n";
#endif
        // x is in [start, end)
        iv.end -= 1;
        size_t start_shift = mid + 1;;
        if (iv.start == iv.end) {
          intervals.erase(intervals.begin() + mid);
          start_shift = mid;
          n_intervals -= 1;
        }
        for (size_t u=start_shift; u<n_intervals; u++) {
          intervals[u].start -= 1;
          intervals[u].end -= 1;
        }
        return;
      }
    }
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
    std::cerr << "SIMP: case 4, low=" << low << "\n";
#endif
    auto iife_first_interval=[&]() -> size_t {
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
      if (low >= intervals.size()) {
        std::cerr << "SIMP: out of range access, low=" << low << " |intervals|=" << intervals.size() << " n_intervals=" << n_intervals << "\n";
        std::cerr << "SIMP: x=" << x << "\n";
        print(std::cerr);
        throw TerminalException{1};
      }
#endif
      if (x < intervals[low].start) {
        return low;
      }
      if (x >= intervals[low].end) {
        return low + 1;
      }
      std::cerr << "SIMP: We should never reach that stage in BlockInterval\n";
      throw TerminalException{1};
    };
    size_t index = iife_first_interval();
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
    std::cerr << "SIMP: case 4, index=" << index << "\n";
#endif
    for (size_t u=index; u<n_intervals; u++) {
      intervals[u].start -= 1;
      intervals[u].end -= 1;
    }
    // Merge the needed intervals if possible.
    if (index > 0) {
      if (intervals[index-1].end == intervals[index].start) {
        // Merging both intervals
        intervals[index-1].end = intervals[index].end;
        intervals.erase(intervals.begin() + index);
      }
    }
#ifdef DEBUG_REMOVE_ENTRY_AND_SHIFT
    std::cerr << "SIMP: case 4, exiting\n";
#endif
  }

  // Calling the function and checking it
  void remove_entry_and_shift(size_t x, [[maybe_unused]] std::ostream& os) {
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<Interval> saved_intervals = intervals;
    std::vector<size_t> indices_bef = list_indices();
#endif
    remove_entry_and_shift_inner(x);
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<size_t> indices_aft_expected;
    for (auto & y : indices_bef) {
      if (y < x) {
        indices_aft_expected.push_back(y);
      }
      if (y > x) {
        indices_aft_expected.push_back(y - 1);
      }
    }
    std::vector<size_t> indices_aft_real = list_indices();
    if (indices_aft_expected != indices_aft_real) {
      std::swap(intervals, saved_intervals);
      os << "SIMP: -------------------------------------------\n";
      os << "SIMP: remove_entry_and_shift, before x=" << x << "\n";
      print(os);
      os << "SIMP: indices_bef: ";
      print_vector_val(indices_bef, os);
      std::swap(intervals, saved_intervals);
      os << "SIMP: remove_entry_and_shift, after\n";
      print(os);
      os << "SIMP: indices_aft_expected: ";
      print_vector_val(indices_aft_expected, os);
      os << "SIMP: indices_aft_real: ";
      print_vector_val(indices_aft_real, os);
      std::cerr << "SIMP: indices_aft_expected does not match indices_aft_real (A)\n";
      throw TerminalException{1};
    }
    test_correctness();
#endif
  }


  // Add one entry at a position and shift.
  // This corresponds to an insertion in the map.
  void insert_entry_and_shift_inner(size_t x) {
    size_t n_intervals = intervals.size();
    if (n_intervals == 0) {
      Interval iv{x, x + 1};
      intervals.push_back(iv);
      return;
    }
    if (x == intervals[n_intervals-1].end) {
      intervals[n_intervals-1].end += 1;
      return;
    }
    if (x > intervals[n_intervals-1].end) {
      Interval iv{x, x+1};
      intervals.push_back(iv);
      return;
    }
    size_t low = 0;
    size_t high = n_intervals;
    //    std::cerr << "SIMP: low=" << low << " high=" << high << "\n";

    while (low < high) {
      size_t mid = low + (high - low) / 2;
      //      std::cerr << "SIMP: low=" << low << " high=" << high << " mid=" << mid << "\n";
      Interval& iv = intervals[mid];

      if (x < iv.start) {
        //        std::cerr << "SIMP: Case 1\n";
        high = mid;
      } else if (x > iv.end) {
        //        std::cerr << "SIMP: Case 2\n";
        low = mid + 1;
      } else {
        //        std::cerr << "SIMP: Case 3\n";
        // x is in [start, end)
        iv.end += 1;
        size_t start_shift = mid + 1;;
        for (size_t u=start_shift; u<n_intervals; u++) {
          intervals[u].start += 1;
          intervals[u].end += 1;
        }
        return;
      }
    }
    //    std::cerr << "SIMP: Case 4 low=" << low << "\n";
    // x should be outside of the intervals.
    auto iife_first_interval=[&]() -> size_t {
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
      if (low >= intervals.size()) {
        std::cerr << "SIMP: out of range access. low=" << low << "\n";
        std::cerr << "SIMP: x=" << x << "\n";
        print(std::cerr);
        throw TerminalException{1};
      }
#endif
      if (x < intervals[low].start) {
        return low;
      }
      if (x >= intervals[low].end) {
        return low + 1;
      }
      std::cerr << "SIMP: We should never reach that stage in BlockInterval\n";
      throw TerminalException{1};
    };
    size_t index = iife_first_interval();
    //    std::cerr << "SIMP: Case 4 index=" << index << "\n";
    Interval new_iv{x, x+1};
    intervals.insert(intervals.begin() + index, new_iv);
    //    std::cerr << "SIMP: Case 4 |intervals|=" << intervals.size() << "\n";
    for (size_t u=index + 1; u<n_intervals + 1; u++) {
      //      std::cerr << "SIMP: Case 4 u=" << u << "\n";
      intervals[u].start += 1;
      intervals[u].end += 1;
    }
  }

  void insert_entry_and_shift(size_t x, [[maybe_unused]] std::ostream& os) {
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<Interval> saved_intervals = intervals;
    std::vector<size_t> indices_bef = list_indices();
#endif
    insert_entry_and_shift_inner(x);
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<size_t> indices_aft_expected;
    indices_aft_expected.push_back(x);
    for (auto & y : indices_bef) {
      if (y < x) {
        indices_aft_expected.push_back(y);
      }
      if (y >= x) {
        indices_aft_expected.push_back(y + 1);
      }
    }
    std::sort(indices_aft_expected.begin(), indices_aft_expected.end());
    std::vector<size_t> indices_aft_real = list_indices();
    if (indices_aft_expected != indices_aft_real) {
      std::swap(intervals, saved_intervals);
      os << "SIMP: -------------------------------------------\n";
      os << "SIMP: insert_entry_and_shift, before x=" << x << "\n";
      print(os);
      os << "SIMP: indices_bef: ";
      print_vector_val(indices_bef, os);
      std::swap(intervals, saved_intervals);
      os << "SIMP: insert_entry_and_shift, after\n";
      print(os);
      os << "SIMP: indices_aft_expected: ";
      print_vector_val(indices_aft_expected, os);
      os << "SIMP: indices_aft_real: ";
      print_vector_val(indices_aft_real, os);
      std::cerr << "SIMP: indices_aft_expected does not match indices_aft_real (B)\n";
      throw TerminalException{1};
    }
    test_correctness();
#endif
  }

  // Shift the index and
  void noinsert_and_shift_inner(size_t x) {
    size_t n_intervals = intervals.size();
    size_t low = 0;
    size_t high = n_intervals;
    if (n_intervals == 0) {
      return;
    }
    //#define DEBUG_NOINSERT_AND_SHIFT
#ifdef DEBUG_NOINSERT_AND_SHIFT
    std::cerr << "SIMP: ------------------- START noinsert_and_shift_inner -------------------\n";
    std::cerr << "SIMP: low=" << low << " high=" << high << "\n";
#endif

    while (low < high) {
      size_t mid = low + (high - low) / 2;
#ifdef DEBUG_NOINSERT_AND_SHIFT
      std::cerr << "SIMP: low=" << low << " high=" << high << " mid=" << mid << "\n";
#endif
      Interval& iv = intervals[mid];

      if (x < iv.start) {
#ifdef DEBUG_NOINSERT_AND_SHIFT
        std::cerr << "SIMP: Case 1\n";
#endif
        high = mid;
      } else if (x >= iv.end) {
#ifdef DEBUG_NOINSERT_AND_SHIFT
        std::cerr << "SIMP: Case 2\n";
#endif
        low = mid + 1;
      } else {
#ifdef DEBUG_NOINSERT_AND_SHIFT
        std::cerr << "SIMP: Case 3\n";
#endif
        // x is in [start, end)
        if (x == iv.start) {
#ifdef DEBUG_NOINSERT_AND_SHIFT
          std::cerr << "SIMP: Case 3.1\n";
#endif
          // First element in the list
          for (size_t u=mid; u<n_intervals; u++) {
            intervals[u].start += 1;
            intervals[u].end += 1;
          }
          return;
        }
        if (x == iv.end) {
#ifdef DEBUG_NOINSERT_AND_SHIFT
          std::cerr << "SIMP: Case 3.2\n";
#endif
          // Last element in the list
          for (size_t u=mid+1; u<n_intervals; u++) {
            intervals[u].start += 1;
            intervals[u].end += 1;
          }
          return;
        }
#ifdef DEBUG_NOINSERT_AND_SHIFT
        std::cerr << "SIMP: Case 3.3, mid=" << mid << "\n";
#endif
        // Break down the interval in two places.
        size_t end = iv.end;
        Interval new_iv{x + 1, end+1};
        iv.end = x;
        intervals.insert(intervals.begin() + mid + 1, new_iv);
#ifdef DEBUG_NOINSERT_AND_SHIFT
        std::cerr << "SIMP: Case 3.3, |intervals|=" << intervals.size() << " n_intervals=" << n_intervals << "\n";
#endif
        for (size_t u=mid + 2; u<n_intervals + 1; u++) {
          intervals[u].start += 1;
          intervals[u].end += 1;
        }
        return;
      }
    }
#ifdef DEBUG_NOINSERT_AND_SHIFT
    std::cerr << "SIMP: Case 4\n";
#endif
    // x should be outside of the intervals.
    auto iife_first_interval=[&]() -> size_t {
      if (x < intervals[low].start) {
        return low;
      }
      if (x >= intervals[low].end) {
        return low + 1;
      }
      std::cerr << "SIMP: We should never reach that stage in BlockInterval\n";
      throw TerminalException{1};
    };
    size_t index = iife_first_interval();
#ifdef DEBUG_NOINSERT_AND_SHIFT
    std::cerr << "SIMP: Case 4, index=" << index << "\n";
#endif
    for (size_t u=index; u<n_intervals; u++) {
      intervals[u].start += 1;
      intervals[u].end += 1;
    }
  }

  void noinsert_and_shift(size_t x, [[maybe_unused]] std::ostream& os) {
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<Interval> saved_intervals = intervals;
    std::vector<size_t> indices_bef = list_indices();
#endif
    noinsert_and_shift_inner(x);
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    std::vector<size_t> indices_aft_expected;
    for (auto & y : indices_bef) {
      if (y < x) {
        indices_aft_expected.push_back(y);
      }
      if (y >= x) {
        indices_aft_expected.push_back(y + 1);
      }
    }
    std::vector<size_t> indices_aft_real = list_indices();
    if (indices_aft_expected != indices_aft_real) {
      std::swap(intervals, saved_intervals);
      os << "SIMP: -------------------------------------------\n";
      os << "SIMP: noinsert_and_shift, before x=" << x << "\n";
      print(os);
      os << "SIMP: indices_bef: ";
      print_vector_val(indices_bef, os);
      os << "SIMP: -------\n";
      std::swap(intervals, saved_intervals);
      os << "SIMP: noinsert_and_shift, after\n";
      print(os);
      os << "SIMP: indices_aft_expected: ";
      print_vector_val(indices_aft_expected, os);
      os << "SIMP: indices_aft_real: ";
      print_vector_val(indices_aft_real, os);
      std::cerr << "SIMP: indices_aft_expected does not match indices_aft_real (C)\n";
      throw TerminalException{1};
    }
    test_correctness();
#endif
  }

};













// The tool for simplifying a list of generators.
// The transformations being applied are Tietze transformations.
// We could add some conjugacy operations like  U V U^{-1}.
// But those are more expensive
//
// The existing design is not good as there are plenty of
// repeated operations. Once we find a new simplification, we redo
// almost everything and that is wildly inefficient.
//
// So, what is a good design?
// * First of all the ordering by the complexity makes full sense.
//   No discussion about that.
// * The use of nonces is also perfectly fine and could be
//   generalized. Or eliminated if we completely avoid having
//   repetitions!
// * What we mean when we treat a new element is that we match it
//   against all the elements. The big problem is what to do when
//   an interruption happens.
// * When we have an interruption, we have a change of the order
//   and we need to adjust immediately.
// * We could have a storing of the data that need to be treated
//   by an interval at first. When changes are made, then the
//   structure is adjusted. Intervals are split and so on.
//   The underlying idea is that yes splitting would occur but
//   they would not have a very detrimental impact.
// * So, the structure could be:
//   struct Interval {
//     size_t start;
//     size_t end;
//   }
//   struct BlockInterval {
//     std::vector<Interval> intervals;
//   }
//   with the following functions:
//   + Removal of one element from the list
//   + Insert an element at a specified position (and shift the intervals).
//   + Iterating: We only need one function:
//     Access to the first element of the list (returns an std::optional<size_t>)
//     and remove it.
template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexityKernel_V2(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, std::ostream& os) {
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  NanosecondTime time_total;
#endif
  using TtypePair = std::pair<Ttype, Ttype>;
  using TcombPair = std::pair<TtypePair, Tnorm>;
  using Tcomb = std::pair<Ttype, Tnorm>;
  auto f_comp=[](TcombPair const& a, TcombPair const& b) -> bool {
    if (a.second < b.second) {
      return true;
    }
    if (a.second > b.second) {
      return false;
    }
    return a.first.first < b.first.first;
  };
  auto get_comb=[&](Ttype const& eM) -> Tcomb {
    Tnorm comp = f_complexity(eM);
    return {eM, comp};
  };
  auto get_comb_pair=[&](Tcomb const& p) -> TcombPair {
    Ttype p_inv = f_invers(p.first);
    TtypePair p_pair{p.first, p_inv};
    return {p_pair, p.second};
  };
  std::map<TcombPair, BlockInterval, decltype(f_comp)> map(f_comp);
  size_t n_matrix = ListM.size();
  for (auto & eM: ListM) {
    Tcomb comb1 = get_comb(eM);
    TcombPair comb2 = get_comb_pair(comb1);
    BlockInterval blk_int;
    map[comb2] = blk_int;
  }
  size_t index = 0;
  for (auto & kv: map) {
    BlockInterval & blk_int = kv.second;
    blk_int.insert_interval(index + 1, n_matrix);
    index += 1;
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  auto get_complexity=[&]() -> Tnorm {
    Tnorm total_complexity(0);
    for (auto & kv: map) {
      total_complexity += kv.first.second;
    }
    return total_complexity;
  };
  os << "SIMP: total_complexity(begin)=" << get_complexity() << "\n";
#endif
  // We need a vector since the map is not random access.
  // That is using iterators and advancing them has complexity O(n) for
  // the map but O(1) for the vector.
  std::vector<TcombPair> vect;
  for (auto & kv: map) {
    vect.push_back(kv.first);
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  auto check_map_vect=[&](std::string const& context) -> void {
    size_t index = 0;
    for (auto & kv: map) {
      if (kv.first != vect[index]) {
        std::cerr << "SIMP: Error in map/vect at index=" << index << " context=" << context << "\n";
        throw TerminalException{1};
      }
      index += 1;
    }
  };
  check_map_vect("initial");
#endif
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  size_t n_get_best_candidate = 0;
#endif
  // Generate the possible ways to simplify the pair of elements.
  // and select the one with the smallest norm
  auto f_get_best_candidate=[&](TcombPair const& a, TcombPair const& b) -> Tcomb {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: f_get_best_candidate, beginning\n";
#endif
    Ttype const& a_dir = a.first.first;
    Ttype const& b_dir = b.first.first;
    Ttype const& a_inv = a.first.second;
    Ttype const& b_inv = b.first.second;
    //
    Ttype prod_gen = f_product(a_dir, b_dir);
    Tcomb pair_ret = get_comb(prod_gen);
    //
    prod_gen = f_product(a_inv, b_dir);
    Tcomb pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
    prod_gen = f_product(a_dir, b_inv);
    pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
    prod_gen = f_product(a_inv, b_inv);
    pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
    n_get_best_candidate += 1;
#endif
    return pair_ret;
  };
  // Iterate the reduction algorithm over pairs of elements.
  // The result of the iteration might be a 0, 1 or 2 new elements.
  auto f_reduce=[&](TcombPair const& a, TcombPair const& b) -> std::pair<size_t, std::vector<TcombPair>> {
    TcombPair a_work = a;
    TcombPair b_work = b;
    size_t n_change = 0;
    while(true) {
      Tcomb cand = f_get_best_candidate(a_work, b_work);
      Tnorm a_norm = a_work.second;
      Tnorm b_norm = b_work.second;
      Tnorm cand_norm = cand.second;
      bool do_something = true;
      if (cand_norm < a_norm && cand_norm < b_norm) {
        if (a_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          a_work = get_comb_pair(cand);
        }
      } else {
        if (cand_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          if (cand_norm < a_norm) {
            a_work = get_comb_pair(cand);
          } else {
            do_something = false;
          }
        }
      }
      if (!do_something) {
        std::vector<TcombPair> npair{a_work, b_work};
        return {n_change, npair};
      } else {
        n_change += 1;
      }
    }
  };
  struct FoundImprov {
    std::vector<TcombPair> list_delete;
    std::vector<TcombPair> list_insert;
  };
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  size_t i_search = 0;
#endif
  auto f_search=[&]() -> std::optional<FoundImprov> {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    size_t n_reduce_calls = 0;
    size_t idx1 = 0;
#endif
    for (auto & kv: map) {
      TcombPair const& x1 = kv.first;
      BlockInterval & blk_int = kv.second;
      while(true) {
        std::optional<size_t> opt = blk_int.get_first();
        if (opt) {
          size_t idx2 = *opt;
          TcombPair const& x2 = vect[idx2];
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
          n_reduce_calls += 1;
#endif
          std::pair<size_t, std::vector<TcombPair>> pair = f_reduce(x1, x2);
          size_t n_changes = pair.first;
          if (n_changes > 0) {
            bool x1_attained = false;
            bool x2_attained = false;
            std::vector<TcombPair> list_delete;
            std::vector<TcombPair> list_insert;
            for (auto & ent : pair.second) {
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
            FoundImprov found_improv{list_delete, list_insert};
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
            os << "SIMP: f_search, i_search=" << i_search
               << " |map|=" << map.size()
               << " n_reduce_calls=" << n_reduce_calls
               << " n_changes=" << n_changes
               << " |list_insert|=" << list_insert.size()
               << " |list_delete|=" << list_delete.size()
               << " x1_att=" << x1_attained << " x2_att=" << x2_attained
               << " idx1=" << idx1 << " idx2=" << idx2 << "\n";
#endif
            return found_improv;
          }
        } else {
          break;
        }
      }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      idx1 += 1;
#endif
    }
    return {};
  };
  //
  // Update the data structure of map/vect
  //
  auto delete_entry=[&](TcombPair const& val) -> void {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val should be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t pos = std::distance(map.begin(), iter);
    map.erase(iter);
    vect.erase(vect.begin() + pos);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    std::string context = "delete_entry_at_" + std::to_string(pos);
    check_map_vect(context);
#endif
    for (auto & kv : map) {
      BlockInterval & blk_int = kv.second;
      blk_int.remove_entry_and_shift(pos, os);
    }
  };
  auto insert_entry=[&](TcombPair const& val) -> void {
    BlockInterval blk_int;
#ifdef SANITY_CHECK_MATRIX_GROUP_SIMPLIFICATION
    auto iter_deb = map.find(val);
    if (iter_deb != map.end()) {
      size_t pos_deb = std::distance(map.begin(), iter_deb);
      std::cerr << "SIMP: insert_entry, entry val already at position " << pos_deb << "\n";
      throw TerminalException{1};
    }
#endif
    auto result = map.insert({val, blk_int});
    if (!result.second) {
      std::cerr << "SIMP: Entry is already present. Unexpected\n";
      throw TerminalException{1};
    }
    auto iter = result.first;
    size_t pos = std::distance(map.begin(), iter);
    vect.insert(vect.begin() + pos, val);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    std::string context = "insert_entry_at_" + std::to_string(pos);
    check_map_vect(context);
#endif
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
  };

  // Now iterating looking for improvements.
  while(true) {
    std::optional<FoundImprov> opt = f_search();
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    i_search += 1;
#endif
    if (opt) {
      FoundImprov found_improv = *opt;
      for (auto & val : found_improv.list_delete) {
        delete_entry(val);
      }
      for (auto & val : found_improv.list_insert) {
        insert_entry(val);
      }
    } else {
      break;
    }
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  os << "SIMP: total_complexity(after)=" << get_complexity() << "\n";
#endif
  std::vector<Ttype> new_list_gens;
  for (auto & kv: map) {
    new_list_gens.push_back(kv.first.first.first);
  }
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  int64_t delta = time_total.eval_int64();
  double delta_d = static_cast<double>(delta);
  double n_get_best_candidate_d = static_cast<double>(n_get_best_candidate);
  double avg_cost_best = delta_d / n_get_best_candidate_d;
  os << "|SIMP: ExhaustiveReductionComplexityKernel, avg(f_get_best_candidate)|=" << avg_cost_best << "\n";
  os << "|SIMP: ExhaustiveReductionComplexityKernel|=" << delta << "\n";
#endif
  return new_list_gens;
}

template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexity_V2(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, std::ostream& os) {
  auto f_total_comp=[&](std::vector<Ttype> const& ListM) -> Tnorm {
    Tnorm Tcomp(0);
    for (auto & eM : ListM) {
      Tcomp += f_complexity(eM);
    }
    return Tcomp;
  };
  Tnorm curr_total_comp = f_total_comp(ListM);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  os << "SIMP: total_complexity(start)=" << curr_total_comp << "\n";
#endif
  std::vector<Ttype> ListMwork = ListM;
  while(true) {
    std::vector<Ttype> ListMwork = ExhaustiveReductionComplexityKernel_V2<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListMwork, f_complexity, f_invers, f_product, os);
    Tnorm new_total_comp = f_total_comp(ListMwork);
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: new_total_comp=" << new_total_comp << " curr_total_comp=" << curr_total_comp << "\n";
#endif
    if (new_total_comp == curr_total_comp) {
      return ListMwork;
    }
    curr_total_comp = new_total_comp;
  }
}





template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexityKernel_V1(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, [[maybe_unused]] std::ostream& os) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  using TtypePair = std::pair<Ttype, Ttype>;
  using TcombPair = std::pair<TtypePair, Tnorm>;
  using Tcomb = std::pair<Ttype, Tnorm>;
  auto f_comp=[](TcombPair const& a, TcombPair const& b) -> bool {
    if (a.second < b.second) {
      return true;
    }
    if (a.second > b.second) {
      return false;
    }
    return a.first.first < b.first.first;
  };
  auto get_comb=[&](Ttype const& eM) -> Tcomb {
    Tnorm comp = f_complexity(eM);
    return {eM, comp};
  };
  auto get_comb_pair=[&](Tcomb const& p) -> TcombPair {
    Ttype p_inv = f_invers(p.first);
    TtypePair p_pair{p.first, p_inv};
    return {p_pair, p.second};
  };
  std::map<TcombPair, size_t, decltype(f_comp)> map(f_comp);
  size_t nonce = 0;
  for (auto & eM: ListM) {
    Tcomb comb1 = get_comb(eM);
    TcombPair comb2 = get_comb_pair(comb1);
    map[comb2] = nonce;
    nonce += 1;
  }
  std::unordered_set<std::pair<size_t,size_t>> set_treated;
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  NanosecondTime time_total;
#endif
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  size_t n_operation = 0;
  size_t i_run = 0;
#endif
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  size_t n_get_best_candidate = 0;
#endif
  // Generate the possible ways to simplify the pair of elements.
  // and select the one with the smallest norm
  auto f_get_best_candidate=[&](TcombPair const& a, TcombPair const& b) -> Tcomb {
    Ttype const& a_dir = a.first.first;
    Ttype const& b_dir = b.first.first;
    Ttype const& a_inv = a.first.second;
    Ttype const& b_inv = b.first.second;
    //
    Ttype prod_gen = f_product(a_dir, b_dir);
    Tcomb pair_ret = get_comb(prod_gen);
    //
    prod_gen = f_product(a_inv, b_dir);
    Tcomb pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
    prod_gen = f_product(a_dir, b_inv);
    pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
    prod_gen = f_product(a_inv, b_inv);
    pair_gen = get_comb(prod_gen);
    if (pair_gen.second < pair_ret.second) {
      pair_ret = pair_gen;
    }
    //
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
    n_get_best_candidate += 1;
#endif
    return pair_ret;
  };
  // Iterate the reduction algorithm over pairs of elements.
  // The result of the iteration might be a 0, 1 or 2 new elements.
  auto f_reduce=[&](TcombPair const& a, TcombPair const& b) -> std::pair<size_t, std::vector<TcombPair>> {
    TcombPair a_work = a;
    TcombPair b_work = b;
    size_t n_change = 0;
    while(true) {
      Tcomb cand = f_get_best_candidate(a_work, b_work);
      Tnorm a_norm = a_work.second;
      Tnorm b_norm = b_work.second;
      Tnorm cand_norm = cand.second;
      bool do_something = true;
      if (cand_norm < a_norm && cand_norm < b_norm) {
        if (a_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          a_work = get_comb_pair(cand);
        }
      } else {
        if (cand_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          if (cand_norm < a_norm) {
            a_work = get_comb_pair(cand);
          } else {
            do_something = false;
          }
        }
      }
      if (!do_something) {
        std::vector<TcombPair> npair{a_work, b_work};
        return {n_change, npair};
      } else {
        n_change += 1;
      }
    }
  };
  auto erase_entry=[&](TcombPair const& val) -> void {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val shoulf be present in order to be removed\n";
      throw TerminalException{1};
    }
    map.erase(iter);
  };
  auto get_pos=[&](TcombPair const& val) -> size_t {
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
      os << "SIMP: Starting increment_uv with u=" << u << " v=" << v << "\n";
#endif
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
    for (auto & kv: map) {
      os << "SIMP: pos=" << pos << " nonce=" << kv.second << " norm=" << kv.first.second << " eM=\n";
      WriteMatrix(os, kv.first.first);
      pos += 1;
    }
#endif
    while(true) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP: " << u << " / " << v << " n_operation=" << n_operation << " |set|=" << map.size() << " |set_treated|=" << set_treated.size() << " i_run=" << i_run << "\n";
#endif
      auto iter1 = map.begin();
      std::advance(iter1, u);
      TcombPair a = iter1->first;
      size_t nonce_a = iter1->second;
      //
      auto iter2 = map.begin();
      std::advance(iter2, v);
      TcombPair b = iter2->first;
      size_t nonce_b = iter2->second;
      std::pair<size_t, size_t> nonce_pair{nonce_a, nonce_b};
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
      os << "SIMP:   Complexities, a.second=" << a.second << " b.second=" << b.second << "\n";
#endif
      //
      bool already_treated = false;
      std::pair<size_t, std::vector<TcombPair>> pair{0,{}};
      if (set_treated.find(nonce_pair) != set_treated.end()) {
        already_treated = true;
      } else {
        pair = f_reduce(a, b);
      }
      size_t n_changes = pair.first;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
      os << "SIMP:   n_changes=" << n_changes << " already_treated=" << already_treated << " |set_treated|=" << set_treated.size() << "\n";
#endif
      if (n_changes > 0) {
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
        std::vector<TcombPair> new_elt;
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
        os << "SIMP:  |new_elt|=" << new_elt.size() << "\n";
#endif
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  Scenario A\n";
#endif
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  Scenario B, |new_elt|=" << new_elt.size() << " min_distance_bis=" << min_distance_bis << "\n";
#endif
          // We have a new generator, adjusting accordingly.
          if (!a_attained) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
            os << "SIMP:  Scenario B, A\n";
#endif
            if (min_distance_bis == map.size() - 1) {
              // We reach the end of what we can do.
              break;
            }
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
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
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    i_run += 1;
#endif
    total_complexity = new_complexity;
  }
  std::vector<Ttype> new_list_gens;
  for (auto & kv: map) {
    new_list_gens.push_back(kv.first.first.first);
  }
#ifdef TIMINGS_MATRIX_GROUP_SIMPLIFICATION
  int64_t delta = time_total.eval_int64();
  double delta_d = static_cast<double>(delta);
  double n_get_best_candidate_d = static_cast<double>(n_get_best_candidate);
  double avg_cost_best = delta_d / n_get_best_candidate_d;
  os << "|SIMP: ExhaustiveReductionComplexityKernel, avg(f_get_best_candidate)|=" << avg_cost_best << "\n";
#endif
  return new_list_gens;
}

template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexityKernel(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, std::ostream& os) {
#ifdef METHOD_COMPARISON_MATRIX_GROUP_SIMPLIFICATION
  MicrosecondTime time;
  std::vector<Ttype> result_V1 = ExhaustiveReductionComplexityKernel_V1<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListM, f_complexity, f_invers, f_product, os);
  os << "|SIMP: ExhaustiveReductionComplexityKernel_V1|=" << time << " |result_V1|=" << result_V1.size() << "\n";
  std::vector<Ttype> result_V2 = ExhaustiveReductionComplexityKernel_V2<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListM, f_complexity, f_invers, f_product, os);
  os << "|SIMP: ExhaustiveReductionComplexityKernel_V2|=" << time << " |result_V2|=" << result_V2.size() << "\n";
  return result_V2;
#else
  return ExhaustiveReductionComplexityKernel_V1<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListM, f_complexity, f_invers, f_product, os);
#endif
}



template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexity(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, std::ostream& os) {
  std::unordered_set<Ttype> SetMred;
  for (auto & eM : ListM) {
    Ttype eM_inv = f_invers(eM);
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
  if (ListMred.size() <= 1) {
    return ListMred;
  }
  return ExhaustiveReductionComplexityKernel<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListMred, f_complexity, f_invers, f_product, os);
}

template<typename T>
std::vector<MyMatrix<T>> ExhaustiveReductionComplexityGroupMatrix(std::vector<MyMatrix<T>> const& ListM, std::ostream& os) {
#ifdef TRACK_INFO_MATRIX_GROUP_SIMPLIFICATION
  write_matrix_group(ListM, "Call_to_ExhaustiveReductionComplexityGroupMatrix");
#endif
  auto f_complexity=[&](MyMatrix<T> const& M) -> T {
    return get_ell1_complexity_measure(M);
  };
  auto f_invers=[](MyMatrix<T> const& M) -> MyMatrix<T> {
    return Inverse(M);
  };
  auto f_product=[&](MyMatrix<T> const& A, MyMatrix<T> const& B) -> MyMatrix<T> {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: f_product, A=";
    WriteMatrix(os, A);
    os << "SIMP: f_product, B=";
    WriteMatrix(os, B);
    MyMatrix<T> result = A * B;
    os << "SIMP: f_product, result=";
    WriteMatrix(os, result);
#endif
    return A * B;
  };
  return ExhaustiveReductionComplexity<T,MyMatrix<T>,decltype(f_complexity),decltype(f_invers),decltype(f_product)>(ListM, f_complexity, f_invers, f_product, os);
}

std::vector<permutalib::SequenceType<false>> ExhaustiveReductionComplexitySequences(std::vector<permutalib::SequenceType<false>> const& ListS, std::ostream& os) {
  using Tseq = permutalib::SequenceType<false>;
  auto f_complexity=[&](Tseq const& x) -> size_t {
    return x.complexity();
  };
  auto f_invers=[&](Tseq const& x) -> Tseq {
    return Inverse(x);
  };
  auto f_product=[&](Tseq const& x, Tseq const& y) -> Tseq {
    return x * y;
  };
  return ExhaustiveReductionComplexity<size_t,Tseq,decltype(f_complexity),decltype(f_invers),decltype(f_product)>(ListS, f_complexity, f_invers, f_product, os);
}

template<typename T, typename Telt>
std::pair<std::vector<MyMatrix<T>>, std::vector<Telt>> ExhaustiveReductionComplexityGroupMatrixPerm(std::vector<MyMatrix<T>> const& ListM, std::vector<Telt> const& ListPerm, std::ostream& os) {
  using Ttype = std::pair<MyMatrix<T>, Telt>;
  auto f_complexity=[&](Ttype const& pair) -> T {
    return get_ell1_complexity_measure(pair.first);
  };
  auto f_invers=[](Ttype const& pair) -> Ttype {
    return {Inverse(pair.first), Inverse(pair.second)};
  };
  auto f_product=[](Ttype const& p1, Ttype const& p2) -> Ttype {
    return {p1.first * p2.first, p1.second * p2.second};
  };
  std::vector<Ttype> ListPair;
  size_t n_gen = ListM.size();
  for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
    Ttype pair{ListM[i_gen], ListPerm[i_gen]};
    ListPair.push_back(pair);
  }
  std::vector<Ttype> RetPair = ExhaustiveReductionComplexity<T,Ttype,decltype(f_complexity),decltype(f_invers),decltype(f_product)>(ListPair, f_complexity, f_invers, f_product, os);
  std::vector<MyMatrix<T>> RetListM;
  std::vector<Telt> RetListPerm;
  size_t n_gen_ret = RetPair.size();
  for (size_t i_gr=0; i_gr<n_gen_ret; i_gr++) {
    RetListM.push_back(RetPair[i_gr].first);
    RetListPerm.push_back(RetPair[i_gr].second);
  }
  return {RetListM, RetListPerm};
}



template<typename T>
std::vector<MyMatrix<T>> Exhaust_get_total_generators(std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_tot;
  for (auto & eMat: list_mat) {
    MyMatrix<T> eMat_inv = Inverse(eMat);
    list_tot.push_back(eMat);
    list_tot.push_back(eMat_inv);
  }
  return list_tot;
}

/*
  Apply a number of exhaustive tricks in order to reduce the size of the vector
 */
template<typename T>
MyVector<T> ExhaustiveVectorSimplificationKernel(MyVector<T> const& V, std::vector<MyMatrix<T>> const& list_mat) {
  int n = V.size();
  auto f_norm=[&](MyVector<T> const& v) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      norm += T_abs(v(i));
    }
    return norm;
  };
  MyVector<T> V_work = V;
  T norm_work = f_norm(V);
  while(true) {
    int n_succ = 0;
    for (auto & mat : list_mat) {
      MyVector<T> V_cand = mat.transpose() * V_work;
      T norm_cand = f_norm(V_cand);
      if (norm_cand < norm_work) {
        V_work = V_cand;
        norm_work = norm_cand;
        n_succ += 1;
      }
    }
    if (n_succ == 0) {
      return V_work;
    }
  }
}

template<typename T>
MyVector<T> ExhaustiveVectorSimplification(MyVector<T> const& V, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  return ExhaustiveVectorSimplificationKernel(V, list_mat_tot);
}


template<typename T>
std::vector<MyVector<T>> ExhaustiveVectorSimplifications(std::vector<MyVector<T>> const& list_V, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  std::vector<MyVector<T>> list_V_red;
  for (auto& eV: list_V) {
    MyVector<T> eV_red = ExhaustiveVectorSimplificationKernel(eV, list_mat_tot);
    list_V_red.push_back(eV_red);
  }
  return list_V_red;
}

/*
  Take a matrix M in argument.
  Look for ways to multiply by elements of list_mat
  so that the L1 norm of M G decreases
 */
template<typename T>
MyMatrix<T> ExhaustiveMatrixCosetSimplificationKernel(MyMatrix<T> const& M, std::vector<MyMatrix<T>> const& list_mat) {
  int n = M.rows();
  auto f_norm=[&](MyMatrix<T> const& Hin) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        norm += T_abs(Hin(i,j));
      }
    }
    return norm;
  };
  MyMatrix<T> M_work = M;
  T norm_work = f_norm(M);
  while(true) {
    int n_succ = 0;
    for (auto & mat : list_mat) {
      MyMatrix<T> M_cand = M_work * mat;
      T norm_cand = f_norm(M_cand);
      if (norm_cand < norm_work) {
        M_work = M_cand;
        norm_work = norm_cand;
        n_succ += 1;
      }
    }
    if (n_succ == 0) {
      return M_work;
    }
  }
}


template<typename T>
std::vector<MyMatrix<T>> ExhaustiveMatrixCosetSimplifications(std::vector<MyMatrix<T>> const& list_cos, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  std::vector<MyMatrix<T>> list_cos_red;
  for (auto& eCos: list_cos) {
    MyMatrix<T> eCos_red = ExhaustiveMatrixCosetSimplificationKernel(eCos, list_mat_tot);
    list_cos_red.push_back(eCos_red);
  }
  return list_cos_red;
}





/*
  The result of the simplification algorithm for double cosets
 */
template<typename T>
struct DoubleCosetSimplification {
  MyMatrix<T> u_red; // The u reduction element
  MyMatrix<T> d_cos_red; // The reduced coset
  MyMatrix<T> v_red; // The v_reduction element
};


// The double coset is U x V
template<typename T>
DoubleCosetSimplification<T> ExhaustiveMatrixDoubleCosetSimplifications(MyMatrix<T> const& d_cos, std::vector<MyMatrix<T>> const& u_gens, std::vector<MyMatrix<T>> const& v_gens, size_t const& max_iter, [[maybe_unused]] std::ostream& os) {
  std::vector<MyMatrix<T>> u_gens_tot = Exhaust_get_total_generators(u_gens);
  std::vector<MyMatrix<T>> v_gens_tot = Exhaust_get_total_generators(v_gens);
  int n_gens_u = u_gens_tot.size();
  int n_gens_v = v_gens_tot.size();
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
  os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications |u_gens_tot|=" << u_gens_tot.size() << " |v_gens_tot|=" << v_gens_tot.size() << "\n";
#endif
  int n = d_cos.rows();
  auto f_norm=[&](MyMatrix<T> const& Hin) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        norm += T_abs(Hin(i,j));
      }
    }
    return norm;
  };
  MyMatrix<T> u_red = IdentityMat<T>(n);
  MyMatrix<T> v_red = IdentityMat<T>(n);
  MyMatrix<T> d_cos_work = d_cos;
  T norm_work = f_norm(d_cos);
  std::vector<int> indices_u;
  for (int i=0; i<n_gens_u; i++) {
    indices_u.push_back(i);
  }
  std::vector<int> indices_v;
  for (int i=0; i<n_gens_v; i++) {
    indices_v.push_back(i);
  }


  auto f_search_uv=[&]() -> bool {
    for (int u=0; u<n_gens_u; u++) {
      for (int v=0; v<n_gens_v; v++) {
        int u2 = indices_u[u];
        int v2 = indices_v[v];
        MyMatrix<T> const& u_gen = u_gens_tot[u2];
        MyMatrix<T> const& v_gen = v_gens_tot[v2];
        MyMatrix<T> d_cos_cand = u_gen * d_cos_work * v_gen;
        T norm_cand = f_norm(d_cos_cand);
        if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
          os << "DCOS_SIMP: Improving with u2=" << u2 << " v2=" << v2 << " XXX u=" << u << " v=" << v << "\n";
#endif
          d_cos_work = d_cos_cand;
          norm_work = norm_cand;
          u_red = u_gen * u_red;
          v_red = v_red * v_gen;
          return true;
        }
      }
    }
    return false;
  };
  auto f_search_u=[&]() -> bool {
    for (int u=0; u<n_gens_u; u++) {
      int u2 = indices_u[u];
      MyMatrix<T> const& u_gen = u_gens_tot[u2];
      MyMatrix<T> d_cos_cand = u_gen * d_cos_work;
      T norm_cand = f_norm(d_cos_cand);
      if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
        os << "DCOS_SIMP: Improving with u2=" << u2 << "\n";
#endif
        d_cos_work = d_cos_cand;
        norm_work = norm_cand;
        u_red = u_gen * u_red;
        return true;
      }
    }
    return false;
  };
  auto f_search_v=[&]() -> bool {
    for (int v=0; v<n_gens_v; v++) {
      int v2 = indices_v[v];
      MyMatrix<T> const& v_gen = v_gens_tot[v2];
      MyMatrix<T> d_cos_cand = d_cos_work * v_gen;
      T norm_cand = f_norm(d_cos_cand);
      if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
        os << "DCOS_SIMP: Improving with v2=" << v2 << "\n";
#endif
        d_cos_work = d_cos_cand;
        norm_work = norm_cand;
        v_red = v_red * v_gen;
        return true;
      }
    }
    return false;
  };
  auto f_search=[&]() -> bool {
    int chosen_method = 3;
    f_random_transpose(indices_u);
    f_random_transpose(indices_v);
    if (chosen_method == 1) {
      return f_search_uv();
    }
    if (chosen_method == 2) {
      bool test_u = f_search_u();
      if (test_u) {
        return true;
      }
      return f_search_v();
    }
    if (chosen_method == 3) {
      bool test_u = f_search_u();
      if (test_u) {
        return true;
      }
      bool test_v = f_search_v();
      if (test_v) {
        return true;
      }
      return f_search_uv();
    }
    return false;
  };
  size_t n_iter = 0;
  while(true) {
    bool test = f_search();
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
    os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_iter=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
    if (!test) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
      os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_final_iter(A)=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
      return {u_red, d_cos_work, v_red};
    }
    if (n_iter == max_iter) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
      os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_final_iter(B)=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
      return {u_red, d_cos_work, v_red};
    }
    n_iter += 1;
  }
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
