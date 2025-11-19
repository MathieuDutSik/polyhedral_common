// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ENUM_SCHEMES_POLY_CELLCOMPLEX_H_
#define SRC_ENUM_SCHEMES_POLY_CELLCOMPLEX_H_

// clang-format off
#include "boost_serialization.h"
#include "basic_datafile.h"
#include "Timings.h"
#include <unordered_map>
#include <utility>
#include <vector>
#include <set>
#include <functional>
#include <string>
// clang-format on

const size_t seed_partition = 10;
const size_t seed_hashmap = 20;


template<typename Tequiv>
struct BoundSerial {
  size_t i_orb;
  Tequiv eq;
};


/*
  We want to implement a fully remplatized system for
  computing the

  Context of usage:
  ---Computing the well-rounded cell-complex for perfect forms.
     Usage: Cohomology computation.
     Gen: Gets the facets of the polyhedral cone.
     Equiv: Take the isobarycenter of the extreme rays.
  ---Computing the full complex for the perfect forms, not just
     Usage: Hecke operators.
     Gen: Take the facets as before.
     Equiv: Use the full dimensional cells for tracking the lower
            dimensional cells.
  ---Computing the cells for Sp2(Z) and others from Dan Yazaki.
     Usage: Hecke operators.
     Gen: Seems to be hardcoded.
     Equiv: Use the full dimensional cells.
  ---Computing the cell-complex of the Delaunay cells.
     Usage: Quantization stuff.
     Gen: Gets the facets of the polytope.

  Considerations:
  ---Using the full-dimensional cells for tracking forces having
     all the stuff in memory. Price to pay, it seems.
     But it does not seem that this impacts the generic code.
  ---The classic usage is for cohomology.
     So, the result of the generation becomes
     boundary: std::vector<Combined<Tequiv>>.
     But we do not always need that.
  ---But sometimes, we just need the list of orbits.

  Design choices:
  ---The computation is done from one level to the next one.
     No doing all the computation in one single go.
  ---The type Tobj is the same on input and output.
  ---We use the database access similar to the one
     for
  ---We use f_spann similarly (which opens the possibility of
     reduction of complexity)
  ---

  Types:
  ---Type Tobj, which is the same on input and output.
  ---Type Tequiv representing the equivalence of objects.
  ---Lambda f_n_treated returns the number of entries
     that have already been treated
     auto f_n_treated() -> size_t;
  ---Lambda f_next_input which gets us the data on input.
     auto f_next_input() -> std::optional<Tobj>;
     When returning None, there is nothing left to treat.
  ---Lambda f_next_output which gets us the data on output.
     auto f_next_output() -> std::optional<Tobj>;
  ---Lambda f_generate which generates the subcells
     auto f_generate(Tobj const& x) -> std::vector<Tobj>;
  ---Lambda f_hash works the same way as the generation stuff.
     auto f_hash(size_t seed, Tobj const& x) -> size_t;
  ---Lambda f_repr for testing equivalence of object.
     auto f_repr(Tobj const& x,Tobj const& y) -> std::optional<Tequiv>
  ---Lambda f_insert_obj for inserting a new object.
     auto f_insert_obj(Tobj const& x) -> bool;
     If returning true then early termination of the enumeration
     is done.
  ---Lambda f_idx_obj is returning the object on position pos from storage.
     auto f_idx_obj(size_t const& pos) -> Tobj;
  ---Lambda f_spann for getting the object
     auto f_spann(Tobj const& x) -> std::pair<Tobj,Tequiv>;
  ---Lambda f_set_boundary_status for setting the boundary and status
     auto f_set_boundary_status(size_t const&idx, std::vector<BoundEquiv<Tequiv>> const& bound) -> void;






 */
template<typename Tobj, typename Tequiv, typename Fn_treated, typename Fnext_input, typename Fnext_output,
  typename Fgenerate, typename Fhash, typename Frepr, typename Finsert_obj, typename Fidx_obj,
  typename Fspann, typename Fset_boundary_status>
bool compute_next_level_serial(int const &max_time_second, bool const& compute_boundary,
                               Fn_treated f_n_treated,
                               Fnext_input f_next_input, Fnext_output f_next_output,
                               Fgenerate f_generate, Fhash f_hash, Frepr f_repr,
                               Finsert_obj f_insert_obj, Fidx_obj f_idx_obj,
                               Fspann f_spann, Fset_boundary_status f_set_boundary_status) {
  SingletonTime start;
  size_t n_obj = 0;
  std::unordered_map<size_t, std::vector<size_t>> indices_by_hash;
  //
  // Loading what as already been done.
  //
  auto insert_load=[&](Tobj const& x) -> void {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    std::vector<size_t> &vect = indices_by_hash[hash_hashmap];
    vect.push_back(n_obj);
    n_obj++;
  };
  auto load_output=[&]() -> void {
    while(true) {
      std::optional<Tobj> opt = f_next_output();
      if (*opt) {
        Tobj const& x = *opt;
        insert_load(x);
      } else {
        break;
      }
    }
  };
  load_output();
  //
  // Accessing entries
  //
  auto f_get_equiv=[&](Tobj const& x) -> BoundSerial<Tequiv> {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    std::vector<size_t> &vect = indices_by_hash[hash_hashmap];
    for (auto &idx : vect) {
      Tobj y = f_idx_obj(idx);
      std::optional<Tequiv> opt = f_repr(y, x);
      if (opt) {
        Tequiv const& eq = *opt;
        BoundSerial<Tequiv> bnd{idx, eq};
        return bnd;
      }
    }
    std::pair<Tobj, Tequiv> pair = f_spann(x);
    f_insert_obj(pair.first);
    BoundSerial<Tequiv> bnd{n_obj, pair.second};
    n_obj += 1;
    return bnd;
  };
  //
  // Treat entries
  //
  auto treat_entry=[&](size_t const& idx, Tobj const& x) -> void {
    std::vector<Tobj> l_elt = f_generate(x);
    std::vector<BoundSerial<Tequiv>> l_bound;
    for (auto & y: l_elt) {
      BoundSerial<Tequiv> bnd = f_get_equiv(y);
      if (compute_boundary) {
        l_bound.push(bnd);
      }
    }
    f_set_boundary_status(idx, l_bound);
  };
  size_t n_treated = f_n_treated();
  while(true) {
    std::optional<Tobj> opt = f_next_input();
    if (*opt) {
      Tobj const& obj = *opt;
      treat_entry(n_treated, obj);
    } else {
      break;
    }
    if (max_time_second > 0 && si(start) > max_time_second) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "CELL_SCH: returning false due to si(start) > max_time_second > 0\n";
#endif
      return false;
    }
  }
  return true;
}


template<typename Tobj, typename Tequiv>
struct ResultNextLevel {
  std::vector<Tobj> l_obj;
  std::vector<std::vector<BoundSerial<Tequiv>> l_bound;
};


template <typename Tdata, typename Fincorrect>
ResultNextLevel<typename Tdata::Tobj, typename Tdata::Tequiv> EnumerateAndStoreLevel_Serial(Tdata &data,
                                                                                            std::vector<Tobj> const& l_input,
                                                                                            Fincorrect f_incorrect,
                                                                                            int const &max_runtime_second,
                                                                                            bool const& compute_boundary) {
  using Tobj = typename Tdata::Tobj;
  using Tequiv = typename Tdata::Tequiv;
  using TadjO = typename Tdata::TadjO;
  using TadjO_work = AdjO_Serial<TadjO>;
  std::ostream &os = data.get_os();
  auto f_hash = [&](size_t const &seed, Tobj const &x) -> size_t {
    return data.f_hash(seed, x);
  };
  auto f_repr = [&](Tobj const &x, Tobj const &y) -> std::optional<Tequiv> {
    return data.f_repr(x, y);
  };
  auto f_spann = [&](Tobj const &x) -> std::pair<Tobj, Tequiv> {
    return data.f_spann(x);
  };
  std::vector<Tobj> l_obj;
  std::vector<std::vector<BoundSerial<Tequiv>> ll_bound;
  auto f_adj = [&](int const &i_orb) -> std::vector<TadjI> {
    Tobj &x = l_obj[i_orb].x;
    return data.f_adj(x);
  };
  auto f_idx_obj = [&](size_t const &idx) -> Tobj {
    return l_obj[idx];
  };
  auto f_n_treated = [&]() -> size_t {
    return 0;
  };
  size_t n_input = l_input.size();
  size_t pos = 0;
  auto f_next_input = [&]() -> std::optional<Tobj> {
    if (pos < n_input) {
      pos += 1;
      return l_input[pos - 1];
    } else {
      return {};
    }
  };
  auto f_next_output = [&]() -> std::optional<Tobj> {
    return {};
  };
  auto f_insert_obj = [&](Tobj const &x) -> bool {
    l_obj.push_back(x);
#ifdef DEBUG_CELL_SCHEME
    os << "ADJ_SCH: EnumerateAndStore_Serial: |l_obj|=" << l_obj.size() << "\n";
#endif
    bool test = f_incorrect(x);
    return test;
  };
  auto f_save_status = [&](size_t const &pos, bool const &val) -> void {
    uint8_t val_i = static_cast<uint8_t>(val);
    if (l_status.size() <= pos) {
      l_status.push_back(val_i);
    } else {
      l_status[pos] = val_i;
    }
  };
  auto f_set_boundary_status = [&](int const &i_orb, std::vector<BoundSerial<Tequiv>> const &l_bound) -> void {
    ll_bound.push_back(l_bound);
  };
  (void)compute_next_level_serial<Tobj, Tequiv, decltype(f_n_treated), decltype(f_next_input),
    decltype(f_next_output), decltype(f_generate), decltype(f_hash), decltype(f_repr),
    decltype(f_insert_obj), decltype(f_idx_obj), decltype(f_spann), decltype(f_set_boundary_status)>
    (max_time_second, compute_boundary,
     f_n_treated,
     f_next_input, f_next_output,
     f_generate, f_hash, f_repr,
     f_insert_obj, f_idx_obj,
     f_spann, f_set_boundary_status);
  return l_obj;
}


// clang-format off
#endif  // SRC_ENUM_SCHEMES_POLY_CELLCOMPLEX_H_
// clang-format on
