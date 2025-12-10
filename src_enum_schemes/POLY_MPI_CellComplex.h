// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ENUM_SCHEMES_POLY_MPI_CELLCOMPLEX_H_
#define SRC_ENUM_SCHEMES_POLY_MPI_CELLCOMPLEX_H_

// clang-format off
#include "boost_serialization.h"
#include "basic_datafile.h"
#include "POLY_CellComplex.h"
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

template <typename Tobj, typename Tequiv, typename Fn_treated,
          typename Fnext_input, typename Fnext_output, typename Fgenerate,
          typename Fhash, typename Frepr, typename Finsert_obj,
          typename Fidx_obj, typename Fspann, typename Fset_boundary_status>
bool compute_next_level_mpi(int const &max_time_second,
                            bool const &compute_boundary,
                            Fn_treated f_n_treated, Fnext_input f_next_input,
                            Fnext_output f_next_output, Fgenerate f_generate,
                            Fhash f_hash, Frepr f_repr,
                            Finsert_obj f_insert_obj, Fidx_obj f_idx_obj,
                            Fspann f_spann,
                            Fset_boundary_status f_set_boundary_status) {
  SingletonTime start;
  size_t n_obj = 0;
  std::unordered_map<size_t, std::vector<size_t>> indices_by_hash;
  //
  // Loading what as already been done.
  //
  auto insert_load = [&](Tobj const &x) -> void {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    std::vector<size_t> &vect = indices_by_hash[hash_hashmap];
    vect.push_back(n_obj);
    n_obj++;
  };
  auto load_output = [&]() -> void {
    while (true) {
      std::optional<Tobj> opt = f_next_output();
      if (*opt) {
        Tobj const &x = *opt;
        insert_load(x);
      } else {
        break;
      }
    }
  };
  load_output();
}

// clang-format off
#endif  // SRC_ENUM_SCHEMES_POLY_MPI_CELLCOMPLEX_H_
// clang-format on
