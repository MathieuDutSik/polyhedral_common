// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ZP_ACTION_H_
#define SRC_LATT_ZP_ACTION_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Stor.h"
// clang-format on

template<typename T>
size_t vector_to_size_t(int dim, MyVector<T> const& v, T const& mod_val) {
}

template<typename T>
void size_t_to_vector(int dim, size_t const& x, T const& mod_val, MyVector<T> const& v_out) {


}


template<typename T>
size_t OnPoints(size_t const& x, ContainerMatrix<T> const& gen) {
  int dim = gen.rows();
  MyVector<T> v = size_t_to_vector(dim, x, gen.mod_val);
  MyVector<T> v_img = gen.elt.transpose() * v;
  return vector_to_size_t(dim, v_img, mod_val);
}



struct ResultModEnumeration {
  std::vector<size_t> vect_orbits;
  std::vector<size_t> orbit_sizes;
}


template<typename T>
ResultModEnumeration get_partition(int dim, std::vector<MyMatrix<T> const& l_gens, T const& mod_val) {
  size_t nbPoint = 1;
  for (int u=0; u<dim; u++) {
    nbPoint *= mod_val;
  }

  IntegerSubsetStorage Vlist = VSLT_InitializeStorage(nbPoint);
  std::vector<size_t> vect_orbits(nbPoint);
  std::vector<size_t> orbit_sizes;
  for (size_t i = 0; i < nbPoint; i++) {
    VSLT_StoreValue(Vlist, i);
  }
  IntegerSubsetStorage set_orbit = VSLT_InitializeStorage(nbPoint);
  std::vector<size_t> gen_orbit;
  gen_orbit.reserve(nbPoint);
  MyVector<T> v1(dim), v2(dim);
  auto f_act=[&](size_t const& x, MyMatrix<T> const& gen) -> size_t {
    size_t_to_vector(dim, x, mod_val, v_out);
    v2 = gen.transpose() * v1;
    return vector_to_size_t(dim, v2, mod_val);
  };
  size_t i_orbit = 0;
  auto treat_point=[&](size_t TheFirst) -> void {
    gen_orbit.clear();
    VSLT_ZeroAssignment(set_orbit);
    auto f_insert=[&](size_t pt) -> void {
      if (!VSLT_IsItInSubset(set_orbit, pt)) {
        VSLT_StoreValue(set_orbit, pt);
        gen_orbit.push_back(pt);
        vect_orbits[pt] = i_orbit;
      }
    };
    f_insert(TheFirst);
    size_t start = 0;
    while(true) {
      size_t len = gen_orbit.size();
      for (size_t u=start; u<len; u++) {
        size_t x = gen_orbit[u];
        for (auto &gen: l_gens) {
          size_t x_img = f_act(x, gen);
          f_insert(x_img);
        }
      }
      start = len;
      if (start == gen_orbit.size()) {
        break;
      }
    }
    orbit_sizes.push_back(gen_orbit.size());
  };
  while (true) {
    if (VSLT_IsEmpty(Vlist)) {
      break;
    }
    size_t TheFirst = VSLT_TheFirstPosition(Vlist);
    treat_point(TheFirst);
    i_orbit += 1;
  }
  return {std::move(vect_orbits), std::move(orbit_sizes)};
}



// clang-format off
#endif  // SRC_LATT_ZP_ACTION_H_
// clang-format on
