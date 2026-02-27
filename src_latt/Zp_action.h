// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ZP_ACTION_H_
#define SRC_LATT_ZP_ACTION_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "MAT_Matrix.h"
// clang-format on

template<typename T>
size_t vector_to_size_t(int dim, MyVector<T> const& v, size_t const& mod_val_s) {
  size_t pos = 0;
  size_t eProd = 1;
  for (int i = 0; i < dim; i++) {
    pos += eProd * v(i);
    eProd *= mod_val_s;
  }
  return pos;
}

template<typename T>
void size_t_to_vector(int dim, size_t const& x, size_t const& mod_val_s, MyVector<T> & v_out) {
  size_t v_work = x;
  for (int u=0; u<dim; u++) {
    size_t res = v_work % mod_val_s;
    v_out[u] = res;
    v_work = (v_work - res) / mod_val_s;
  }
}

struct ResultModEnumeration {
  std::vector<size_t> vect_orbits;
  std::vector<size_t> orbit_sizes;
};


template<typename T>
ResultModEnumeration get_partition_section1(int dim, std::vector<MyMatrix<T>> const& l_gens, T const& mod_val) {
  size_t mod_val_s = UniversalScalarConversion<int32_t,T>(mod_val);
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
    size_t_to_vector(dim, x, mod_val_s, v1);
    v2 = gen.transpose() * v1;
    for (int u=0; u<dim; u++) {
      T val = ResInt(v2[u], mod_val);
      v2[u] = val;
    }
    return vector_to_size_t(dim, v2, mod_val_s);
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

template<typename T, typename Twork>
ResultModEnumeration get_partition_section2(int dim, std::vector<MyMatrix<T>> const& l_gens, T const& mod_val) {
  std::vector<MyMatrix<Twork>> l_gens_b;
  for (auto & gen: l_gens) {
    MyMatrix<Twork> M(dim, dim);
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
        T val1 = gen(i,j);
        T val2 = ResInt(val1, mod_val);
        Twork val3 = UniversalScalarConversion<Twork,T>(val2);
        M(i,j) = val3;
      }
    }
    l_gens_b.push_back(M);
  }
  Twork mod_val_b = UniversalScalarConversion<Twork,T>(mod_val);
  return get_partition_section1(dim, l_gens_b, mod_val_b);
}



template<typename T>
ResultModEnumeration get_partition(int dim, std::vector<MyMatrix<T>> const& l_gens, T const& mod_val) {
  T max_coeff = mod_val * mod_val * dim;
  int8_t max_val_i8 = std::numeric_limits<int8_t>::max();
  int16_t max_val_i16 = std::numeric_limits<int16_t>::max();
  int32_t max_val_i32 = std::numeric_limits<int32_t>::max();
  int64_t max_val_i64 = std::numeric_limits<int64_t>::max();
  T max_val_i8_T = UniversalScalarConversion<T,int8_t>(max_val_i8);
  T max_val_i16_T = UniversalScalarConversion<T,int16_t>(max_val_i16);
  T max_val_i32_T = UniversalScalarConversion<T,int32_t>(max_val_i32);
  T max_val_i64_T = UniversalScalarConversion<T,int64_t>(max_val_i64);
  if (max_coeff < max_val_i8_T) {
    return get_partition_section2<T,int8_t>(dim, l_gens, mod_val);
  }
  if (max_coeff < max_val_i16_T) {
    return get_partition_section2<T,int16_t>(dim, l_gens, mod_val);
  }
  if (max_coeff < max_val_i32_T) {
    return get_partition_section2<T,int32_t>(dim, l_gens, mod_val);
  }
  if (max_coeff < max_val_i64_T) {
    return get_partition_section2<T,int64_t>(dim, l_gens, mod_val);
  }
  std::cerr << "Failed to find a matching entry. Very unlikely to happen. Probably a bug\n";
  throw TerminalException{1};
}




// clang-format off
#endif  // SRC_LATT_ZP_ACTION_H_
// clang-format on
