// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_PARALLEL_CLASSES_TYPES_H_
#define SRC_POLY_PARALLEL_CLASSES_TYPES_H_

#include <set>
#include <string>
#include <utility>
#include <vector>

//
// Fundamental struct data types
//

template <typename Tequiv> struct EquivInfo {
  int iEntry;
  Tequiv TheEquiv;
  int nbEntryRelevant;
};

template <typename T> struct invariant_info {};

template <typename T> struct equiv_info {};

template <typename T> struct DataBank_ResultQuery {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  bool test;
  T eEnt;
  Tequiv TheEquiv;
};

template <typename T> struct FctsDataBank {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  std::function<std::optional<Tequiv>(T const &, T const &)> FctEquiv;
  std::function<int(T const &)> FctSize;
};

template <typename T> struct PairT_Tinv {
  typedef typename invariant_info<T>::invariant_type Tinv;
  T x;
  Tinv xInv;
};

// clang-format off
#endif  // SRC_POLY_PARALLEL_CLASSES_TYPES_H_
// clang-format on
