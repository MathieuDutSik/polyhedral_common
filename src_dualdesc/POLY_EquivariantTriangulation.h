// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
#define SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_

// clang-format off
#include "POLY_RecursiveDualDesc.h"
// clang-format on


template<typename T>
struct EquivariantTriangulation {
  MyMatrix<T> RelVectors;
  vectface vf_trig;
};

template<typename T>
struct ExtendibleEquivariantTriangulation {
  std::vector<MyVector<T>> ListVert;
  std::unordered_map<MyVector<T>, size_t> MapVert;
  std::vector<std::vector<size_t>> ListTrig;
};














// clang-format off
#endif  // SRC_DUALDESC_POLY_EQUIVARIANTTRIANGULATION_H_
// clang-format on
