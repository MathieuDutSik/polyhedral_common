// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_H_

// clang-format off
#include "POLY_LinearProgramming.h"
#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include <string>
#include <vector>
// clang-format on


template<typename Tint>
struct SingleEquiv {
  Face incd;
  MyMatrix<Tint> P;
  int i_orb;
};

template<typename Tint, typename Tgroup>
struct SingleDelaunay {
  MyMatrix<Tint> EXT;
  Tgroup GRPlatt;
  std::vector<SingleEquiv<Tint>> l_adj;
};

template<typename Tint, typename Tgroup>
struct DelaunayTesselation {
  std::vector<SingleDelaunay<Tint,Tgroup>> l_dels;
};







// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on
