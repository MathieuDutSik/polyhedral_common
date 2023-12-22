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

/*
  Code for the L-type domains.

  Two main use case:
  ---Lattice case: Then the Tvert is actually a Tint and can be mpz_class, int32_t, etc.
  ---Periodic structure case: Then the coordinates are no longer integral.
    Also the equivalence are no longer integral. Sure the matrix transformation is
    integral, but the translation vector is not necessarily so.
 */

template<typename Tvert>
struct SingleEquiv {
  Face incd;
  MyMatrix<Tvert> P;
  int i_orb;
};

template<typename Tvert, typename Tgroup>
struct SingleDelaunay {
  MyMatrix<Tvert> EXT;
  Tgroup GRPlatt;
  std::vector<SingleEquiv<Tvert>> l_adj;
};

template<typename Tvert, typename Tgroup>
struct DelaunayTesselation {
  std::vector<SingleDelaunay<Tvert,Tgroup>> l_dels;
};







// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_H_
// clang-format on
