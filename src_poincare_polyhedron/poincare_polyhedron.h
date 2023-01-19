// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
#define SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_


#include "POLY_DirectDualDesc.h"
#include "Namelist.h"

// The initial data for the Poincare Polyhedron Theorem
// ---a point x
// ---a list of group element which is of course assumed to generate the group
template<typename T>
struct DataPoincare {
  MyVector<T> x;
  std::vector<MyMatrix<T>> ListGroupElt;
};


// The data structure for keeîng track of the group elements:
// ---Positive value (so element 0 correspond to 1, X to X+1, ...) are the elements themselves.
// ---Negative values correspond to their inverse (that is -1 correspond to inverse of generator 0)
// ---0 should never show up in the list.
// DI stands for "Direct or Inverse"
struct TrackGroup {
  std::vector<int> ListDI;
};

TrackGroup ProductTrack(TrackGroup const& tg1, TrackGroup const& tg2) {
  std::vector<int> ListDI = tg1.ListDI;
  ListDI.insert(ListDI.end(), tg2.ListDI.begin(), tg2.ListDI.end());
  return {ListDI};
}

TrackGroup Inverse(TrackGroup const& tg) {
  std::vector<int> ListDI_ret;
  size_t len = tg.ListDI.size();
  for (size_t u=0; u<len; u++) {
    size_t v = len - 1 - u;
    ListDI_ret.push_back( - tg.ListDI[v]);
  }
  return {ListDI_ret};
}

template<typename T>
struct StepEnum {
  MyMatrix<T> L
}











// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
// clang-format on
