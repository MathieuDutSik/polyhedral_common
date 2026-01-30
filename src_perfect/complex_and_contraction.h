// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_
#define SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_

template<typename Tint>
struct BoundEntry {
  int iOrb;
  int sign;
  MyMatrix<Tint> M;
};

template<typename Tint>
struct ListBoundEntry {
  std::vector<BoundEntry<Tint>> l_bound;
};

template<typename Tint>
struct FullBoundary {
  std::vector<ListBoundEntry<Tint>> ll_bound;
};



// clang-format off
#endif  // SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_
// clang-format on

