// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_APPROXDUALDESC_H_
#define SRC_RANKIN_APPROXDUALDESC_H_

// clang-format off
#include <unordered_set>
#include <vector>
// clang-format on

/*
  We want to deal with the dual descriptions computed over the
  double number and we want to reduce them.
  ---
  Several phenomena can occur:
  * An artificial facet shows up that does not actually correspond to a real
  facet.
  * A facet that actually merged into a bigger one.
  The first phenomenon is detected by the tolThr and the second by the tolFacet.

  */
template <typename T>
vectface ReduceVectfaceApproximate(MyMatrix<T> const &EXT, vectface const &vf,
                                   T const &tolThr, T const &tolFacet) {
  int n_row = EXT.rows();
  int dim = EXT.cols();
  std::unordered_set<Face> set;
  for (auto &eFace : vf) {
    std::vector<size_t> eList = FaceToVector<size_t>(eFace);
    SelectionRowCol<T> eSelect =
        TMat_SelectRowColMaxPivot_subset(EXT, eList, tolThr);
    if (eSelect == dim - 1) {
      Face FullFace(n_row);
      for (int i_row = 0; i_row < n_row; i_row++) {
        T eSum = 0;
        for (int i = 0; i < dim; i++) {
          eSum += EXT(i_row, i) * eSelect.NSP(0, i);
        }
        T absSum = T_abs(eSum);
        if (absSum < tolFacet) {
          FullFacet[i_row] = 1;
        }
      }
      set.insert(FullFace);
    }
  }
  vectface vf_ret(m_row);
  for (auto &eFace : set) {
    vf_ret.push_back(eFace);
  }
  return vf_ret;
}

// clang-format off
#endif  // SRC_RANKIN_APPROXDUALDESC_H_
// clang-format on
