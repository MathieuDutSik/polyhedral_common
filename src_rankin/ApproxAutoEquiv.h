// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_APPROXAUTOEQUIV_H_
#define SRC_RANKIN_APPROXAUTOEQUIV_H_

#include "Shvec_exact.h"


template<typename T, typename Tint>
WeightMatrix<true, T, uint32_t> GetWeightMatrix(MyMatrix<T> const& eG, MyMatrix<Tint> const& SHV, T const& tol) {
  int n = eG.rows();
  int n_row = SHV.rows();
  // We have a linear algorithm, but we could use a std::set<T> to get a linear algorithm.
  std::vector<T> ListVal;
  auto get_idx=[&](T const& val) -> size_t {
    for (size_t u=0; u<ListVal.size(); u++) {
      if (T_abs(ListVal[u] - val) < tol) {
        return u;
      }
    }
    size_t pos = ListVal.size();
    ListVal.push_back(val);
    return pos;
  };
  MyMatrix<size_t> M(n_row, n_row);
  std::vector<MyVector<T>> ListV;
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<Tint> eV_i = GetMatrixRow(SHV, i_row);
    MyVector<T> eV = UniversalVectorConversion<T,Tint>(eV_i);
    ListV.push_back(eV);
  }
  for (int i_row=0; i_row<n_row; i_row++) {
    MyVector<T> V = eG * ListV[i_row];
    for (int j_row=0; j_row<n_row; j_row++) {
      T scal = V.dot(ListV[j_row]);
      size_t pos = get_idx(scal);
      M(i_row, j_row) = pos;
    }
  }

}



// clang-format off
#endif  // SRC_RANKIN_APPROXAUTOEQUIV_H_
// clang-format on
