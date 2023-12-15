// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_PERFECTION_H_
#define SRC_RANKIN_PERFECTION_H_

template<typename T, typename Tint>
MyMatrix<T> GetProjectorMatrix(MyMatrix<Tint> const& P, MyMatrix<T> const& eG) {
  MyMatrix<T> P_T = UniversalMatrixConversion<T,Tint>(P);
  MyMatrix<T> eG_red = P_T * eG * P_T.transpose();
  MyMatrix<T> eGinv = Inverse(eG_red);
  MyMatrix<T> ProjMat = P_T.transpose() * eGinv * P_T;
  return ProjMat;
}




// clang-format off
#endif  // SRC_RANKIN_PERFECTION_H_
// clang-format on
