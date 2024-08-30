// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHORTESTUNIVERSAL_H_
#define SRC_LATT_SHORTESTUNIVERSAL_H_

// clang-format off
#include "Shvec_exact.h"
#include <string>
// clang-format on

template <typename T, typename Tint>
resultCVP<T, Tint>
CVPVallentinProgram_choice(MyMatrix<T> const &GramMat, MyVector<T> const &eV,
                           std::string const &NameMeth, std::ostream &os) {
  //
  if (NameMeth == "SVexact") {
    return CVPVallentinProgram_exact<T, Tint>(GramMat, eV, os);
  }
  //
  std::cerr << "No matching method found\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
resultCVP<T, Tint> CVPVallentinProgram(MyMatrix<T> const &GramMat,
                                       MyVector<T> const &eV,
                                       std::ostream &os) {
  return CVPVallentinProgram_exact<T, Tint>(GramMat, eV, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector(MyMatrix<T> const &GramMat, T const &MaxNorm,
                             std::ostream &os) {
  return T_ShortVector_exact<T, Tint>(GramMat, MaxNorm, os);
}

template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVector(MyMatrix<T> const &eMat, std::ostream &os) {
  T MinNorm = MinimumDiagonal(eMat);
  MyMatrix<Tint> TheSHVall = T_ShortVector<T, Tint>(eMat, MinNorm, os);
  return SelectShortestVector(eMat, TheSHVall);
}

// clang-format off
#endif  // SRC_LATT_SHORTESTUNIVERSAL_H_
// clang-format on
