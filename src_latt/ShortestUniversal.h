// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHORTESTUNIVERSAL_H_
#define SRC_LATT_SHORTESTUNIVERSAL_H_

// clang-format off
#include "Shvec_exact.h"
#ifdef USE_LIBSHORT
#include "Shvec_double.h"
#endif
#include <string>
// clang-format on

template <typename T, typename Tint>
resultCVP<T, Tint>
CVPVallentinProgram_choice(MyMatrix<T> const &GramMat, MyVector<T> const &eV,
                           std::string const &NameMeth, std::ostream &os) {
  /*
  bool DoCheck = false;
  if (DoCheck) {
    resultCVP<T, Tint> res1 = CVPVallentinProgram_exact<T, Tint>(GramMat, eV,
  os); resultCVP<T, Tint> res2 = CVP_N23_24A1<T, Tint>(eV); if (res1 != res2) {
      std::cerr << "res1.TheNorm=" << res1.TheNorm << "\n";
      std::cerr << "res2.TheNorm=" << res2.TheNorm << "\n";
      std::cerr << "|res1.ListVect|=" << res1.ListVect.rows()
                << " |res2.ListVect|=" << res2.ListVect.rows() << "\n";
      std::cerr << "res1.ListVect=\n";
      WriteMatrixGAP(std::cerr, res1.ListVect);
      std::cerr << "res2.ListVect=\n";
      WriteMatrixGAP(std::cerr, res2.ListVect);
      std::cerr << "Clear error in the code\n";
      throw TerminalException{1};
    }
    return res1;
  }
  */
  //
  if (NameMeth == "SVexact")
    return CVPVallentinProgram_exact<T, Tint>(GramMat, eV, os);
    //
#ifdef USE_LIBSHORT
  if (NameMeth == "SVdouble")
    return CVPVallentinProgram_double<T, Tint>(GramMat, eV);
#endif
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
