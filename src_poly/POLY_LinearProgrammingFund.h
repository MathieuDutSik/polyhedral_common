// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_
#define SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
// clang-format on

template <typename T>
struct LpSolution {
  std::string method;
  bool PrimalDefined = false;
  bool DualDefined = false;
  //
  MyVector<T> DualSolution;
  //
  T OptimalValue;
  //
  MyVector<T> DirectSolution;
  MyVector<T> DirectSolutionExt;
  Face eFace;
  int rankDirectSol = -1;
  std::string Answer;
};

template <typename T>
void PrintLpSolution(LpSolution<T> const &eSol, std::ostream &os) {
  os << "method=" << eSol.method << "\n";
  os << "PrimalDefined=" << eSol.PrimalDefined << "\n";
  os << "DualDefined=" << eSol.DualDefined << "\n";
  os << "DualSolution=" << StringVector(eSol.DualSolution) << "\n";
  os << "OptimalValue=" << eSol.OptimalValue << "\n";
  os << "DirectSolution=" << StringVector(eSol.DirectSolution) << "\n";
  os << "DirectSolutionExt=" << StringVector(eSol.DirectSolutionExt) << "\n";
  os << "rankDirectSol=" << eSol.rankDirectSol << "\n";
  os << "Answer=" << eSol.Answer << "\n";
}


// clang-format off
#endif  // SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_
// clang-format on
