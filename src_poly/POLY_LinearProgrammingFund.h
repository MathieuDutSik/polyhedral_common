// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_
#define SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include <string>
// clang-format on

#ifdef DEBUG
#define DEBUG_LINEAR_PROGRAMMING_FUND
#endif

#ifdef TIMINGS
#define TIMINGS_LINEAR_PROGRAMMING_FUND
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LINEAR_PROGRAMMING_FUND
#endif

template <typename T> struct LpSolutionSimple {
  bool PrimalDefined;
  T OptimalValue;
  int nbRow;
  int nbCol;
  MyVector<T> DirectSolution;
  MyVector<T> DirectSolutionExt;
  // Value 0 for not assigned.
  // Value 1 for "B"
  // Value 2 for "NF"
  // Value 3 for "NL"
  MyVector<int> RowStatus;
  MyVector<int> ColumnStatus;
};

template <typename T> struct LpSolution {
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
  os << "Answer=" << eSol.Answer << "\n";
}

template<typename T>
Face ComputeFaceLpSolution(MyMatrix<T> const& EXT, LpSolution<T> const& eSol) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  Face eFace(nbRow);
#ifdef SANITY_CHECK_LINEAR_PROGRAMMING_FUND
  if (!eSol.DualDefined || !eSol.PrimalDefined) {
    std::cerr << "We should have DualDefined and PrimalDefined for the computation to make sense\n";
    throw TerminalException{1};
  }
#endif
  // The comparison with values makes sense only of the dual program is
  // defined. Otherwise, what we may get is actually a primal_direction and
  // it would just not make sense with negative values for eSum.
#ifdef TIMINGS_LINEAR_PROGRAMMING_FUND
  MicrosecondTime time;
#endif
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eSum(0);
    for (int iCol = 0; iCol < nbCol; iCol++) {
      eSum += eSol.DirectSolutionExt(iCol) * EXT(iRow, iCol);
    }
#ifdef SANITY_CHECK_LINEAR_PROGRAMMING_FUND
    if (eSum < 0) {
      std::cerr << "CDD_LinearProgramming Error iRow=" << iRow
                << " eSum=" << eSum << "\n";
      std::cerr << "DualDefined=" << eSol.DualDefined
                << " PrimalDefined=" << eSol.PrimalDefined << "\n";
      std::cerr << "DirectSolutionExt =";
      for (int iCol = 0; iCol < nbCol; iCol++)
        std::cerr << " " << eSol.DirectSolutionExt(iCol);
      std::cerr << "\n";
      std::cerr << "EXT=\n";
      WriteMatrix(std::cerr, EXT);
      std::cerr << "eVect=\n";
      WriteVectorNoDim(std::cerr, eVect);
      std::cerr << "Obtained vertex solution is not valid\n";
      std::cerr << "Please debug. Before calling TerminalEnding\n";
      throw TerminalException{1};
    }
#endif
    if (eSum == 0) {
      eFace[iRow] = 1;
    }
  }
#ifdef TIMINGS_LINEAR_PROGRAMMING_FUND
  std::cerr << "|eFace|=" << time << "\n";
#endif
  return eFace;
}

// clang-format off
#endif  // SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_
// clang-format on
