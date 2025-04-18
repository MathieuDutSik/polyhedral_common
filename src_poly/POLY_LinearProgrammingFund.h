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
  bool PrimalDefined = false;
  bool DualDefined = false;
  T OptimalValue;
  MyVector<T> DualSolution;
  MyVector<T> DirectSolution;
};

template <typename T>
void PrintLpSolution(LpSolution<T> const &eSol, std::ostream &os) {
  os << "PrimalDefined=" << eSol.PrimalDefined << "\n";
  os << "DualDefined=" << eSol.DualDefined << "\n";
  os << "DualSolution=" << StringVector(eSol.DualSolution) << "\n";
  os << "OptimalValue=" << eSol.OptimalValue << "\n";
  os << "DirectSolution=" << StringVector(eSol.DirectSolution) << "\n";
}

template <typename T>
Face ComputeFaceLpSolution(MyMatrix<T> const &EXT, LpSolution<T> const &eSol) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  Face eFace(nbRow);
#ifdef SANITY_CHECK_LINEAR_PROGRAMMING_FUND
  if (!eSol.DualDefined || !eSol.PrimalDefined) {
    std::cerr << "We should have DualDefined and PrimalDefined for the "
                 "computation to make sense\n";
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
    T eSum = EXT(iRow, 0);
    for (int iCol = 0; iCol < nbCol - 1; iCol++) {
      eSum += eSol.DirectSolution(iCol) * EXT(iRow, iCol + 1);
    }
#ifdef SANITY_CHECK_LINEAR_PROGRAMMING_FUND
    if (eSum < 0) {
      std::cerr << "CDD_LinearProgramming Error iRow=" << iRow
                << " eSum=" << eSum << "\n";
      std::cerr << "DualDefined=" << eSol.DualDefined
                << " PrimalDefined=" << eSol.PrimalDefined << "\n";
      std::cerr << "DirectSolution =" << StringVector(eSol.DirectSolution)
                << "\n";
      std::cerr << "EXT=\n";
      WriteMatrix(std::cerr, EXT);
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

// If eVect = C + sum_i alpha_i v_i
// if v_i.x >= 0 and alpha_i >= 0 then
template <typename T>
bool CheckDualSolutionGetOptimal(MyMatrix<T> const &EXT,
                                 MyVector<T> const &eVect,
                                 LpSolution<T> const &eSol,
                                 [[maybe_unused]] std::ostream& os) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyVector<T> V(nbCol - 1);
  T objDual = eVect(0);
  for (int iCol = 0; iCol < nbCol - 1; iCol++) {
    V(iCol) = eVect(iCol + 1);
  }
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T scal = eSol.DualSolution(iRow);
    for (int iCol = 0; iCol < nbCol - 1; iCol++) {
      V(iCol) += scal * EXT(iRow, iCol + 1);
    }
    objDual += scal * EXT(iRow, 0);
  }
  for (int iCol = 0; iCol < nbCol - 1; iCol++) {
    if (V(iCol) != 0) {
#ifdef DEBUG_LINEAR_PROGRAMMING_FUND
      os << "LPFUND: Returning false as the vector V is not zero at coordinate iCol=" << iCol << "\n";
#endif
      return false;
    }
  }
  if (objDual != eSol.OptimalValue) {
#ifdef DEBUG_LINEAR_PROGRAMMING_FUND
    os << "LPFUND: Returning false as objDual=" << objDual << " eSol.OptimnalValue=" << eSol.OptimalValue << "\n";
#endif
    return false;
  }
#ifdef DEBUG_LINEAR_PROGRAMMING_FUND
  os << "LPFUND: Returning true in CheckDualSolutionGetOptimal\n";
#endif
  return true;
}

template <typename T>
MyVector<T> GetDirectSolutionExt(LpSolution<T> const &eSol) {
  int siz = eSol.DirectSolution.size();
  MyVector<T> V(siz + 1);
  V(0) = 1;
  for (int i = 0; i < siz; i++) {
    V(i + 1) = eSol.DirectSolution(i);
  }
  return V;
}

// clang-format off
#endif  // SRC_POLY_POLY_LINEARPROGRAMMINGFUND_H_
// clang-format on
