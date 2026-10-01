// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_ILP_FUNDAMENTAL_H_
#define SRC_MILP_ILP_FUNDAMENTAL_H_

// clang-format off
#include "MAT_Matrix.h"
// clang-format on

/*
  The types shared by the integer programming methods, see
  integer_linear_programming.h for the conventions.
 */

enum class IlpStatus { Optimal, Infeasible, Unbounded };

template <typename T> struct IlpSolution {
  IlpStatus status;
  MyVector<T> solution; // Only meaningful when status == Optimal
  T OptimalValue;       // Only meaningful when status == Optimal
};

template <typename T>
T ILP_EvaluateRow(MyVector<T> const &row, MyVector<T> const &x) {
  T val = row(0);
  int d = x.size();
  for (int j = 0; j < d; j++) {
    val += row(j + 1) * x(j);
  }
  return val;
}

template <typename T>
bool ILP_IsFeasible(MyMatrix<T> const &ListIneq, MyMatrix<T> const &ListEqua,
                    MyVector<T> const &x) {
  for (int i = 0; i < ListIneq.rows(); i++) {
    MyVector<T> row = GetMatrixRow(ListIneq, i);
    if (ILP_EvaluateRow(row, x) < 0) {
      return false;
    }
  }
  for (int i = 0; i < ListEqua.rows(); i++) {
    MyVector<T> row = GetMatrixRow(ListEqua, i);
    if (ILP_EvaluateRow(row, x) != 0) {
      return false;
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_MILP_ILP_FUNDAMENTAL_H_
// clang-format on
