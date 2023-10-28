// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHVEC_EXACT_POLYTOPE_H_
#define SRC_LATT_SHVEC_EXACT_POLYTOPE_H_

// clang-format off
#include "Shvec_exact.h"
#include "POLY_LinearProgramming.h"
// clang-format on


template <typename T, typename Tint, typename Finsert>
int computeIt_polytope(const T_shvec_request<T> &request, const T &bound,
                       const MyMatrix<T> &FAC, Finsert f_insert) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n_rows = FAC.rows();
  int n_col = FAC.cols();
  if (n_col != request.dim + 1) {
    std::cerr << "Error in the size of FAC\n";
    throw TerminalException{1};
  }
  std::cerr << "Beginning of computeIt_polytope\n";
  auto f_set_bound = [&](const T &eQuot, const T &eSum,
                         [[maybe_unused]] const MyMatrix<T> &q,
                         const MyVector<Tint> &x, const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Infinitesimal_Floor<T, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<T, Tint>(eQuot, eSum);
    int len = 2 + i;
    MyMatrix<T> FACwork(n_rows, len);
    for (int i_row = 0; i_row < n_rows; i_row++) {
      for (int i_col = 0; i_col < len; i_col++)
        FACwork(i_row, i_col) = FAC(i_row, i_col);
      for (int i_col = len; i_col < n_col; i_col++) {
        FACwork(i_row, 0) += x(i_col - 1) * FAC(i_row, i_col);
      }
    }
    LpSolution<T> eSol;
    MyVector<T> Vminimize = ZeroVector<T>(len);
    //
    Vminimize(1 + i) = 1;
    eSol = CDD_LinearProgramming(FACwork, Vminimize);
    if (eSol.DualDefined && eSol.PrimalDefined) {
      // Well defined so we get a potential lower bound
      Tint eLow = UniversalCeilScalarInteger<Tint, T>(eSol.OptimalValue);
      if (eLow > lower)
        lower = eLow;
    }
    if (!eSol.DualDefined && eSol.PrimalDefined) {
      // Infinite direction. Therefore no better bound available
    }
    if (!eSol.PrimalDefined) {
      // No feasible solution. Therefore not feasible.
      // This will lead to a backtrack operation
      upper = lower - 1;
      return;
    }
    //
    Vminimize(1 + i) = -1;
    eSol = CDD_LinearProgramming(FACwork, Vminimize);
    if (eSol.DualDefined && eSol.PrimalDefined) {
      // Well defined so we get a potential upper bound
      Tint eUpp = UniversalFloorScalarInteger<Tint, T>(-eSol.OptimalValue);
      if (eUpp < upper)
        upper = eUpp;
    }
    if (!eSol.DualDefined &&
        eSol.PrimalDefined) { // Infinite direction. Therefore no bound
                              // available
    }
    if (!eSol.PrimalDefined) { // No feasible solution. Therefore not feasible.
      upper = lower - 1;       // This will lead to a backtrack operation
      return;
    }
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, bound, f_insert, f_set_bound);
}

// clang-format off
#endif  // SRC_LATT_SHVEC_EXACT_POLYTOPE_H_
// clang-format on
