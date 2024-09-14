// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LAMINATIONS_H_
#define SRC_POLY_POLY_LAMINATIONS_H_

// clang-format off
#include "COMB_Combinatorics.h"
#include "MAT_Matrix.h"
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LAMINATIONS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LAMINATIONS
#endif

template <typename T, typename F>
void compute_two_laminations_f(MyMatrix<T> const &M, F f) {
#ifdef SANITY_CHECK_LAMINATIONS
  if (RankMat(M) != M.cols()) {
    std::cerr << "M should be full dimensional\n";
    throw TerminalException{1};
  }
#endif
  int dim = M.cols();
  int nbRow = M.rows();
  SelectionRowCol<T> src = TMat_SelectRowCol(M);
#ifdef DEBUG_LAMINATIONS
  std::cerr << "LAM: We have src\n";
#endif
  MyMatrix<T> Mred(dim, dim);
  for (int iRow = 0; iRow < dim; iRow++) {
    int eCol = src.ListRowSelect[iRow];
    for (int iCol = 0; iCol < dim; iCol++) {
      Mred(iRow, iCol) = M(eCol, iCol);
    }
  }
#ifdef DEBUG_LAMINATIONS
  std::cerr << "LAM: We have Mred\n";
#endif
  MyMatrix<T> MredInv = Inverse(Mred);
  Face f_select(nbRow);
  for (auto &eVal : src.ListRowSelect) {
    f_select[eVal] = 1;
  }
#ifdef DEBUG_LAMINATIONS
  std::cerr << "LAM: We have f_select\n";
#endif
  std::vector<int> Voth;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (f_select[iRow] == 0) {
      Voth.push_back(iRow);
    }
  }
  MyMatrix<T> M2 = M * MredInv;
#ifdef DEBUG_LAMINATIONS
  std::cerr << "LAM: We have M2\n";
#endif
  BlockCppIterator bci(dim - 1, 2);
  auto is_non_zero = [&](std::vector<int> const &V) -> bool {
    for (auto &val : V) {
      if (val > 0) {
        return true;
      }
    }
    return false;
  };
  MyVector<T> vect = ZeroVector<T>(dim);
  Face f_ass(nbRow);
  auto is_correct = [&](std::vector<int> const &V) -> bool {
    if (!is_non_zero(V)) {
      return false;
    }
    for (int i = 0; i < dim - 1; i++) {
      vect(i + 1) = UniversalScalarConversion<T, int>(V[i]);
    }
#ifdef DEBUG_LAMINATIONS
    std::cerr << "LAM: We have vect\n";
#endif
    for (auto &eRow : Voth) {
      T sum(0);
      for (int iCol = 0; iCol < dim; iCol++) {
        sum += vect(iCol) * M2(eRow, iCol);
      }
#ifdef DEBUG_LAMINATIONS
      std::cerr << "LAM: eRow=" << eRow << " sum=" << sum << "\n";
#endif
      if (sum != 0 && sum != 1) {
        return false;
      }
    }
    return true;
  };
  auto set_f_ass = [&]() -> void {
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T sum(0);
      for (int iCol = 0; iCol < dim; iCol++) {
        sum += vect(iCol) * M2(iRow, iCol);
      }
#ifdef DEBUG_LAMINATIONS
      if (sum != 0 && sum != 1) {
        std::cerr << "sum should be 0 or 1\n";
        throw TerminalException{1};
      }
#endif
      if (sum == 0) {
        f_ass[iRow] = 0;
      }
      if (sum == 1) {
        f_ass[iRow] = 1;
      }
    }
#ifdef DEBUG_LAMINATIONS
    std::cerr << "LAM: f_ass done\n";
#endif
  };
  for (auto &eV : bci) {
    if (is_correct(eV)) {
      set_f_ass();
      bool test = f(f_ass);
      if (test) {
        return;
      }
    }
  }
}

template <typename T>
vectface compute_all_two_laminations(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  vectface vf(nbRow);
  auto f = [&](Face const &f_ass) -> bool {
    vf.push_back(f_ass);
    return false;
  };
  compute_two_laminations_f(M, f);
  return vf;
}

template <typename T>
std::optional<Face> compute_one_two_laminations(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  vectface vf(nbRow);
  std::optional<Face> opt;
  auto f = [&](Face const &f_ass) -> bool {
    opt = f_ass;
    return true;
  };
  compute_two_laminations_f(M, f);
  return opt;
}

// clang-format off
#endif  // SRC_POLY_POLY_LAMINATIONS_H_
// clang-format on
