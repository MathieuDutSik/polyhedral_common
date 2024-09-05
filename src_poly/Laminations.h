// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LAMINATIONS_H_
#define SRC_POLY_POLY_LAMINATIONS_H_

#include "MAT_Matrix.h"
#include "COMB_Combinatorics.h"

#ifdef DEBUG
#define DEBUG_LAMINATIONS
#endif



template<typename T, typename F>
void compute_two_laminations_f(MyMatrix<T> const& M, F f) {
#ifdef DEBUG_LAMINATIONS
  if (RankMat(M) != M.cols()) {
    std::cerr << "M should be full dimensional\n";
    throw TerminalException{1};
  }
#endif
  int dim = M.cols();
  int nbRow = M.rows();
  SelectionRowCol<T> src = TMat_SelectRowCol(M);
  MyMatrix<T> Mred(dim, dim);
  for (int iRow=0; iRow<dim; iRow++) {
    int eCol = src.ListRowSelect[iRow];
    for (int iCol=0; iCol<dim; iCol++) {
      Mred(iRow, iCol) = M(eCol, iCol);
    }
  }
  MyMatrix<T> MredInv = Inverse(Mred);
  Face f_select(nbRow);
  for (auto & eVal : src.ListRowSelect) {
    f_select[eVal] = 1;
  }
  std::vector<int> Voth;
  for (int iRow=0; iRow<nbRow; iRow++) {
    if (f_select[iRow] == 0) {
      Voth.push_back(iRow);
    }
  }
  MyMatrix<T> M2 = MredInv * M;
  BlockCppIterator bci(dim-1, 2);
  auto is_non_zero=[&](std::vector<int> const& V) -> bool {
    for (auto & val: V) {
      if (val > 0) {
        return true;
      }
    }
    return false;
  };
  MyVector<T> vect = ZeroVector<T>(dim);
  Face f_ass(nbRow);
  auto is_correct=[&](std::vector<int> const& V) -> bool {
    if (!is_non_zero(V)) {
      return false;
    }
    for (int i=0; i<dim; i++) {
      vect(i+1) = UniversalScalarConversion<T,int>(V[i]);
    }
    for (auto & eRow : Voth) {
      T sum(0);
      for (int iCol=0; iCol<dim; iCol++) {
        sum += vect(iCol) * M2(eRow, iCol);
      }
      if (sum != 0 && sum != 1) {
        return false;
      }
    }
    return true;
  };
  auto set_f_ass=[&]() -> void {
    for (int iRow=0; iRow<nbRow; iRow++) {
      T sum(0);
      for (int iCol=0; iCol<dim; iCol++) {
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
  };
  for (auto & eV : bci) {
    if (is_correct(eV)) {
      set_f_ass();
      bool test = f(f_ass);
      if (test) {
        return;
      }
    }
  }
}

template<typename T>
vectface compute_all_two_laminations(MyMatrix<T> const& M) {
  int nbRow = M.rows();
  vectface vf(nbRow);
  auto f=[&](Face const& f_ass) -> bool {
    vf.push_back(f_ass);
    return false;
  };
  compute_two_laminations_f(M, f);
  return vf;
}

template<typename T>
std::optional<Face> compute_one_two_laminations(MyMatrix<T> const& M) {
  int nbRow = M.rows();
  vectface vf(nbRow);
  std::optional<Face> opt;
  auto f=[&](Face const& f_ass) -> bool {
    opt = f_ass;
    return true;
  };
  compute_two_laminations_f(M, f);
  return opt;
}


// clang-format off
#endif  // SRC_POLY_POLY_LAMINATIONS_H_
// clang-format on

