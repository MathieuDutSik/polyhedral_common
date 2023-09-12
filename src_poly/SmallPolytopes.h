// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_SMALL_POLYTOPES_H_
#define SRC_POLY_SMALL_POLYTOPES_H_

// clang-format off
#include "POLY_PolytopeFct.h"
// clang-format on

// For simplicial and near simplicial polytopes,
// there is exact solution and so no need to use
// lrs / cdd and the like that would be more
// expensive.

template<typename T>
vectface Simplicial_Incidence(MyMatrix<T> const& EXT) {
  int n = EXT.cols();
  Face f_full(n);
  for (int i=0; i<n; i++)
    f[i] = 1;
  vectface vf(n);
  for (int i=0; i<n; i++) {
    Face f = f_full;
    f[i] = 0;
    vf.push_back(f);
  }
  return vf;
}

template<typename T>
vectface NearSimplicial_Incidence(MyMatrix<T> const& EXT) {
  int n = EXT.cols();
  int n_ext = EXT.rows();
  MyMatrix<T> NSP = NullspaceMat(EXT);
#ifdef DEBUG_SMALL_POLYTOPE
  if (NSP.rows() != 1) {
    std::cerr << "The rank is incorrect\n";
    throw TerminalException{1};
  }
#endif
  std::vector<int> V_p, V_m, V_z;
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    T val = NSP(0,i_ext);
    if (val > 0)
      V_p.push_back(i_ext);
    if (val < 0)
      V_m.push_back(i_ext);
    if (val == 0)
      V_z.push_back(i_ext);
  }
  Face f_full(n_ext);
  for (int i=0; i<n_ext; i++)
    f[i] = 1;
  vectface vf(n_ext);
  for (auto & x_p : V_p) {
    for (auto & x_m : V_m) {
      Face f = f_full;
      f[x_p] = 0;
      f[x_m] = 0;
      vf.push_back(f);
    }
  }
  for (auto & x_z : V_z) {
    Face f = f_full;
    f[x_z] = 0;
    vf.push_back(f);
  }
  return vf;
}

template<typename T>
vectface SmallPolytope_Incidence(MyMatrix<T> const& EXT) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  if (n_row == n_col)
    return Simplicial_Incidence(EXT);
  if (n_row == n_col+1)
    return NearSimplicial_Incidence(EXT);
  std::cerr << "Unfortunately, there is no simple strategy for those higher incidence cases\n";
  throw TerminalException{1};
}


template<typename T>
MyMatrix<T> SmallPolytope_Ineq(MyMatrix<T> const& EXT) {
  using Tint = typename underlying_ring<T>::ring_type;
  vectface vf = Simplicial_Incidence(EXT);
  MyMatrix<Tint> EXT_int = Get_EXT_int(EXT);
  SubsetRankOneSolver<T> solver(EXT_int);
  size_t n_fac = EXT.rows();
  int n_col = EXT.cols();
  MyMatrix<T> FAC(n_fac, n_col);
  size_t i_fac = 0;
  for (auto & eFace : vf) {
    MyVector<Tint> Vint = solver.GetPositiveKernelVector(eFace);
    for (int i_col=0; i_col<n_col; i_col++) {
      FAC(i_fac,i_col) = UniversalScalarConversion<T,Tint>(Vint(i_col));
    }
    i_fac++;
  }
  return FAC;
}

template<typename T, typename Fprocess>
void SmallPolytope_FaceIneq(MyMatrix<T> const& EXT, Fprocess f_process) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  using Tint = typename underlying_ring<T>::ring_type;
  vectface vf = Simplicial_Incidence(EXT);
  MyMatrix<Tint> EXT_int = Get_EXT_int(EXT);
  SubsetRankOneSolver<T> solver(EXT_int);
  std::pair<Face, MyVector<T>> pair{Face(n_row), MyVector<T>(n_col)};
  for (auto & eFace : vf) {
    MyVector<Tint> Vint = solver.GetPositiveKernelVector(eFace);
    for (int i_col=0; i_col<n_col; i_col++) {
      pair.second(i_col) = UniversalScalarConversion<T,Tint>(Vint(i_col));
    }
    pair.first = eFace;
    f_process(pair);
  }
}

// clang-format off
#endif  // SRC_POLY_SMALL_POLYTOPES_H_
// clang-format on
