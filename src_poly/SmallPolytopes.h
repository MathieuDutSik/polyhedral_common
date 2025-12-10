// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_SMALLPOLYTOPES_H_
#define SRC_POLY_SMALLPOLYTOPES_H_

// clang-format off
#include "POLY_Fundamental.h"
#include <vector>
#include <utility>
// clang-format on

// For simplicial and near simplicial polytopes,
// there is exact solution and so no need to use
// lrs / cdd and the like that would be more
// expensive.

#ifdef DEBUG
#define DEBUG_SMALL_POLYTOPE
#endif

#ifdef DISABLE_DEBUG_SMALL_POLYTOPE
#undef DEBUG_SMALL_POLYTOPE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SMALL_POLYTOPE
#endif

template <typename T> vectface Simplicial_Incidence(MyMatrix<T> const &EXT) {
  int n = EXT.cols();
  Face f_full(n);
  for (int i = 0; i < n; i++)
    f_full[i] = 1;
  vectface vf(n);
  for (int i = 0; i < n; i++) {
    Face f = f_full;
    f[i] = 0;
    vf.push_back(f);
  }
  return vf;
}

template <typename T>
vectface NearSimplicial_Incidence(MyMatrix<T> const &EXT,
                                  [[maybe_unused]] std::ostream &os) {
  int n_ext = EXT.rows();
  MyMatrix<T> NSP = NullspaceMat(EXT);
#ifdef SANITY_CHECK_SMALL_POLYTOPE
  if (NSP.rows() != 1) {
    std::cerr << "The rank is incorrect\n";
    throw TerminalException{1};
  }
#endif
  std::vector<int> V_p, V_m, V_z;
  for (int i_ext = 0; i_ext < n_ext; i_ext++) {
    T val = NSP(0, i_ext);
    if (val > 0)
      V_p.push_back(i_ext);
    if (val < 0)
      V_m.push_back(i_ext);
    if (val == 0)
      V_z.push_back(i_ext);
  }
#ifdef DEBUG_SMALL_POLYTOPE
  os << "n_ext=" << n_ext << "\n";
  os << "|V_p|=" << V_p.size() << " |V_m|=" << V_m.size()
     << " |V_z|=" << V_z.size() << "\n";
#endif
#ifdef SANITY_CHECK_SMALL_POLYTOPE
  auto check = [&](Face const &f) -> void {
    CheckFacetInequality(EXT, f, "SmallPolytopes");
  };
#endif
  Face f_full(n_ext);
  for (int i = 0; i < n_ext; i++)
    f_full[i] = 1;
  vectface vf(n_ext);
  for (auto &x_p : V_p) {
    for (auto &x_m : V_m) {
      Face f = f_full;
      f[x_p] = 0;
      f[x_m] = 0;
#ifdef DEBUG_SMALL_POLYTOPE
      os << "x_p=" << x_p << " x_m=" << x_m << " |f|=" << f.size() << " / "
         << f.count() << "\n";
#endif
#ifdef SANITY_CHECK_SMALL_POLYTOPE
      check(f);
#endif
      vf.push_back(f);
    }
  }
  for (auto &x_z : V_z) {
    Face f = f_full;
    f[x_z] = 0;
#ifdef DEBUG_SMALL_POLYTOPE
    os << "x_z=" << x_z << " |f|=" << f.size() << " / " << f.count() << "\n";
#endif
#ifdef SANITY_CHECK_SMALL_POLYTOPE
    check(f);
#endif
    vf.push_back(f);
  }
  return vf;
}

template <typename T>
vectface SmallPolytope_Incidence(MyMatrix<T> const &EXT, std::ostream &os) {
#ifdef DEBUG_SMALL_POLYTOPE
  os << "SmallPolytope_Incidence, begin\n";
#endif
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  if (n_row == n_col)
    return Simplicial_Incidence(EXT);
  if (n_row == n_col + 1)
    return NearSimplicial_Incidence(EXT, os);
  std::cerr << "Unfortunately, there is no simple strategy for those higher "
               "incidence cases\n";
  throw TerminalException{1};
}

template <typename T>
MyMatrix<T> SmallPolytope_Ineq(MyMatrix<T> const &EXT, std::ostream &os) {
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  vectface vf = SmallPolytope_Incidence(EXT, os);
  MyMatrix<Tint> EXT_int = Get_EXT_int(EXT);
  SubsetRankOneSolver<T> solver(EXT_int);
  size_t n_fac = EXT.rows();
  int n_col = EXT.cols();
  MyMatrix<T> FAC(n_fac, n_col);
  size_t i_fac = 0;
  for (auto &eFace : vf) {
    MyVector<Tint> Vint = solver.GetPositiveKernelVector(eFace);
    for (int i_col = 0; i_col < n_col; i_col++) {
      FAC(i_fac, i_col) = UniversalScalarConversion<T, Tint>(Vint(i_col));
    }
    i_fac++;
  }
  return FAC;
}

template <typename T, typename Fprocess>
void SmallPolytope_FaceIneq(MyMatrix<T> const &EXT, Fprocess f_process,
                            std::ostream &os) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
#ifdef DEBUG_SMALL_POLYTOPE
  os << "SmallPolytope_FaceIneq, n_row=" << n_row << " n_col=" << n_col << "\n";
#endif
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  vectface vf = SmallPolytope_Incidence(EXT, os);
  MyMatrix<Tint> EXT_int = Get_EXT_int(EXT);
  SubsetRankOneSolver<T> solver(EXT_int);
  std::pair<Face, MyVector<T>> pair{Face(n_row), MyVector<T>(n_col)};
  for (auto &eFace : vf) {
    MyVector<Tint> Vint = solver.GetPositiveKernelVector(eFace);
    for (int i_col = 0; i_col < n_col; i_col++) {
      pair.second(i_col) = UniversalScalarConversion<T, Tint>(Vint(i_col));
    }
    pair.first = eFace;
#ifdef SANITY_CHECK_SMALL_POLYTOPE
    for (int i_row = 0; i_row < n_row; i_row++) {
      T scal(0);
      for (int i_col = 0; i_col < n_col; i_col++) {
        scal += pair.second(i_col) * EXT(i_row, i_col);
      }
      if (scal < 0) {
        std::cerr << "Negative scalar product at i_row=" << i_row
                  << " scal=" << scal << "\n";
        throw TerminalException{1};
      }
      bool test1 = false;
      if (scal == 0) {
        test1 = true;
      }
      bool test2 = false;
      if (eFace[i_row] == 1) {
        test2 = true;
      }
      if (test1 != test2) {
        std::cerr << "At i_row=" << i_row << " We have test1=" << test1
                  << " test2=" << test2 << "\n";
        throw TerminalException{1};
      }
    }
#endif
    f_process(pair);
  }
}

// clang-format off
#endif  // SRC_POLY_SMALLPOLYTOPES_H_
// clang-format on
