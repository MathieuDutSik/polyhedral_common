// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_
#define SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_


#include "POLY_LinearProgramming.h"
#include "POLY_PolytopeFct.h"
#include "POLY_lrslib.h"


template<typename T>
MyVector<T> AbsoluteRescaleVector(MyVector<T> const& V) {
  int dim = V.size();
  T sum = 0;
  for (int i=0; i<dim; i++)
    sum += T_abs(V(i));
  if (sum == 0)
    return ZeroVector<T>(dim);
  return V / sum;
}


template<typename T, typename Fdual>
MyMatrix<T> POLY_DualDescription_PrimalDual_Kernel(MyMatrix<T> const& FAC, Fdual f_dual, std::ostream &os) {
  std::set<MyVector<T>> SetFAC;
  int n_rows_fac = FAC.rows();
  for (int i_row=0; i_row<n_rows_fac; i_row++) {
    MyVector<T> eRow = GetMatrixRow(FAC, i_row);
    MyVector<T> eRowRed = AbsoluteRescaleVector(eRow);
    SetFAC.insert(eRowRed);
  }
  std::vector<MyVector<T>> ListEXT;
  std::set<MyVector<T>> SetEXT;
  auto f_insert=[&](MyVector<T> const& eEXT) -> void {
    MyVector<T> eEXTred = AbsoluteRescaleVector(eEXT);
    if (SetEXT.count(eEXTred) == 0) {
      SetEXT.insert(eEXTred);
      ListEXT.push_back(eEXTred);
    }
  };
  int dim = FAC.cols();
  int nb = 100;
  int iter = 0;
  MyMatrix<T> EXT;
  while(true) {
    vectface vf = FindVertices(FAC, nb);
    for (auto & face : vf) {
      MyVector<T> eEXT = FindFacetInequality(FAC, face);
      f_insert(eEXT);
    }
    EXT = MatrixFromVectorFamily(ListEXT);
    int rnk = RankMat(EXT);
    os << "iter=" << iter << " rnk=" << rnk << " |EXT|=" << EXT.rows() << "\n";
    if (rnk == dim) {
      os << "Exiting the loop\n";
      break;
    }
  }
  os << "Now going into the second iteration\n";
  iter = 0;
  while (true) {
    vectface vf_myfac = f_dual(EXT);
    int n_before = SetEXT.size();
    os << "------------------------\n";
    os << "iter=" << iter << " |vf_myfac|=" << vf_myfac.size() << " n_before=" << n_before << "\n";
    std::vector<MyVector<T>> ListNew;
    for (auto & face : vf_myfac) {
      MyVector<T> eFAC = FindFacetInequality(EXT, face);
      MyVector<T> eFACred = AbsoluteRescaleVector(eFAC);
      if (SetFAC.count(eFACred) == 0) {
        Face face = FindViolatedFace(FAC, eFACred);
        MyVector<T> eEXT = FindFacetInequality(FAC, face);
        f_insert(eEXT);
      }
    }
    int n_after = SetEXT.size();
    os << "n_after=" << n_after << "\n";
    EXT = MatrixFromVectorFamily(ListEXT);
    if (n_before == n_after) {
      os << "Exiting the second loop\n";
      break;
    }
  }
  return EXT;
}

template<typename T>
vectface POLY_DualDescription_PrimalDual(MyMatrix<T> const& FAC, std::ostream &os) {
  int dim = FAC.cols();
  auto f_dual=[&](MyMatrix<T> const& FACin) -> vectface {
    return lrs::DualDescription_incd(FACin);
  };
  MyMatrix<T> EXT = POLY_DualDescription_PrimalDual_Kernel(FAC, f_dual, os);
  int nbRow = FAC.rows();
  vectface vf(nbRow);
  for (int i_ext=0; i_ext<EXT.rows(); i_ext++) {
    Face f(nbRow);
    for (int i_row=0; i_row<nbRow; i_row++) {
      T sum = 0;
      for (int i=0; i<dim; i++)
        sum += FAC(i_row,i) * EXT(i_ext,i);
      if (sum == 0)
        f[i_row] = 1;
    }
    vf.push_back(f);
  }
  return vf;
}




// clang-format off
#endif  // SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_
// clang-format on


