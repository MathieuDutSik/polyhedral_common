// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_
#define SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_


#include "POLY_LinearProgramming.h"
#include "POLY_SamplingFacet.h"
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
  os << "Before GetFullRankFacetSet\n";
  vectface vf = GetFullRankFacetSet(FAC, os);
  for (auto & face : vf) {
    MyVector<T> eEXT = FindFacetInequality(FAC, face);
    (void)f_insert(eEXT);
  }
  os << "Now going into the second iteration\n";
  int iter = 0;
  MyMatrix<T> EXT = MatrixFromVectorFamily(ListEXT);
  while (true) {
    vectface vf_myfac = f_dual(EXT);
    int n_ext = SetEXT.size();
    os << "iter=" << iter << " |vf_myfac|=" << vf_myfac.size() << " n_ext=" << n_ext << "\n";
    std::vector<MyVector<T>> ListNewEXT;
    auto has_violating_facet=[&](MyVector<T> const& eFAC) -> bool {
      for (auto & eNewEXT : ListNewEXT) {
        T scal = eNewEXT.dot(eFAC);
        if (scal < 0) {
          return true;
        }
      }
      return false;
    };
    for (auto & face : vf_myfac) {
      MyVector<T> eFAC = FindFacetInequality(EXT, face);
      MyVector<T> eFACred = AbsoluteRescaleVector(eFAC);
      if (SetFAC.count(eFACred) == 0) { // Missing so operation is needed
        if (!has_violating_facet(eFACred)) { // Check if we already had something matching
          Face face = FindViolatedFaceFast(FAC, eFACred);
          MyVector<T> eEXT = FindFacetInequality(FAC, face);
          f_insert(eEXT);
          ListNewEXT.push_back(eEXT);
        }
      }
    }
    if (ListNewEXT.size() == 0) {
      os << "Exiting the infinite loop\n";
      break;
    }
    EXT = MatrixFromVectorFamily(ListEXT);
    iter++;
  }
  os << "Returning EXT\n";
  return EXT;
}

template<typename T>
MyMatrix<T> POLY_DualDescription_PrimalDualInequalities(MyMatrix<T> const& FAC, std::ostream &os) {
  auto f_dual=[&](MyMatrix<T> const& FACin) -> vectface {
    return lrs::DualDescription_incd(FACin);
  };
  return POLY_DualDescription_PrimalDual_Kernel(FAC, f_dual, os);
}

template<typename T>
vectface POLY_DualDescription_PrimalDualIncidence(MyMatrix<T> const& FAC, std::ostream &os) {
  MyMatrix<T> EXT = POLY_DualDescription_PrimalDualInequalities(FAC, os);
  int dim = FAC.cols();
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
  os << "|vf|=" << vf.size() << "\n";
  return vf;
}

template<typename T, typename Fprocess>
void POLY_DualDescription_PrimalDualFaceIneq(MyMatrix<T> const& FAC, Fprocess f_process, std::ostream &os) {
  MyMatrix<T> EXT = POLY_DualDescription_PrimalDualInequalities(FAC, os);
  int dim = FAC.cols();
  int nbRow = FAC.rows();
  vectface vf(nbRow);
  for (int i_ext=0; i_ext<EXT.rows(); i_ext++) {
    MyVector<T> eEXT = GetMatrixRow(EXT, i_ext);
    Face f(nbRow);
    for (int i_row=0; i_row<nbRow; i_row++) {
      T sum = 0;
      for (int i=0; i<dim; i++)
        sum += FAC(i_row,i) * eEXT(i);
      if (sum == 0)
        f[i_row] = 1;
    }
    std::pair<Face,MyVector<T>> pair{f, eEXT};
    f_process(pair);
  }
}


// clang-format off
#endif  // SRC_POLY_POLY_DUALDESCRIPTION_PRIMALDUAL_H_
// clang-format on


