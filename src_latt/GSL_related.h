// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_GSL_RELATED_H_
#define SRC_LATT_GSL_RELATED_H_

// clang-format off
#include "MAT_Matrix.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
// clang-format on

template <typename T>
MyVector<T> T_FindNegativeVector(MyMatrix<T> const &eMat) {
  gsl_matrix *Gram;
  gsl_matrix *EigenVectors;
  gsl_vector *eigenvalues;
  gsl_vector *FundamentalLevel;
  gsl_eigen_symmv_workspace *workspace;
  int dimension;
  int i, j, iColSel, ScalarMult;
  double eEigSel, TheEig;
  T eVal, eSum;
  dimension = eMat.rows();
  eigenvalues = gsl_vector_alloc(dimension);
  Gram = gsl_matrix_alloc(dimension, dimension);
  EigenVectors = gsl_matrix_alloc(dimension, dimension);
  workspace = gsl_eigen_symmv_alloc(dimension);
  for (i = 0; i < dimension; i++)
    for (j = 0; j < dimension; j++) {
      eVal = eMat(i, j);
      double eVal_d = UniversalScalarConversion<double, T>(eVal);
      gsl_matrix_set(Gram, i, j, eVal_d);
    }
  gsl_eigen_symmv(Gram, eigenvalues, EigenVectors, workspace);
  iColSel = -1;
  eEigSel = 0;
  for (i = 0; i < dimension; i++) {
    TheEig = static_cast<double>(gsl_vector_get(eigenvalues, i));
    if (TheEig < eEigSel) {
      eEigSel = TheEig;
      iColSel = i;
    }
  }
  /*  gsl_eigen_symmv_sort(eigenvalues, EigenVectors, GSL_EIGEN_SORT_ABS_ASC);*/
  /* lowest=gsl_vector_get(eigenvalues, 0);*/
  MyVector<T> eVec = MyVector<T>(dimension);
  if (eEigSel < 0) {
    FundamentalLevel = gsl_vector_alloc(dimension);
    gsl_matrix_get_col(FundamentalLevel, EigenVectors, iColSel);
    ScalarMult = 1;
    while (true) {
      for (i = 0; i < dimension; i++) {
        double f1 = static_cast<double>(gsl_vector_get(FundamentalLevel, i));
        double f2 = static_cast<double>(ScalarMult);
        eVec(i) = static_cast<int>(rint(f1 * f2));
      }
      eSum = 0;
      for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++) {
          eVal = eMat(i, j);
          eSum = eSum + eVal * eVec(i) * eVec(j);
        }
      if (eSum < 0)
        break;
      ScalarMult++;
    }
    eVec->n = dimension;
    gsl_vector_free(FundamentalLevel);
  } else {
    std::cerr << "Matrix does seem positive semidefinite in double precision\n";
    std::cerr << "Most likely you will need to reprogram\n";
    throw TerminalException{1};
  }
  gsl_vector_free(eigenvalues);
  gsl_matrix_free(Gram);
  gsl_matrix_free(EigenVectors);
  gsl_eigen_symmv_free(workspace);
  return eVec;
}

// clang-format off
#endif  // SRC_LATT_GSL_RELATED_H_
// clang-format on
