#ifndef GSL_RELATED_INCLUDE
#define GSL_RELATED_INCLUDE

#include "MAT_Matrix.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>


template<typename T>
MyVector<T> T_FindNegativeVector(MyMatrix<T> const & eMat)
{
  gsl_matrix * Gram;
  gsl_matrix * EigenVectors;
  gsl_vector * eigenvalues;
  gsl_vector * FundamentalLevel;
  gsl_eigen_symmv_workspace * workspace;
  int dimension;
  int i, j, iColSel, ScalarMult;
  double eEigSel, TheEig, eVal_d;
  int *CURRENT;
  T eVal, eSum;
  dimension=eMat.rows();
  eigenvalues = gsl_vector_alloc (dimension);
  Gram = gsl_matrix_alloc (dimension, dimension);
  EigenVectors = gsl_matrix_alloc (dimension, dimension);
  workspace = gsl_eigen_symmv_alloc (dimension);
  if ((CURRENT = (int*)malloc(dimension*sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  for (i=0; i<dimension; i++)
    for (j=0; j<dimension; j++) {
      eVal=eMat(i, j);
      GET_DOUBLE(eVal, eVal_d);
      gsl_matrix_set(Gram, i, j, eVal_d);
    }
  gsl_eigen_symmv (Gram, eigenvalues, EigenVectors, workspace);
  iColSel=-1;
  eEigSel=0;
  for (i=0; i<dimension; i++) {
    TheEig=(double)gsl_vector_get(eigenvalues, i);
    if (TheEig < eEigSel) {
      eEigSel=TheEig;
      iColSel=i;
    }
  }
  /*  gsl_eigen_symmv_sort(eigenvalues, EigenVectors, GSL_EIGEN_SORT_ABS_ASC);*/
  /* lowest=gsl_vector_get(eigenvalues, 0);*/
  MyVector<T> eVec=MyVector<T>(dimension);
  if (eEigSel<0) {
    FundamentalLevel = gsl_vector_alloc (dimension);
    gsl_matrix_get_col(FundamentalLevel, EigenVectors, iColSel);
    ScalarMult=1;
    while(true) {
      for (i=0; i<dimension; i++)
	CURRENT[i]=(int)rint((double)gsl_vector_get(FundamentalLevel, i)*((double)ScalarMult));
      eSum=0;
      for (i=0; i<dimension; i++)
	for (j=0; j<dimension; j++) {
	  eVal=eMat(i, j);
	  eSum=eSum + eVal*CURRENT[i]*CURRENT[j];
	}
      if (eSum<0)
	break;
      ScalarMult++;
    }
    eVec->n=dimension;
    for (i=0; i<dimension; i++)
      TVec_Assign(eVec, i, CURRENT[i]);
    gsl_vector_free(FundamentalLevel);
  }
  else {
    std::cerr << "Matrix does seem positive semidefinite in double precision\n";
    std::cerr << "Most likely you will need to reprogram\n";
    throw TerminalException{1};
  }
  gsl_vector_free(eigenvalues);
  gsl_matrix_free(Gram);
  gsl_matrix_free(EigenVectors);
  gsl_eigen_symmv_free(workspace);
  free(CURRENT);
  return eVec;
}



#endif
