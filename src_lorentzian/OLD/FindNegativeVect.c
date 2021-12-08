#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/* This program is for choosing a vector of negative norm
   that is hopefully small */

int main(int argc, char *argv[])
{
  int i, j, dimension;
  gsl_matrix * Gram;
  gsl_matrix * EigenVectors;
  gsl_vector * eigenvalues;
  gsl_vector * FundamentalLevel;
  gsl_eigen_symmv_workspace * workspace;
  int eRet;
  int* CURRENT;
  int** GramInt;
  int Sum, ScalarMult;
  FILE *FileMat;
  FILE *FileRes;
  double TheEig,  eEigSel;
  int iColSel;

  if (argc !=3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "FindNegativeVect [MatFile] [ResFile]\n");
      return -1;
    }
  FileMat=fopen(argv[1], "r");
  if (FileMat == NULL)
    {
      fprintf(stderr,"graph file %s was not found\n",argv[1]);
      return -1;
    }
  eRet=fscanf(FileMat, "%d", &dimension);

  if ((GramInt = (int**)malloc(dimension*sizeof(int*))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<dimension; i++)
    if ((GramInt[i] = (int*)malloc(dimension*sizeof(int))) == 0)
      exit (EXIT_FAILURE);
  if ((CURRENT = (int*)malloc(dimension*sizeof(int))) == 0)
    exit (EXIT_FAILURE);

  eigenvalues = gsl_vector_alloc (dimension);
  Gram = gsl_matrix_alloc (dimension, dimension);
  EigenVectors = gsl_matrix_alloc (dimension, dimension);
  workspace = gsl_eigen_symmv_alloc (dimension);

  gsl_matrix_fscanf (FileMat, Gram);
  for (i=0; i<dimension; i++)
    for (j=0; j<dimension; j++)
      GramInt[i][j]=(int) gsl_matrix_get(Gram, i, j);
  fclose (FileMat);

  gsl_eigen_symmv (Gram, eigenvalues, EigenVectors, workspace);

  iColSel=-1;
  eEigSel=0;
  for (i=0; i<dimension; i++)
    {
      TheEig=(double)gsl_vector_get(eigenvalues, i);
      /* fprintf(stderr,"i=%d eig=%f\n", i, TheEig);*/
      if (TheEig < eEigSel)
	{
	  eEigSel=TheEig;
	  iColSel=i;
	}
    }
  /*  gsl_eigen_symmv_sort(eigenvalues, EigenVectors, GSL_EIGEN_SORT_ABS_ASC);*/
  /* lowest=gsl_vector_get(eigenvalues, 0);*/
  FileRes=fopen(argv[2], "w");
  /*  fprintf(stderr, "eEigSel=%f\n", eEigSel);*/
  if (eEigSel<0)
    {
      FundamentalLevel = gsl_vector_alloc (dimension);
      gsl_matrix_get_col(FundamentalLevel, EigenVectors, iColSel);
      ScalarMult=1;
      while(1)
	{
	  for (i=0; i<dimension; i++)
	    CURRENT[i]=(int)rint((double)gsl_vector_get(FundamentalLevel, i)*((double)ScalarMult));
	  Sum=0;
	  for (i=0; i<dimension; i++)
	    for (j=0; j<dimension; j++)
	      Sum+=GramInt[i][j]*CURRENT[i]*CURRENT[j];
	  if (Sum<0)
	      break;
	  /*	  fprintf(stderr, "ScalarMult=%d Sum=%d\n", ScalarMult, Sum);*/
	  ScalarMult++;
	}
      fprintf(FileRes, "return rec(pos_semidef:=false, eVect:=[");
      for (i=0; i<dimension; i++)
	{
	  fprintf(FileRes, "%d", CURRENT[i]);
	  if (i<dimension-1)
	    fprintf(FileRes, ",");
	}
      fprintf(FileRes, "]);\n");
      gsl_vector_free(FundamentalLevel);
    }
  else
    {
      fprintf(FileRes, "return rec(pos_semidef:=true);\n");
    }
  fclose(FileRes);
  
  gsl_vector_free(eigenvalues);
  gsl_matrix_free(Gram);
  gsl_matrix_free(EigenVectors);
  gsl_eigen_symmv_free(workspace);
  free(CURRENT);
  for (i=0; i<dimension; i++)
    free(GramInt[i]);
  free(GramInt);
  return 0;
}
