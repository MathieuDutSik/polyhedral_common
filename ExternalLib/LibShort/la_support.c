/* la_support.c  low-level functions for vectors and matrices */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

#include "la_support.h"
#include "panic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define LU_EPSILON 0.00000001 /* if  fabs(a) < LU_EPSILON then a = 0.0 */


/* ======================
 * = internal functions =
 * ======================
 */

static void *laMalloc(size_t size)
{
  void
    *ptr;

  ptr = (void *) calloc(1, size);
  if (ptr == NULL) {
    PANIC;
  }
  return ptr;
}


static void laFree(void *ptr)
{
  if (ptr == NULL) {
    PANIC;
  }
  free(ptr);
}


static void pivotSearch(int dim,
			double **mat,
			int row,
			int *pivot_row)
{
  int
    i;
  double
    pivot;

  *pivot_row = row;

  pivot = fabs(mat[row][row]);
  for (i = row + 1; i < dim; i++)
    if (fabs(mat[i][row]) > pivot) {
      pivot = fabs(mat[i][row]);
      *pivot_row = i;
    }
}


static void solveRecursively(int dim,
			     double **u_mat,
			     double *b,
			     double *x)
{
  int
    i,
    j;

  for (i = dim - 1; i >= 0; i = i - 1) {
    x[i] = b[i];
    for (j = i + 1; j < dim; j++)
      x[i] = x[i] - u_mat[i][j] * x[j];
    x[i] = x[i] / u_mat[i][i];
  }
}

/* ==============================
 * = number theoretic functions =
 * ==============================
 */

int factorial(int n)
{
  int
    i,
    f;

  f = 1;
  for (i = 1; i <= n; i++)
    f = f*i;
  return f;
}

/* =================================
 * = functions for double-matrices =
 * =================================
 */

void doubleMakeSquareMatrix(int dim,
			    double ***matrix)
{
  int
    i;

  *matrix = (double **) laMalloc((dim+1) * sizeof(double *));
  (*matrix)[0] = (double *) laMalloc(dim * dim * sizeof(double));
  for (i = 1; i < dim; i++)
    (*matrix)[i] = (*matrix)[0] + i * dim;
  (*matrix)[dim] = (*matrix)[0];
}


void doubleDestroySquareMatrix(int dim,
			       double ***matrix)
{
  laFree((*matrix)[dim]);
  laFree(*matrix);
  *matrix = NULL;
}


void doubleCopySquareMatrix(int dim,
			    double **from,
			    double **to)
{
  int
    i;

  for (i = 0; i < dim; i++)
    memcpy(to[i], from[i], dim * sizeof(double));
}


void doubleAddSquareMatrix(int dim,
			   double factor,
			   double **from,
			   double **to)
{
  int
    i,
    j;
  
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      to[i][j] = to[i][j] + factor*from[i][j];
}


void doubleSetIdentitySquareMatrix(int dim,
				   double **matrix)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      if (i == j) matrix[i][j] = 1;
      else matrix[i][j] = 0.0;
}


void doubleMultiplySquareMatrixVector(int dim,
				  double **matrix,
				  double *vec1,
				  double *vec2)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++) {
    vec2[i] = 0.0;
    for (j = 0; j < dim; j++) 
      vec2[i] += matrix[i][j] * vec1[j];
  }
}


void doubleRowAddSquareMatrix(int dim,
			      double **matrix,
			      double factor,
			      int from,
			      int to)
{
  int
    i;

  for (i = 0; i < dim; i++)
    matrix[to][i] = matrix[to][i] + factor * matrix[from][i];
}


void doubleRowChangeSquareMatrix(int dim,
				 double **matrix,
				 int index1,
				 int index2)
{
  int
    i;
  double
    t;

  for (i = 0; i < dim; i++) {
    t = matrix[index1][i];
    matrix[index1][i] = matrix[index2][i];
    matrix[index2][i] = t;
  }
}


int doubleComputeLUSquareMatrix(int dim,
				double **mat,
				double **l_mat,
				double **u_mat,
				double **p_mat)
{
  int
    i,
    j,
    k,
    pivot_row;
  double
    temp;

  /* 
   * initialization of l_mat, u_mat and p_mat
   */

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      u_mat[i][j] = mat[i][j];
      if (i == j) {
	l_mat[i][j] = 1;
	p_mat[i][j] = 1;
      }
      else {
	l_mat[i][j] = 0.0;
	p_mat[i][j] = 0.0;
      }
    }

  /*
   * determine l_mat and u_mat row by row
   */

  for (i = 0; i < dim; i++) {

    pivotSearch(dim, u_mat, i, &pivot_row);
    
    if (fabs(u_mat[pivot_row][i]) < LU_EPSILON) {
      return 0;
    }

    /* 
     * permute the i-th row and the pivot_row-th row, and save this in p_mat
     */

    for (j = 0; j < dim; j++) {
      temp = u_mat[i][j];
      u_mat[i][j] = u_mat[pivot_row][j];
      u_mat[pivot_row][j] = temp;
      temp = p_mat[i][j];
      p_mat[i][j] = p_mat[pivot_row][j];
      p_mat[pivot_row][j] = temp;
    }

    for (j = 0; j < i; j++) {
      temp = l_mat[i][j];
      l_mat[i][j] = l_mat[pivot_row][j];
      l_mat[pivot_row][j] = temp;
    }

    /*
     * define l_mat
     */

    for (j = i + 1; j < dim; j++)
      l_mat[j][i] = u_mat[j][i] / u_mat[i][i];

    /*
     * define u_mat
     */

    for (j = i + 1; j < dim; j++) {
      for (k = i + 1; k < dim; k++)
	u_mat[j][k] = u_mat[j][k] - 
	  u_mat[j][i] / u_mat[i][i] * u_mat[i][k];
      u_mat[j][i] = 0.0;
    }
  } /* for i */

  return 1;
}


void doubleTransposeSquareMatrix(int dim,
				 double **mat,
				 double **t_mat)
{
  int
    i,
    j;
  double
    temp;

  for (i = 0; i < dim; i++)
    t_mat[i][i] = mat[i][i];

  for (i = 0; i < dim; i++)
    for (j = i; j < dim; j++) {
      temp = mat[i][j];
      t_mat[i][j] = mat[j][i];
      t_mat[j][i] = temp;
    }
}


void doubleSolveLUSquareMatrix(int dim,
			       double **l_mat,
			       double **u_mat,
			       double **p_mat,
			       double *b,
			       double *x)
{
  double
    *c,
    *d;
  double
    **inv_l_mat;

  doubleMakeVector(dim, &c);
  doubleMakeVector(dim, &d);
  doubleMakeSquareMatrix(dim, &inv_l_mat);
  
  /* Ax = b  <=>  PAx = Pb  <=>  LUx = Pb  <=>  Ux = L'Pb */
  doubleMultiplySquareMatrixVector(dim, p_mat, b, c);

  doubleInverseLowerSquareMatrix(dim, l_mat, inv_l_mat);
  doubleMultiplySquareMatrixVector(dim, inv_l_mat, c, d);
  solveRecursively(dim, u_mat, d, x);

  doubleDestroySquareMatrix(dim, &inv_l_mat);

  doubleDestroyVector(dim, &d);
  doubleDestroyVector(dim, &c);
}


int doubleInverseSquareMatrix(int dim,
			      double **mat,
			      double **inv_mat)
{
  int
    i,
    col,
    ok;
  double
    **l_mat,
    **u_mat,
    **p_mat;
  double
    *x_vec,
    *e_vec;

  /*
   * 1. allocate memory
   */

  doubleMakeVector(dim, &x_vec);
  doubleMakeVector(dim, &e_vec);
  doubleMakeSquareMatrix(dim, &l_mat);
  doubleMakeSquareMatrix(dim, &u_mat);
  doubleMakeSquareMatrix(dim, &p_mat);

  /*
   * 2. LU-factorization of mat
   */

  ok = doubleComputeLUSquareMatrix(dim, mat, l_mat, u_mat, p_mat);

  /*
   * 3. compute the inverse column by column
   */
  
  if (ok) {

    for (col = 0; col < dim; col++) {
    
      /* e_vec = col-th vector of the standard basis */
      for(i = 0; i < dim; i++) {
	if (i == col) e_vec[i] = 1;
	else e_vec[i] = 0;
      }
    
      /* determine col-th column of the inverse */
      doubleSolveLUSquareMatrix(dim, l_mat, u_mat, p_mat, e_vec, x_vec);
  
      for (i = 0; i < dim; i++) {
	inv_mat[i][col] = x_vec[i];
      }
    }
  }  

  /*
   * 4. free memory
   */

  doubleDestroySquareMatrix(dim, &u_mat);
  doubleDestroySquareMatrix(dim, &l_mat);
  doubleDestroySquareMatrix(dim, &p_mat);
  doubleDestroyVector(dim, &e_vec);
  doubleDestroyVector(dim, &x_vec);

  return ok;
}


void doubleInverseLowerSquareMatrix(int dim,
				    double **l_mat,
				    double **inv_l_mat)
{
  double
    *e_vec;
  int
    i,
    row,
    col;

  doubleMakeVector(dim, &e_vec);
  
  for (col = 0; col < dim; col++) {
    for (i = 0; i < dim; i++)
      if (col == i) e_vec[i] = 1;
      else e_vec[i] = 0.0;

    for (row = 0; row < dim; row++) {
      inv_l_mat[row][col] = e_vec[row];
      for (i = 0; i < row; i++)
	inv_l_mat[row][col] = inv_l_mat[row][col]
	  - l_mat[row][i] * inv_l_mat[i][col];
      inv_l_mat[row][col] = inv_l_mat[row][col] / l_mat[row][row];
    }
  }
  doubleDestroyVector(dim, &e_vec);
}


void doubleInverseUpperSquareMatrix(int dim,
			double **u_mat,
			double **inv_u_mat)
{
  doubleTransposeSquareMatrix(dim, u_mat, u_mat);
  doubleInverseLowerSquareMatrix(dim, u_mat, inv_u_mat);
  doubleTransposeSquareMatrix(dim, inv_u_mat, inv_u_mat);
  doubleTransposeSquareMatrix(dim, u_mat, u_mat);
}


void doubleAbsoluteDeterminantSquareMatrix(int dim,
					   double **l_mat,
					   double **u_mat,
					   double *det)
{
  int
    i;
  double
    d_det;

  d_det = 1.0;
  for (i = 0; i < dim; i++)
    d_det = d_det * l_mat[i][i];
  for (i = 0; i < dim; i++)
    d_det = d_det * u_mat[i][i];
  *det = fabs(d_det);
}


int doubleComputeCholeskySquareMatrix(int dim,
				      double **mat,
				      double **u_mat)
{
  int
    i,
    j,
    k,
    l;

  for (i = 0; i < dim; i++)
    for (j = i; j < dim; j++)
      u_mat[i][j] = mat[i][j];

  for (i = 0; i < dim - 1; i++) {
    for (j = i + 1; j < dim; j++) {
      u_mat[j][i] = u_mat[i][j];
      u_mat[i][j] = u_mat[i][j] / u_mat[i][i];
    }
    for (k = i + 1; k < dim; k++)
      for (l = k; l < dim; l++)
	u_mat[k][l] = u_mat[k][l] - u_mat[k][i] * u_mat[i][l];
  }

  for (i = 0; i < dim; i++)
    for (j = 0; j < i; j++)
      u_mat[i][j] = 0;

  for (i = 0; i < dim; i++)
    if (u_mat[i][i] < 0)
      return 0;
    else
      u_mat[i][i] = sqrt(u_mat[i][i]);

  for (i = 0; i < dim; i++)
    for (j = i + 1; j < dim; j++)
      u_mat[i][j] = u_mat[i][i] * u_mat[i][j];

  return 1;
}


void doubleComputeSymmetrizationSquareMatrix(int dim,
					     double **mat,
					     double **sym_mat)
{
  int
    i,
    j,
    k;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      sym_mat[i][j] = 0.0;
      for (k = 0; k < dim; k++)
	sym_mat[i][j] = sym_mat[i][j] + mat[k][i] * mat[k][j];
    }
}


void doubleMultiplySquareMatrix(int dim,
				double **a,
				double **b,
				double **c)
{
  int
    i,
    j,
    k;
  double*
    t;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      c[i][j] = 0.0;
      for (k = 0; k < dim; k++) {
	t = a[i];
	c[i][j] += t[k] * b[k][j];
      }
    }

}


void doubleComputeRowNormSquareMatrix(int dim,
				      double **mat,
				      int row,
				      double *row_norm)
{
  int
    i;
  
  *row_norm = 0.0;
  for (i = 0; i < dim; i++)
    *row_norm = *row_norm +  mat[row][i] * mat[row][i];
}


void doublePermuteColumnsSquareMatrix(int dim,
				      double **mat,
				      int index1,
				      int index2)
{
  int
    i;
  double
    t;

  for (i = 0; i < dim; i++) {
    t = mat[i][index1];
    mat[i][index1] = mat[i][index2];
    mat[i][index2] = t;
  }
}


void doublePermuteRowsSquareMatrix(int dim,
				   double **mat,
				   int index1,
				   int index2)
{
  int
    i;
  double
    t;

  for (i = 0; i < dim; i++) {
    t = mat[index1][i];
    mat[index1][i] = mat[index2][i];
    mat[index2][i] = t;
  }
}


void double2intSquareMatrix(int dim,
			     double **d_mat,
			     int **l_mat)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      l_mat[i][j] = (int) floor(d_mat[i][j] + 0.5);      
    }
}


void doublePrintSquareMatrix(int dim,
			     double **mat)
{
  int
    i,
    j;

  printf("[");
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (j < dim - 1)
	printf("%3.4f, ", mat[i][j]);
      if (j == dim - 1 && i != dim - 1)
	printf("%3.4f;\\\n ", mat[i][j]);
    }
  }
  printf("%3.4f]\n", mat[dim - 1][dim - 1]);
}


void doublePrintGramMatrix(int dim,
			   double **gram_matrix)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j <= i; j++) {
      printf("% 3.4f", gram_matrix[i][j]);
    }
    printf("\n");
  }
}


/*
 * ================================
 * = functions for double-vectors =
 * ================================
 */

void doubleMakeVector(int dim,
		      double **vector)
{
  *vector = (double *) laMalloc(dim * sizeof(double));
}


void doubleDestroyVector(int dim,
			 double **vector)
{
  laFree(*vector);
  *vector = NULL;
}


void doubleCopyVector(int dim,
		      double *from,
		      double *to)
{
  memcpy(to, from, dim * sizeof(double));
}


void doubleMultiplyVector(int dim,
			  double factor,
			  double *from,
			  double *to)
{
  int
    i;
  
  for (i = 0; i < dim; i++)
    to[i] = factor * from[i];
}

void doubleAddVector(int dim,
		     double factor,
		     double *from,
		     double *to)
{
  int
    i;

  for (i = 0; i < dim; i++)
    to[i] = to[i] + factor * from[i];
}


void doubleInnerProduct(int dim,
			double *vector1,
			double *vector2,
			double *product)
{
  int
    i;

  *product = 0.0;

  for (i = 0; i < dim; i++) {
    *product = *product + vector1[i] * vector2[i];
  }
}


void doubleVectorNorm(int dim,
		      double **gram_matrix,
		      double *vector,
		      double *norm)
{
  int
    i,
    j;
  double
    temp;

  *norm = 0.0;

  for (i = 0; i < dim; i++) {
    temp = 0.0;
    for (j = 0; j < dim; j++)
      temp = temp + gram_matrix[i][j] * vector[j];
    *norm = *norm + temp * vector[i];
  }
}


void doublePrintVector(int dim,
		       double *vector)
{
  int
    i;
  
  for (i = 0; i < dim; i++)
    printf("% 3.10f", vector[i]);
  printf("\n");
}


int doubleIsMultipleVector(int dim,
			   double *vector1,
			   double *vector2)
{
  int
    i,
    j,
    k;

  i = 0;
  while (vector1[i] == 0.0 && i < dim)
    i++;

  j = 0;
  while (vector2[j] == 0.0 && j < dim)
    j++;
  
  if (i != j) return 0;

  k = 0;
  while (fabs(vector2[k]*vector1[i] - vector1[k]*vector2[i]) < LU_EPSILON && 
	 k < dim)
    k++;

  if (k != dim) return 0;
  else return 1;

}


/* =============================
 * = functions for int-vectors =
 * =============================
 */


void intMakeVector(int dim,
		   int **vector)
{
  *vector = (int *) laMalloc(dim * sizeof(int));
}


void intDestroyVector(int dim,
		      int **vector)
{
  laFree(*vector);
  *vector = NULL;
}


void intCopyVector(int dim,
		   int *from,
		   int *to)
{
  memcpy(to, from, dim * sizeof(int));
}


void intPrintVector(int dim,
		    int *vector)
{
  int
    i;
  
  for (i = 0; i < dim; i++)
    printf(" %2d", vector[i]);
  printf("\n");
}


void intInnerProductVector(int dim,
			   int *vector1,
			   int *vector2,
			   int *product)
{
  int
    i;

  *product = 0;

  for (i = 0; i < dim; i++) {
    *product = *product + vector1[i] * vector2[i];
  }
}


void intVectorNorm(int dim,
		   double **gram_matrix,
		   int *vector,
		   double *norm)
{
  int
    i,
    j;
  double
    temp;

  *norm = 0.0;

  for (i = 0; i < dim; i++) {
    temp = 0.0;
    for (j = 0; j < dim; j++)
      temp = temp + gram_matrix[i][j] * (double) vector[j];
    *norm = *norm + temp * (double) vector[i];
  }
}



/* ==============================
 * = functions for int-matrices =
 * ==============================
 */

void intMakeSquareMatrix(int dim,
			 int ***matrix)
{
  int
    i;

  *matrix = (int **) laMalloc(dim * sizeof(int *));
  (*matrix)[0] = (int *) laMalloc(dim * dim * sizeof(int));
  for (i = 1; i < dim; i++)
    (*matrix)[i] = (*matrix)[0] + i * dim;
}


void intDestroySquareMatrix(int dim,
			    int ***matrix)
{
  laFree((*matrix)[0]);
  laFree(*matrix);
  *matrix = NULL;
}


void intCopySquareMatrix(int dim,
			 int **matrix1,
			 int **matrix2)
{
  int
    i;

  for (i = 0; i < dim; i++)
    memcpy(matrix2[i], matrix1[i], dim * sizeof(int));
}


void intMultiplySquareMatrixVector(int dim,
				   int **matrix,
				   int *vec1,
				   int *vec2)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++) {
    vec2[i] = 0;
    for (j = 0; j < dim; j++) 
      vec2[i] = vec2[i] + matrix[i][j] * vec1[j];
  }
}


void intPrintGramMatrix(int dim,
			int **gram_matrix)
{
  int
    i,
    j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j <= i; j++) {
      printf("% 3d", gram_matrix[i][j]);
    }
    printf("\n");
  }
}


void intPrintSquareMatrix(int dim,
			  int **matrix)
{
  int
    i,
    j;

  printf("\n");
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      printf(" %2d ", matrix[i][j]);
    printf("\n");
  }
}
