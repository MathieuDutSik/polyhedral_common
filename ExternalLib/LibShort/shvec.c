/* shvec.c  computing short and close vectors in lattices     */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

#include "shvec.h"
#include "la_support.h"
#include "lll_basis.h"
#include "panic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*
 * ======================
 * = constants & macros =
 * ======================
 */

#define SHVEC_MAX 1000 /* allocate always SHVEC_MAX vectors */
#define STOP_COMPUTATION 666
#define ROUND(x) (int) floor(x + 0.5)
#define SHVEC_TEST_EPSILON 0.01



/*
 * ======================
 * = internal functions =
 * ======================
 */


void die_shvec(char* last_words)
{
  printf("shvec.c: ");
  printf("last_words=%s\n", last_words);
  exit(EXIT_FAILURE);
}



/*
 * internal memory management
 * ==========================
 * SHVEC_malloc, SHVEC_realloc, SHVEC_free
 */

static void *SHVEC_malloc(size_t size)
{
  void
    *ptr;

  ptr = (void *) calloc(1, size);
  if (ptr == NULL) {
    PANIC;
  }
  return ptr;
}

static void *SHVEC_realloc(void *ptr, size_t size)
{
  void
    *new_ptr;
  
  new_ptr = (void *) realloc(ptr, size);
  if (new_ptr == NULL) {
    PANIC;
  }
  return new_ptr;
}

static void SHVEC_free(void *ptr)
{
  if (ptr == NULL) {
    PANIC;
  }
  free(ptr);
}


/*
 * insertBound
 * ===========
 * don't declare static
 */

int insertBound(shvec_info info,
		int *vector,
		int *c,
		int coset,
		double norm)
{
  int *sigma, number, i, vector_norm=-1, shell_bound, dim;

  dim = info->request.ints.dim;
  number = info->short_vectors_number;

  if (info->request.ints.integral && !coset) {
    vector_norm = ROUND(norm);
    norm = (int) vector_norm;
  }

  sigma = info->request.ints.sigma;

  /*
   * increase arrays if necessary
   */

  if (number % SHVEC_MAX == 0) {
    info->short_vectors =
      SHVEC_realloc(info->short_vectors,
		    SHVEC_MAX * (number/SHVEC_MAX + 1) * sizeof(int *));
    info->short_vectors[number] =
      SHVEC_malloc(dim * sizeof(int) * SHVEC_MAX);
    for (i = 1; i < SHVEC_MAX; i++) {
      info->short_vectors[i + number] =
	info->short_vectors[number] + i * dim;
    }
    info->short_vectors_norm =
      SHVEC_realloc(info->short_vectors_norm,
		    SHVEC_MAX * (number/SHVEC_MAX + 1) * sizeof(double));
  }

  /* 
   * theta_series and shells only for integral lattices without coset option
   */
  /*  fprintf(stderr, "ComputeTheta=%d\n", info->request.ints.ComputeTheta); */
  if (info->request.ints.integral && info->request.ints.ComputeTheta && !coset) {
    fprintf(stderr, "theta_series_bound=%f vector_norm=%d\n", info->theta_series_bound, vector_norm);
    if (info->theta_series_bound < vector_norm) {
      fprintf(stderr, "reallocating\n");
      info->theta_series =
	SHVEC_realloc(info->theta_series, (vector_norm + 1) * sizeof(int));
      info->shells =
	SHVEC_realloc(info->shells, (vector_norm + 1) * sizeof(int *));
      for (i = info->theta_series_bound + 1; i <= vector_norm; i++) {
	info->theta_series[i] = 0;
	info->shells[i] = NULL;
      }
      info->theta_series_bound = vector_norm;
    }
    
    shell_bound = info->theta_series[vector_norm]/2;
    fprintf(stderr, "Just after assignation shell_bound=%d\n", shell_bound);
    if (shell_bound % SHVEC_MAX == 0) {
      info->shells[vector_norm] =
	SHVEC_realloc(info->shells[vector_norm],
		      SHVEC_MAX * (shell_bound / SHVEC_MAX + 1) * sizeof(int));
    }
  }

  /*
   * insert values
   */

  for (i = 0; i < dim; i++)
    (info->request.ints.dummy_vector)[i] = vector[sigma[i]] - c[sigma[i]];

  intMultiplySquareMatrixVector(dim,
				info->request.ints.transformation,
				info->request.ints.dummy_vector,
				info->short_vectors[number]);
  /*
   * coset
   */

  info->short_vectors_norm[number] = norm;
  
  if (info->request.ints.integral && info->request.ints.ComputeTheta && !coset) {
    info->theta_series[vector_norm] = info->theta_series[vector_norm] + 2;
    fprintf(stderr, "vector_norm=%d  shell_bound=%d\n", vector_norm, shell_bound);
    info->shells[vector_norm][shell_bound] = number;
    fprintf(stderr, "After assignation\n");
  }

  /*
   * counters
   */

  info->short_vectors_count++;
  info->short_vectors_number++;

  /*
   * minimum
   */
  
  if (norm < info->minimum || info->minimum == -1.0)
    info->minimum = norm;


  if (info->short_vectors_number == info->request.number)
    return STOP_COMPUTATION;
  else
    return 0;
}


int insertThetaSeries(shvec_info info,
		      int *vector,
		      int *c,
		      int coset,
		      double norm)
{
  int vector_norm, i;
  info->short_vectors_count++;
  vector_norm = ROUND(norm);
  info->short_vectors_number++;
  if (info->theta_series_bound < (double) vector_norm) {
    info->theta_series =
      SHVEC_realloc(info->theta_series,
		    (vector_norm + 1) * sizeof(int));
    for (i = info->theta_series_bound + 1; i <= vector_norm; i++)
      info->theta_series[i] = 0;
    info->theta_series_bound = (double) vector_norm;
  }
  info->theta_series[vector_norm] = info->theta_series[vector_norm] + 2;
  if (norm < info->minimum || info->minimum == -1.0) {
    info->minimum = norm;
  }
  return 0;
}



/*
 * insertStop
 * ==========
 * don't declare static
 */

int insertStop(shvec_info info,
		      int *vector,
		      int *c,
		      int coset,
		      double norm)
{
  info->minimum = norm;
  return STOP_COMPUTATION;
}



/*
 * =================================
 * = algorithm of Fincke and Pohst =
 * =================================
 */

int computeIt(shvec_info info,
	      int (*insert)())
{
  int sigma, coset, result, not_finished, needs_new_bound, *x, *c, i, j, dim;
  double *double_vec, *C, sum, *X, bound, *T, *U, *L, **q, Z;

  result = 0;
  dim = info->request.ints.dim;
  bound = info->request.bound;

  doubleMakeVector(dim, &T);
  doubleMakeVector(dim, &U);
  doubleMakeVector(dim, &L);
  doubleMakeVector(dim, &X);
  doubleMakeVector(dim, &C);
  doubleMakeVector(dim, &double_vec);
  intMakeVector(dim, &x);
  intMakeVector(dim, &c);
  
  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info->request.coset[i] != 0.0);
    i++;
  }
  
  /* 
   * Normalize the entries of the coset-vectors to N = [-0.5, 0.5),
   * then info->request.coset[i] = C[i] + c[i] and C[i] lies in N.
   * Take care of sigma.
   */
    
  doubleMultiplySquareMatrixVector(dim, info->request.ints.inv_transformation,
				   info->request.coset, double_vec);
  
  for (i = 0; i < dim; i++) {
    sigma = info->request.ints.sigma[i];
    c[sigma] = (int) double_vec[i];
    C[sigma] = double_vec[i] - c[sigma];
    if (C[sigma] >= 0.5) {
      C[sigma] = C[sigma] - 1.0;
      c[sigma]++;
    }
    else if (C[sigma] < -0.5) {
      C[sigma] = C[sigma] + 1.0;
      c[sigma]--;
    }
  }

  /*
   * Computation!
   */

  not_finished = 1;
  needs_new_bound = 1;
  i = dim - 1;
  T[i] = bound;
  U[i] = 0.0;
  q = info->request.ints.quadratic_form;
  
  while (not_finished) {

    if (needs_new_bound) {
      Z = sqrt((T[i] + SHVEC_ELLIPSOID_EPSILON) / q[i][i]);
      L[i] = floor(Z - U[i] - C[i]);
      X[i] = ceil(-Z - U[i] - C[i]) - 1.0;
      x[i] = (int) X[i];
      needs_new_bound = 0;
    }

    x[i]++;
    X[i] = X[i] + 1.0;
    
    if (X[i] <= L[i]) {
	
      if (i == 0) {

	if (!coset) {
	  j = dim - 1;
	  not_finished = 0;
	  while (j >= 0 && !not_finished) {
	    not_finished = (x[j] != 0);
	    j--;
	  }
	  if (!not_finished) goto urgent_stop;
	}
	result = insert(info, x, c, coset,
			bound - T[0] + q[0][0] * (X[0] + C[0] + U[0]) *
			(X[0] + C[0] + U[0]));
	if (result == STOP_COMPUTATION) goto urgent_stop;
      }
      else {
	i--;
	sum = 0.0;
	for (j = i + 1; j < dim; j++)
	  sum += q[i][j] * (X[j] + C[j]);
	U[i] = sum;
	T[i] = 
	  T[i+1] - q[i+1][i+1] * (X[i+1] + C[i+1] + U[i+1]) * (X[i+1] + C[i+1] + U[i+1]);
	needs_new_bound = 1;
      }
    }
    else {
      i++;
      if (i >= dim) not_finished = 0;
    }
  } /* while (not_finished) */
  
urgent_stop:
    
  /*
   * free memory
   */

  doubleDestroyVector(dim, &T);
  doubleDestroyVector(dim, &U);
  doubleDestroyVector(dim, &L);
  doubleDestroyVector(dim, &X);
  doubleDestroyVector(dim, &C);
  doubleDestroyVector(dim, &double_vec);
  intDestroyVector(dim, &c);
  intDestroyVector(dim, &x);

  return result;
}

static int computeMinimum(shvec_info info)
{
  int result, finished, dim, coset, i;
  double *double_vec, *C, step_size;
  dim = info->request.ints.dim;

  doubleMakeVector(dim, &C);
  doubleMakeVector(dim, &double_vec);

  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info->request.coset[i] != 0.0);
    i++;
  }


  /* 
   * Normalize the entries of the coset-vectors to N = [-0.5, 0.5),
   * then info->request.coset[i] = C[i] + c[i] and C[i] lies in N.
   */

  for (i = 0; i < dim; i++) {
    C[i] = info->request.coset[i] - (int) info->request.coset[i];
    if (C[i] >= 0.5)
      C[i] = C[i] - 1.0;
    else if (C[i] < -0.5)
      C[i] = C[i] + 1.0;
  }

  if (coset)
    doubleVectorNorm(dim, info->request.ints.gram_matrix,
		     C, &(info->minimum));
  else {
    info->minimum = info->request.ints.gram_matrix[0][0];
    for (i = 1; i < dim; i++)
      if (info->minimum > info->request.ints.gram_matrix[i][i])
	info->minimum = info->request.ints.gram_matrix[i][i];
  }

  doubleDestroyVector(dim, &C);
  doubleDestroyVector(dim, &double_vec);
  
  if (coset || !info->request.ints.integral) step_size = SHVEC_STEP_EPSILON;
  else if (info->request.ints.even) step_size = 2.0;
  else step_size = 1.0;

  finished = 0;
  while (!finished) {
    info->request.bound = info->minimum - step_size;
    result = computeIt(info, insertStop);
    if (result != STOP_COMPUTATION) finished = 1;
  }  
  info->minimum = info->request.bound + step_size;
  
  return 0;
}



/*
 * copyShvecReq
 * ============
 * copy a busy request to an empty request
 */

static void copyShvecReq(shvec_request req_from,
			 shvec_request req_to)
{
  int dim;
  dim = req_from->ints.dim;

  req_to->mode = req_from->mode;
  req_to->bound = req_from->bound;
  req_to->number = req_from->number;
  doubleMakeVector(dim, &(req_to->coset));
  doubleCopyVector(dim, req_from->coset, req_to->coset);
  req_to->ints.dim = dim;
  doubleMakeSquareMatrix(dim, &(req_to->ints.gram_matrix));
  doubleCopySquareMatrix(dim, req_from->ints.gram_matrix, 
			 req_to->ints.gram_matrix);
  req_to->ints.integral = req_from->ints.integral;
  req_to->ints.ComputeTheta = req_from->ints.ComputeTheta;
  req_to->ints.cholesky_decomposition = req_from->ints.cholesky_decomposition;
  doubleMakeSquareMatrix(dim, &(req_to->ints.quadratic_form));
  doubleCopySquareMatrix(dim, req_from->ints.quadratic_form,
			 req_to->ints.quadratic_form);
  intMakeSquareMatrix(dim, &(req_to->ints.transformation));
  intCopySquareMatrix(dim, req_from->ints.transformation,
		      req_to->ints.transformation);
  intMakeVector(dim, &(req_to->ints.sigma));
  doubleMakeSquareMatrix(dim, &(req_to->ints.inv_transformation));
  doubleCopySquareMatrix(dim, req_from->ints.inv_transformation,
			 req_to->ints.inv_transformation);
  intCopyVector(dim, req_from->ints.sigma, req_to->ints.sigma);
  intMakeVector(dim, &(req_to->ints.dummy_vector));
  req_to->ints.even = req_from->ints.even;
  req_to->ints.check = req_from->ints.check;
}  



/*
 * ======================
 * = exported functions =
 * ======================
 */


void initShvecReq(int dim,
		  double **gram_matrix,
		  int (*cholesky_decomposition)(),
		  int check,
		  shvec_request request)
{
  double
    **up_trig_mat,
    **inv_up_trig_mat,
    **double_transformation,
    **reduced_mat,
    **inv_reduced_mat,
    **q,
    *row_norm,
    temp;
  int
    success,
    i,
    j,
    t;

  if (dim < 2 || gram_matrix == NULL)
    {
      printf("shvec: (initShvecReq) wrong input!\n");
      printf("dimension too low or unassigned matrix\n");
      exit(EXIT_FAILURE);
    }

  doubleMakeVector(dim, &(request->coset));

  request->ints.dim = dim;
  request->ints.check = check;

  /*
   * test whether gram matrix is integral and
   * if yes, test if it is even
   */

  request->ints.integral = 1;
  i = 0; j = 0;
  while (i < dim && request->ints.integral) {
    while (j < dim && request->ints.integral) {
      request->ints.integral = 
	(gram_matrix[i][j] - (int) gram_matrix[i][j] == 0.0);
      j++;
    }
    i++;
  }

  if (request->ints.integral) {
    request->ints.even = 1;
    for (i = 0; i < dim; i++)
      if (ROUND(gram_matrix[i][i]) % 2 != 0)
	request->ints.even = 0;
  }
  else
    request->ints.even = 0;

  
  /* 
   * copy gram matrix to shvec_request and test if it is even
   */

  doubleMakeSquareMatrix(dim, &(request->ints.gram_matrix));
  doubleCopySquareMatrix(dim, gram_matrix, request->ints.gram_matrix);

  
  /*
   * compute the other fields
   */

  request->mode = 0;
  request->bound = 0.0;
  request->number = 0;

  if (cholesky_decomposition == NULL)
    request->ints.cholesky_decomposition = doubleComputeCholeskySquareMatrix;
  else
    request->ints.cholesky_decomposition = cholesky_decomposition;

  
  /*
   * ... some precomputations ...
   */

  doubleMakeSquareMatrix(dim, &(request->ints.quadratic_form));
  intMakeSquareMatrix(dim, &(request->ints.transformation));
  doubleMakeSquareMatrix(dim, &(request->ints.inv_transformation));
  intMakeVector(dim, &(request->ints.sigma));
  intMakeVector(dim, &(request->ints.dummy_vector));

  doubleMakeSquareMatrix(dim, &up_trig_mat);
  doubleMakeSquareMatrix(dim, &inv_up_trig_mat);
  doubleMakeSquareMatrix(dim, &double_transformation);
  doubleMakeSquareMatrix(dim, &reduced_mat);
  doubleMakeSquareMatrix(dim, &inv_reduced_mat);
  doubleMakeVector(dim, &row_norm);
  
  success = request->ints.cholesky_decomposition(dim, gram_matrix,
						 up_trig_mat);
  if (!success) 
    die_shvec("(initShvecReq) cholesky went wrong!\n");

  doubleInverseUpperSquareMatrix(dim, up_trig_mat, inv_up_trig_mat);
  computeLLLReducedBasis(dim, inv_up_trig_mat, inv_reduced_mat, check);
  doubleInverseSquareMatrix(dim, inv_reduced_mat, reduced_mat);

  doubleMultiplySquareMatrix(dim, inv_up_trig_mat, reduced_mat, 
			 double_transformation);
  double2intSquareMatrix(dim, double_transformation,
			 request->ints.transformation);
  doubleInverseSquareMatrix(dim, double_transformation, 
			    request->ints.inv_transformation);

  /*
   * sort rows and columns
   */
  
  for (i = 0; i < dim; i++) {
    doubleComputeRowNormSquareMatrix(dim, inv_reduced_mat, i, &(row_norm[i]));
    request->ints.sigma[i] = i;
  }

  for (i = 0; i < dim; i++)
    for (j = i + 1; j < dim; j++) {
      if (row_norm[i] < row_norm[j]) {
	t = request->ints.sigma[i];
	request->ints.sigma[i] = request->ints.sigma[j];
	request->ints.sigma[j] = t;
	temp = row_norm[i];
	row_norm[i] = row_norm[j];
	row_norm[j] = temp;
      }
    }

  /* DIRTY: use up_trig_mat as temporary space */

  doubleCopySquareMatrix(dim, reduced_mat, up_trig_mat);
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      reduced_mat[j][request->ints.sigma[i]] = up_trig_mat[j][i];

  doubleComputeSymmetrizationSquareMatrix(dim, reduced_mat, up_trig_mat);

  q = request->ints.quadratic_form;

  success = request->ints.cholesky_decomposition(dim, up_trig_mat, q);

  if (!success)
    die_shvec("(initShvecReq) cholesky went wrong!\n");
  
  for (i = 0; i < dim; i++)
    q[i][i] = q[i][i] * q[i][i];

  for (i = 0; i < dim; i++)
    for (j = i + 1; j < dim; j++)
      q[i][j] = q[i][j] / sqrt(q[i][i]);

  for (i = 0; i < dim; i++)
    for (j = 0; j < i; j++)
      q[i][j] = q[j][i];

  doubleDestroyVector(dim, &row_norm);
  doubleDestroySquareMatrix(dim, &up_trig_mat);
  doubleDestroySquareMatrix(dim, &inv_up_trig_mat);
  doubleDestroySquareMatrix(dim, &double_transformation);
  doubleDestroySquareMatrix(dim, &reduced_mat);
  doubleDestroySquareMatrix(dim, &inv_reduced_mat);

}



/*
 * destroyShvecReq
 * ===============
 */

void destroyShvecReq(shvec_request request)
{
  int
    dim;

  if (request == NULL)
    die_shvec("(destroyShvecReq) wrong input!\n");

  dim = request->ints.dim;

  request->mode = 0;
  request->bound = 0.0;
  request->number = 0;
  if (request->coset != NULL)
    doubleDestroyVector(dim, &(request->coset));
  request->ints.cholesky_decomposition = NULL;

  request->ints.dim = 0;
  if (request->ints.gram_matrix != NULL)
    doubleDestroySquareMatrix(dim, &(request->ints.gram_matrix));
  request->ints.integral = 0;
  if (request->ints.quadratic_form != NULL)
    doubleDestroySquareMatrix(dim, &(request->ints.quadratic_form));
  if (request->ints.transformation != NULL)
    intDestroySquareMatrix(dim, &(request->ints.transformation));
  if (request->ints.inv_transformation != NULL)
    doubleDestroySquareMatrix(dim, &(request->ints.inv_transformation));
  if (request->ints.sigma != NULL)
    intDestroyVector(dim, &(request->ints.sigma));
  if (request->ints.dummy_vector != NULL)
    intDestroyVector(dim, &(request->ints.dummy_vector));
  request->ints.even = 0;
  request->ints.check = 0;
}

void initShvecInfo(shvec_info info)
{

  info->request.mode = 0;
  info->request.bound = 0.0;
  info->request.number = 0;
  info->request.coset = NULL;

  info->request.ints.dim = 0;
  info->request.ints.gram_matrix = NULL;
  info->request.ints.integral = 0;
  info->request.ints.cholesky_decomposition = NULL;
  info->request.ints.quadratic_form = NULL;
  info->request.ints.transformation = NULL;
  info->request.ints.inv_transformation = NULL;
  info->request.ints.sigma = NULL;
  info->request.ints.dummy_vector = NULL;
  info->request.ints.even = 0;
  info->request.ints.check = 0;

  info->short_vectors_count = 0;
  info->short_vectors_number = 0;
  info->short_vectors = NULL;
  info->short_vectors_norm = NULL;
  info->minimum = -1.0;
  info->theta_series_bound = 0.0;
  info->theta_series = NULL;
  info->shells = NULL;
 
}



/*
 * destroyShvecInfo
 * ================
 */

void destroyShvecInfo(shvec_info info)
{
  int i;
  if (info == NULL)
    die_shvec("(destroyShvecInfo) wrong input!");

  info->short_vectors_count = 0;
  if (info->short_vectors != NULL) {
    int res=info->short_vectors_number % SHVEC_MAX;
    int nbBlock;
    if (res == 0) {
      nbBlock=info->short_vectors_number / SHVEC_MAX;
    }
    else {
      nbBlock=info->short_vectors_number / SHVEC_MAX + 1;
    }
    for (i = 0; i < nbBlock; i++) {
      SHVEC_free(info->short_vectors[i * SHVEC_MAX]);
    }
    SHVEC_free(info->short_vectors);
  }
  info->short_vectors = NULL;

  if (info->short_vectors_norm != NULL)
    SHVEC_free(info->short_vectors_norm);
  info->short_vectors_norm = NULL;
  
  info->short_vectors_number = 0;
  info->minimum = -1.0;

  if (info->theta_series != NULL)
    SHVEC_free(info->theta_series);
  info->theta_series = NULL;

  if (info->shells != NULL) {
    for (i = 0; i <= info->theta_series_bound; i++)
      if (info->shells[i] != NULL)
	SHVEC_free(info->shells[i]);
    SHVEC_free(info->shells);
  }
  info->shells = NULL;
  
  info->theta_series_bound = 0.0;

  destroyShvecReq(&info->request);
}


void computeShvec(shvec_request request, shvec_info info)
{
  int dim, result, coset, i;

  if (request == NULL || info == NULL) 
    die_shvec("(computeShvec) wrong input!");
  /* Not sure about that one below */
  /*  destroyShvecInfo(info); */
  
  /*
   * test the request
   */

  dim = (request->ints).dim;

  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (request->coset[i] != 0.0);
    i++;
  }
  
  switch (request->mode)
    {
    case SHVEC_MODE_BOUND:
      if (request->bound <= 0.0)
	{
	  fprintf(stderr, "bound=%f\n", request->bound);
	  die_shvec("shvec.c (computeShvec): MODE_BOUND request->bound !\n");
	}
      if (request->number < 0)
	die_shvec("shvec.c (computeShvec): MODE_BOUND request->number!\n");
      break;
    case SHVEC_MODE_SHORTEST_VECTORS:
      if (request->bound != 0.0 || request->number < 0)
	die_shvec("shvec.c (computeShvec): wrong options MODE_SHORTEST_VECTORS!\n");
      break;
    case SHVEC_MODE_MINIMUM:
      if (request->bound != 0.0 || request->number != 0)
	die_shvec("shvec.c (computeShvec): wrong options MODE_MINIMUM!\n");
      break;
    case SHVEC_MODE_THETA_SERIES:
      if (request->bound <= 0.0 || request->number != 0 || 
	  request->ints.integral == 0 || coset)
	die_shvec("shvec.c (computeShvec): wrong options MODE_THETA_SERIES!\n");
      break;
    case SHVEC_MODE_VINBERG:
      break;
    default:
      die_shvec("shvec.c (computeShvec): wrong options (default)!\n");
    }

  /*
   * copy request
   */
  
  copyShvecReq(request, &(info->request));

  if (request->mode == SHVEC_MODE_BOUND)
    result = computeIt(info, insertBound);
  else if (request->mode == SHVEC_MODE_SHORTEST_VECTORS) {
    result = computeMinimum(info);
    info->request.bound = info->minimum;
    result = computeIt(info, insertBound);
  }
  else if (request->mode == SHVEC_MODE_MINIMUM)
    {
      result = computeMinimum(info);
    }
  else if (request->mode == SHVEC_MODE_THETA_SERIES)
    {
      result = computeIt(info, insertThetaSeries);
    }
  else if (request->mode == SHVEC_MODE_VINBERG)
    {
      result = computeIt(info, insertBound);
    }
  fprintf(stderr, "result=%d\n", result);
  /*
   * 0-th coefficient of theta-series is always 1
   */
  
  if (info->theta_series != NULL) info->theta_series[0] = 1;
  if (info->shells != NULL) info->shells[0] = NULL;

  return;
}
