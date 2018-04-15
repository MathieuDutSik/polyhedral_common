/* finddel.c                                                      */
/* Version September 20, 2005                                     */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com     */


#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

#include "shvec.h"
#include "la_support.h"
#include "lrslib.h"

void die_finddel(char* last_words)
{
  printf("finddel.c: ");
  printf(last_words);
  exit(EXIT_FAILURE);
}


void solve_lp(int dim,
	      int* objective_function,
	      int* inequalities,
	      int inequalities_count,
	      double* vertex)
{
  int
    row,
    col;
  long 
    num[dim+1],
    den[dim+1];
  double
    det;
  lrs_dic
    *P;
  lrs_dat
    *Q;

  if (!lrs_init ("*solvelp. "))
    die_finddel("!lrs_init\n");

  Q = lrs_alloc_dat("LRS globals");
  if (Q == NULL)
    die_finddel("Q == NULL.\n");

  Q->m = inequalities_count;
  Q->n = dim+1;
  Q->verbose = FALSE;
  Q->lponly = TRUE;
  P = lrs_alloc_dic(Q);
  if (P == NULL) 
    die_finddel("P == NULL.\n");

  for (col = 0; col < dim+1; col++)
    den[col] = 1L;

  for (row = 1; row <= inequalities_count; row++) {
    for (col = 0; col < dim + 1; col++) {
      num[col] = (long) inequalities[(row - 1)*(dim + 1) + col];
    }
    lrs_set_row(P, Q, row, num, den, GE);
  }

  for(col = 0; col < dim+1; col++)
    num[col] = (long) objective_function[col];
  lrs_set_obj(P, Q, num, den, MAXIMIZE);

  itomp(ONE, Q->Lcm[0]); /* circumvent a bug in lrs */

  lrs_solve_lp(P,Q);
  
  lrs_mp_vector output = lrs_alloc_mp_vector(Q->n);
  for (col = 0; col < Q->n; col++)
    lrs_getsolution(P, Q, output, col);
  
  det = mptodouble(output[0]);
  for (col = 1; col < dim + 1; col++) {
    vertex[col - 1] = mptodouble(output[col]);
    vertex[col - 1] = vertex[col - 1] / det;
  }

  lrs_clear_mp_vector(output, Q->n);

  lrs_free_dic(P, Q);
  lrs_free_dat(Q);
  lrs_close("finished. ");
} 



void find_random_delone_polytope(int dim,
				 shvec_request request,
				 shvec_info info,
				 double **gram,
				 int print_radii) 
{
  int
    i,
    j,
    k,
    is_finished,
    finished,
    *rand_direction,
    inequalities_count,
    *inequalities;
  double 
    *vertex1,
    *vertex2,
    **inv_gram,
    norm,
    radius;

  finished = 0;

  intMakeVector(dim+1, &rand_direction);
  doubleMakeVector(dim, &vertex1);
  doubleMakeVector(dim, &vertex2);

  doubleMakeSquareMatrix(dim, &inv_gram);
  doubleInverseSquareMatrix(dim, gram, inv_gram);

  /* initialize list_facets with standard basis vectors */

  inequalities_count = 2*dim;  
  inequalities = (int *) malloc((dim + 1) * inequalities_count * sizeof(int));
  if (inequalities == NULL)
    die_finddel("inequalities == NULL.\n");
  
  for (i = 0; i < dim; i++) {
    inequalities[(dim+1)*i] = (int) gram[i][i];
    for (j = 1; j < dim+1; j++)
      if (i == j - 1) inequalities[(dim+1)*i + j] = 1;
    else inequalities[(dim+1)*i + j] = 0;
  }

  for (i = 0; i < dim; i++) {
    inequalities[(dim+1)*(dim+i)] = (int) gram[i][i];
    for (j = 1; j < dim+1; j++)
      if (i == j - 1) inequalities[(dim+1)*(dim+i) + j] = -1;
    else inequalities[(dim+1)*(dim+i) + j] = 0;
  }


  /* choose a random direction */

  for (i = 1; i < dim+1; i++) {
    rand_direction[i] = mrand48();
  }

  do {
    solve_lp(dim, rand_direction, inequalities, inequalities_count, vertex1);
    doubleMultiplySquareMatrixVector(dim, inv_gram, vertex1, vertex2);
    doubleAddVector(dim, -0.5, vertex2, vertex2);
    
    doubleVectorNorm(dim, gram, vertex2, &radius);
    
    request->mode = SHVEC_MODE_SHORTEST_VECTORS;
    request->bound = 0;
    request->number = 0;
  
    for (i = 0; i < dim; i++) {
      vertex2[i] = -vertex2[i];    
      request->coset[i] = vertex2[i];
    }
    
    computeShvec(request, info);

    for (k = 0; k < info->short_vectors_number; k++) {
      inequalities = (int *) realloc(inequalities, (dim+1) * (inequalities_count+2) * sizeof(int));
      if (inequalities == NULL)
	die_finddel("inequalities == NULL.\n");
      
      /* recycling: vertex2 */
      for (j = 0; j < dim; j++) {
	inequalities[(dim + 1)*inequalities_count + j + 1] = 
	  info->short_vectors[k][j];
	inequalities[(dim + 1)*(inequalities_count + 1) + j + 1] = 
	  -info->short_vectors[k][j];
	vertex2[j] = (double) info->short_vectors[k][j];
      }
      doubleVectorNorm(dim, gram, vertex2, &norm);
      inequalities[(dim + 1)*inequalities_count] = (int) norm;
      inequalities[(dim + 1)*(inequalities_count + 1)] = (int) norm;
      inequalities_count = inequalities_count + 2;

      /* we are finished if 0 is among the short_vectors */
      is_finished = 1;
      for (j = 0; j < dim; j++)
	if (info->short_vectors[k][j] != 0)
	  is_finished = 0;
      if (!finished) finished = is_finished;
    }
  } while (!finished);

  printf("%d\n", info->short_vectors_number);
  for (k = 0; k < info->short_vectors_number; k++) {
    for (j = 0; j < dim; j++) {
      printf("%d ", info->short_vectors[k][j]);
    }
    printf("\n");
  }
  if (print_radii)
    printf("squared radius = %3.10f\n", radius);

  free(inequalities);
  intDestroyVector(dim+1, &rand_direction);
  doubleDestroyVector(dim, &vertex1);
  doubleDestroyVector(dim, &vertex2);
  doubleDestroySquareMatrix(dim, &inv_gram);
}


int main(int argc, char **argv)
{

  double
    **gram_matrix;
  int
    tries,
    seed,
    dim,
    print_radii,
    c,
    i, 
    j,
    t;
  
  print_radii = 0;
  tries = 1;
  seed = 1;

  while ((c = getopt (argc, argv, "hs:t:r")) != -1)
    switch (c)
      {
      case 'h':
	printf("Usage: finddel [options] <file\n");
	printf("Options:\n");
	printf(" -h    show this help\n");
	printf(" -s N  set seed to N\n");
	printf(" -t N  try N times\n");
	printf(" -r    print squared circumradii\n");
	return 0;
	break;
      case 's':
	seed = atoi(optarg);
	break;
      case 't':
	tries = atoi(optarg);
	break;
      case 'r':
	print_radii = 1;
	break;
      default:
	printf("finddel: invalid option\n");
	printf("Try 'finddel -h' for more information.\n");
	return 1;
      }

  shvec_request request;
  shvec_info info;

  srand48(seed);
  
  /* input gram matrix */

  scanf("%d", &dim);

  doubleMakeSquareMatrix(dim, &gram_matrix);

  for (i = 0; i < dim; i++)
    for (j = 0; j <= i; j++) {
      scanf("%lf", &gram_matrix[i][j]);
      gram_matrix[j][i] = gram_matrix[i][j];
    }


  /* allocate shvec_structs */
  
  request = (shvec_request) malloc(sizeof(struct shvec_request_struct));
  if (request == NULL)
    die_finddel("request == NULL.\n");
  
  info = (shvec_info) malloc(sizeof(struct shvec_info_struct));
  if (info == NULL)
    die_finddel("info == NULL.\n");
  
  initShvecReq(dim, gram_matrix, NULL, 0, request);
  initShvecInfo(info);

  for (t = 0; t < tries; t++) {
    find_random_delone_polytope(dim, request, info, gram_matrix, print_radii);
    printf("\n");
  }
   
  /* clean */
    
  destroyShvecReq(request);
  destroyShvecInfo(info);
  doubleDestroySquareMatrix(dim, &gram_matrix);
  return 0;
}

