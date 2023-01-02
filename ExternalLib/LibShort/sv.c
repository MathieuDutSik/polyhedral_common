/* sv.c  simple driver for shvec                              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */
/* Further correction by Mathieu Dutour Sikiric mathieu.dutour@gmail.com */

#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "la_support.h"
#include "shvec.h"

void die_sv(char *last_words) {
  printf("sv.c: ");
  printf(last_words);
  exit(EXIT_FAILURE);
}


void print_vectors(int dim, shvec_info info) {
  int i;
  printf("%d\n", info->short_vectors_number);
  for (i = 0; i < info->short_vectors_number; i++) {
    printf("[%f]: ", info->short_vectors_norm[i]);
    intPrintVector(dim, info->short_vectors[i]);
  }
}


int main(int argc, char *argv[]) {
  int check, gram, i, j, mode, number, dim, coset;
  double bound, **gram_matrix;
  shvec_request request;
  shvec_info info;

  /* default values */

  bound = 0.0;
  number = 0;
  mode = SHVEC_MODE_BOUND;
  check = 0;
  coset = 0;
  gram = 0;

  /*  user options */

  int c;

  while ((c = getopt(argc, argv, "hb:s:v:t:mMcge")) != -1)
    switch (c) {
    case 'h':
      printf("Usage: sv [options] <file\n");
      printf("-h  show this help\n");
      printf("-bN compute vectors v with (v, v) <= N\n");
      printf("-sN print maximal N vectors\n");
      printf("-vN determine the vectors with (v-c, v-c) <= N\n");
      printf("-tN compute the first N coefficients of the theta-series\n");
      printf("-m  determine the minimum\n");
      printf("-M  compute minimal vectors\n");
      printf("-c  find shortest vectors in a coset\n");
      printf("-g  print gram matrix\n");
      printf("-e  do some additional checks\n");
      return 0;
      break;
    case 'b':
      mode = SHVEC_MODE_BOUND;
      bound = atof(optarg);
      break;
    case 's':
      number = atoi(optarg);
      break;
    case 'M':
      mode = SHVEC_MODE_SHORTEST_VECTORS;
      break;
    case 'm':
      mode = SHVEC_MODE_MINIMUM;
      break;
    case 't':
      mode = SHVEC_MODE_THETA_SERIES;
      bound = atof(optarg);
      break;
    case 'v':
      mode = SHVEC_MODE_VINBERG;
      coset = 1;
      bound = atof(optarg);
      break;
    case 'c':
      coset = 1;
      break;
    case 'g':
      gram = 1;
      break;
    case 'e':
      check = 1;
      break;
    default:
      die_sv("invalid option\nTry 'sv -h' for more information.\n");
      return 1;
    }

  /* input gram matrix */

  scanf("%d", &dim);
  doubleMakeSquareMatrix(dim, &gram_matrix);
  for (i = 0; i < dim; i++)
    for (j = 0; j <= i; j++) {
      scanf("%lf", &gram_matrix[i][j]);
      gram_matrix[j][i] = gram_matrix[i][j];
    }
  if (gram)
    doublePrintGramMatrix(dim, gram_matrix);
  /* initialization of request and info */
  request = (shvec_request)malloc(sizeof(struct shvec_request_struct));
  if (request == NULL)
    die_sv("request == NULL.\n");
  if (mode == SHVEC_MODE_THETA_SERIES)
    request->ints.ComputeTheta = 1;
  else
    request->ints.ComputeTheta = 0;

  initShvecReq(dim, gram_matrix, NULL, check, request);

  request->bound = bound;
  request->mode = mode;
  request->number = number;

  if (coset) {
    for (i = 0; i < dim; i++)
      scanf("%lf", &(request->coset[i]));
  }

  info = (shvec_info)malloc(sizeof(struct shvec_info_struct));
  if (info == NULL)
    die_sv("info == NULL.\n");

  initShvecInfo(info);

  /* compute short vectors */

  computeShvec(request, info);

  /* output results */

  if (mode == SHVEC_MODE_BOUND) {
    print_vectors(dim, info);
  }

  else if (mode == SHVEC_MODE_SHORTEST_VECTORS) {
    print_vectors(dim, info);
  }

  else if (mode == SHVEC_MODE_MINIMUM) {
    printf("Minimum: %f\n", info->minimum);
  }

  else if (mode == SHVEC_MODE_THETA_SERIES) {
    for (i = 0; i < info->theta_series_bound; i++) {
      if (i == 0)
        printf("%d + ", info->theta_series[i]);
      else if (i % 2 == 0 && info->theta_series[i] != 0)
        printf("%d q^{%d} + ", info->theta_series[i], i / 2);
      else if (i % 2 == 1 && info->theta_series[i] != 0)
        printf("%d q^{%d.5} + ", info->theta_series[i], i / 2);
    }
    i = info->theta_series_bound;
    if (i % 2 == 0)
      printf("%d q^{%d}\n", info->theta_series[i], i / 2);
    else
      printf("%d q^{%d.5}\n", info->theta_series[i], i / 2);
  } else if (mode == SHVEC_MODE_VINBERG) {
    printf("%d\n", info->short_vectors_number);
    for (i = 0; i < info->short_vectors_number; i++) {
      printf("[%f]: ", info->short_vectors_norm[i]);
      intPrintVector(dim, info->short_vectors[i]);
    }
  }

  /* clean up */

  destroyShvecInfo(info);
  destroyShvecReq(request);
  doubleDestroySquareMatrix(dim, &gram_matrix);
  free(request);
  free(info);

  return 0;
}
