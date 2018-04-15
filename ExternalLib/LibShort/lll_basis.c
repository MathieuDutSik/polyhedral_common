/* lll_basis.c  LLL reduction of a lattice basis              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

#include "lll_basis.h"
#include "la_support.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define LLL_EPSILON 0.00000000001
#define LOVASZ_CONSTANT 0.99999


/*
 * ======================
 * = internal functions =
 * ======================
 */

/*
 * red
 * ===
 */

static void red(long dim,
		double **reduced_basis,
		double **mu,
		long k,
		long l)
{
  double
    q;
  long
    i;

  if (fabs(mu[k][l]) > 0.5) {
    q = floor(0.5 + mu[k][l]);
    doubleRowAddSquareMatrix(dim, reduced_basis, -q, l, k);
    mu[k][l] = mu[k][l] - q;
    for (i = 0; i <= l - 1; i++)
      mu[k][i] = mu[k][i] - q * mu[l][i];
  } 
}


/*
 * swap
 * ====
 */

static void swap(long dim,
		 double **reduced_basis,
		 double **o_basis,
		 double **mu,
		 double *B,
		 double *b0,
		 long k,
		 long k_max)
{
  long
    i,
    j;
  double
    mu0,
    B0,
    t;

  doubleRowChangeSquareMatrix(dim, reduced_basis, k, k - 1);

  for (j = 0; j <= k - 2; j++) {
    t = mu[k][j];
    mu[k][j] = mu[k - 1][j];
    mu[k - 1][j] = t;
  }

  mu0 = mu[k][k - 1];
  B0 = B[k] + mu0 * mu0 * B[k - 1];
  mu[k][k - 1] = mu0 * B[k - 1] / B0;

  doubleCopyVector(dim, o_basis[k - 1], b0);
  
  doubleCopyVector(dim, o_basis[k], o_basis[k - 1]);
  doubleAddVector(dim, mu0, b0, o_basis[k - 1]);

  doubleRowAddSquareMatrix(dim, o_basis, -mu[k][k - 1] - 1, k, k);
  doubleAddVector(dim, B[k]/B0, b0, o_basis[k]);

  B[k] = B[k - 1] * B[k] / B0;
  B[k - 1] = B0;
  for (i = k + 1; i <= k_max; i++) {
    t = mu[i][k];
    mu[i][k] = mu[i][k - 1] - mu0 * t;
    mu[i][k - 1] = t + mu[k][k - 1] * mu[i][k];
  }

}
  


/* ======================
 * = exported functions =
 * ======================
 */

/*
 * computeLLLReducedBasis
 * ======================
 */

void computeLLLReducedBasis(long dim,
			    double **basis,
			    double **reduced_basis,
			    int check)
{
  long
    i,
    j,
    k,
    k_max,
    l;
  double
    prod,
    prod1,
    prod2,
    *b0,
    *B,
    **mu,
    **o_basis;

  doubleMakeSquareMatrix(dim, &mu);
  doubleMakeSquareMatrix(dim, &o_basis);
  doubleMakeVector(dim, &B);
  doubleMakeVector(dim, &b0);
  
  doubleCopySquareMatrix(dim, basis, reduced_basis);

  /*
   * 1. initialize
   */

  k = 1;
  k_max = 0;

  doubleCopyVector(dim, reduced_basis[0], o_basis[0]);

  doubleInnerProduct(dim, reduced_basis[0], reduced_basis[0], &(B[0]));

  do {

    /*
     * 2. incremental Gram-Schmidt
     */

    if (k > k_max) {
      k_max = k;
      doubleCopyVector(dim, reduced_basis[k], o_basis[k]);

      for (j = 0; j <= k - 1; j++) {
	doubleInnerProduct(dim, reduced_basis[k], o_basis[j], &prod); 
	mu[k][j] = 1/B[j] * prod;
	doubleRowAddSquareMatrix(dim, o_basis, -mu[k][j], j, k);
      }
      
      doubleInnerProduct(dim, o_basis[k], o_basis[k], &(B[k]));

      if (B[k] < LLL_EPSILON) {
	fprintf(stderr, "computeLLLReducedBasis: linear dependent vectors!\n");
	abort();
      }
    }

    /*
     * 3. test LLL condition
     */

    red(dim, reduced_basis, mu, k, k - 1);
    if (B[k] < (LOVASZ_CONSTANT - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1]) {
      swap(dim, reduced_basis, o_basis, mu, B, b0, k, k_max);
      if (k > 1) k --;
    }
    else {
      for (l = k - 2; l >= 0; l--) {
	red(dim, reduced_basis, mu, k, l);
      }
      k++;
    }

  } while (k < dim);
    
  if (check) {
    for (j = 0; j < dim; j++)
      for (i = j + 1; i < dim; i++)
	if (fabs(mu[i][j]) > 0.5) {
	  fprintf(stderr, "lll_basis: size condidition not valid!\n");
	  abort();
	}
    for (i = 1; i < dim; i++) {
      doubleInnerProduct(dim, o_basis[i], o_basis[i], &prod1);
      doubleInnerProduct(dim, o_basis[i - 1], o_basis[i - 1], &prod2);
      if (prod1 < (LOVASZ_CONSTANT - mu[i][i - 1] * mu[i][i - 1]) *
	  prod2)
	{
	  fprintf(stderr, "lll_basis.c: lovasz condition not valid!\n");
	  abort();
	}
    }
  }

  doubleDestroyVector(dim, &b0);
  doubleDestroyVector(dim, &B);
  doubleDestroySquareMatrix(dim, &o_basis);
  doubleDestroySquareMatrix(dim, &mu);
}
