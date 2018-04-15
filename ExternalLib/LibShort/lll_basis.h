#ifndef LLL_BASIS_H
#define LLL_BASIS_H

#ifdef __cplusplus
extern "C" {
#endif

/* lll_basis.h  LLL reduction of a lattice basis              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

/*
 * computeLLLReducedBasis
 * ======================
 *
 * description:
 * ------------
 * computes a LLL-reduced basis using the standard inner product
 *
 * input:
 * ------
 * long dim               - dimension of the lattice
 * double **basis         - the vectors of the basis as a dim x dim-matrix
 * double **reduced_basis - allocated dim x dim-matrix
 * int check              - if check == 1, then the LLL-reduced conditions are
 *                          checked a posteriori
 *
 * output:
 * -------
 * double **reduced_basis - a LLL-reduced basis of the lattice
 */

extern void computeLLLReducedBasis(long dim,
				   double **basis,
				   double **reduced_basis,
				   int check);
#ifdef __cplusplus
}
#endif

#endif
