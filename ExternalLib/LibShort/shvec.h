#ifndef SHVEC_H
#define SHVEC_H

#ifdef __cplusplus
extern "C" {
#endif

/* shvec.h  computing short and close vectors in lattices     */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

/*
 * description:
 * ============
 * The function "computeShortVectors" computes short vectors of a lattice L
 * which is given by a Gram-matrix. There are 4 different modes.
 * - (SHVEC_MODE_BOUND)
 *   compute all vectors of L which length is bounded by a given constant
 * - (SHVEC_MODE_SHORTEST_VECTORS)
 *   compute the minimal vectors of L
 * - (SHVEC_MODE_MINIMUM)
 *   compute only the minimum of L
 * - (SHVEC_MODE_THETA_SERIES)
 *   count only the number of vectors of L which length is bounded by a
 *   given constant
 * Additionally there are many options.
 * 
 * The algorithm used is based on the algorithm of Fincke and Pohst.
 * (see Henri Cohen --  A Course in Computational Algebraic Number Theory)
 */


/*
 * internal epsilon
 * ================
 */

#define SHVEC_ELLIPSOID_EPSILON 0.0000001
#define SHVEC_STEP_EPSILON 0.000001



/*
 * modes
 * =====
 */

#define SHVEC_MODE_BOUND             0
#define SHVEC_MODE_SHORTEST_VECTORS  1
#define SHVEC_MODE_MINIMUM           2
#define SHVEC_MODE_THETA_SERIES      3
#define SHVEC_MODE_VINBERG           4



/*
 * structures
 * ==========
 */

struct internals_struct {
  int dim;
  double **gram_matrix;
  int integral;
  int ComputeTheta;
  int (*cholesky_decomposition)();
  double **quadratic_form;
  int **transformation;
  double **inv_transformation;
  int *sigma;
  int *dummy_vector;
  int even;
  int check;
};

typedef struct internals_struct internals;


struct shvec_request_struct {
  int mode;
  double bound;
  int number;
  double *coset;
  internals ints;
};

typedef struct shvec_request_struct *shvec_request;


struct shvec_info_struct {
  struct shvec_request_struct request;
  int short_vectors_count;
  int short_vectors_number;
  int **short_vectors;
  double *short_vectors_norm;
  double minimum;
  double theta_series_bound;
  int *theta_series;
  int **shells;
};

typedef struct shvec_info_struct *shvec_info;



/*
 * functions
 * ==========
 */

extern void initShvecReq(int dim,
			 double **gram_matrix,
			 int (*cholesky_decomposition)(),
			 int check,
			 shvec_request request);

extern void destroyShvecReq(shvec_request request);

extern void initShvecInfo(shvec_info info);

extern void destroyShvecInfo(shvec_info info);

extern void computeShvec(shvec_request request,
			 shvec_info info);

#ifdef __cplusplus
}
#endif

#endif
