#ifndef LA_SUPPORT_H
#define LA_SUPPORT_H

#ifdef __cplusplus
extern "C" {
#endif

/* la_support.h  low-level functions for vectors and matrices */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

/* ==============================
 * = number theoretic functions =
 * ==============================
 */

extern int factorial(int n);

/* =================================
 * = functions for double-matrices =
 * =================================
 */

extern void doubleMakeSquareMatrix(int dim,
				   double ***matrix);
extern void doubleDestroySquareMatrix(int dim,
				      double ***matrix);
extern void doubleCopySquareMatrix(int dim,
				   double **from,
				   double **to);
extern void doubleAddSquareMatrix(int dim,
				  double factor,
				  double **from,
				  double **to);
extern void doubleSetIdentitySquareMatrix(int dim,
			      double **matrix);
extern void doubleMultiplySquareMatrixVector(int dim,
					 double **matrix,
					 double *vec1,
					 double *vec2);
extern void doubleRowAddSquareMatrix(int dim,
			 double **matrix,
			 double factor,
			 int from,
			 int to);
extern void doubleRowChangeSquareMatrix(int dim,
					double **matrix,
					int index1,
					int index2);
extern int doubleComputeLUSquareMatrix(int dim,
				       double **mat,
				       double **l_mat,
				       double **u_mat,
				       double **p_mat);
extern void doubleTransposeSquareMatrix(int dim,
					double **mat,
					double **t_mat);
extern void doubleSolveLUSquareMatrix(int dim,
		    double **l_mat,
		    double **u_mat,
		    double **p_mat,
		    double *b,
		    double *x);
extern int doubleInverseSquareMatrix(int dim,
				     double **mat,
				     double **inv_mat);
extern void doubleInverseLowerSquareMatrix(int dim,
			       double **l_mat,
			       double **inv_l_mat);
extern void doubleInverseUpperSquareMatrix(int dim,
			       double **u_mat,
			       double **inv_u_mat);
extern void doubleAbsoluteDeterminantSquareMatrix(int dim,
			    double **l_mat,
			    double **u_mat,
			    double *det);
extern int doubleComputeCholeskySquareMatrix(int dim,
					     double **mat,
					     double **u_mat);
extern void doubleComputeSymmetrizationSquareMatrix(int dim,
						    double **mat,
						    double **sym_mat);
extern void doubleMultiplySquareMatrix(int dim,
				   double **a,
				   double **b,
				   double **c);
extern void doubleComputeRowNormSquareMatrix(int dim,
		    double **mat,
		    int row,
		    double *row_norm);
extern void doublePermuteColumnsSquareMatrix(int dim,
					     double **mat,
					     int index1,
					     int index2);
extern void doublePermuteRowsSquareMatrix(int dim,
		    double **mat,
		    int index1,
		    int index2);
extern void doublePrintSquareMatrix(int dim,
				    double **mat);
extern void doublePrintGramMatrix(int dim,
				  double **gram_matrix);
extern void double2intSquareMatrix(int dim,
				    double **d_mat,
				    int **l_mat);



/* ================================
 * = functions for double-vectors =
 * ================================
 */

extern void doubleMakeVector(int dim,
			     double **vector);
extern void doubleDestroyVector(int dim,
				double **vector);
extern void doubleCopyVector(int dim,
			     double *from,
			     double *to);
extern void doubleMultiplyVector(int dim,
			    double factor,
			    double *from,
			    double *to);
extern void doubleAddVector(int dim,
			    double factor,
			    double *from,
			    double *to);
extern void doubleInnerProduct(int dim,
			       double *vector1,
			       double *vector2,
			       double *product);
extern void doubleVectorNorm(int dim,
			     double **gram_matrix,
			     double *vector,
			     double *norm);
extern void doublePrintVector(int dim,
			      double *vector);
extern int doubleIsMultipleVector(int dim,
				  double *vector1,
				  double *vector2);


/* =============================
 * = functions for int-vectors =
 * =============================
 */

extern void intMakeVector(int dim,
			  int **vector);
extern void intDestroyVector(int dim,
			     int **vector);
extern void intCopyVector(int dim,
			  int *from,
			  int *to);
extern void intPrintVector(int dim,
			   int *vector);
extern void intInnerProductVector(int dim,
				  int *vector1,
				  int *vector2,
				  int *product);
extern void intVectorNorm(int dim,
			  double **gram_matrix,
			  int *vector,
			  double *norm);




/* ==============================
 * = functions for int-matrices =
 * ==============================
 */

extern void intMakeSquareMatrix(int dim,
				int ***matrix);
extern void intDestroySquareMatrix(int dim,
				   int ***matrix);
extern void intCopySquareMatrix(int dim,
				int **matrix1,
				int **matrix2);
extern void intMultiplySquareMatrixVector(int dim,
				      int **matrix,
				      int *vec1,
				      int *vec2);
extern void intPrintGramMatrix(int dim,
			       int **gram_matrix);
extern void intPrintSquareMatrix(int dim,
				 int **matrix);
#ifdef __cplusplus
}
#endif

#endif
