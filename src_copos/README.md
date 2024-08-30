Copositive
==========

This C++11 code allows to compute the copositive minimum of a non-necessarily
positive definite matrix A. IT relies on decomposing cone into smaller ones
and it should be faster than existing matlab codes.

Overview
--------

A symmetric nxn-matrix A is called copositive if for all vector x in R_+^n
we have A[x] >= 0. It is called strictly copositive is for all non-zero vectors
x in R_+^n we have A[x] > 0.

It is called completely positive if there exist vectors v_1, ....., v_M such
that A = sum_i p(v_i) with p(v_i) the corresponding rank 1 matrix to v_i.

The copositive minimum of a strictly copositive matrix A is the minimum over
all non-zero vectors v in Z_+^n of A[x]. The set of vectors matching the bound
is denoted min_{CP}(A).

An integral completely positive expression of A is m integral vectors
v_1, ..., v_M such that A = sum_i p(v_i).


Algorithms
----------

The fundamental step is to determine the copositive set min_{CP}(A). The algorithm used is a
combination of following two:
  * The splitting of the cone {x in R^n, s.t. x_i >= 0} into a number of simplicial cone such that
  on each one of them the determination of the minimum is easy. By easy we mean that if a cone C is spanned by v_1, ..., v_n such that v_i A v_j >= 0 then the minimum over this cone is easy to determine.
  * On each of the cone, we apply the Fincke-Pohst algorithm suitably adapted for the enumeration.

Programs
--------

The following programs are available
 * **CP_ComputeCopositiveMin**: It tests is a matrix is copositive.
 * **CP_TestCompletePositivity**: We test if a matrix is completely positive.
 * **CP_CopositiveMin**: Get the list of positive vectors of small norms

References
----------

For applications of the closest vector problem the user is referred to following papers and books:
  * Berman, Abraham and Shaked-Monderer, Naomi, Completely positive matrices, World Scientific
  * Peter J. C. Dickinson, Ph.{D}. thesis: The Copositive Cone, the Completely Positive Cone and their Generalisations, University of Groningen
  * Dur Miriam, Copositive programming — a survey, Recent Advances in Optimization and its Applications in Engineering, Springer, 2010
  * Fincke, U. and Pohst, M., Improved methods for calculating vectors of short length in a lattice, including a complexity analysis, Mathematics of Computation, 44, 1985, 463--471

This work was used for the following two articles:
  * Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin, A simplex algorithm for rational CP-factorization, preprint at arxiv:1807.01382, Mathematical Programming 187 (2020) 25--45
  * Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin, Rational Factorizations of Completely Positive Matrices, preprint at arxiv:1701.03148, Linear Algebra and its Applications 523 (2017) 46--51

Dependencies
------------

Following dependencies are needed for compiling the code:

  * Eigen: http://eigen.tuxfamily.org/
  * GNU MultiPrecision Library (GMP): https://gmplib.org/
