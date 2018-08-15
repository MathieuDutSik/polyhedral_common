SV exact
========

This C++11 code allows to compute the closest vector problem and shortest
vector problems using arbitrary precision arithmetic. This is based on the
code by Frank Vallentin which uses double arithmetic and is so faster.
The code is template based, so double precision arithmetic can be obtained
easily. No LLL reduction is provided which so has to be done separately.

Overview
--------

The Closest vector problem and the Shortest Vector Problems are fundamental
problems in Geometry of Numbers:
  * The Closest vector problem ask given a positive definite quadratic form Q a vector c in R^n to find the vectors v in Z^n that minimize Q[v - c].
  * The Shortest vector problem ask given a positive definite quadratic form Q to find the non-zero vectors v in Z^n minimizing Q.

The algorithm used in this work is the one by
  * Fincke, U. and Pohst, M., Improved methods for calculating vectors of short length in a lattice, including a complexity analysis, Mathematics of Computation, 44, 1985, 463--471

Example
-------

The Examples directory contain the example of the E8 lattice for which we compute the closest vector. On the first line is the dimension, then the matrix with the lower triangular part and then the vector in question. In this example, only the vector is expressed as fractions, but actually everything can be done with fractions.

References
----------

For applications of the closest vector problem the user is referred to following papers and books:
  * A. Schürmann, Computational geometry of positive definite quadratic forms. Polyhedral reduction theories, algorithms and applications. American Mathematical Society, xvi+162 pp.
  * M. Dutour Sikirić, A. Schürmann, F. Vallentin, Complexity and algorithms for computing Voronoi cells of lattices, Mathematics of Computation 78 (2009) 1713--1731, preprint at arxiv:0804.0036
  * M. Dutour Sikirić, A. Schürmann, F. Vallentin, Classification of eight dimensional perfect forms, Electronic Research Announcements of the American Mathematical Society 13 (2007) 21--32, preprint at arxiv:0609388

Dependencies
------------

Following dependencies are needed for compiling the code:

  * Eigen: http://eigen.tuxfamily.org/
  * Boost: http://www.boost.org/
  * GNU MultiPrecision Library (GMP): https://gmplib.org/
