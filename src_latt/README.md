Lattice based computation
=========================

The functions in this subdirectory allow computing with lattices

General references of the lattice:
  * A. Schürmann, Computational geometry of positive definite quadratic forms. Polyhedral reduction theories, algorithms and applications. American Mathematical Society, xvi+162 pp.

LLL algorithm
-------------

The LLL algorithm allow to reduce a quadratic form into one with smaller
coefficients.

The relevant program is:
  * **LATT_lll** This is for computing the LLL of a positive definite quadratic form.


Closest Vector Problems
-----------------------

The Closest vector problem and the Shortest Vector Problems are fundamental
problems in Geometry of Numbers:
  * The Closest vector problem ask given a positive definite quadratic form Q a vector c in R^n to find the vectors v in Z^n that minimize Q[v - c].
  * The Shortest vector problem ask given a positive definite quadratic form Q to find the non-zero vectors v in Z^n minimizing Q.

The relevant programs are:
  * **LATT_near** It allows to do those operations using LLL.
  * **sv_exact** the same with different APIs.

The algorithm used in this work is the one by
  * Fincke, U. and Pohst, M., Improved methods for calculating vectors of short length in a lattice, including a complexity analysis, Mathematics of Computation, 44, 1985, 463--471

Stabilizer/Equivalence of positive definite quadratic form
----------------------------------------------------------

There is a set of function for testing stabilizer/equivalence for positive definite quadratic forms.

Programs:
  * **LATT_GenerateCharacteristicVectorSet** Computing the set of characteristic vectors of a positive definite form.
  * **LATT_Automorphism** The automorphism group of a positive definite form
  * **LATT_Isomorphism** The equivalence of two positive definite forms.

Reference:
  * W. Plesken, B. Souvignier, Computing isometries of lattices, J. Symbolic Computation 24 (1997) 327--334

Canonical form of positive definite forms
-----------------------------------------

For a positive definite quadratic form, we want to find the canonical form of a form as this is very
useful for enumeration purposes.

Programs:
  * **LATT_canonicalize** Canonicalize a positive definite form
  * **LATT_canonicalizeMultiple** Canonicalize several positive forms with the first one being positive definite.
  * **LATT_canonicalizeSymplectic** Finding the canonical form with the equivalence being a symplectic matrix.

Reference:
  * Mathieu Dutour Sikirić, Anna Haensch, John Voight, Wessel Van Woerden, A canonical form for positive definite matrices, Proceedings of the Fourteenth Algorithmic Number Theory Symposium (ANTS-XIV), edited by Steven Galbraith, Open Book Series 4, Mathematical Sciences Publishers, Berkeley, 2020.

T-space of quadratic forms
--------------------------

A T-space is a vector space of quadratic forms which contains one positive definite form.

Programs:
  * **TSPACE_Equivalence** Computing equivalence within the T-space.
  * **TSPACE_Stabilizer** Computing stabilizer within the T-space
  * **TSPACE_FileFormatConversion** Converting format of T-spaces

Reference:
  * Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, A generalization of Voronoi's reduction theory and its application, preprint at arxiv:0601084, Duke Mathematical Journal 142 (2008) 127--164
  * Herbert Gangl, Paul Gunnells, Jonathan Hanke, Achill Schürmann, Mathieu Dutour Sikirić, Dan Yasaki, On the cohomology of linear groups over imaginary quadratic fields, preprint at arxiv:1307.1165, Journal of Pure and Applied Algebra 220-7 (2016) 2564--2589
  * Achill Schuermann, Enumerating perfect forms, arXiv:0901.1587


Lattice Delaunay
----------------

The enumeration of Delaunay polytopes is an important geometric problem.

References for computing the Delaunay of a lattice:
  * Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin, Complexity and algorithms for computing Voronoi cells of lattices, preprint at arxiv:0804.0036, Mathematics of Computation 78 (2009) 1713--1731


Iso-Delaunay domains
--------------------

For a T-space, we want to enumerate the iso-Delaunay domains inside of it.

References:
  * Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, A generalization of Voronoi's reduction theory and its application, preprint at arxiv:0601084, Duke Mathematical Journal 142 (2008) 127--164
