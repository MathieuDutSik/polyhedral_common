Polyhedral: Computational Tools for Polytopes, Lattices, and Quadratic Forms
============================================================================

This is the set of functionality for dealing with polytopes,
quadratic forms and lattices.

The general approach is to use polytope and polyhedral
structures and groups are used to make everything faster
as a general rule.

The goal is to get extreme speed in order to solve record
problems.


General organization of the code
--------------------------------

The program are in several independent subdirectory. The software of each
directory can be compiled independently of the others:
  * *src_group*: for computting the groups of a polytope.
  * *src_poly*: for polyhedral computations.
  * *src_delaunay*: for computing Delaunay polytopes and space of Delaunay tesselations.
  * *src_dualdesc*: for computing dual description on serial computers.
  * *src_copos*: for copositivity / strict copositivity functionalities.
  * *src_short*: for short vector related computations.
  * *src_latt*: for fundamental lattice related computations (canonicalization of positive definite matrices and shortest vector mostly).
  * *src_ctype*: for computing C-types. This was done for computing all the C-types in dimension 6.
  * *src_perfect*: for perfect form related computations.
  * *src_indefinite*: for indefinite form reduction.
  * *src_lorentzian*: for using the Vinberg/Edgewalk algorithm of hyperbolic forms. Also computes perfect hyperbolic forms.
  * *src_sparse_solver*: Sparse solver for linear equations.
  * *src_isotropy*: Testing if quadratic forms are isotropic and finding a zero.

Works in Progress:
  * *src_rankin*: for computing rankin constants.
  * *src_single_delaunay*: About computing space for a single Delaunay.
  * *src_poincare_polyhedron*: Applying the Poincare Polyhedron Theorem to some tiling spaces.
  * *src_robust_covering*: Finding robust covering density.
  * *src_k_coverings*: Finding k-coverings of lattices.


Arithmetics
-----------

Many of the programs of this repository allow to take different arithmetic.
This is both from a functional aspect and a speed aspect. There is no systematic
way this is done but the following rules are applied. The types are always
template parameter.

Functional types:
  * **Qsqrt2**, **Qsqrt3** and **Qsqrt5**. Those types are used for the fields **Q(sqrt(2))** and similar. For each field we need a type. So, if you need say **Q(sqrt(6))**, you need to modify the code.
  * The function **RealAlgebraic=FileDesc** is for working with real algebraic numbers with the description in **FileDesc**. That file should contain the continous fraction approximant up to some chosen precision. If insufficient, a clean error of failure will be reported.

For both those types, the coefficient are written in entries like **1**, **1+x**, **x/4**, **-3-3x/4**, **3+(3*4)*x^3** or such. No space in the entry.

Speed types:
  * The **mpq_class** is from gmp and is the standard rational type being used.
  * The **mpz_class** is from gmp and is the standard integer type being used.
  * The **boost::multiprecision::cpp_rational** is the pure boost based rational type. Not as fast as gmp, but purely heaer based.
  * The **boost::multiprecision::cpp_int** same as above but for the integers.
  * The **boost::multiprecision::mpq_rational** is boost based encapsulation of gmp library.
  * The **boost::multiprecision::mpq_int** same as above but for the integers.
  * The **Rational<T>** For the rational implementation from an integer type.
  * The **SafeInt64** for integer computation in int64. If the computation go over the limit, then an exception is raised and typically the computation stops. But no wrong result is reported.
  * The **Rational<SafeInt64>** the same as above but for rational fields.

The choice is usually for **mpq_class** which has the least issues.


Usage
-----

The primary use of the programs is via the CLI. But this can be relatively impractical.

Therefore, we plan to have interfacing with common Computer Algebra Systems:
* GAP (Group Algebra Programming) which is an independent system.
* Oscar which is based on Julia.
* Sage which is based on Python.
