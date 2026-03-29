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


Usage
-----

The primary use of the programs is via the CLI. But this can be relatively impractical.

Therefore, we plan to have interfacing with common Computer Algebra Systems:
* GAP (Group Algebra Programming) which is an independent system.
* Oscar which is based on Julia.
* Sage which is based on Python.
