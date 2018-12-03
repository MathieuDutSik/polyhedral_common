Polytopes, lattices and quadratic forms programs
================================================

This is the set of functionality for dealing with polytopes.

Since this repository uses submodules, the cloning command is

```sh
$ git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive
```

In order to update the submodule the command is
```sh
$ git submodule update --remote
```



There is a number of programs for polytopes, lattices and quadratic
which are in a number of independent directories:
  * *src_copos*: for copositivity / strict copositivity functionalities.
  * *src_poly*: for polyhedral computations.
  * *src_short*: for short vector related computations.
  * *src_latt*: for lattice related computations
  * *src_perfect*: for perfect form related computations.


Copositivity
------------

In *src_copos* there is a number of functionality for working with copositive forms.
The subprograms are:
 * **CP_ComputeCopositiveMin**: It tests is a matrix is copositive.
 * **CP_TestCompletePositivity**: We test if a matrix is completely positive.
 * **CP_CopositiveMin**: Get the list of positive vectors of small norms


Shortest vectors
----------------

The *src_latt* directory contains the **sv_exact** program for computing the set of
closest points of a positive definite form. The code is general and has
two template types:
 * an integer type **Tint** for the integer coefficients.
 * a coefficient type **T** for the coefficient of the Gram matrix.


LLL computation
---------------

The *src_latt* directory contains the **LATT_lll** programs for
computing the LLL reduced form of a positive definite form.
The code uses **Tint**/**T** template types.


Canonical form
--------------

The *src_latt* directory contains the **LATT_canonicalize** for computing
the canonical form of a positive definite form.
The code uses **Tint**/**T** template types.


Short vector configurations
---------------------------

The directory *src_short* contains a set of functions for dealing with
configuration of vectors that can occur as set of minimum vectors
of positive definite forms.


Sparse solver
-------------

The directory *src_sparse_solver* contains the code for a solver of the
equations **Ax = b** by finding sparse solution **x** for sparse matrices
**A** and **b**. The algorithmic method used is Generalized Approximate
Message Passing.

Compilation
-----------

The compilation of the software is relatively complex. However in
script_docker/Dockerfile a dockerfile is given that should make everything
clear.
