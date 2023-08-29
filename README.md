Polytopes, lattices and quadratic forms programs
================================================

This is the set of functionality for dealing with polytopes,
quadratic forms and lattices.

The general approach is to use polytope and polyhedral
structures and groups are used to make everything faster
as a general rule.

The goal is to get extreme speed in order to solve record
problems.


Access to the source code
-------------------------

Since this repository uses submodules, the cloning command is

```sh
$ git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive
```

In order to update the submodule the command is
```sh
$ git submodule update --remote
```


Compilation
-----------

The compilation of the software is relatively complex. However in
**script_docker/Dockerfile** a dockerfile is given that should make everything
clear. It should also allow any user to install the system on their
computer and use it fairly easily.


General organization of the code
--------------------------------

The program are in several independent subdirectory. The software of each
directory can be compiled independently of the others:
  * *src_copos*: for copositivity / strict copositivity functionalities.
  * *src_poly*: for polyhedral computations.
  * *src_short*: for short vector related computations.
  * *src_latt*: for lattice related computations (canonicalization of positive definite matrices and shortest vector mostly)
  * *src_dualdesc*: for computing dual description on serial computers.
  * *src_ctype_mpi*: for computing C-types. This was done for computing all the C-types in dimension 6.
  * *src_perfect*: for perfect form related computations.
  * *src_vinberg*: for using the Vinberg algorithm of hyperbolic forms.


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

Compilation options
-------------------

There are several environment variables that can be used during runtime.
* `DEBUG` for making some print statements and making some checks. Those checks have to be fast.
* `KEY_VALUE` for printing some `KEY(....) VALUE=(....)` that can be used for optimization of the 
* `CHECK` for making checks that can be expensive to do.
* `TIMINGS` for printing some runtime information.
* `SANITY_CHECK` for doing some checks and stopping if incoherence are detected.

If we want more modular checking, then something like `DEBUG_LINEAR_PROGRAM`
can be used.

Tests
-----

The directory `Exmpl_Bench` contains a bunch of test and development.
They are typically work in progress and not necessarily in a state of being
finished.

The directory `CI_tests` contains some tests that are run in CI on GitHub.
Their runtime should be short, from 5 minutes to 1 hour. Together they
should cover as much as possible of the functionality of the code. If the
test is working like normally, then it should be scheduled in the cron to
run once per month. We do not want to overflow the credit that we have.

List of CI tests:
* `Canonicalization`: This is for running the canonicalization of Gram
matrices.
* `ConeIntersection`: This is checking different methods for computing the
intersection of polyhedral cones.
* `Copositivity`: This is for checking the copositivity code.
* `Reflective`: This is for using the edgewalk algorithm of Allcock for
building the polyheral cone.
* `SimpleDualDesc`: This is for the code for computing the dual description
of polyhedral cones.
* `WythoffH4`: This is for computing the dual description of the facets of
the orbit of `x W(H4)` for x a random vector (there are 4 orbits).
* `Redundant`: We test 3 methods for reducing by redundancy.

More needs to be added:
* Enumerating integral points in polytopes (the 8-dim perfect Delaunay polytopes
are good examples).
* Computing the orbits of K-dim faces of polytopes (the G6 and G7 are good examples).
* Running of the various orbit splitting algorithm.
* Computing automorphism group of polytope with LinPolytope.
* Computing automorphism group of skeletton of a polytope.
* Checking for pointedness of cones.
* The sampling of facets (using the Various CUT polytopes)
* Short vector enumeration.
