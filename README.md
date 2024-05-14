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
  * *src_group*: for computting the groups of a polytope.
  * *src_poly*: for polyhedral computations.
  * *src_dualdesc*: for computing dual description on serial computers.
  * *src_copos*: for copositivity / strict copositivity functionalities.
  * *src_short*: for short vector related computations.
  * *src_latt*: for lattice related computations (canonicalization of positive definite matrices and shortest vector mostly)
  * *src_ctype*: for computing C-types. This was done for computing all the C-types in dimension 6.
  * *src_perfect*: for perfect form related computations.
  * *src_indefinite*: for indefinite form reduction.
  * *src_lorentzian*: for using the Vinberg/Edgewalk algorithm of hyperbolic forms.

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

Compilation options related to debug
------------------------------------

There are several environment variables that can be used during runtime.
* `DEBUG` for making some print statements and making some checks. Those checks have to be fast.
* `KEY_VALUE` for printing some `KEY(....) VALUE=(....)` that can be used for optimization of the 
* `CHECK` for making checks that can be expensive to do.
* `TIMINGS` for printing some runtime information.
* `SANITY_CHECK` for doing some checks and stopping if incoherence are detected.

If we want more modular checking, then something like `DEBUG_LINEAR_PROGRAM`
can be used.

A printout to `std::cerr` should occur if an error has been identified and the program
will terminate with a call to `TerminalException`. Other print statement should be
encapsulated in `std::ostream & os` that should be passed by reference from the initial
case. So typically for serial output we pass the `std::cerr` while for the parallel runs
we pass a stream to the output of that process. That way we avoid mixing between
different sources.

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

Usage
-----

The primary use of the programs is via the CLI. But this can be relatively impractical.

Therefore, we plan to have interfacing with common Computer Algebra Systems:
* GAP (Group Algebra Programming) which is an independent system.
* Oscar which is based on Julia.
* Sage which is based on Python.

Dependencies
------------

Following dependencies are needed for compiling the code:

  * Eigen: http://eigen.tuxfamily.org/
  * Boost: http://www.boost.org/
  * GNU MultiPrecision Library (GMP): https://gmplib.org/
  * nauty : https://pallini.di.uniroma1.it/
