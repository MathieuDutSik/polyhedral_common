Polyhedral: Computational Tools for Polytopes, Lattices, and Quadratic Forms
============================================================================

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
$ git clone https://github.com/MathieuDutSik/polyhedral_common --recursive
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


Compilation options related to debug
------------------------------------

There are several environment variables that can be used during runtime.
* `DEBUG` for making some print statements and making some checks. Those checks have to be fast.
* `KEY_VALUE` for printing some `KEY(....) VALUE=(....)` that can be used for postprocessing of the options and heuristic optimization.
* `CHECK` for making checks that can be expensive to do.
* `TIMINGS` for printing some runtime information.
* `SANITY_CHECK` for doing some checks and stopping if incoherence are detected.

The options `TIMINGS` and `DEBUG` enable all the timings and debugging statements.
For a more granular debugging, stuff like `DEBUG_LINEAR_PROGRAM` can be used. See
the top of the header files.

A printout to `std::cerr` should occur if an error has been identified and the program
will terminate with a call to `TerminalException`. Other print statement should be
encapsulated in `std::ostream & os` that should be passed by reference from the initial
call. So typically for serial output we pass the `std::cerr` while for the parallel runs
we pass a stream to the output of that process. That way we avoid mixing between
different sources.


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
