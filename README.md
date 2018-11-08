Polytopes, lattices and quadratic forms programs
================================================

This is the set of functionality for dealing with polytopes.

Since this repository uses submodules, the cloning command is

```sh
$ git clone git@github.com:MathieuDutSik/polyhedral_common.git --recursive
```

There is a number of programs for polytopes, lattices and quadratic
which are in a number of independent directories:
  * src_copos: for copositivity / strict copositivity functionalities.
  * src_poly: for polyhedral computations.
  * src_short: for short vector related computations.
  * src_latt: for lattice related computations
  * src_perfect: for perfect form related computations.


Copositivity
------------

in src_copos


Shortest vectors
----------------

the src_latt directory contains the sv_exact program for 



LLL computation
---------------

The src_latt directory contains the LATT_lll programs for
computing the LLL reduced form of a positive definite form.


Canonical form
--------------

The src_latt directory contains the LATT_canonicalize for computing
the canonical form of a positive definite form.