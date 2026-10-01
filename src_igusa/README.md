Igusa polyhedron of a T-space
=============================

For a T-space T of symmetric n x n matrices, let L be the lattice of the
integral valued matrices of T (X[v] in Z for v in Z^n) and I the positive
definite elements of L. The Igusa polyhedron is P = conv(I). For the full
space of symmetric matrices this is the central cone of Igusa. The C++ code
is a T-space version of `Tspace_Igusa.g` of the MyPolyhedral GAP package.

`IGUSA_EnumerateVertices` enumerates the vertices of P up to the arithmetic
group of the T-space, and for each vertex A gives its stabilizer, its edges
(neighbors and infinite edges) and its facets. It also gives the orbits of
facets of P, as inequalities tr(F X) >= rhs. The facets X[v] >= 1 are
those with F of rank 1, the other ones are specific to P. The T-space is
given by the standard `&TSPACE` block (Classic, ImagQuad for the Gaussian
and Eisenstein integers, InvGroup, Raw, File).

The method
----------

All the computations come down to minimizing a linear function over I. This
is done by cutting planes X[v] >= 1, first in linear programming, then in
integer programming over L. The cuts are valid on all of I and are kept in a
pool shared by the whole enumeration.

For a vertex A the local cone C_A = cone{ D in L : A + D positive definite }
is computed as follows:

* Candidate rays: the D with A + D positive definite and
  tr((A^{-1} D)^2) <= bound, the bound being doubled until they span the
  space.
* Random facets of the cone of the candidates are checked by integer
  programming. A facet f is valid iff there is no X in I with f(X) < f(A).
  The violating X found is the closest to A for tr(A^{-1} X), which is
  positive on the positive semidefinite matrices, so the integer programs
  stay bounded.
* Norm closure: all the candidates up to the largest norm of the extreme
  rays found are added.
* Only then the dual description is computed and all its orbits of facets
  are checked. With the previous steps, the list of rays is normally correct
  at that point, so a single dual description is done. The number done is
  reported as `nbDualDescription` for each vertex.

An extreme ray D, taken primitive in L, gives the edge [A, A + t D] with t
the largest integer such that A + t D is positive definite, or an infinite
edge if D is positive semidefinite.

Namelist
--------

```
&DATA
 arithmetic = "flint"    ! gmp, or flint (faster) with ENABLE_FLINT_SUPPORT
 IlpMethod = "default"   ! default, scip, exact_bb
 OutFormat = "GAP"       ! GAP, PYTHON
 OutFile = "result.g"    ! or stderr, stdout
 NormBound = 2           ! initial bound for the candidate rays
 NbRandomFacet = 10      ! random facets tested per round
 NbCleanRound = 2        ! clean rounds before the dual description
/

&TSPACE
 TypeTspace = "ImagQuad"
 SuperMatMethod = "NotNeeded"
 ListComm = "Use_realimag"
 PtGroupMethod = "Trivial"
 FileListSubspaces = "unset"
 RealImagDim = 3
 RealImagSum = 0         ! 0, 1 for the Gaussian integers
 RealImagProd = 1        ! 1, 1 for the Eisenstein integers
/
```

Output
------

`rec(ListVertex:=[...], ListFacetOrbit:=[...])`. Each vertex has its Gram
matrix, the order of its stabilizer, its number of extreme rays, the number
of dual descriptions done, its facets up to the stabilizer as
`rec(F, rhs)` and its infinite edges, and the neighbors with the matrix
mapping them to their representative. Each orbit of facets has F, rhs, the
rank of F and the list of [vertex, facet] pairs where it appears.

Integer programming
-------------------

The integer programs use SCIP when compiled with it
(`make ENABLE_SCIP_SUPPORT=1`, or the `ENABLE_SCIP_SUPPORT` CMake option)
and otherwise an exact branch and bound, see
`src_milp/integer_linear_programming.h`. The SCIP specific code is in
`src_milp/scip_support.h`, included under `#ifdef ENABLE_SCIP_SUPPORT`.
SCIP works in floating point: the problems are passed with integral
coefficients and the solutions are checked exactly. With `SANITY_CHECK`
each SCIP result is compared with the exact branch and bound.
