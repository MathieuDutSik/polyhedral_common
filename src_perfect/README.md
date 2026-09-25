Perfect form enumeration
========================

This is the same business as for IsoDelaunay domains,
IsoEdge domains and others.

However, due to the particular importance of perfect forms,
we ought to think carefully about what we want to do.

What we want from this is the following:
* Work in a T-space. That is super important. So that means
the space of the perfect domains and the Ryshkov space of the
perfect forms. This is like that in Opgenorth paper.
* The Ryshkov stuff allows to do the enumeration quite
efficiently. On the other hand, the perfect domain are better
for computing cohomology.
* The finding of initial form and the flipping are quite clear
from existing code. The equivalence and stabilizer can also
be done from the Tspace equivalence code. That part is clear.
* How to build the dual cone? Possible sources:
  --- The TechVor.pdf unpublished manuscript.
  --- The GL4(Z[i]) paper.
  --- The GAP code, there are several.
* What we want from the dual cone:
  --- That it embeds into S^n_{>0} so that the existing tools
  for stab/equi can be used.
  --- That the linear form of evaluation correspond to a vector
  in the dual.
* What to do:
  --- It is a little bit unclear whether we can have a general theory
  --- It would make sense to have the dual space put as argument. If
  there is no universal theory, then so be it.
  ---




Hecke operators
===============

`hecke_operators.h` computes the Hecke operators on the homology of the
quotient of the perfect form complex by a finite index subgroup. It is
driven by the `&HECKE` block of `PERF_SerialPerfectComputation`:

* `FileHeckeMatrix`: the rational matrix x of the double coset Gamma x Gamma
  (it must preserve the T-space).
* `FileHeckeChainMap`: output of step A, the chains representing the images
  of the cells under the cosets of G x G.
* `FileHeckeHomology`: output of step B, the Hecke matrices on the homology
  of the quotient by the subgroup, level by level, in a basis of harmonic
  representatives.
* `SubgroupType` / `SubgroupLevel`: `Full`, `Principal`, `Gamma0` (last row
  congruent to (0,...,0,*)) or `Gamma1` (last row congruent to (0,...,0,1)).
* `OnlyWellRoundedHomology`: restrict to the well rounded cells, which
  gives the cohomology of the subgroup by duality; otherwise the homology of
  the full complex with its cells at infinity.

The complex must be computed with `OnlyWellRounded = F`, `ComputeBoundary = T`
and `ComputeContractingHomotopy = T`: the construction starts from the
vertices (mapped to vertices) and goes up by contracting homotopies.
The GAP access point is `PERFCOMP_hecke_operators` in `CI_tests/access_points.g`.
