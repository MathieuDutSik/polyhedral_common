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
* `HeckeMethod`: how step A is built. `Dual` (induction from the perfect
  cones in the well rounded complex), `Vertex` (induction from the vertices
  in the full complex, needed for the homology of the full complex) or
  `Sharbly` (Ash-Rudolph reduction of the modular symbols at the top degree
  and filling of the level below: the top two degrees only). `Sharbly` also
  works for the Hermitian T-spaces (`ImagQuad`) over a norm-Euclidean ring
  O_K (d = -1, -2, -3, -7, -11): the symbols are then m lines over O_K and
  the top two degrees are 5 and 6 for GL_3(O_K).
* `FileHeckeChainMap`: output of step A, the chains representing the images
  of the cells under the cosets of G x G.
* `FileHeckeChainMapCache`: step A does not depend on the subgroup. If this
  file holds the chain map of the same x on the same complex it is read
  and step A is skipped (whatever `HeckeMethod` says); otherwise step A is
  computed and written there. A run reading the cache of a `Sharbly` or
  `Dual` chain map can use `ComputePartialContractingHomotopy = F`, which
  saves the memory of the coverings.
* `FileHeckeHomology`: output of step B, the Hecke matrices on the homology
  of the quotient by the subgroup, level by level, in a basis of harmonic
  representatives.
* `SubgroupType` / `SubgroupLevel`: `Full`, `Principal`, `Gamma0` (last row
  congruent to (0,...,0,*)) or `Gamma1` (last row congruent to (0,...,0,1)),
  and their intersections with SL_n: `SL`, `PrincipalSL`, `Gamma0SL`,
  `Gamma1SL`. The elements of determinant -1 can reverse the orientation of
  the cells and kill them in the quotient complex: the cohomology of the
  literature (Ash, Gunnells, McConnell for SL_4) is the one of the SL
  variants. For a Hermitian T-space the subgroups are those of GL_m(O_K):
  the congruences are on the entries over O_K modulo N O_K (N an integer)
  and the SL condition is that the determinant over O_K is 1.
* `OnlyWellRoundedHomology`: restrict to the well rounded cells, which
  gives the cohomology of the subgroup by duality; otherwise the homology of
  the full complex with its cells at infinity.

The complex must be computed with `OnlyWellRounded = F`, `ComputeBoundary = T`
and `ComputePartialContractingHomotopy = T`: the construction starts from the
vertices (mapped to vertices) and goes up by the partial contracting homotopy.
The GAP access point is `PERFCOMP_hecke_operators` in `CI_tests/access_points.g`.
