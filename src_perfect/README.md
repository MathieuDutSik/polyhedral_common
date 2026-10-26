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
quotient of the perfect form complex by a finite index subgroup. The code
is split in `hecke_common.h` (failures, budgets, the ambient group, the
Hecke double cosets), `hecke_subgroups.h` (the congruence subgroups),
`hecke_chains.h` (chains, walk through the perfect cones),
`hecke_chain_map.h` (step A and the `Vertex` method), `hecke_dual.h` (the
`Dual` method), `hecke_symbols.h` and `hecke_sharbly.h` (the `Sharbly`
method), `hecke_homology.h` (step B), the sparse exact solutions being in
`src_sparse_solver/SparseChainSolver.h`. It is driven by the `&HECKE`
block of `PERF_SerialPerfectComputation`:

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
* `HeckeSharblyDepth`: for `Sharbly`, the number of levels below the
  symbol level that are computed (default 1, the top two degrees). The
  levels further down are filled from the image of the boundary by the
  map of the level above (no reduction of the k-sharblies); depth 2 gives
  H^4 for GL_3(O_K). A cached chain map with fewer levels is extended and
  written back.
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
* `SubgroupLevelIdeal`: for a Hermitian T-space, the level as an ideal n
  of O_K given by generators, `"a,b"` or `"a1,b1;a2,b2"` for a + b omega
  (any ideal, also for a class number > 1); it replaces `SubgroupLevel`.
  `"none"` (the default) for the level `SubgroupLevel` O_K.
* `HeckeHomologyAllLevels`: also give the dimension of the homology at the
  levels where the chain map is not computed (without the Hecke matrix).
* `HeckeBoundaryMap`: the boundary map of the Satake compactification on
  the homology modulo p (needs `OnlyWellRounded = F` in DATA, so that the
  cells at infinity, whose vectors do not span, are in the complex, and the
  Sharbly or Vertex method). The class of a cycle z goes to the class of
  the part of its boundary on the cells at infinity, in the homology of the
  subcomplex at infinity; for GL_2 over Q this is the boundary map of the
  modular symbols to the cusps. Its kernel (`dim_ker_boundary`, with the
  Hecke matrix `HeckeMatrixKerBoundary` on it) contains the interior
  cohomology. For GL_3 over O_K it detects the Eisenstein classes coming
  from the parabolic subgroup fixing a line (eigenvalue q^2 + a for a
  Bianchi form), but not those coming from the one fixing a plane
  (1 + q a), whose stratum in the Satake boundary is a point. Those are
  seen by the dual Satake compactification (the forms replaced by their
  inverses): with theta(g) = (g^T)^{-1}, the composition of theta and of
  the conjugation by diag(N, ..., N, 1) (which maps Gamma_0(N) to theta(
  Gamma_0(N))) is an involution of the cohomology exchanging T(p, k) and
  T(p, n - k), and it maps the kernel of the boundary map to the kernel of
  the dual one. So an eigensystem (a(p, 1), a(p, 2)) is in the dual kernel
  iff (a(p, 2), a(p, 1)) is in the kernel, and the interior cohomology is
  in the eigensystems whose swap is also in the kernel: for GL_3 over
  Q(sqrt(-7)), Gamma_0(11), the kernel is the class (-3, 2) (1 + q a,
  q^2 + a for a = -2), whose swap is not there, while at Gamma_0(13) the
  pair (lambda, conj(lambda)) stays. Note that theta(Gamma_0(N)), being
  conjugate to Gamma_0(N) in GL_n(Q), gives the same boundary map.
* `OnlyWellRoundedHomology`: restrict to the well rounded cells, which
  gives the cohomology of the subgroup by duality; otherwise the homology of
  the full complex with its cells at infinity.

The complex must be computed with `OnlyWellRounded = F`, `ComputeBoundary = T`
and `ComputePartialContractingHomotopy = T`: the construction starts from the
vertices (mapped to vertices) and goes up by the partial contracting homotopy.
The GAP access point is `PERFCOMP_hecke_operators` in `CI_tests/access_points.g`.
