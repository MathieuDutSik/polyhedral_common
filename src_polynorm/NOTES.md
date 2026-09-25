# src_polynorm — design notes and status

Packing and covering of `Z^n` by the translates of a polytope `P`, that is
the first minimum of the difference body and the covering radius of the
polyhedral gauge `||.||_P`. `README.md` describes the programs; this file
records the design choices, what was validated and what is missing.

## Setting

`P` is given by its vertices, recentered at their isobarycenter so that `0`
is interior and the symmetries of the problem are linear: the finite group
`G = {A in GL_n(Z) : PA = P}`, computed as the integral stabilizer of the
vertex set (`LinPolytopeIntegral_Automorphism`). Both parameters are
translation invariant, so the recentering loses nothing, and it is
required: the group of the uncentered polytope could be smaller (its affine
symmetries have non-integral translation parts).

Facets are stored as linear forms `ell_j` with `P = {x : ell_j(x) <= 1}`,
so the gauge is `max_j ell_j(x)`. A symmetry `A` with facet permutation
`sigma` satisfies `ell_{sigma(j)}(xA) = ell_j(x)`, which is what makes the
invariants of the covering search work.

## Packing

`alpha(P) = min_{z != 0} ||z||_{P-P}`: the translates `alpha P + v` and
`alpha P + w` have disjoint interiors exactly when `w - v` is outside the
interior of `alpha (P - P)`. The facets of `P - P` come from the dual
description of the differences of vertices (at most `|V|^2` points, which
is fine in the dimensions where the covering is feasible at all), the
gauge of `e_1` bounds `alpha` from above and the lattice points of the
dilated difference body are enumerated with `GetListIntegralPoint`. The
group is not needed here; the contact vectors are returned and split into
orbits only by the caller if wanted.

## Covering: why a branch and bound and not a Delaunay style enumeration

The Delaunay code enumerates the Delaunay polytopes up to the symmetry
group by walking across their facets, and the covering radius is the
largest circumradius. The analogue here would be the empty translates
`p - mu P` pinned by lattice points, whose centers are the local maxima
of `f(p) = min_z ||p - z||_P`. That walk is not available for polyhedral
gauges:

* the Voronoi cells of the gauge are star shaped but not convex, and the
  bisectors of two points can contain full-dimensional pieces (the
  well-known degeneracy of polyhedral norm bisectors, see Icking, Klein,
  Lê, Ma, Santos and the Deza–Dutour Sikirić paper on Voronoi polytopes for
  polyhedral norms in the literature folder), so there is no polyhedral
  complex of Voronoi vertices and edges to walk on;
* the local maxima of `f` are not connected by a natural adjacency: the
  edges of the arrangement leaving a local maximum go through vertices that
  are not local maxima, and the whole arrangement would have to be
  enumerated.

The function `f` is however a minimum of the convex piecewise linear
functions `g_z(p) = ||p - z||_P`, so `{(p, mu) : mu <= f(p)}` is an
intersection over `z` of unions over the facets `i` of the halfspaces
`mu <= ell_i(p - z)`. This is a disjunctive linear program, and a branch and
bound over the disjunctions is exact and finite:

* a node is a set `C` of pairs `(z, i)`, the constraints being that `p` is
  in the cone of the facet `i` at `z` (`ell_i(p - z) >= ell_j(p - z)` for
  all `j`, so that `g_z(p) = ell_i(p - z)`) and `mu <= ell_i(p - z)`. The
  cone constraints make the children of a branching a partition of space
  around `z'`, which is what keeps the relaxations tight;
* the linear program `max mu` of a node bounds it, the value `f(p*)` at its
  optimum is a lower bound on `mu(P)` (an incumbent), and a node is solved
  when the two agree. Otherwise the lattice point minimizing `g_z(p*)` is
  branched on; it is found among the lattice points of the translate
  `p* - mu* P`, enumerated at each node. This per-node enumeration replaced
  a precomputed set `Z^n cap mu0 (P - P)`: the latter grows like `mu0^n`
  and made every node scan tens of thousands of points when `mu0` was loose
  (the 4-dimensional simplex went from 54 s to 0.04 s);
* `0` is required to be a nearest lattice point of `p` (`p in mu P`), which
  is a valid normalization by translation invariance, so the roots are the
  facets of `P`;
* the linear programs are the exact fraction-free simplex
  (`SIMPLEX_LinearProgramming`) in `n+1` variables with `1 + m (|C| + 1)`
  rows; about 0.2 ms each in dimension 4.

Use of the group:

* the roots are reduced to one per `G`-orbit of facets
  (`OrbitSplittingSet` on the facet incidences, needs only generators);
* when `|G| <= 50000` the elements are listed as matrices with their facet
  permutations, and every node is checked against the processed ones under
  the affine group `G x Z^n`. A cheap invariant is computed for every node:
  per pair `a`, the orbit of its facet and, sorted over the other pairs
  `b`, the relative positions `ell_{i_a}(c_{i_b})` of the two facets
  (`c_i` the facet isobarycenter) and the values of `ell_{i_a}`, `ell_{i_b}`
  and the gauge on `+-(z_b - z_a)`. Nodes are only compared within an
  invariant class, by an explicit equivalence test: the pair 0 of the one
  must go to a pair of the other, which leaves the `|G|/m` elements
  mapping the one facet to the other, and the translation is then fixed.
  Two earlier versions were bottlenecks: a canonical form (minimum over all
  `|G|` elements) for every child cost more than it saved on the
  4-dimensional cross polytope (3.2 s against 1.2 s with no pruning at
  all, although it halves the node count), and the invariant without the
  facet positions put hundreds of nodes in one class on cross polytopes
  (9 s). With the final version the equivalence tests are exactly the
  true duplicates on every test case.

The upper bound `mu0` only bounds the root linear programs. It is obtained
from lattice polytopes `Q` of known covering radius fitted into `tP` by a
linear program in the translation (`mu(P) <= t mu(Q)` by monotonicity and
translation invariance): the simplices `conv(0, +-e_1, ..., +-e_n)`
(`mu = n`), the cross polytope (`n/2`) and the cube `[-1,1]^n` (`1/2`).
This is exact on the simplices, cross polytopes and cubes themselves.

## Validation

`POLYNORM_TestCovering` and `CI_tests/25B_PolyNorm/run_tests.sh`: known
values in dimensions 2 to 5 (cubes, cross polytopes, simplices, the
hexagon, the 24-cell) and the brute force enumeration of the paper (all
systems of `n+1` pairs, one at `z = 0`, solved exactly) on every
2-dimensional case and on the 3-dimensional cube.

Timings (Apple M-series, exact arithmetic, sanity checks on): every case
of dimension at most 4 in under 1.5 s; the 5-dimensional cube in 0.5 s;
the 5-dimensional cross polytope (`m = 32`, `|G| = 3840`, 798 nodes and
22309 linear programs) in 36 s, down from 559 s with the first canonical
form. A random 3-dimensional polytope with 12 vertices and trivial group
takes 0.5 s (128 nodes).

## Known issues / future work

* **Repeated vertices hang the upstream automorphism code.** With two
  equal rows `LinPolytopeIntegral_Automorphism` does not return. The input
  is deduplicated here, but the upstream function should reject or handle
  repeated rows.
* **Large facet counts.** The branching factor is the number of facets
  `m`, and the linear programs dominate the 5-dimensional cross polytope
  (`m = 32`, 36 s). The children that are infeasible or dominated could be
  filtered before their linear program, and a dual simplex warm start from
  the parent would cut the per-child cost.
* **Large groups.** Above 50000 elements only the root reduction is used.
  The stabilizer of a node could be computed from the permutation group
  instead of by filtering the element list.
* **Upper bound.** Only the standard unimodular simplices and the
  coordinate cross polytope and cube are fitted. Fitting the images of
  these under a reduced basis of `Z^n` (LLL with respect to a quadratic
  form adapted to `P`) would help for elongated polytopes.
* **Other arithmetics.** Only `mpq_class` is wired in the drivers; the
  headers are templated on `T` and `Tint` and the other rational types of
  the package should just work.
* **Lattices other than `Z^n`** are handled by a change of basis by the
  caller; a `[FileLattice]` argument would be a convenience.
* **Non-lattice (periodic) point sets** are not handled.
