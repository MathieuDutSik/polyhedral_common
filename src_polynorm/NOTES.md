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

## Arithmetic: gmp against flint

The drivers take `gmp` (`mpq_class`) or, with `make ENABLE_FLINT_SUPPORT=1`,
`flint` (`fmpq_class`). Same results everywhere; the non-debug timings of
`POLYNORM_TestCovering` (Apple M-series):

| case                         | gmp     | flint  |
|------------------------------|---------|--------|
| random 3D, 12 vertices       | 0.51 s  | 0.10 s |
| cross polytope, dimension 4  | 0.94 s  | 0.27 s |
| 24-cell                      | 0.84 s  | 0.29 s |
| cube, dimension 5            | 0.50 s  | 0.14 s |
| cross polytope, dimension 5  | 36.0 s  | 9.8 s  |
| brute force, 3D cube, 4M systems | 7.8 s | 2.4 s |

So flint is the arithmetic to use; gmp stays the default of the CI test
since the runners do not install flint.

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
* **Other arithmetics.** `mpq_class` and `fmpq_class` are wired in the
  drivers; the headers are templated on `T` and `Tint` and the other
  rational types of the package should just work.
* **Lattices other than `Z^n`** are handled by a change of basis by the
  caller; a `[FileLattice]` argument would be a convenience.
* **Non-lattice (periodic) point sets** are not handled.

# Optimal lattices (PolyNorm_LatticeOpt.h, PolyNorm_Packing2D.h)

The programs above compute the packing and covering radii of `P` for a
fixed lattice. The harder problem is to choose the lattice: the densest
lattice packing `delta_L(P) = max vol(P) / det L` over the lattices `L`
with `P + L` a packing, and the thinnest lattice covering
`theta_L(P) = min vol(P) / det L` over the lattices with `P + L = R^n`.
The literature is in `Papers_polynorm_works/` (see its `INDEX.md`).

## What is known (the state of the art the code is tested against)

* Packing, dimension 2: exact for every convex polygon (Minkowski; Mount
  and Silverman give an O(m) algorithm). Quadrilaterals: closed formula of
  Sriamorn (arXiv:1411.4713).
* Packing, dimension 3: Betke and Henk (math/9909172) turned Minkowski's
  1904 characterization of critical lattices into an algorithm: some
  critical lattice of `D = P - P` has a "test set" of 6 or 7 lattice
  points `(1,0,0)_B, ..., (1,1,1)_B` on the boundary of `D`, so one
  enumerates the choices of 6 or 7 facets of `D` and minimizes `det` on
  the affine space of bases putting the points on those facets (a cubic
  polynomial in 2 or 3 variables, whose local minima they find with
  resultants). They computed all Platonic and Archimedean solids:
  tetrahedron 18/49, octahedron 18/19, cuboctahedron 45/49, truncated
  tetrahedron 207/304, ... Some values are irrational (icosahedron,
  snub cube).
* Packing, dimension 4 and more: no algorithm is known. Minkowski's
  method needs the classification of lattice cross-polytopes (known in
  dimensions 4 and 5) and a large case analysis that has never been done;
  Betke and Henk consider it impractical without new ideas. The general
  structural result is Swinnerton-Dyer's: a locally optimal (extreme)
  lattice has at least `n(n+1)/2` pairs of contacts.
* Covering, dimension 2: Fejes Toth (centrally symmetric: largest
  inscribed hexagon), Fary (triangle 3/2, the worst case), Sriamorn
  (quadrilaterals, closed formula).
* Covering, dimension 3: only the ball (Bambah), the space fillers, and
  now the octahedron, 9/8 with the body-centered cubic lattice (Lian and
  Xue, arXiv:2609.23210, after Dougherty and Faber proved strict local
  optimality). For the tetrahedron the best known lattice covering is
  125/63, a strict local optimum conjectured optimal; the best lower bound
  is 25/18 (Fu, Xue, Zong, arXiv:2202.06025).

So there is no general exact algorithm beyond dimension 3 for packing and
beyond dimension 2 for covering, and the practical approach is local
optimization with exact verification. This is what is implemented.

## The formulation

The lattice is `Z^n W`, with the basis vectors the rows of `W`. A lattice
vector is `zW` and `ell_j(zW) = z W ell_j^T` is linear in `W`; the
packing and covering code for `Z^n` is reused with the linear forms
`Lmat W^T`.

* Packing: minimize `|det W|` subject to `||zW||_{P-P} >= 1` for all
  `z != 0`. Each constraint is a disjunction over the facets `d_j` of
  `P - P` of linear constraints `z W d_j^T >= 1`.
* Covering: a hole (an empty translate `p - mu P` pinned by `n+1` pairs
  `(z_a, i_a)`, `ell_{i_a}(p - z_a W) = mu`) has, for fixed pairs, the
  radius `mu_S(W) = <G_S, W>`, LINEAR in `W`: the matrix of the system in
  `(p, mu)`, with rows `(ell_{i_a}, -1)`, does not depend on `W`. So the
  covering radius is a piecewise linear function of `W`, and the problem
  is to maximize `|det W|` subject to `mu(W) <= 1`.

Both are "maximize or minimize det over a region bounded by linear pieces",
the analogue for polyhedral norms of Voronoi's theory of perfect forms
(packing) and of the L-type domains (covering). The difference is that the
parameter is the basis `W` and not the Gram matrix: there is no convexity,
and `det` is a cubic in dimension 3.

## The search

Sequential linear programming with a trust region:

* at the current `W` the near-active constraints are collected: the
  contacts within a ratio of the packing radius, each with one of its
  active facets (chosen at random among ties), or the holes within a
  ratio of the covering radius with all the selections of one facet per
  lattice point;
* `det(W + X)` is replaced by `det(W) (1 + tr(W^{-1} X))` and the linear
  program over `X` (exact, `SIMPLEX_LinearProgramming`) gives the step;
* the trial point is rounded to a dyadic grid finer than the trust region
  and evaluated exactly with the packing or covering code; it is accepted
  when the density improves, and the trust region is adapted. Since every
  trial is evaluated exactly, the linear model only guides the search;
* the basis is LLL reduced from time to time (the problem is invariant
  under `GL_n(Z)`), which keeps the enumerations cheap;
* when the steps stall, a pattern search along random directions with
  exact evaluations is tried, and the linear programming resumes if it
  improves. This was needed for covering: at the optimum the deep holes
  are often degenerate (the hexagon tiling: each hole is at a vertex of
  three translates) and the linear model only sees the facet selections
  that are exactly tight at the current point. Two more measures make it
  see the degeneracy appear: the pairs that are tight up to a tolerance
  proportional to the trust region are included, and each system enters
  with the bound `max(mu, mu_S(W))` so that `X = 0` stays feasible.

## The polishing and the certificate

At the end the near-active constraints are imposed as equalities, and the
critical point of `det` on that affine space is found by Newton's method
in floating point and recognized as a rational matrix by continued
fractions. Its gradient is then checked to vanish exactly. Then, all in
exact arithmetic:

* the radius of the polished lattice is recomputed;
* first order: `s W^{-T} = sum_k lambda_k h_k` with `lambda >= 0`, by a
  linear program that also maximizes the support of `lambda`;
* second order: the Hessian of `det` in direction `X` is
  `det(W) ((tr Y)^2 - tr(Y^2))`, `Y = W^{-1} X`, and `tr Y = 0` on the
  tangent space `T` of the constraints with `lambda > 0`. If `-s tr(Y^2)`
  is positive definite on `T`, done. Otherwise the critical cone
  `C = T cap {<h_k, X> >= 0 : lambda_k = 0}` is used: `C = L + C'` with `L`
  its lineality space and `C'` pointed. The form must be positive definite
  on `L`, and its Schur complement (the minimum over `L`) must be strictly
  copositive on the extreme rays of `C'`. Copositivity is decided exactly
  by enumerating the faces of the simplex and their critical points.

For packing (`s = 1`), the positive definiteness of `-tr(Y^2)` needs
`dim T <= n(n-1)/2`, which is Swinnerton-Dyer's bound again. Near `W` the
feasible region is exactly the union, over the choices of one active facet
per contact, of polyhedra, so the certificate is checked on every choice
("branch"). It proves a strict local optimum. On the tetrahedron and the
cuboctahedron (7 contacts, 4 positive multipliers) the lineality space is
Betke and Henk's affine space, and the copositivity part is needed.

For covering (`s = -1`), choosing one active facet per lattice point of a
hole gives a lower bound `mu_S(W) <= mu(W)` near `W`, valid when the hole
has exactly `n+1` lattice points (each with any number of active facets).
The region `{mu <= 1}` then lies inside the polyhedron of these
constraints and the certificate proves a strict local optimum. A hole with
more than `n+1` lattice points makes the local region a union and is
reported as not certified. The holes are collected from the solved nodes
of the covering branch and bound (`collect_ratio`); two holes of the same
depth in one node give one record, so a hole can be missed. That cannot
make the certificate wrong: fewer constraints give a larger polyhedron,
on which a strict local maximum is harder, not easier, to certify.

When the critical point is not rational (it is algebraic, as for the
icosahedron in Betke and Henk; all the 4-dimensional examples below), the
floating point critical value is reported as `refined_floating_point`,
and the exact density of the search's last lattice, a rigorous bound, as
`density`.

The 2-dimensional packing is solved exactly and globally
(`PolyNorm_Packing2D.h`): Minkowski's theorem gives a critical lattice
with `b1, b2, b1 + b2` on the boundary of `P - P`, and conversely such a
lattice is admissible. For each triple of edges the feasible set is a
segment along which `det` is quadratic, so the minimum is rational and
found exactly.

## Validation

Release build (`CHOICE_COMPILATION` empty), Apple M-series, `gmp`. The
exact values are those of the polished lattices, the certificates are
exact. "Starts" counts how many random starts reach the best value.

| polytope | problem | found | known | certified | starts | time |
|---|---|---|---|---|---|---|
| triangle | packing | 2/3 | 2/3 (Fary) | yes, also exact 2D | 4/4 | 0.3 s |
| quadrilateral `asym2d` | packing | 3/4 | 3/4 (Sriamorn) | yes, also exact 2D | 4/4 | 0.4 s |
| lattice pentagon | packing | 14/15 | (exact 2D: 14/15) | yes | 4/4 | 0.3 s |
| triangle | covering | 3/2 | 3/2 (Fary) | yes | 4/4 | 5 s |
| quadrilateral `asym2d` | covering | 4/3 | 4/3 (Sriamorn) | yes | 4/4 | 8 s |
| hexagon | covering | 1 | 1 (tiling) | yes | 4/4 | 11 s |
| lattice pentagon | covering | 7/6 | | yes | 4/4 | 9 s |
| octahedron | covering | 9/8 (BCC) | 9/8 (Lian, Xue) | yes | 5/6 | 767 s |
| cube | covering | 1 | 1 (tiling) | tiling | 6/6 | 409 s |
| tetrahedron | covering | 125/63 | 125/63 (best known, conjectured optimal) | yes | 1/6 | 2317 s |
| tetrahedron | packing | 18/49 | 18/49 (Hoylman) | yes | 8/8 | 2 s |
| octahedron | packing | 18/19 | 18/19 (Minkowski) | yes | 8/8 | 1 s |
| cuboctahedron | packing | 45/49 | 45/49 (Hoylman) | yes | 8/8 | 2 s |
| truncated tetrahedron | packing | 207/304 | 207/304 (Betke, Henk) | yes | 6/8 | 2 s |
| cube, truncated octahedron | packing | 1 | 1 (tiling) | tiling | 8/8 | |
| 24-cell | packing | 1 | 1 (tiling) | tiling | 14/16 | 21 s |
| 4-cross-polytope | packing | 0.82485791730145 | unknown | numerical | 8/32 | 31 s |
| 4-simplex | packing | 0.173812554167293 | unknown | numerical | 5/32 | 45 s |

`POLYNORM_LatticeCertify` certifies directly the lattices of Betke and
Henk for the tetrahedron (18/49) and the octahedron (18/19), and the
body-centered cubic lattice for the covering by the octahedron (9/8): an
exact computer proof of the strict local optimality that Dougherty and
Faber proved by hand. Its deep holes have 4 lattice points each, every
one on an edge of its translate (two active facets), so the critical cone
and the copositivity test are needed (4 positive multipliers, tangent
space of dimension 5). The other local optimum found for the octahedron
covering is 18/13.

The tetrahedron covering 125/63 of Fiduccia, Forcade and Zito and of
Dougherty and Faber is found and certified (13 holes, each with exactly 4
lattice points, 6 positive multipliers, tangent space of dimension 3), but
only from 1 of 6 random starts: the other starts end at the local optima
32/15 and 69/32 or stall without an exact polish between 2.05 and 2.24.
The covering landscape is much more rugged than the packing one.

## Limits and future work

* Global optimality is not proved, except for tilings (density 1) and for
  the exact 2-dimensional packing. The multistart only makes it likely.
  The global algorithm in dimension 3 is Betke and Henk's enumeration of
  facet choices for Minkowski's test sets. The pieces are here: the
  feasibility of a choice is a linear program, and the minimization of
  `det` on its affine space is what `PolyNorm_Polish` does, though it
  needs all the local minima and not one.
* Covering in dimension 2 could also be solved exactly by Fejes Toth's
  characterization (the largest centrally symmetric hexagon contained in
  `P`, a nonconvex quadratic program with 8 variables).
* The certificate for covering does not cover holes with more than `n+1`
  lattice points.
* Irrational optima: the critical point is algebraic. An exact
  certificate would need the minimal polynomial (resultants, as in Betke
  and Henk) or interval arithmetic around the floating point solution.
* Symmetry is not used: the starts are random and the covering evaluations
  use the trivial group. Starting from lattices invariant under subgroups
  of the symmetry group of `P` would reduce the search space, as in the
  search for perfect forms.
* The cost is dominated, for covering, by the exact covering radius at
  each trial point (a few hundred evaluations per start). The `flint`
  arithmetic is about three times faster.
