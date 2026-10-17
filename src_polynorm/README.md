Packing and covering with polyhedral norms
==========================================

A polytope `P` in `R^n` with `0` in its interior defines the gauge
`||x||_P = min{t >= 0 : x in tP}`, a norm when `P` is centrally symmetric
and an asymmetric distance function otherwise. This directory computes the
two classical lattice parameters of that gauge for the lattice `Z^n`:

  * **Packing**: the largest `alpha` such that the translates `alpha P + z`,
    `z in Z^n`, have pairwise disjoint interiors.
  * **Covering**: the smallest `mu` such that the translates `mu P + z`,
    `z in Z^n`, cover `R^n`. This is the covering radius `mu(P)` of
    Cslovjecsek, Malikiosis, Naszódi and Schymura, *Computing the covering
    radius of a polytope with an application to lonely runners*,
    Combinatorica 42 (2022).

Both parameters are invariant under translations of `P`, so `P` is first
recentered at the isobarycenter of its vertices. The symmetries of the
problem are then the matrices `A` in `GL_n(Z)` with `PA = P` (vectors are
row vectors and matrices act on the right), and the covering computation
uses them as the Delaunay code does: the search is done up to that group.

A lattice `L` with basis `B` other than `Z^n` is handled by applying
`B^{-1}` to the vertices of `P`, which maps `L` to `Z^n`.

Programs
--------

All programs take the vertices of the polytope in a file in the usual
matrix format, one vertex per row in homogeneous coordinates `(1, v)`, for
example the square `[-1/2, 1/2]^2`:

```
4 3
1 -1/2 -1/2
1 1/2 -1/2
1 -1/2 1/2
1 1/2 1/2
```

Repeated rows and rows that are not vertices are dropped. The polytope
must be full dimensional. The arithmetic `arith` is `gmp` (`mpq_class`)
or, when built with `make ENABLE_FLINT_SUPPORT=1`, `flint` (`fmpq_class`,
faster); every result is exact.

  * **POLYNORM_Packing** `[arith] [FileEXT] [OutFormat] [OutFile]` returns the
    packing scalar `alpha` and the contact vectors, the lattice vectors `z`
    for which `alpha P` and `alpha P + z` touch.
  * **POLYNORM_Covering** `[arith] [FileEXT] [OutFormat] [OutFile]` returns the
    covering radius `mu`, a last covered point `p` (in the coordinates of
    the input) and the lattice points on the boundary of the empty translate
    `p - mu P` together with the facets of `P` they touch. This translate is
    the analogue of a Delaunay polytope: its interior is lattice free and
    the lattice points on its boundary pin it.
  * **POLYNORM_TestCovering** `[arith] [FileEXT] [alpha] [mu] [max_brute_force]`
    runs both computations, compares them with the expected values (or
    `none`) and, when the enumeration has at most `max_brute_force`
    systems, with the brute force algorithm of the paper above.
  * **POLYNORM_OptimalLattice** `[arith] [FileEXT] [problem] [n_start] [seed] [FileStart] [OutFile]`
    searches the lattice `L` giving the densest lattice packing
    (`problem=packing`) or the thinnest lattice covering
    (`problem=covering`) of `R^n` by translates of `P`, from `n_start`
    random lattices (and the basis in `FileStart`, or `none`). Each result
    is polished to an exact rational lattice when possible and its local
    optimality is then certified exactly. In dimension 2 the packing problem
    is also solved exactly and globally. See "Optimal lattices" below.
  * **POLYNORM_BetkeHenk** `[arith] [FileEXT] [OutFile] [n_thread]` computes
    the densest lattice packing of a 3-polytope globally, by the algorithm
    of Betke and Henk (Minkowski's test sets, enumeration of facet choices,
    critical points of `det` on the resulting affine spaces). `arith` is
    `gmp`, `Qsqrt2`, `Qsqrt5` or `RealAlgebraic=File`, so that polytopes
    with irrational coordinates are handled exactly. It reproduces their
    table of the Platonic and Archimedean solids; see `NOTES.md`.
  * **POLYNORM_LatticeCertify** `[arith] [FileEXT] [FileLattice] [OutFile]`
    returns the packing and covering densities of `P` for the lattice whose
    basis is the rows of `FileLattice`, and whether that lattice is an
    exactly certified strict local optimum of each problem.
  * **POLYNORM_MaxKissing** `[arith] [FileEXT] [OutFormat] [OutFile]` computes
    the maximal kissing number of `P`: the maximum, over the lattices `L`
    with `P + L` a packing, of the number of translates `P + v`, `v` in `L`,
    touching `P`, with a witness lattice and its contact vectors. The value
    is scale invariant. See "Maximal kissing number" below.

`OutFormat` is `text` (default) or `GAP`, and `OutFile` is a file name or
`stderr` / `stdout`.

Some values, for the standard polytopes in dimension `n`:

| polytope                        | alpha | mu    |
|---------------------------------|-------|-------|
| cube `[0,1]^n`                  | 1     | 1     |
| cross polytope `conv(+-e_i)`    | 1/2   | n/2   |
| simplex `conv(0, e_1, ..., e_n)`| 1     | n     |
| hexagon `conv(+-(1,0), +-(0,1), +-(1,1))` | 1/2 | 2/3 |

Algorithms
----------

**Packing.** The translates `alpha P + v` and `alpha P + w` have disjoint
interiors exactly when `w - v` is not in the interior of `alpha (P - P)`.
So `alpha` is the first minimum of the difference body `D = P - P` for the
lattice `Z^n`: the minimum of `||z||_D` over the nonzero lattice vectors.
The facets of `D` are computed from the differences of vertices, an easy
upper bound is the gauge of a basis vector, and the lattice points of the
correspondingly dilated `D` are enumerated.

**Covering.** With `f(p) = min_z ||p - z||_P`, the covering radius is the
maximum of `f`, a piecewise linear function that is a minimum of convex
functions. The point `p` is not covered by `mu P + z` when some facet `i`
of `P` has `ell_i(p - z) >= mu`, so the feasible region
`{(p, mu) : mu <= f(p)}` is an intersection, over the lattice points `z`,
of unions of halfspaces, one per facet. This disjunctive program is solved
by an exact branch and bound:

  * a node is a set `C` of pairs `(z, i)` meaning that `p` lies in the cone
    of the facet `i` at `z` and that `mu <= ell_i(p - z)`; its linear
    program maximizes `mu`, which bounds the node from above;
  * at the optimum `(p*, mu*)` the value `f(p*)` is a lower bound and the
    node is solved when it equals `mu*`; otherwise a lattice point `z'` with
    `||p* - z'||_P < mu*` is found among the lattice points of `p* - mu* P`
    and the node is split into one child per facet cone of `z'`;
  * translations are used by requiring `0` to be a nearest lattice point
    of `p`, so the root nodes are one per facet of `P`, and the group
    `{A : PA = P}` reduces them to one per orbit of facets and prunes every
    node equivalent to a processed one under the affine group, through an
    invariant and a canonical form;
  * the root linear programs are bounded by an upper bound `mu0` on the
    covering radius, obtained by fitting translates of the unit simplices,
    cross polytope and cube into `tP` by linear programming.

The brute force of the paper (all systems of `n+1` pairs `(z, i)` with
affinely independent normals) is kept in `PolyNorm_CoveringBruteForce` as a
cross-check; it is exponential in the number of relevant lattice points and
only runs for small cases.

See `NOTES.md` for the design discussion and the open issues.

Optimal lattices
----------------

`POLYNORM_OptimalLattice` looks for the lattice `L = Z^n W` (basis vectors
the rows of `W`) of densest packing `vol(P) / det L` (with `P + L` a
packing) or of thinnest covering (with `P + L = R^n`):

  * the constraints are linear in `W`: `z W d_j^T >= 1` for a contact `z`
    and a facet `d_j` of `P - P`, and for covering the radius of a hole
    pinned by `n+1` pairs is `<G_S, W>`, linear in `W`;
  * a trust region sequential linear programming method on `W`, with every
    trial lattice evaluated exactly by the packing and covering code, and
    a pattern search fallback;
  * the result is polished to the exact critical point of `det` on the
    active constraints and certified exactly: Karush-Kuhn-Tucker
    multipliers, and second order positivity on the critical cone
    (copositivity on its extreme rays);
  * in dimension 2 the packing problem is solved exactly and globally from
    Minkowski's characterization of critical lattices.

This finds and certifies the known optima: Fary's triangle (packing 2/3,
covering 3/2), Sriamorn's quadrilaterals, the Betke-Henk values for the
tetrahedron (18/49), octahedron (18/19), cuboctahedron (45/49) and
truncated tetrahedron (207/304), the octahedron covering 9/8 of Lian and
Xue and the tetrahedron covering 125/63 (best known). Global optimality is not proved by this search beyond dimension 2
(except for tilings); in dimension 3, `POLYNORM_BetkeHenk` gives the
global optimum of the packing problem. See `NOTES.md`, section "Optimal lattices", for the method, the
validation table and the limits.

Maximal kissing number
----------------------

The translates `P` and `P + v` touch exactly when `v` is on the boundary of
the difference body `D = P - P`, so the maximal kissing number is the
maximum number of points of an admissible lattice (no nonzero point in the
interior of `D`) on the boundary of `D`. `POLYNORM_MaxKissing` enumerates
the contact configurations as the Betke-Henk style packing enumeration
does, but without the criticality (the maximum is attained by a lattice
whose contact set spans `R^n`, which is all that the cross-polytope and
determinant conditions need), maximizing the number of contacts instead of
minimizing the determinant. A maximal configuration is decided exactly: a
linear program puts every contact on its facet and in `D`, and a branch
and bound certifies that some lattice of the configuration is admissible,
by the exact enumeration of its points in `D` (a violating point is pushed
out of the interior one facet at a time). The witness lattice, its contact
vectors and `vol(P)/det` are returned, and a configuration exhausting the
branch and bound budget is reported, turning the result into an interval.

The configurations are pruned by the closure of the contact set: a lattice
point `(a +- b) / 2` of two contacts is on the boundary (in `D` by
convexity, not interior by admissibility), so it is forced to be a contact;
a lattice point `(a +- b) / t` with `t >= 3`, or a candidate divisible in
the lattice, is interior and forbids the configuration. The search is
seeded with `Z^n` scaled to its first minimum and stops at the bound
`3^n - 1` of Hadwiger.

Some values: the square and the cube reach `3^n - 1` (8 and 26, the
parallelohedra of maximal kissing), the triangle and the hexagon reach 6,
the octahedron reaches 18, by the lattice `D_3 = {x : sum x_i even}` whose
18 minimal vectors lie on the boundary of `2 P`, and the simplex of
dimension 3 reaches 18, by a lattice of the index two class: the value 18
of Talata for the translative kissing number of the tetrahedron is
attained by a lattice packing. The 2D cases and the cube take seconds; the
exactness proofs for the octahedron and the tetrahedron take minutes.

Files
-----

  * `PolyNorm_Basic.h`: the recentered polytope, its facets as linear forms,
    the gauge, the difference body, the symmetry group as permutations of
    the vertices and as integral matrices, the orbits of facets.
  * `PolyNorm_Packing.h`: the packing scalar and the contact vectors.
  * `PolyNorm_Covering.h`: the upper bounds, the branch and bound and the
    brute force.
  * `PolyNorm_LatticeOpt.h`: the search of optimal lattices, the polishing
    to exact critical points and the certificates of local optimality.
  * `PolyNorm_Packing2D.h`: the exact densest lattice packing of a polygon.
  * `PolyNorm_BetkeHenk.h`: the densest lattice packing of a 3-polytope by
    the algorithm of Betke and Henk.
  * `PolyNorm_Kissing.h`: the maximal kissing number by the contact
    configuration enumeration of `PolyNorm_PackingContacts.h` and the exact
    admissibility branch and bound.
  * The CI test is `CI_tests/25B_PolyNorm/run_tests.sh`.
