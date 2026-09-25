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

Files
-----

  * `PolyNorm_Basic.h`: the recentered polytope, its facets as linear forms,
    the gauge, the difference body, the symmetry group as permutations of
    the vertices and as integral matrices, the orbits of facets.
  * `PolyNorm_Packing.h`: the packing scalar and the contact vectors.
  * `PolyNorm_Covering.h`: the upper bounds, the branch and bound and the
    brute force.
  * The CI test is `CI_tests/25B_PolyNorm/run_tests.sh`.
