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
  * **POLYNORM_Covering** `[arith] [FileEXT] [OutFormat] [OutFile] [method]`
    returns the covering radius `mu`, a last covered point `p` (in the coordinates of
    the input) and the lattice points on the boundary of the empty translate
    `p - mu P` together with the facets of `P` they touch. This translate is
    the analogue of a Delaunay polytope: its interior is lattice free and
    the lattice points on its boundary pin it. The `method` is
    `branch_and_bound` (default) or `brute_force`, the enumeration of all
    the candidate systems of the paper above, for small cases only, which
    gives `mu` and `p` only.
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

  * **POLYNORM_MinimalCoveringPolytope** `[arith] [FileEXT] [OutFormat] [OutFile] [method] [param ...]`
    takes a covering polytope `P` (`P + Z^n = R^n`) and returns a covering
    polytope `Q` contained in `P`, minimal when possible: no convex body
    strictly contained in `Q` is covering (Codenotti, Freyer, Krivokuca,
    *Minimal covering bodies and Brunn-Minkowski type inequalities for the
    covering radius*). `P` is first shrunk to `mu(P) P`. The `method` is
    `grid [N] [n_round]` (default), which keeps the vertices of a
    polytope with vertices on a grid, or `reduction [MaxStep]`, which cuts
    the vertices; the output says whether `Q` is proved minimal. See
    "Minimal covering polytopes" below.

  * **POLYNORM_IsMinimalCovering** `[arith] [FileEXT] [OutFormat] [OutFile]`
    decides whether `P` is minimal covering and, when it is covering, lists
    the vertices that can be cut keeping the covering. It is the local
    criterion of Corollary 2.5 of the paper above, much cheaper than the
    search of a minimal polytope. See "Minimal covering polytopes" below.

  * **POLYNORM_CoveringPath** `[arith] [FileP1] [FileP2] [OutFormat] [OutFile]`
    computes exactly the function `lambda -> 1 / mu((1 - lambda) P1 + lambda
    P2)` on `[0,1]`, which is continuous and piecewise affine: its affine
    pieces and break points, every piece being proved. See "The covering
    radius along a Minkowski path" below.

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

Minimal covering polytopes
--------------------------

A covering polytope `K` is minimal exactly when, for every vertex `v` and
every `u` exposing it, the truncations `K cap {<x,u> <= <v,u> - eps}` are
not covering (Lemma 2.2 of the paper above). Removing vertices is not
enough: the rectangle `[0,3/2] x [0,1]` contains the unit square, but the
convex hull of three of its vertices has area 3/4. So the vertices are cut.

The deepest covering cut `t* = min{t : K cap {<x,u> <= t} is covering}` is
`sup_x h(x)` with `h(x) = min{<x - z, u> : x - z in K}`, computed exactly by
a branch and bound like the covering radius one: after translating, `x` is
in `K` with `<x,u> >= s`, and for every other lattice point `z`, either
`x - z` is in `K` with `<x - z, u> >= s`, or `x - z` is outside `K` through
some facet `j`, in the cone of `j`. A node is solved when no lattice point
of `x* - K` violates its bound. Since `h` is not continuous, only the nodes
whose region has a nonempty interior count (a second linear program).

Whether a covering polytope `K` is minimal is decided locally (Corollary
2.5 of the paper): at each vertex `v`, some direction `d` entering `K` at
`v` must leave `K` at every point `v - z` with `z != 0` a lattice point and
`v - z` in `K`, that is `ell_i(d) < 0` for the facets `i` through `v` and,
for each `z`, `ell_j(d) > 0` for some facet `j` through `v - z`. A depth
first search over the choices of `j`, each partial system checked by a
linear program, decides it; when some `v - z` is interior to `K` the vertex
can be cut at once. This is what `POLYNORM_IsMinimalCovering` does, and
what certifies the results of the grid method.

The method `reduction` applies at each step the reduction of smallest
volume among the removal of a vertex, the deepest translate of a facet and
the deepest cut of a vertex. Each reduction is exact, but in dimension 3
they can go on without end, shallower and shallower (octahedron,
tetrahedron).

The method `grid` (default) works on a finite set `S` of candidate
points: the vertices of `mu(P) P` and its points in `(1/D) Z^n`, with `D`
equal to `N` times the lcm of the denominators of these vertices. A vertex
of `conv(S)` is dropped from `S` when the hull of the rest still covers
(a covering test that stops at the first uncovered point), until none
can be dropped. This takes at most `|S|` tests: a point that cannot be
dropped never can later, as `S` only shrinks. The result is then proved
minimal when its volume is 1 (a lattice tile) or when no vertex can be cut;
otherwise the grids `(1/2D) Z^n`, `(1/3D) Z^n`, ... are tried, for at most
`n_round` rounds (default 3).

Some results of the grid method, all proved minimal (times with the
sanity checks on, Apple M-series):

| `P` | `Q` | time |
|---|---|---|
| `[0,3/2] x [0,1]` | `[1/2,3/2] x [0,1]` | seconds |
| `conv(0, 3e_1, 3e_2)` | a lattice parallelogram | seconds |
| `conv(+-e_1, +-e_2, +-(e_1 + e_2))` | a centrally symmetric hexagon, area 1 | seconds |
| `conv(0, 3e_1, 3e_2, 3e_3)` | the cube `[0,1]^3` | seconds |
| `(3/2) conv(+-e_i)`, `n = 3` | 12 vertices, 14 facets, volume 5/4, not a tile | 20 s |
| `conv(0, 2e_i, (2,2,1), (1,2,2), (2,1,2))` | 15 vertices, 23 facets, volume 7/6, not a tile | 110 s |

The covering radius along a Minkowski path
------------------------------------------

For `K_lambda = (1 - lambda) P1 + lambda P2`, the facet normals `c_i` of
`P1 + P2` describe `K_lambda` for every `lambda` in `[0,1]`, with the right
hand sides `h_i(lambda)` affine in `lambda` (support functions add). With
`q = p / mu` and `s = 1 / mu`, a point `p` uncovered at level `mu` is
`for all z there is i with c_i.q - s c_i.z >= h_i(lambda)`, which is linear
in `(q, s, lambda)`, and `g(lambda) = 1 / mu(K_lambda)` is the minimum of
`s` over this disjunctive set `G`. A hole pinned by `n + 1` pairs gives `(q,
s)` affine in `lambda`, hence an affine piece of `g`, valid on an interval
computed exactly. So `g` is continuous and piecewise affine, while `mu` is
piecewise of the form `1 / (a + b lambda)`.

`POLYNORM_CoveringPath` computes the pieces of the deepest holes at 0 and 1
(by a branch and bound over `(q, s, lambda)` with `lambda` fixed), takes
their minimum `U`, an upper bound on `g`, and certifies each linear stretch
`s_H` of `U` by a branch and bound over `(q, s, lambda)` with `lambda` in
the stretch, minimizing `s - s_H(lambda)` over `G`: no point of negative
value proves `g = s_H` there, and a point of negative value is a hole whose
pieces lower `U`. Finitely many hole systems exist, so this terminates; all
the numbers are rational and exact. The holes are taken at vertices of the
linear programs, so that their tight pairs pin them.

On the parallelograms of Section 3 of Codenotti, Freyer and Krivokuca,
`P1 = {a (1,0) + b (1/2,1)}` and `P2 = {a (1,1) + b (-3/2,-1/2)}`, `0 <= a, b
<= 1`, both lattice tiles:

| `lambda` | `1 / mu` | `mu` |
|---|---|---|
| `[0, 7/24]` | `1 + lambda/14` | `1` to `48/49` |
| `[7/24, 1/2]` | `7/6 - lambda/2` | `48/49` to `12/11` |
| `[1/2, 17/24]` | `2/3 + lambda/2` | `12/11` to `48/49` |
| `[17/24, 1]` | `15/14 - lambda/14` | `48/49` to `1` |

the two local maxima and the local minimum of `1 / mu` that the paper
describes (6 s). The facet cones of the covering radius code are bilinear
in `(q, lambda)` and are not used, so the children of a branching overlap:
two boxes in dimension 3 take 166 s.

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
  * `PolyNorm_MinimalCovering.h`: the deepest covering cut and the
    reduction to a minimal covering polytope.
  * `PolyNorm_CoveringPath.h`: the covering radius of `(1 - lambda) P1 +
    lambda P2` as an exact piecewise affine function of `lambda`.
  * The CI test is `CI_tests/25B_PolyNorm/TestPolyNorm.g`.
