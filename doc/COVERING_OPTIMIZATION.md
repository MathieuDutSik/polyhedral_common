# Optimizing the sphere covering density

This documents the pipeline that computes the best sphere covering density
attainable by a lattice, or by a periodic point set with a prescribed coset
structure, in a given dimension.

The mathematics is Section 6.2 of

> Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin,
> *A generalization of Voronoi's reduction theory and its application*,
> Duke Mathematical Journal 142 (2008) 127--164, arXiv:math/0601084

and the implementation is `src_delaunay/CoveringMaxdet.h`, whose header
comment is the detailed reference.

## The idea in one paragraph

The covering density `Theta(Q) = kappa_d mu(Q)^d / sqrt(det Q)` is invariant
under scaling of `Q`, so minimizing it means normalizing `mu(Q) <= 1` and
maximizing `det Q`. Inside one iso-Delaunay domain the combinatorics of the
Delone subdivision is fixed, and there `mu(Q) <= 1` is a *linear matrix
inequality*: for a simplex `conv{0, v_1, ..., v_d}` spanned by vertices of a
Delone polytope, its circumradius is at most 1 exactly when

```
        /  4         (v_1,v_1)  ...  (v_d,v_d)  \
        | (v_1,v_1)  (v_1,v_1)  ...  (v_1,v_d)  |
BR(Q) = |    .           .       .       .      |  >= 0,   (y,z) = y^T Q z,
        \ (v_d,v_d)  (v_d,v_1)  ...  (v_d,v_d)  /
```

every entry being linear in `Q`. One block per orbit of Delone polytopes,
plus one `1 x 1` block per facet of the L-type cone, turns the covering
optimization over the domain into a determinant maximization problem

```
minimize  -log det Q(x)   subject to   F(x) >= 0.
```

Minimizing over *all* domains is then: enumerate the iso-Delaunay domains,
solve this over each, take the best.

Nothing in the argument needs the point set to be a lattice — only that the
vertices of the Delone polytopes are known and that `Q` is the only unknown.
For the periodic point sets of `src_delaunay/PeriodicStructures.h` the cosets
are fixed rational vectors, so the very same formulation applies, and this is
what makes the search for periodic coverings beating the lattice ones
possible.

## The solver

The solver is a barrier path-following method inside `CoveringMaxdet.h`. It
depends on nothing outside this repository, so the covering optimization is
available in a plain build.

The solve is **numerical**. The problem data (simplices, LMI blocks, cone
facets) is built in exact rational arithmetic, but the optimum is a floating
point Gram matrix and the density evaluated at it is an estimate, not a
certificate. Turning a promising domain into a proven record needs the
optimum rounded to a rational `Q` and re-checked exactly;
`MaxSquaredCircumRadius` and `IsInCone` are type-generic for exactly that and
are the entry points for such a follow-up.

## The lattice workflow

1. Enumerate the iso-Delaunay domains, dumping each of them:

   ```
   &DATA
    arithmetic = "gmp"
    FileDualDescription = "unset"
    CommonGramMat = "unset"
    PrefixIsoDelaunayDomains = "/path/dom_"
    CVPmethod = "SVexact"
   /
   &TSPACE
    TypeTspace = "Classic"
    ClassicDim = 5
   /
   ```

   ```sh
   src_delaunay/LATT_SerialLattice_IsoDelaunayDomain enum.nml
   ```

   This writes `dom_0`, `dom_1`, ... — one boost archive per domain — and the
   count into the `OutFile`.

2. Optimize over each domain:

   ```
   &SYSTEM
    OutFile = "/path/ana_0.out"
   /
   &DATA
    arithmetic = "gmp"
    FileIsoDelaunay = "/path/dom_0"
   /
   &QUERIES
    FileCoveringOptimum = "/path/cov_0.g"
   /
   ```

   ```sh
   src_delaunay/LATT_AnalysisIsoDelaunay ana_0.nml
   ```

   `cov_0.g` is a GAP record with `covering_density`, `covering_radius_sq`,
   `det`, `gap_bound` and the optimal `GramMat`; the `OutFile` repeats the
   headline numbers as `covering_density=...` lines.

3. The answer for the dimension is the smallest `covering_density` over the
   domains. This reproduces the classical values:

   | dimension | domains | best density | attained by |
   | --------- | ------- | ------------ | ----------- |
   | 3         | 1       | 1.4635030689668180 = 5 pi sqrt(5) / 24 | `A_3^*` |
   | 4         | 3       | 1.7655285081493524 | `A_4^*` |
   | 5         | 222     | 2.1242859089916246 | `A_5^*` |

   Beware of a capped `max_runtime_second` in step 1: the enumeration then
   stops early, reports a *partial* list of domains and still terminates
   normally (the log says "returning true (valid partial result)"). The
   minimum over a partial list is an over-estimate of the best density, so a
   run meant to answer the question has to be given an uncapped budget.

## The periodic workflow

The same, with `LATT_SerialPeriodic_IsoDelaunayDomain` in step 1. Its `DATA`
block additionally takes

* `FileCosets` — the rational coset matrix of the periodic point set, one
  coset per row, the zero coset included;
* `PrefixIsoDelaunayDomains` — the per-domain dump, as above;
* `FileLinSpaceOut` — the T-space actually used, written out for step 2.

Step 2 then has to be told both of those, because the defining inequalities
of the domains live in the coordinates of that T-space, and because the
covering density of a periodic point set differs from the lattice one:

```
&DATA
 arithmetic = "gmp"
 FileIsoDelaunay = "/path/perdom_0"
 FileLinSpace = "/path/linspa.txt"
 FileCosets = "/path/cosets.txt"
/
```

With `m` cosets of common denominator `N`, the density is the lattice formula
multiplied by `m / N^d`. `LATT_AnalysisIsoDelaunay` recovers that factor from
the coset file and reports it as `point_density` in the output record.
**Omitting `FileCosets` on a periodic domain silently computes the lattice
density instead**, which is wrong by exactly that factor.

The coset structure is an input, not something the enumeration searches over:
each choice of `(m, N, cosets)` is a separate run. Searching for a periodic
covering better than the lattice record therefore means sweeping over coset
structures, and for each of them over its iso-Delaunay domains.

Two constraints narrow which coset structures can be given:

* The cosets must not form a group modulo `Z^d` — such a point set is a
  lattice and has to be treated as one. Both periodic programs reject it with
  a message saying so.
* Every generator of the pointwise stabilizer of the T-space has to admit an
  affine extension preserving the point set, which the enumeration needs for
  its equivalences to be the right ones. For the classic `GL_d(Z)` T-space
  that stabilizer is generated by `-I`, so the coset set has to satisfy
  `-C = C` modulo `Z^d` up to a translation. `Z^3 + {0, (1/3,0,0)}` passes
  (`-1/3 = 2/3` and the set is `{0, 1/3}` translated), whereas
  `Z^3 + {0, (1/3,0,0), (0,1/3,0)}` does not and is rejected; adding the
  opposites, `{0, (1/3,0,0), (2/3,0,0), (0,1/3,0), (0,2/3,0)}`, restores the
  symmetry.

## Testing

`CI_tests/27B_CoveringMaxdet` runs the whole pipeline in dimensions 3, 4 and 5
for lattices and on `Z^3 + {0, (1/3,1/3,1/3)}` for the periodic case; see its
`README.md`.
