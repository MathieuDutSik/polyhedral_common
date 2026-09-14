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

## Searching for a record periodic covering

Enumerating every domain answers the question outright, and that is what the
two workflows above do. It only scales as far as the domain count does: 222
in dimension 5 for a lattice, but far more for a periodic point set, and the
count grows with the dimension and with the coset structure. Beyond that,
the search has to walk instead of enumerate.

`PERIODIC_LookForRecordCovering` does the walk. Given a coset configuration
it descends on the per-domain covering optimum — moving to the best adjacent
domain, jumping randomly out of local minima — until a domain beats a record
or the runtime budget runs out. It is the strategy of
`LATT_LookForFullRankRayDomain`, which is how the dimension 6 record
covering was found.

```
&SYSTEM
 max_runtime_second = 7200
 Prefix = "/path/hits/"
 OutFile = "/path/result.g"
/
&DATA
 arithmetic = "gmp"
 FileDualDescription = "unset"
 FileCosets = "/path/cosets.txt"
/
&SEARCH
 RecordToBeat = "auto"
 n_walk_steps = 20
/
&TSPACE
 TypeTspace = "Classic"
 ClassicDim = 5
/
```

`RecordToBeat = "auto"` takes the least dense known lattice covering of the
dimension: `A_n^*` up to dimension 5, from its closed form
`kappa_n sqrt(n+1) (n(n+2)/(12(n+1)))^{n/2}`, and the `L^c_n` of Table 2 of
the reference from dimension 6 on — in particular Vallentin's `L^c_6` at
2.464801 in dimension 6. Any other value can be given explicitly.

`OutFile` is a GAP record with `found_record`, `best_density`, the Gram
matrix attaining it, and the counters of the walk. The Gram matrix of every
improvement is also written under `Prefix` as it is found, so a run that is
killed still leaves its best find behind.

**Not finding a record is the expected outcome.** In dimensions 3 to 5 the
best covering is conjectured to be the lattice one, so the program reports
the best density it saw and terminates normally. Running out of budget,
landing on a domain with no flippable wall, and a domain too thin to be
optimized in floating point are all reported rather than thrown.

`max_runtime_second` is honoured at the granularity of one step of the walk,
and a step computes every domain adjacent to the current one before anything
can be checked. In a dimension where a domain has many walls the budget can
therefore be overrun by the cost of one such step, so it bounds the search
loosely rather than exactly.

### Sweeping over configurations

The coset structure is an input, so a search is a sweep over configurations.
`PERIODIC_RandomCosets` draws a random admissible one — checking with the
enumeration's own predicates that the cosets are not a group modulo `Z^d`
and that `-I` preserves the point set — and
`scripts/periodic_covering_sweep.sh` loops the two:

```sh
scripts/periodic_covering_sweep.sh -n 5 -o /tmp/sweep5 -k 20 -t 600
```

It stops on the first record and otherwise reports the best density of every
configuration it tried. Exit status 1 means "swept everything, found no
record", which is the normal outcome, not a failure.

### The drift of the walk

A flip re-expresses the tessellation in the same lattice basis, and nothing
brings it back, so a long walk drifts into ever more skewed representatives
of the domains it visits. On `Z^3 + {0, (1/2,1/2,1/2), (3/4,3/4,3/4)}`,
whose 46 domains have Gram matrices of largest entry at most 108, a
twenty-minute walk was reaching representatives of largest entry 1e15, with
a median of 1e12. Those are the same domains in worse coordinates — the
covering optimum is an invariant and does not change — but the forms
interior to them become anisotropic enough that the optimization cannot be
carried out in double, and 95% of the evaluations of that run were lost.

The enumeration never sees this, mapping every domain it reaches to a
canonical representative; a walk has no such step. Nor can the
representative simply be reduced: a unimodular reduction of the form moves
the cosets with it and would give a different point set, unless taken in the
subgroup preserving the one at hand.

So the walk restarts. `max_gram_entry` in the `SEARCH` block (default 1e6)
is the largest Gram entry tolerated; above it, or when an optimization
fails, the current domain is dropped for a fresh one, whose coordinates are
small by construction. The best found so far is kept across restarts, so
only the position is lost. On the configuration above this takes the
evaluations lost from 95% to none, and the walk then reaches the same
2.0881700162 that enumerating all 46 domains gives.

`n_domain_failed` and `n_restart` in the result record are what to watch: a
run with a large `n_domain_failed` reports the best over a small part of
what it visited, and the sweep script prints both.

### Two traps that turn a search into a wrong answer

**A time-capped enumeration is partial, and used to say so only in a debug
build.** `EnumerateAndStore_Serial` returns the same thing whether it
finished or ran out of time, and nothing in the returned database
distinguishes the two, so a truncated run reports its partial orbit count as
if it were the total. It now warns unconditionally, in every build. Either
read that warning or run enumerations with `max_runtime_second = 0`.

**`srand` does not seed `random()` outside glibc.** The walks, the tie
breaking and `PERIODIC_RandomCosets` all draw from `random()`, whose state
on the BSD derived platforms is separate from `rand()`'s and was left at its
default -- so every run replayed one fixed sequence, the "random" coset
draws returned the same configuration every time, and repeated searches
explored the same trajectory. The programs now call `srandom` beside
`srand`. The shared helper `srand_random_set` of
`basic_common_cpp/src_basic/Basic_random.h` has the same defect and should
be fixed upstream. `PERIODIC_RandomCosets` also takes an optional explicit
seed, which is what a script drawing several configurations in a row needs:
the default seed is the clock plus the pid, and two calls in the same second
are close enough that `rand()` returns the same first values.

### What a record would and would not prove

The optimization is numerical, so a density that lands within a few units of
the last digit of the record proves nothing, all the more so from dimension
6 on where the record itself is a rounded decimal. A hit has to be confirmed
by rounding the reported Gram matrix to a rational one and re-checking it
exactly, which `MaxSquaredCircumRadius` and `IsInCone` are type-generic for.

Note also that a walk that fails proves nothing at all: it visited some
domains, not every one.

## What has been searched in dimension 5

Recorded so the search is not repeated from scratch. All of it is for
periodic point sets `Z^5 + {c_1, ..., c_m}` over the classic `GL_5(Z)`
T-space, against the record `Theta(A_5^*) = 2.1242859090`.

The best found is **2.160060765053**, for `Z^5 + {0, (3/4,1/4,1/2,1/4,0)}` —
1.684% above the record, so **no record was found**. The optimum is
algebraic: every entry of the optimal form lies in `Q(sqrt2)`,

```
(9+2V)/80   (-9+6V)/80   0        V/10        9/80
(-9+6V)/80  (27+18V)/80  -9/80    3V/10      -9/80          V = sqrt2
0           -9/80        9/40     9/80        9/80
V/10        3V/10        9/80    (9+16V)/40   0
9/80        -9/80        9/80     0           9/40
```

and the exact density there agrees with the numerical walk to 4e-16.
Replacing `sqrt2` by the convergent `3363/2378` and clearing denominators
gives an integral form whose tessellation was checked in exact arithmetic:
two orbits of 5-simplices, squared circumradii `113097681/1189` and
`110813083/1189`, `det Q = 176394205620236639052`, point density `1/512`,
density `2.160060765053`.

For `m = 2` the family is one equivalence class per denominator `N` (see
above), so the table below is the whole of it:

| N | 3 | 4 | 5 | 6 | 7 | 8 |
| - | - | - | - | - | - | - |
| best | 2.691 | **2.16006** | 2.507 | 2.200 | 3.23 | 3.35 |

The last two entries had 40 minutes each and are starved rather than
measured: `N=5` fell from 2.676 to 2.507 and `N=6` from 2.282 to 2.200 once
the run length was matched to the `N^5` size of the coset space. Larger `N`
with run length scaled accordingly is the one direction left untested, and
is where an irrational optimal translation would show up.

For `m >= 3`, a screen of 136 configurations was monotone and much worse:
3.16, 3.18, 3.73, 4.61, 4.17, 5.10 for `m = 3..8`.

### Why the picture looks like that

`A_5^*` has **exactly one orbit** of Delone cells, of exact squared
circumradius `35/12`, which reproduces 2.124286. Every hole is equally
deep, so a *local* move gains nothing: a periodic set with `m` points per
fundamental cell multiplies the density by `m` and so needs `mu` to fall by
`m^(1/5)`, and killing a few of the 120 equally deep holes per lattice
point does not come close. That explains the monotone table and the
stalls of the seeded descents.

It is a local argument, not a global obstruction. An optimal covering can
perfectly well have only a subset of its Delone cells at the covering
radius -- the dimension 6 record `L^c_6` does, and beats the balanced
`A_6^*` -- and indeed the best periodic configuration found here, at
2.160060765053, is a local optimum of the continuum problem with one of
its two orbits slack at 0.98. What supports "no periodic covering of
dimension 5 beats `A_5^*`" is therefore only the accumulated search
evidence, all of it converging to 2.160060765053 with nothing below.

### The continuum: freeing the coset position

The searches above fix the cosets and walk over forms. The complementary
question — perturb `c` continuously and descend jointly in `(Q, c)` — is
answered by the numerical tooling of `PythonScript/joint_covering/` (see its
README for the method and the caveats). What it established:

* The certified 2.160060765053 configuration is a **local optimum of the
  full continuum problem**: 576 of its 864 translation classes of Delone
  cells are simultaneously active, and no joint direction improves it. The
  optimal coset is exactly `(3/4,1/4,1/2,1/4,0)` even though nothing
  constrains it to be rational.
* Freeing `c` at the other families' walk optima descends genuinely — the
  N=5 values drop by 0.2–0.27, so there the denominators were binding — but
  the descents drain into a small set of continuum basins (2.2411, 2.3149,
  ...), and the deepest reachable one is again 2.160060765053: the N=6 walk
  optimum at 2.1997 flows exactly there, with `c` moving 0.25 through
  irrational positions.
* Random multistarts of the joint descent land only in shallower basins
  (2.5–2.7): the deep attractors have small basins from generic forms.

Together with the per-family walks, this is strong evidence that
2.160060765053 is the optimum of the whole 2-point periodic family in
dimension 5 — still 1.684% above `Theta(A_5^*)`.

### The double-precision C++ search tool

`src_delaunay/PERIODIC_JointCoveringSearch` (built by
`Makefile_joint_double`, requiring Eigen and qhull) is the fast successor
of the Python prototype: covering density evaluation, joint local descent
on `(Q, C)` and multistart, entirely in double precision, for any dimension
and any number of cosets. The tessellation is qhull's Delaunay of the
Q-ball of points enumerated through an LLL-reduced basis, each retained
cell certified by its circumsphere lying inside the ball; the optimization
minimizes the smooth surrogate
`(n/2) log softmax_beta(R^2_S) - (1/2) log det Q` by L-BFGS with analytic
gradients (finite-difference-checked to 6e-9), annealing beta and
re-tessellating between stages, every stage verified against a fresh
tessellation and the input configuration held as the first incumbent so a
descent can never return something worse than its seed.

Validated digit for digit against the exact values 2.1600607650528 (m=2)
and 2.3765773479972 (m=3). Candidates found by this tool are numerical and
must be re-verified by the exact pipeline.

**Two descent engines.** The original `descend`/`multistart` use a soft-max
surrogate minimized by L-BFGS with analytic gradients. The faster
`descend-alt`/`multistart-alt` exploit the block structure of the problem
directly, following the observation that

* fixing the cosets `c`, the constraint `r(Delta_i)^2 <= h` is an LMI linear
  in `(Q, h)` (Delone-Dolbilin-Ryshkov-Stogrin), so the `Q`-block is an exact
  convex SDP;
* fixing `Q`, minimizing `max_i r(Delta_i)^2(c)` is a smooth minimax whose
  nonconvexity is confined to `c`.

`descend-alt` alternates an interior-point `Q`-step (the exact per-`(c,T)`
optimum, `QStep`, gradient-checked to 1e-7 and holding at the certified
optimum) with a `c`-step (`CStepMinimax`: steepest descent for the maximum,
the direction being the min-norm element of the active gradients' convex
hull found by Frank-Wolfe, with an Armijo line search on the true maximum),
re-tessellating once per round. It removes the beta-annealing of the
soft-max engine entirely and, from a well-rounded seed, converges in ~10
rounds at a few seconds each.

**The real bottleneck is the tessellation, not the optimizer.** A Delaunay
tessellation of a well-rounded `(Q,c)` costs a few seconds, but an
anisotropic form or near-cocircular cosets send the qhull-plus-ball-growth
loop into tens of seconds or worse. Both engines pay this, and it dominates
the run time. The mitigations in place -- LLL-reduced enumeration, an
anisotropy guard on the `Q`-step output (condition number `< 300`), a
per-tessellation deadline and a per-descent wall-clock budget -- keep the
search robust (no infinite stalls) but do not remove the cost; a genuinely
fast search would need incremental Delaunay updates or a covering-radius
evaluation that avoids a full tessellation. `multistart-alt` seeds from
well-rounded forms (`Q = R^T R`, `R = I + small`), where covering-optimal
forms live and tessellations stay cheap.

Both engines confirm, as expected, that random multistart lands in the
shallow basins (2.5-3.8 for `m=2`) and only rarely reaches the deep 2.16006
basin: the difficulty is the basin structure, which no change of optimizer
removes.

**On the optimizer.** Two structured alternatives were tried and both are
regressions against the soft-max L-BFGS: an alternation (exact SDP `Q`-step +
a coset step) and a joint first-order minimax step. The alternation jams --
alternating minimization stalls at non-stationary points of a non-smooth
objective, and the covering radius `max_i r_i^2` is non-smooth, so a joint
descent direction that must move `Q` and `c` together is invisible to it
(demonstrated: it stalls at 2.2435 where L-BFGS reaches 2.2301 and escapes
its stall point). The joint first-order step avoids the jam in principle but
loses to L-BFGS in practice, because L-BFGS's quasi-Newton curvature handles
the very different scales of the `Q` and `c` variables automatically while a
plain steepest-descent does not. The conclusion is that the soft-max L-BFGS
was already the right method -- joint, smooth, quasi-Newton -- and the effort
belongs on the tessellation, not the optimizer.

**Incremental re-tessellation.** The dominant cost is the Delaunay
tessellation, and most of it is wasted: after a small `(Q,c)` move the
previous cell list is usually still the Delaunay triangulation. `TryReuseCells`
checks this cheaply -- no simplex degenerated, and no point of the set lies
strictly inside any simplex's circumsphere (the empty-sphere property, which
a flip would violate) -- and reuses the cells when it holds, falling back to
a full qhull recompute only when a flip actually occurred. It is exact (the
reused list is genuinely the Delaunay triangulation) and gives bit-identical
results; on an anisotropic seed where tessellation is expensive it cut a
descent from 375 s to 151 s (2.5x), reusing 5 of 12 tessellations, and it
helps most exactly where the tessellation is slowest. It is wired into the
L-BFGS `descend`/`multistart`.

## Testing

`CI_tests/27B_CoveringMaxdet` runs the whole pipeline in dimensions 3, 4 and 5
for lattices and on `Z^3 + {0, (1/3,1/3,1/3)}` for the periodic case, and
checks that the random-walk search recovers on that same point set the
optimum the full enumeration gives; see its `README.md`.
