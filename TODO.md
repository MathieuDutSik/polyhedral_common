# TODO

## Dual description: speed up the double description (the workhorse)

The double description (`src_poly/POLY_DualDesc_double_description.h`) is what
copes with the degenerate cut / metric / TSP / perfect-form cones that dominate
this codebase (lrs reverse-search explodes on them). It is COMBINATORICS-bound:
the Fukuda-Prodon adjacency bitset scan dominates, not arithmetic (measured -
mpq and int64 kernels are within noise; a TryInt64 attempt/fallback was tried
and reverted as neutral). So the levers are all combinatorial/structural, taken
from reading normaliz `source/libnormaliz/full_cone.cpp` (which computes support
hyperplanes with the polar dual of our kernel). Ranked by expected payoff:

1. **Parallelization** -- the single biggest practical lever, and our
   architecture already suits it. The sign-classification loop and the
   positive/negative pair matching are embarrassingly parallel (per-thread
   `new_rays` buffers merged after the loop). The symmetry Adjacency
   Decomposition is also embarrassingly parallel across orbit representatives.
   For our scale this dwarfs any single-core tweak (cf. mplrs near-linear MPI
   scaling; normaliz parallel pyramid evaluation).
2. **Simplicial-ray fast path.** Rays whose tight set has exactly d-1 rows
   ("simplicial") admit direct enumeration of their (d-2)-subsets; matching
   positive/negative simplicial pairs through a hash map replaces the quadratic
   pair loop on the simplicial part (normaliz: `Pos_Simp`/`Neg_Simp` deques in
   `find_new_facets`). `zero_count` already classifies rays for free.
3. **Adaptive adjacency test.** Complement the Fukuda-Prodon comparison test
   with a fraction-free rank test (rank of the common tight rows must be d-2)
   and choose per-pair by a cost estimate. The rank test wins in the degenerate
   regime (large tight sets, long `row_tight` candidate lists, CUT polytopes);
   comparison wins when the answer is mostly "not adjacent" (normaliz
   full_cone.cpp ~line 2611).
4. **Pyramid decomposition** (Bruns & Ichim, J. Symbolic Computation 74, 2016).
   Past a pair-count threshold, replace global matching: for each negative
   ray r, the new rays adjacent to r are the extreme rays of the subcone cut
   out by the inequalities tight on r plus the inserted inequality -- a small
   recursive dual description. Independent subproblems, parallelizes well.
   Largest restructuring, so last; but it is THE algorithmic win for degenerate
   cones (what lets normaliz survive the pos x neg explosion).
5. **Cheap experiment: insertion order.** We insert in lexicographic-minimum
   order; normaliz sorts generators by degree / zero counts
   (`sort_gens_by_degree`). Intermediate ray counts depend strongly on the
   order; benchmark a few orders on the corpus below.

Re-run the benchmark trio after each step. Evidence (2026-08-04, single-thread,
normaliz 3.11.1 no-OpenMP; s; ContactE8 240x9/19440, Example3 48x21/11432,
CUT_K7 64x22/116764):

| method                   | ContactE8 | Example3 | CUT_K7 |
|--------------------------|-----------|----------|--------|
| internal cdd (mpq_class) | 0.43      | 0.49     | 51     |
| normaliz                 | 0.13      | 1.10     | 36     |

normaliz wins where intermediate cones are simplicial-rich (ContactE8, 3x) and
on degenerate CUT_K7 (1.4x), loses on Example3 (2.2x): its edge is algorithmic
and bounded, worth importing but not a rewrite-level gap.

## Dual description: other algorithm families (literature survey 2026-08)

Honest framing: the vertex/facet enumeration ("representation conversion")
problem has NO known algorithm polynomial in input+output size (open problem).
There is no untapped serial family that would revolutionize facet enumeration
of degenerate symmetric polytopes -- the field's answer for exactly our inputs
IS double description + symmetry decomposition + pyramid decomposition +
parallelism, which we have or are planning above. Fourier-Motzkin is correctly
excluded (intermediate-redundancy explosion). Ranked verdict on the rest:

- **Primal-dual driven by our internal DD as the oracle** (Bremner-Fukuda-
  Marzetta, DCG 1998). We have pd_lrs (primal-dual with lrs as checker); the
  distinct idea is to use our fast ring-arithmetic DD as the LP/adjacency
  oracle instead. Payoff: avoids the intermediate non-facets that pure
  insertion generates -- output-sensitive-ish when facets << intermediate ray
  count. Moderate novelty, low effort (we own the pieces). Worth a prototype
  in POLY_DirectDualDesc.
- **Lawrence signed decomposition / volume-by-descent** (Normaliz volume
  paper, 2021). Only if exact VOLUME becomes a goal: volume as a signed sum
  over simplices from one apex, no stored triangulation. Not facet enumeration.
- **Skip (assessed, not worth it for our workload):**
  - Normaliz project-and-lift / dual mode -- solve lattice points / Hilbert
    bases, a DIFFERENT problem, not dual description.
  - BDD-based 0/1 enumeration (Behle-Eisenbrand azove) -- enumerates 0/1
    VERTICES given inequalities; for cut polytopes the vertices are already
    known, the hard direction is facets. Niche.
  - Optimization-oracle combinatorial generation (Merino-Mutze, SIAM J Comp
    2023, O(t_LP log n) delay) -- elegant and exact but VERTICES only, needs
    combinatorial structure; does not touch facets.
  - Gift-wrapping / output-sensitive oracle hulls (Emiris-Fisikopoulos) --
    floating-point / low-dimensional / projection-oriented, subsumed.

Refs: Bremner-Dutour Sikiric-Schurmann (representation conversion up to
symmetry, arXiv:math/0702239); Bruns-Ichim (pyramid decomposition, JSC 74);
Bremner-Fukuda-Marzetta (primal-dual, DCG 20:333, 1998); Avis-Jordan (mplrs,
MPC 2017; comparative results, arXiv:1510.02545); Behle-Eisenbrand (BDD 0/1,
ALENEX 2007); Merino-Mutze (arXiv:2304.08567); Seidel (Handbook chap. 26).

## Integer arithmetic of the dual-description kernels (normaliz / lrs / bb / dd)

Background: the kernels currently run their machine-integer fast path over the
deferred-overflow TryCarryInt64 (exact per-operation detection, exact-ring
fallback). Micro-benchmarks (Apple Silicon) put the per-operation checking
cost at 1.6x-3x over raw int64_t on the hot loop shapes (scalar product, FM
combination, Bareiss row); end-to-end the arithmetic share is smaller, so the
observed whole-run overhead against a raw-int64 implementation is ~1.1x-1.4x.
The reference Normaliz long-long mode is NOT rigorous (unchecked hand-unrolled
scalar products; range checks applied to results after products that may
already have wrapped, which is UB); we do not want to copy that. The items
below are the rigorous ways to close the gap.

* **Hadamard-certified raw int64 for dimension <= 15.** For a submatrix with
  entries bounded by G, every Bareiss intermediate is a product of two
  consecutive minors, so Had_k(G) * Had_{k+1}(G) < 2^63 with
  Had_k(G) = (G*sqrt(k))^k certifies the whole elimination a priori --
  no runtime checks at all. With G = 1 this holds up to dimension ~15
  (fails at 16: ~1.4e20). Same regime covers scalar products and FM
  combinations (facet coefficients bounded by Had_{d-1}(G)). Applicable to
  all four backends on low-dimensional (sub)problems; useless at d ~ 98,
  where the Hadamard bound is ~1e96 while the measured coefficients on the
  P_9_6 hard family are <= 96.

* **Precondition-gated raw int64_t for high dimension.** Runtime-gated, but
  rigorous (contrast Normaliz): verify the precondition BEFORE running a
  block raw, never check results after possibly-wrapped operations.
  Invariants: generator entries <= G (checked at input conversion), stored
  facet coefficients <= B (checked at the number_hyperplane choke point,
  after gcd reduction), with B = sqrt(2^62 / (dim * G)); then every
  intermediate of a scalar product (<= dim*G*B) and of an FM combination
  (<= 2*dim*G*B^2) is provably wrap-free. For Bareiss eliminations track the
  max workspace entry E per step (one compare per write) and require
  2*E^2 < 2^63 before each row update; violation falls back per call
  (rank: to the comparison test), per object, or per run (to TryCarryInt64,
  then the exact ring). Measured headroom on the P_9_6 hard family:
  facet coefficients <= 96 against a budget of ~2.2e8 (the
  NMZ_COEFF_STATS instrumentation in POLY_DualDesc_normaliz.h measures
  this; the max-Bareiss-entry counterpart is still to be measured).
  Expected gain: the arithmetic share only, est. 5-15% on the hard family.

* **NEON widening multiplies for the classification loop.** AArch64 NEON has
  no 64x64 vector multiply (and AVX2 none either; vpmullq is AVX-512DQ), so
  int64 loops cannot vectorize on current targets -- this is why
  TrySimdInt64 cannot beat TryCarryInt64 there (measured: it is 10-30%
  slower; kept selectable via POLY_NORMALIZ_TRY_SIMD / POLY_LRS_TRY_SIMD for
  future ISAs). But 32x32->64 widening multiplies (smull/smlal) DO exist:
  when the tracked bounds certify coefficients < 2^31 (P_9_6 family: <= 96),
  storing facet normals as int32 makes the scalar products of the
  classification loop vectorizable 2-4 wide with 64-bit accumulation.
  To be pursued after the precondition-gated int64 mode exists, since it
  reuses the same bound bookkeeping.

## Method-selection reference (for the heuristics)

On CI_tests/23B_SimpleDualDesc: cdd is fastest; lrs is far behind on the
degenerate cut/metric/TSP polytopes (reverse search explodes -- lrs on CUT_K7/K8
does not finish); bb sits between (Example6 METP6: bb 0.18s vs lrs 5.4s). The
heuristics should keep steering degenerate / high-incidence cases away from lrs,
to cdd (or bb as fallback).

## Benchmark corpus

scratch bench/ (regenerate from CI_tests + the session's python generator):
24cell, ContactE8, cross8, cube5/6, CUT_K7, CUT_K8, cyclic4_12, cyclic6_12,
Example3, rand01_10_40, and CI_tests/23B_SimpleDualDesc Example4-7. CUT_K8 is
the ~90s gold-standard cdd case; lrs corpus is the _lrs_-marked SimpleDualDesc
entries (Example6/7 overflow int64) plus ContactE8 and Example3.
