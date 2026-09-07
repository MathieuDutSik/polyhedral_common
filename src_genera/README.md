# src_genera: enumeration of a genus of even positive definite lattices

`GENUS_Enumerate` lists the isometry classes of a genus by Kneser
p-neighbours, deduplicating by canonical form and stopping when the
Smith-Minkowski-Siegel mass is reached.

    GENUS_Enumerate [arith] [Genus] [Lattice] [Mass] [OutFormat] [OutFile]

* `Genus`   -- key/value lines `rank`, `det`, and optionally `prime`
               (0 or absent selects the smallest prime not dividing det).
* `Lattice` -- one or more Gram matrices in the ListMatrix format, the seeds.
* `Mass`    -- the mass of the genus as `num den`, or a single integer.
* `OutFormat` -- `GAP` (default), `CPP`, or `Summary`.

## Why this is not an adjacency scheme

The enumerations in `src_perfect`, `src_delaunay` and friends are built on
an adjacency relation: orbits are joined by flipping across a facet, and the
enumeration is over when the orbit list is closed under that relation. Their
`DataFunc` machinery encodes exactly that, and it does not apply here. In a
genus every class is a p-neighbour of some other after enough steps, there is
no facet, and closure of the neighbour graph is not what certifies the answer.

What certifies it is arithmetic:

    mass(G) = sum over the classes L of the genus of 1 / |Aut(L)| .

The enumeration accumulates that sum and stops on equality. Since the mass is
supplied independently, this is a proof of completeness rather than a
termination heuristic: an error in the neighbour construction cannot return a
short list silently, it shows up as a mass that never closes.

## Spinor genera

The p-neighbour graph is connected on a SPINOR genus, not on a genus. A genus
with several spinor genera therefore needs one seed lattice per spinor genus,
which is why the `Lattice` file takes a list. Given a single seed on such a
genus the enumeration converges to a proper subset and the mass falls short;
the program then reports `complete = false` and says why rather than
presenting a partial list as an answer.

## Getting the input

The genus symbol, a representative and the mass are all cheap to obtain from
Hecke (or Sage, or Magma) and are the theory-heavy part; the expensive part is
the enumeration, which is what this program does. Timings for
`representative(G)` and `mass(G)` in Hecke on rank 14: 0.02 to 0.8 s per genus.

## Performance notes

* The neighbour is LLL-reduced before anything else touches it. The basis
  coming from `GetZbasis` divided by p is arbitrary and can be very skewed,
  and both the canonical form and the automorphism group start by enumerating
  short vectors.
* The projective points are cut into orbits under the automorphism group
  before neighbours are formed. In rank 14 at p = 2 that is 16383 points, and
  on the determinant-243 seed they fall into 8 orbits.
* Deduplication is by canonical form (`ComputeCanonicalFormFullRank`), so
  recognising an already known class is a dictionary lookup and not a test
  against every class found so far; the pairwise scheme is quadratic in the
  class number, which is what makes it the bottleneck of a large enumeration.
* Both the canonical form and the automorphism group work from the FULL RANK
  invariant vector family rather than from one spanning Z^n. On the rank-14
  lattices of interest those families have 2358 and 24702 members, the gap
  coming from a single index-2 obstruction (see `TestData/SlowCanonic`), and
  it decides whether the computation is feasible: the extraction alone drops
  from 1.3 s to 13 ms.
* `|Aut|` is read off the permutation action of the generators returned by
  `GetIntAutomorphism_ListMat_Vdiag`, which restricts to the transformations
  preserving the LATTICE. Reading it off the raw permutation group of a full
  rank family instead counts rational isometries that do not preserve the
  lattice: on one class of the determinant-243 genus that gave 713451110400
  instead of 356725555200.

## Validation

`TestData/` holds three genera generated with Hecke. The determinant-243
rank-14 genus is checked class by class against a certified complete list
(21 classes, every `|Aut|` known), which is how both bugs listed above were
found.
