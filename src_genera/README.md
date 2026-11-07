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
* Deduplication is by canonical form (`ComputeCanonicalFormWeylChamber`,
  the characteristic vectors in one Weyl chamber and the simple roots), so
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

## The mass of a genus of lattices with an isometry

`GENUS_IsometryMass` computes, for an even positive definite L and an
isometry g of finite order (g^T A g = A), the mass

    mass(L, g) = sum over the classes of pairs (L', g) in the genus of the
                 pair of 1 / |C_{O(L')}(g)| ,

the certificate an enumeration of the g-stable lattices needs.

    GENUS_IsometryMass mass  [Gram] [Isometry] [HermMass] [OutFile]
    GENUS_IsometryMass check [Gram] [Isometry] [HermMass] [OutFile]

With L_d = L cap ker Phi_d(g) and L_0 = (+)_d L_d, the stabiliser of L in
the adelic centraliser is a subgroup of finite index of that of L_0, so

    mass(L, g) = prod_d mass(L_d, g) * prod_q |G_q . H_q| ,

H_q the glue L / L_0 at q in the discriminant group of L_0 and G_q the
image of the local stabiliser of L_0 (see `IsometryMass.h`). The masses of
the orthogonal components (d = 1, 2) come from `GenusMass`; those of the
Hermitian components (d >= 3, lattices over Z[zeta_d]) are given in
`HermMass` as lines `d num den`, from Hecke for instance (`none` if there
are none). The glue factor is computed exhaustively: the isometries of
each primary discriminant module commuting with g are listed by
backtracking, and the orbit of the glue is closed under them. For d = 1, 2
these are the image of the local group (Nikulin, Theorem 1.9.5); for d >= 3
the code takes all of them, and says so (`image_assumed`).

Limits: the exhaustive glue computation stops (with an error) beyond
2^20 elements in a discriminant module, 2 10^6 isometries or 2 10^6
subgroups in an orbit. A large 2-elementary glue, as for an involution
with two large eigenspaces in rank 19, is beyond them; it needs the orders
of the finite orthogonal groups and of their parabolic subgroups rather
than a listing.

Validation (`check` mode): the elements tau of O(L) with the
characteristic polynomial of g are listed, grouped by the mass of (L, tau)
and the genera of the components, and for a genus of class number one the
count of each group divided by |O(L)| times the mass must be the number of
genera of pairs in it. On A1^3, A3, A2+A2, D4, D5 and E6 with every
product of orthogonal reflections, and on A2, D4 and E6 with elements of
order 3, 4, 5 and 6 whose Hermitian components have rank one (mass 1/6,
1/4, 1/10 over Z[zeta_3], Z[i], Z[zeta_5]), every ratio is exactly 1,
with glue factors up to 6, including cases where `image_assumed` is set.

## Enumeration of the lattices with an isometry

`GENUS_EnumerateIsometry` enumerates the classes of pairs (L, g) in the
genus of a pair, g an isometry of finite order m, and stops when the mass
computed by `GENUS_IsometryMass` is reached.

    GENUS_EnumerateIsometry [Gram] [Isometry] [Mass] [OutFormat] [OutFile]
                            [prime] [max_prime]

The code is in `GenusEnumerationGroup.h` and reuses `GenusEnumeration.h`
without changing its cost: `GenusNeighbor` is split into
`GenusNeighborKernel`, which also returns the basis of the neighbour, and
the reduction.

* The neighbours are the sigma-stable ones: at a good prime p (p not
  dividing m det), N = L_v + Z v/p for v an eigenvector modulo p^2 of an
  eigenvalue in F_p, obtained by the projector R(P) / R(lambda') with
  lambda' the Teichmueller lift and R = (x^m - 1) / (x - lambda'). For
  lambda != +-1 every eigenvector is isotropic; for lambda = +-1 the
  isotropic lines are corrected modulo p^2 inside the eigenlattice.
* The lines are cut into orbits of the centralizer C_{O(L)}(g), eigenspace
  by eigenspace, the small ones first.
* A pair is identified by the canonical form of (A, A g^T) and its weight
  is 1 / |C_{O(L)}(g)|, the automorphism group of that pair of forms.
* The primes: 2 when it is good (m and det odd: the fixed space over F_2 is
  the cheapest projective space), then the good primes p = 1 mod m, at
  which all the eigenvalues are in F_p. A single prime can be imposed.

The seeds (one per genus of pairs and per component of the neighbour
graph) are input. `CI_tests/15B_GenusEnumeration/RunEnumerateIsometryTests.sh`
checks the even unimodular genus of rank 16 with -1, a root reflection and
an element of order 3, and the genus of rank 8 and determinant 25 with -1
and a root reflection.
