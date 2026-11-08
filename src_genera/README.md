# src_genera: enumeration of a genus of even positive definite lattices

`GENUS_Enumerate` lists the isometry classes of a genus by Kneser
p-neighbours, deduplicating by canonical form and stopping when the
Smith-Minkowski-Siegel mass is reached.

    GENUS_Enumerate [file.nml]

The namelist has four blocks; run the program without argument for the
full list of entries, their defaults and their documentation.

    &DATA
     arithmetic = "gmp"
     FileLattice = "E8_lattice.txt"
     FileMass = "E8_mass.txt"
    /

    &METHOD
     strategy = "iterate"
     prime = 0
     recognition = "canonical"
     canonical_method = "weyl_chamber"
     isomorphism_method = "plesken_souvignier"
     automorphism_method = "weyl_chamber"
     hecke = F
     order = "mixed"
    /

    &RANDOM_WALK
     seed = 1
     hits_factor = 0.25
     failed_draws_factor = 16
     walk_prime_min = 2
     walk_prime_max = 2
    /

    &SYSTEM
     OutFormat = "Summary"
     OutFile = "stdout"
    /

* `DATA` -- the seeds, one or more Gram matrices in the ListMatrix format
  (one per spinor genus), and the mass of the genus as `num den` or a single
  integer. The seeds must have the same rank and determinant.
* `METHOD` -- the `strategy`: `out_prime` (the default, one prime, `prime`
  or the smallest one not dividing det), `iterate` or `random_walk` (see
  `NeighborStrategy` in `GenusEnumeration.h`); the `recognition` of a known
  class, by its canonical form (`canonical`) or by the isometry test within
  the buckets of a cheap invariant (`isomorphism`), with the
  `canonical_method`, `isomorphism_method` and `automorphism_method` of
  `src_latt/LatticeMethods.h`, the same methods as `LATT_Canonicalize`,
  `LATT_Isomorphism` and `LATT_Automorphism` (see `GenusMethods`); `hecke`
  to skip the expansions whose neighbours are all expanded classes (see
  `GenusEnumerationState`); and the `order` of the expansions, `creation`,
  `deficit`, `aut` or `mixed` with its proportions `mixed_deficit` and
  `mixed_aut` (see `GenusExpansionOrder`).
* `RANDOM_WALK` -- the parameters of the `random_walk` strategy.
* `SYSTEM` -- `OutFormat`, `GAP` (default), `CPP` or `Summary`, and
  `OutFile`, a file name, `stderr` (default) or `stdout`.

A block that is omitted takes its defaults. The file names are relative to
the directory the program runs in: `TestData/E8.nml` is run from
`TestData`.

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

* The neighbour is LLL-reduced once before anything else touches it. The
  basis coming from `GetZbasis` divided by p is arbitrary and can be very
  skewed, and both the canonical form and the automorphism group start by
  enumerating short vectors. The reduction is done by the dispatchers of
  `src_latt/LatticeMethods.h`, which reduce their input once and give the
  result back in the basis of the input, for every method.
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
    GENUS_IsometryMass group [Gram] [ListIsometry] [OutFile]
    GENUS_IsometryMass groupcheck [Gram] [Isometry] [OutFile]

The group modes are for a finite abelian group G generated by commuting
isometries: the components are then the L cap V_chi for the Galois orbits
of characters chi, on which G acts through one element g_chi, and the mass
is computed as for g_chi on each.

With L_d = L cap ker Phi_d(g) and L_0 = (+)_d L_d, the stabiliser of L in
the adelic centraliser is a subgroup of finite index of that of L_0, so

    mass(L, g) = prod_d mass(L_d, g) * prod_q |G_q . H_q| ,

H_q the glue L / L_0 at q in the discriminant group of L_0 and G_q the
image of the local stabiliser of L_0 (see `IsometryMass.h`). The masses of
the orthogonal components (d = 1, 2) come from `GenusMass`. Those of the
Hermitian components (d >= 3, lattices over Z[zeta_d]) of rank one are
computed (`CyclotomicHermitianMass.h`: the mass is 2^{1-n} |L(0,
chi_{E/K})|, halved when d is a prime power, that is h^- / w, the
L-values being products of generalized Bernoulli numbers taken exactly
through norms in Q(zeta_e)). In higher rank the mass is the standard
mass times the local factors of Hecke.jl at the bad primes, from the
Hermitian Jordan decompositions read off the trace form and the action of
zeta_d: at the primes of K unramified in E (all of them when d is not a
prime power) through the idempotents of Z_q[x]/Phi_d and the invariants of
the discriminant module under 1 - zeta^{d'}, and at the ramified prime
(d = q^a, q odd) by a Jordan splitting of the Hermitian form computed in E,
with the norm classes of the blocks of even scale. At the dyadic ramified
prime (4 | d) the factor is that of Hecke.jl: the formula for a maximal
lattice when L is maximal there, and otherwise the generic method, a chain
of overlattices of index 2 up to a maximal one with, at each step, the
numbers of sublattices and superlattices of index 2 locally isometric to
the neighbours in the chain, the local isometry being decided by the genus
symbol of Kirschmer. All the masses over Q(zeta_d) are thus computed; a
mass given in `HermMass` as a line `d num den` still takes precedence.
Checked against Hecke.jl on 78 Hermitian lattices of rank 2 and 3 over
Q(zeta_d), d = 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16.
`GENUS_IsometryMass hermitian d m` prints the standard mass of rank m over
Q(zeta_d) and, for m = 1, the mass. The glue factor is computed exhaustively: the isometries of
each primary discriminant module commuting with g are listed by
backtracking, and the orbit of the glue is closed under them. For d = 1, 2
these are the image of the local group (Nikulin, Theorem 1.9.5); for d >= 3
the code takes all of them, and says so (`image_assumed`).

When the groups are too large to be listed, the orbit length is |G| /
|Stab_G(H)| with both orders taken from graph automorphism groups: a
bijection of a module preserving the nondegenerate b is additive, so the
isometries (commuting with g) are the automorphisms of the complete graph
on the elements with the colours q(x), b(x, y) and y = g(x), and adding one
vertex per element of the glue, joined to its components, gives the
stabilizer. This works up to some thousands of vertices: an involution of
the determinant-30 genus with eigenspaces of rank 9 and 10 (glue (Z/2)^9,
modules of 512 and 1024 elements, groups of order up to 9 10^10) takes
20 s; E8+E8 with the exchange of the factors gives |O+_8(2)| = 348364800
as glue factor, as it must.

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
* The primes: the good ones (p not dividing m det) in the order of the
  number of points of the projective spaces of their eigenspaces over F_p.
  At p != 1 mod m only some eigenvalues are in F_p, which is often enough
  and much cheaper. For g = +-1 every prime not dividing det is good and
  the lattice alone gives the canonical form and the automorphism group.
  A single prime can be imposed.
* Several commuting isometries per seed (k per Gram matrix in the Isometry
  file) give the classes of (L, G) for the abelian group they generate:
  the neighbours of the common eigenvectors, the classes by the forms
  (A, A g_1^T, ..., A g_k^T).

The seeds (one per genus of pairs and per component of the neighbour
graph) are input, or come from gluing:

    GENUS_GlueIsometry [glue|enumerate] [Gram1] [Isom1] [Gram2] [Isom2] [det]
                       [OutFormat] [OutFile]

gives the classes of the g-stable primitive extensions of (L_1, g_1) (+)
(L_2, g_2), g = g_1 (+) g_2, of the given determinant or of any, from the
g-stable isotropic subgroups of D(L_1) (+) D(L_2) meeting neither summand,
taken up to the centralisers of the g_i (`IsometryGlue.h`). In the
enumerate mode each genus of pairs met is enumerated by its mass; since the
genus of a glued pair depends only on the genera of the pieces and the
local glue, this gives all the classes of pairs whose pieces lie in the
genera of the (L_i, g_i). `RunGlueIsometryTests.sh` checks it. `CI_tests/15B_GenusEnumeration/RunEnumerateIsometryTests.sh`
checks the even unimodular genus of rank 16 with -1, a root reflection and
an element of order 3, and the genus of rank 8 and determinant 25 with -1
and a root reflection.
