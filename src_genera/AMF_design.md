# Algebraic modular forms on top of src_genera: design note

Status (2026-10-28): steps a, b, c, e and f done for the split and inert
primes and the trivial weight; see "State" at the end.

## Goal

Spaces of algebraic modular forms (Gross) for the group G = Aut(L, h) of a
definite lattice, with level structure and Hecke operators, in one code that
covers:

- orthogonal groups over Z (the present `src_genera`, d = 1),
- unitary groups over an imaginary quadratic order O = Z[w], w^2 = s w - p
  (the check of the Picard computations: definite U(3) over Z[i]),
- later quaternionic Hermitian lattices (same mechanism, O non commutative).

The first use cases are the two independent checks of the Hecke papers:

1. U(3) over Z[i]: reproduce the PARI prototype (`amf_gauss.gp`,
   `amf_inert.gp`) at the eight levels of `su21_hecke.tex`.
2. O(5) over Z: the genus of quinary lattices of discriminant 61 gives the
   paramodular forms of level 61 (Ladd, Rama-Tornaria); recover f61 with
   lambda_2 = -7, lambda_3 = -3 of `sp4_hecke.tex`.

## 1. One representation for all rings

An O-lattice of O-rank m is stored as Z^(d m) (d = rank_Z O), in the realified
layout of `src_projective_systems` (a_1..a_m, b_1..b_m), with

- the Gram matrix B = Tr_{K/Q}(h) of the trace form (positive definite),
- the matrix W of the multiplication by w (Z-linear, d m x d m).

For d = 1, W is absent and everything reduces to the current code.

The ring part of `ProjectiveRing` (d, m, s, p, W, units, the trace form
`get_real_form`, the residue maps) moves to one shared header, used by both
`src_projective_systems` and `src_genera`, rather than being copied.

## 2. Isometries, automorphisms, canonical forms: the family {B, B W}

A Z-linear g preserving B and B W (the bilinear form (x, y) -> B(x, w y))
commutes with W: B(gx, W g y) = B(x, W y) = B(gx, g W y) and B is
non-degenerate. So the O-linear isometries of h are exactly the integral
automorphisms of the family {B, B W}, and isometry classes are the classes of
that family. This uses the existing machinery unchanged:

- `ArithmeticAutomorphismGroupMultiple`, `ArithmeticEquivalenceMultiple`:
  the computation on the whole family goes through `FCT_ListMat_Vdiag`, which
  detects a non-symmetric family and builds a directed weight matrix; the
  absolute trick on the antipodal pairs did not (fixed, 4f6bec41).
- The invariant vector family of B is O-stable (the units preserve B), so the
  shells / V_cv families are the right ones.

**A bug to fix first.** `ComputeCanonicalFormMultiple` reaches
`CanonicallyReorder_SHV_kernel`, which hard-codes `is_symm = true`
(LatticeStabEquiCan.h:153): only the upper triangle of B W is read, so for a
non-symmetric family the result depends on the input order of the vectors and
is not a canonical form. `ComputeCanonicalFormSymplectic` (LATT_Canonicalize
Symplectic) is affected in the same way. The fix: dispatch on
`is_family_symmmetric` to `WeightMatrix<false, ...>`, and check the antipodal
(AbsTrick) path for the same assumption. Test: canonical forms of random
GL_n(Z) conjugates of a symplectic and of a Hermitian family must coincide.

The canonical form is what makes the Hecke operators cheap: the canonical
basis of a neighbour L' and that stored for its class give the isometry
L' -> L_j directly (g = B_j^{-1} B'), with no equivalence test. That is
needed for the level structure (section 4).

## 3. Neighbours over O

At a prime P of O over p, coprime to the level and to disc(L):

- **Split** p = P Pbar (unitary) or p in Z (orthogonal): residue field F_p.
  For a line x of L/PL, isotropic for the induced form, adjusted so that
  p | h(x, x): M = {y in L : h(y, x) in Pbar}, L' = M + O x / varpi. This is
  the construction of `amf_gauss.gp`; for d = 1 it is `GenusNeighbor`.
- **Inert** p: residue field F_{p^2}, x isotropic in L/pL for the Hermitian
  form over F_{p^2}, h(x, x) adjusted to 0 mod p^2, L' = M + O x / p
  (`amf_inert.gp`).
- Ramified primes: later, not needed for the first checks.

All of it is Z-linear algebra: M is the kernel of a Z-linear map
L -> O/Pbar, the generators of O x / varpi are x W^k / varpi, then an HNF.

Needed new pieces: the residue map O -> O/P (F_p or F_{p^2}) and a small
F_{p^2} arithmetic (there is none in basic_common_cpp); the projective space
P^(m-1)(F_q) with the index bijection of `ProjectiveLineOfIndex`, generalised
from F_p to F_q, so that `ProjectiveLineOrbitRepresentatives` (no
materialisation, one bit per point) works for both; the isotropic-line
enumeration.

The p = 2 valuation correction of `GenusNeighbor` stays as it is; the Z
path keeps its own specialised code where it matters for speed.

## 4. Level structure and Hecke matrices

Level K(n) for n coprime to the Hecke primes. The space is the functions on
the pairs (L_i, orbit of Aut(L_i) on X(L_i)), X(L) = lines of L/nL (isotropic
ones for inert n), as in the PARI prototype; trivial weight first, the weight
representations V_lambda later.

T(P) acts by: for each class L_i, each isotropic line x of L_i/PL_i (all of
them, or the orbits of Aut(L_i) weighted by their lengths when the level is
trivial), the neighbour L', its class j and isometry g: L' -> L_j from the
canonical form; a level datum l of L_i is transported to L'/nL' = L_i/nL_i
(n coprime to P) and then by g, giving the column of (L_i, l).

The matrix is checked by: commutation of T(P) for two primes, T(P) and
T(Pbar) adjoint for the mass inner product (weights 1/|Stab|), the row sums
equal to the number of neighbours, and the Ramanujan bound on the
non-Eisenstein eigenvalues.

## 5. Enumeration of the genus

`GenusEnumeration` is generic in the neighbour construction, the canonical
form and the automorphism data; it gets a policy (Z-lattice / O-lattice)
rather than a second copy. Termination:

- over Z: the mass, as now;
- over O: no Hermitian mass formula in the tree yet. The neighbour graph at a
  split prime connects the special genus; the class number of the genus is
  then cross-checked with Hecke.jl (installed: `~/.julia/packages/Hecke`),
  which enumerates Hermitian genera, and the Hermitian mass formula
  (Gan-Hanke-Yu) is a later addition giving the same certification as over
  Z.

## 6. Program

`GENUS_AMF` (or `AMF_HeckeOperators`): input the description file of
`src_projective_systems` (ring, form), the level n, the primes; output the
Hecke matrices, their characteristic polynomials factored, and the
commutation checks, in the format of the `PROJ_HeckeOperators` outputs so that
the lifts and identifications of the papers apply unchanged.

## 7. Work plan

| Step | Content | Estimate |
|------|---------|----------|
| a | Non-symmetric canonical form in LatticeStabEquiCan + test (also fixes the symplectic one) | 0.5 day |
| b | Shared ring header extracted from projective_complex.h | 0.5 day |
| c | O-neighbours (split, inert), F_{p^2}, P^(m-1)(F_q) orbits | 1 day |
| d | GenusEnumeration with the lattice policy | 0.5-1 day |
| e | Hecke matrices with level | 1 day |
| f | Validation: U(3)/Z[i] vs PARI at 8 levels; quinary disc 61 vs f61; Hecke.jl class numbers over Q(sqrt -3), Q(sqrt -7) | 0.5-1 day |

Runs follow the usual 10G / 30 min caps through `scripts/run_capped.sh`.

## State

- `RealifiedRing.h`: the ring, the form, the description head and the
  `a+bw` parser, shared with `src_projective_systems`.
- `OrderLatticeNeighbor.h`: the primes of O (split of degree one, inert;
  over Z), the lines of L / P L through an adapted basis, the neighbours.
- `AlgebraicModularForms.h`, `GENUS_AlgebraicModularForms`: the classes in
  canonical form, the level lines mod N and their orbits, the Hecke
  matrices; the genus is the closure under the Hecke primes, the mass
  sum 1/|Aut| is printed.
- Validated: the 8 levels of su21_hecke.tex against the PARI prototype
  (identical characteristic polynomials, T(2+i) and T(3+2i) commuting at
  7+2i); the quinary genus of half discriminant 61 (8 classes, mass 31/96
  equal to GENUS_LocalDensity genus_mass, f61 with lambda = -7, -3 at 2, 3);
  the class numbers, automorphism orders and masses of O_K^m, m = 2, 3, 4,
  for D = -4, -3, -7, -8, -11, equal to those of Hecke.jl.
- Not yet: the orbits of Aut(L) on the Hecke lines (each neighbour is
  canonicalized), the levels over Z, the ramified primes, the weights, the
  Hermitian mass formula (the closure is certified only by comparison).
