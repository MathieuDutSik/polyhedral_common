# Hecke operators for Sp_4(Z): the projective systems

A port of the GAP code of `MyPolyhedral/lib/ProjectiveSystem.g` and
`Symplectic.g`, for the functionality used by the driver
`HeckeOperatorSp4/TheBigComputation.g`: the Hecke operators on the homology
of the subgroups `Gamma0_A(N)` (the Siegel parabolic subgroup modulo N) and
`Gamma0_B(N)` (the Klingen parabolic subgroup) of Sp_4(Z).

## The complex

The cells are the configurations of vectors of Z^4, up to the signs, whose
cones of rank one forms meet the Siegel space: the intersection of the
Voronoi decomposition of the positive definite forms of rank 4 with the
Siegel space, of dimension 6, which has 7 levels. The top cells are the
perfect cones A_4 (10 vectors) and D_4 (12 vectors); the cells of the other
levels are the subsets of their vectors equivalent under Sp_4(Z) to one of
the listed configurations, down to a single vector. The cells whose
vectors span Z^4 are the well rounded ones, the others are at infinity.

`sp4_complex.h` builds it (`compute_sp4_cell_complex`):

* the stabilizers in Sp_4(Z) of the top configurations (automorphisms of
  the vector configuration with the integral symplectic ones kept),
* the embeddings of the configurations into the top cells (for the
  configurations at infinity, which span Lagrangian subspaces, through the
  primitive isotropic sublattices),
* the neighbors of the top cells across their ridges and the orbits of top
  cells, with generators of Sp_4(Z),
* the orbits of cells of each level as triples (top cell, subset of its
  vectors, matrix), with their stabilizers,
* the boundaries with the signs from the faces of codimension two
  (`untangle_all_signs`) and the orientation characters of the stabilizers.

The levels are indexed as in `src_perfect`: the top cells at the index 0 and
the faces of a cell of the index k at the index k + 1 (the level L of GAP
is the index 7 - L here).

## The Hecke operators

`sp4_hecke.h` computes the chain map of a double coset Sp_4(Z) x Sp_4(Z)
on the complex (`compute_sp4_hecke_chain_map`), in the format of the chain
maps of `src_perfect` (`HeckeChainMap`):

* a single vector goes to the cell of its image, found by a walk through
  the top cells across the facets of the cones of rank one forms,
* a cell of a lower index goes to a filling of the image of its boundary,
  found on a covering of top cells by the exact sparse solver of
  `src_sparse_solver`, averaged over the stabilizer of the coset and
  propagated to the orbit of the coset.

Step B is the one of `src_perfect` (`compute_hecke_homology`), through the
view `Sp4ComplexView` of the complex: the cosets of the subgroup, the
orbits of the cells of the quotient, the homology (rational or modulo a
prime) and the Hecke matrices. The Klingen subgroup is the type `Gamma0` of
`src_perfect` (the last row congruent to (0, 0, 0, *)), the Siegel one the
type `Gamma0Siegel`.

## The program

```
SP4_HeckeOperators [choice] [p] [eDimMax] [subgroups] [levels] [prime] [prefix]
```

* `choice`: `A`, `B` or `C`, the Hecke matrices of
  `SYMPL_GetHeckeOperatorCandidate_A/B/C` (`B` = diag(p, p, 1, 1) is the
  case of the driver),
* `p`: their prime,
* `eDimMax`: the levels of GAP up to which the chain map is computed (4 in
  the driver),
* `subgroups`: `Gamma0_A`, `Gamma0_B` or both, comma separated,
* `levels`: the levels N, comma separated,
* `prime`: 0 for the rational homology, otherwise the homology modulo it,
* `prefix`: the results go to `prefix_[subgroup]_[N].g` (GAP readable).

The homology is that of the well rounded cells (the cells at infinity
eliminated, `EliminateInfinity = true` in the driver).

## Differences with the GAP code

* The cells of the quotient by the subgroup are obtained from the cosets
  of the subgroup (the orbit splitting of `src_perfect`), not from a
  separate enumeration of the top cells of the subgroup.
* The Hecke matrices on the homology of the stored GAP results
  (`DATA/TheResult_*_Hecke_B2_2`) are not integral (for instance 4/3 on a
  space of dimension one), which is not possible for Hecke eigenvalues; the
  matrices here are to be compared with the integrality and with the
  dimensions of the homology instead.

## Results

The driver's range (`B`, p = 2, eDimMax = 4), computed modulo 2^31 - 1. The
whole range takes about 5 minutes and 2.4G; the complex alone takes 2.5 s.
Below are the factors of the characteristic polynomial of T(2) on the
homology of the GAP level 4 (exponents give multiplicities). For Gamma0_B(17)
and (19) the polynomial was lifted by the Chinese remainder theorem from
three primes; the rational homology does not finish in 20 minutes.

| N  | Gamma0_A (Siegel)                  | Gamma0_B (Klingen)                                          |
|----|------------------------------------|-------------------------------------------------------------|
| 2  | (x-2)                              | (x-10)(x-3)                                                 |
| 3  | (x-15)                             | (x-15)^2                                                    |
| 5  | (x-15)                             | (x-15)^2 (x-2)                                              |
| 7  | (x-15)                             | (x-15)^2 (x-5)                                              |
| 11 | (x-15)                             | (x-15)^2 (x^2-14x+46) (x^2+6x+6)^2                          |
| 13 | (x-15) (x-1)                       | (x-15)^2 (x-1)^4 (x^2-13x+38)                               |
| 17 | (x-15) (x-3) (x+5)^4               | (x-15)^2 (x-3)^4 (x+5)^2 (x^3-19x^2+96x-76) (x^3+5x^2-16x-12)^2 |
| 19 | (x-15) (x-3) (x+3)^4               | (x-15)^2 (x-3)^4 (x+2)^4 (x+3)^2 (x^3-21x^2+126x-178) (x^3-3x^2-18x+38)^2 |

The eigenvalue 15 = 1 + 2 + 4 + 8 is the degree of T(2), the Eisenstein
class.
