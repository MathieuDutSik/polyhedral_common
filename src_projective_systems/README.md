# Hecke operators for projective systems

A port of the GAP code of `MyPolyhedral/lib/ProjectiveSystem.g` and
`Symplectic.g` (the functionality used by the driver
`HeckeOperatorSp4/TheBigComputation.g`), made independent of the group:
the group and its cell complex are given by a description file, of which
`data/Sp4.txt` is the one of Sp_4(Z).

## The complex

The group G is the group of the matrices of GL_m(R) preserving a form h
on R^m (g h g^* = h, acting on the row vectors on the right), for R = Z
or R = O_K = Z[omega], omega^2 = s omega - p, K imaginary quadratic. The
cells are configurations of vectors of R^m taken up to the units of R:
the description file lists representative configurations of each level,
the top cells first, and the cells of a level are the subsets of the
vectors of the top cells that are equivalent under G to one of the
representatives of the level, down to the single vectors.

`projective_complex.h` builds the complex (`compute_projective_cell_complex`):

* the stabilizers in G of the top configurations (automorphisms of the
  unit multiples of the vectors, the elements of G kept),
* the embeddings of the configurations into the top cells: an isometry
  (for the form) between a basis of the span of the configuration and
  vectors of the top cell; for the configurations that do not span R^m,
  an isometry of the saturated sublattices of their spans, assumed to
  extend to G (true for the primitive isotropic sublattices of Sp_4(Z)),
* the neighbors of the top cells across their ridges and the orbits of top
  cells, with generators of G,
* the orbits of cells of each level as triples (top cell, subset of its
  vectors, matrix), with their stabilizers,
* the boundaries with the signs from the faces of codimension two
  (`untangle_all_signs`) and the orientation characters of the stabilizers.

The levels are indexed as in `src_perfect`: the top cells at the index 0 and
the faces of a cell of the index k at the index k + 1 (the level L of GAP
is the index n_levels - L here).

Over O_K the vectors and matrices are realified as in `src_perfect`
(`FiniteIndexSubgroup`): x_j = a_j + b_j omega gives (a_1, ..., a_m, b_1,
..., b_m), and the form is replaced by the trace form Tr_{K/Q} h(u, v).

## The description file

Tokens separated by blanks, `#` starts a comment:

```
ring Z                      # or: ring QuadraticImag s p
dimension m
form <matrix>               # over O_K: two matrices, h = h_a + h_b omega
at_infinity rank            # rank: the cells not spanning R^m are at infinity;
                            # single: only the single vectors
walk voronoi                # how the Hecke operators locate a vector
levels L
level k <matrix> ... <matrix>   # k configurations, top cells first
...
```

A matrix is its numbers of rows and columns followed by its entries; the
configurations have one vector per row (realified over O_K).

## The Hecke operators

`projective_hecke.h` computes the chain map of a double coset G x G on the
complex (`compute_projective_hecke_chain_map`), in the format of the chain
maps of `src_perfect` (`HeckeChainMap`):

* a single vector goes to the cell of the primitive vector of its image,
  found by a walk through the top cells across the facets of the cones of
  the rank one forms (`walk voronoi`),
* a cell of a lower index goes to a filling of the image of its boundary,
  found on a covering of top cells by the exact sparse solver of
  `src_sparse_solver`, averaged over the stabilizer of the coset and
  propagated to the orbit of the coset.

Step B is the one of `src_perfect` (`compute_hecke_homology`), through the
view `ProjectiveComplexView` of the complex: the cosets of the subgroup,
the orbits of the cells of the quotient, the homology (rational or modulo a
prime) and the Hecke matrices. For Sp_4(Z) the Klingen subgroup is the type
`Gamma0` of `src_perfect` (the last row congruent to (0, 0, 0, *)), the
Siegel one the type `Gamma0Siegel`.

## The program

```
PROJ_HeckeOperators [description] [hecke] [eDimMax] [subgroups] [levels] [prime] [prefix]
```

* `description`: the description file (`data/Sp4.txt`),
* `hecke`: the Hecke matrix, realified, as a matrix of rationals
  (`data/Sp4_T2.txt` = diag(2, 2, 1, 1), `data/Sp4_T3.txt`),
* `eDimMax`: the number of levels, from the single vectors, on which the
  chain map is computed (4 in the Sp_4 driver),
* `subgroups`: subgroup types of `src_perfect`, comma separated (for
  Sp_4(Z): `Gamma0Siegel`, `Gamma0`),
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

Sp_4(Z), the driver's range (`data/Sp4.txt`, `data/Sp4_T2.txt`,
eDimMax = 4), computed modulo 2^31 - 1. The prime levels up to 73, with
the lifts identified, are in `Paper_hecke_operators/sp4_hecke.tex`.
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
