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

## The Picard group U(2,1; Z[i])

`data/Picard_Gaussian.txt` is the ideal cell complex dual to the spine of
D. Yasaki (*Explicit reduction for SU(2,1; Z[i])*): its cells are his nine
strongly admissible sets of isotropic vectors (the 4-cells J^5, J^8 on top,
the cusps at infinity, `at_infinity single`). The stabilizers (4 times
those of Yasaki, U(2,1) containing the scalar units) and the incidences
agree with his tables, and the ranks of H^1, H^2, H^3 of Gamma_0(n) with
his Table 3 (n = 2, 3, 4, 5, 2+i, 3+2i, 4+i, 5+2i).

The cells are not cones of a linear fan (J^8 spans a space of rank 6), so
the walk across the facets does not apply. `walk cusp_reduction` uses the
reduction theory of Yasaki instead (`picard_reduction.h`): a Euclidean
algorithm on the cusps, translations of the integral Heisenberg group and
the element w, gives for a cusp v an element of the group mapping it to
e_1 and a path of edges equivalent to {e_1, e_3} from e_1 to v (the
analogue of the continued fractions of Manin). The top cells containing
these edges give the cell of a vertex and the coverings of the fillings.

The Hecke matrices (realified, row vectors) are `data/Picard_Gaussian_T2+i.txt`
= diag(5, 2+i, 1), `T2-i`, `T1+i` = diag(2, 1+i, 1) and `T3` = diag(9, 3, 1).
The levels are given as `a+bw` (w = i), for instance

```
PROJ_HeckeOperators data/Picard_Gaussian.txt data/Picard_Gaussian_T2+i.txt 3 Gamma0 3,3+2w,4+w 2147483647,2147483629 out
```

computes T(2+i) on H^2 (index 2) and H^3 (index 3) of Gamma_0(n) modulo two
primes. `identity` as the Hecke matrix gives the homology only.

## The paramodular groups and Sp_2g(Z): PROJ_SiegelVoronoiCells

```
PROJ_SiegelVoronoiCells [g] [N] [mode] [n_start] [seed] [output]
```

for the group Sp(J, Z), J = [[0, D], [-D, 0]], D = diag(1, ..., 1, N) on
Z^{2g}: Sp_2g(Z) for N = 1, the paramodular group K(N) for g = 2. Its
cells are the cells of the Siegel space of J (the positive forms A with
(A J)^2 = -lambda I) cut out by the Voronoi cones of the perfect forms in
2g variables, as MacPherson and McConnell did for Sp_4(Z)
(`siegel_voronoi.h`). The modes:

* `full` (g = 2): the description file of the complex, all the levels;
* `top`: the top cells only, with the orders of their stabilizers;
* `fast`: as `top`, the facets tested with the margin 0.01 / |X| only (the
  same cells for Sp_4(Z), K(2), K(3), about ten times faster);
* `rays`: as `top`, the ridges found by shooting rays from an interior
  point instead of testing all the facets (fast, but misses ridges: 20 or
  21 top cells instead of 22 for K(3)).

* The perfect cones of GL_2g(Z) and their adjacencies come from
  `src_perfect` (`get_full_complex_enumeration_kernel`); the first top cell
  is the cone containing a random rational point of the Siegel space
  (`ConeWalker::find_containing_cone_form`), the others are reached across
  the facets meeting the Siegel space in dimension 5, up to Sp(J_N, Z)
  (`list_embedding`).
* Whether a face C(X) meets the Siegel space, and the dimension of the
  intersection, is decided numerically (`SiegelOracle`): a point with all
  weights at least 0.01 / |X| is searched by Levenberg-Marquardt
  (`n_start` random starts), the dimension is the kernel of the
  differential of the equations on the span of C(X). Points found only
  with smaller weights are reported as doubtful and not counted.
* The top cells and the ridges crossed by the walk are certified exactly
  (`certify_top_cell`, `certify_ridge`): the Siegel space has the rational
  parametrization A(X, W) = S^{-1} [[W, W X], [X W, X W X + W^{-1}]] S^{-1},
  S = diag(I, D). A top cell gets a rational point A(X, W) in the open
  cone. A ridge gets a rational curve (X, W)(t) from one side to the
  other, on which, by Sturm sequences, the functional of the facet has one
  simple root in (0, 1) and the other facet functionals of the two cones
  none: the facet meets the Siegel space in dimension g (g + 1) - 1. Only
  the certified ridges are crossed. In the modes other than `full`, 200
  random rational points of the Siegel space are located in the perfect
  cones as an independent check that every top cell met is known.
* The faces at infinity are the faces spanning an isotropic line or
  Lagrangian plane.
* A configuration Y of the next level contained in X is a face only if the
  cell of Y is in the closure of the cell of X (`SiegelOracle::in_closure`);
  the others are written in the optional section

  ```
  non_faces K
  index iConf <matrix>       # K times: the vectors of Y, rows of l_conf[index][iConf]
  ```

  which `compute_projective_cell_complex` uses to drop these incidences.

For N = 1 the output is the complex of `data/Sp4.txt` (the fourteen sets,
no non faces, the same homology and Hecke operators). The top cells, all
certified: 2 for Sp_4(Z), 6 for K(2), 22 for K(3). For N = 2 the orbits
of cells of the dimensions 6, ..., 0 are 6, 25, 41, 33, 12, 4, 2, but the
cells of dimension 1 end at interior points of the cells at infinity of
two vectors (Lagrangian pairs), so the complex is not regular at infinity
and the orientations of the low dimensional cells are not obtained yet.

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
