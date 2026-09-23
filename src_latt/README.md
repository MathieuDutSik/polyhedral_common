Lattice based computation
=========================

The functions in this subdirectory allow computing with lattices. But a more restrictive description is that
we do the following in this directory:
  * Computing equivalence, stabilizer of positive definite quadratic forms.
  * Computing set of shortest vectors.
  * Computing in T-spaces of quadratic forms.
  * Computing Delaunay polytopes and Iso-Delaunay domains.

General references of the lattice:
  * A. Schürmann, Computational geometry of positive definite quadratic forms. Polyhedral reduction theories, algorithms and applications. American Mathematical Society, xvi+162 pp.

LLL algorithm
-------------

The LLL algorithm allow to reduce a quadratic form into one with smaller
coefficients.

The relevant program is:
  * **LATT_lll** This is for computing the LLL of a positive definite quadratic form.

Reference:
  * [LLL basis reduction algorithm](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm)


Closest Vector Problems
-----------------------

The Closest vector problem and the Shortest Vector Problems are fundamental
problems in Geometry of Numbers:
  * The Closest vector problem ask given a positive definite quadratic form Q a vector c in R^n to find the vectors v in Z^n that minimize Q[v - c].
  * The Shortest vector problem ask given a positive definite quadratic form Q to find the non-zero vectors v in Z^n minimizing Q.

The relevant programs are:
  * **LATT_near** It allows to do those operations using LLL.

The algorithm used in this work is the one by
  * Fincke, U. and Pohst, M., Improved methods for calculating vectors of short length in a lattice, including a complexity analysis, Mathematics of Computation, 44, 1985, 463--471

Stabilizer/Equivalence of positive definite quadratic form
----------------------------------------------------------

There is a set of function for testing stabilizer/equivalence for positive definite quadratic forms.

Programs:
  * **LATT_GenerateCharacteristicVectorSet** Computing the set of characteristic vectors of a positive definite form.
  * **LATT_Automorphism** The automorphism group of a positive definite form
  * **LATT_Isomorphism** The equivalence of two positive definite forms.

Reference:
  * W. Plesken, B. Souvignier, Computing isometries of lattices, J. Symbolic Computation 24 (1997) 327--334

Canonical form of positive definite forms
-----------------------------------------

For a positive definite quadratic form, we want to find the canonical form of a form as this is very
useful for enumeration purposes.

Programs:
  * **LATT_Canonicalize** Canonicalize a positive definite form
  * **LATT_CanonicalizeMultiple** Canonicalize several positive forms with the first one being positive definite.
  * **LATT_CanonicalizeSymplectic** Finding the canonical form with the equivalence being a symplectic matrix.

Reference:
  * Mathieu Dutour Sikirić, Anna Haensch, John Voight, Wessel Van Woerden, A canonical form for positive definite matrices, [preprint at arxiv:2004.14022](https://arxiv.org/abs/2004.14022), Proceedings of the Fourteenth Algorithmic Number Theory Symposium (ANTS-XIV), edited by Steven Galbraith, Open Book Series 4, Mathematical Sciences Publishers, Berkeley, 2020.

T-space of quadratic forms
--------------------------

A T-space is a vector space of quadratic forms which contains one positive definite form.

Programs:
  * **TSPACE_Equivalence** Computing equivalence within the T-space.
  * **TSPACE_Stabilizer** Computing stabilizer within the T-space
  * **TSPACE_FileFormatConversion** Converting format of T-spaces

Reference:
  * Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, A generalization of Voronoi's reduction theory and its application, [preprint at arxiv:math/0601084](https://arxiv.org/abs/math/0601084), Duke Mathematical Journal 142 (2008) 127--164
  * Herbert Gangl, Paul Gunnells, Jonathan Hanke, Achill Schürmann, Mathieu Dutour Sikirić, Dan Yasaki, On the cohomology of linear groups over imaginary quadratic fields, [preprint at arxiv:1307.1165](https://arxiv.orb/abs/1307.1165), Journal of Pure and Applied Algebra 220-7 (2016) 2564--2589
  * Achill Schuermann, Enumerating perfect forms, [preprint at arXiv:0901.1587](https://arxiv.org/abs/0901.1587)


Lattice Delaunay
----------------

The enumeration of Delaunay polytopes is an important geometric problem for various stuff.

Program:
  * **LATT_MPI_ComputeDelaunay** for computing the Delaunay polytopes in MPI framework.

References for computing the Delaunay of a lattice:
  * Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin, Complexity and algorithms for computing Voronoi cells of lattices, [preprint at arxiv:0804.0036](https://arxiv.org/abs/0804.0036), Mathematics of Computation 78 (2009) 1713--1731


Iso-Delaunay domains
--------------------

For a T-space, we want to enumerate the iso-Delaunay domains inside of it.

Program:
* **LATT_MPI_Lattice_IsoDelaunayDomain** for running the enumeration by using MPI parallelism.

References:
  * Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, A generalization of Voronoi's reduction theory and its application, [preprint at arxiv:math/0601084](https://arxiv.org/abs/math/0601084), Duke Mathematical Journal 142 (2008) 127--164


Lattice reduction benchmark
---------------------------

A harness for comparing basis-reduction algorithms on the hidden-good-basis
experiment: a lattice is presented by a Gram matrix known to be good, a random
unimodular matrix destroys that presentation, and the algorithm sees only the
destroyed one. It is kept separate from the reducers themselves so that a
disappointing experiment can be attributed to the representation, the
objective, the move set or the search strategy rather than to all four at once.

Programs:
  * **TEST_ReductionBenchmark** `[dim] [n_iter] [seed]`, comparing the
    available reducers over the families Zn, An, Dn, E8 and a low-symmetry
    random family, at three strengths of destruction.
  * **TEST_SlideReduction** `[dim] [n_iter] [seed]`, validating the slide
    reduction of `SlideReduction.h`. It rebuilds every block from the output
    and re-tests both families of conditions, the dual half -- which goes
    through the reversed dual and back -- being the part worth testing hardest.
  * **TEST_BKZ** `[dim] [n_iter] [seed]`, validating the BKZ reduction of
    `BKZ.h`. Its decisive check is that every `b_j^*` really is a shortest
    vector of its block, the projected blocks being rebuilt from the output and
    re-enumerated by a routine sharing no state with the descent. It also
    checks that a BKZ-beta reduced basis passes the test at every smaller block
    size, which the definition implies and an implementation can break.
  * **TEST_DeepLLL** `[dim] [n_iter] [seed]`, validating the Schnorr-Euchner
    deep insertion of `src_isotropy/DeepLLL.h`. Its decisive check is that the
    output really satisfies the deep condition, the Gram-Schmidt data being
    recomputed from scratch by a routine sharing no state with the descent.
  * **TEST_SeysenReduction** `[dim] [n_iter] [seed]`, validating the Seysen
    reduction of `src_isotropy/SeysenReduction.h` on the same instances. Its
    decisive check is that the output is a local minimum of Seysen's measure,
    established by trying every transvection with a small coefficient rather
    than by trusting the closed form for the optimal one.

The two questions the harness keeps apart are *recovery*, whether the original
presentation is found back up to signed permutation, and *reduction*, whether
the output is at least as good as the hidden presentation. Only the second is
the lattice-reduction problem; recovery is strictly stronger and is not
geometrically meaningful, since every signed permutation of a good basis is
equally good.

Everything is measured on the Gram side in exact arithmetic. The orthogonal
factor of the Iwasawa decomposition never appears, being killed by
`G = B B^T`, and the remaining Iwasawa data is read off the Bareiss
decomposition that `FullGramInfo` in `Shvec_exact.h` already computes, through
`a_i^2 = d(i)/d(i-1)` and `mu_{i,j} = Nmat(i,j)/d(i)`.


BKZ
---

`BKZ.h` implements block Korkine-Zolotarev reduction. With `pi_j` the
projection orthogonal to the span of the first `j-1` vectors and
`k = min(j+beta-1, n)`, the condition is that `b_j^*` be a **shortest vector**
of the lattice generated by `pi_j(b_j), ..., pi_j(b_k)`. `beta = 2` is the
Lovasz condition, `beta = n` is Hermite-Korkine-Zolotarev, and in between
`beta` is one dial between a cheap local condition and an expensive global one.

What changes from LLL, deep insertion and Seysen is the **oracle**: the
question at each index is no longer a comparison of quantities already at hand
but a shortest-vector problem solved exactly in dimension `beta`. That routine
is `T_ShortestVector` of `Shvec_exact.h`, which is why BKZ lives here rather
than in `src_isotropy`: the dependency runs `src_latt -> src_isotropy` and not
the other way.

The block is extracted **incrementally and over the integers**. The Gram matrix
of the projected block is a Schur complement, rational in general, but scaled
by the leading minor `d_j` it is exactly what `j` steps of the fraction-free
Bareiss elimination leave in the trailing block. Advancing `j` by one is one
more elimination step, so a whole tour extracts all its blocks in `O(n^3)`
rather than recomputing each. A shortest vector is unchanged by a positive
rescaling, so the enumerator runs on the scaled matrix directly, and the test
against `delta |b_j^*|^2` is a comparison of two integers.

Termination is the deep-insertion argument, not LLL's: an insertion at `j`
leaves `D_1..D_{j-1}` fixed and strictly decreases `D_j` while saying nothing
about the minors between `j` and `k`, so `(D_1, ..., D_{n-1})` decreases
lexicographically. The intervening LLL passes never increase any `D_i`.

Measured on the benchmark at dimension 20, geometric mean of the LLL potential
relative to LLL, with total time over the 30 cases:

| method | potential vs LLL | time | failures |
|---|---|---|---|
| LLL | 1.00x | 60 ms | 3/30 |
| deep insertion | 1.92x | 108 ms | 1/30 |
| BKZ-4 | 1.43x | 147 ms | 2/30 |
| BKZ-8 | **2.06x** | 337 ms | 2/30 |
| BKZ-12 | 2.01x | 657 ms | 1/30 |

BKZ-12 buys nothing over BKZ-8 at twice the price: on these instances the dial
has already saturated by `beta = 8`. At dimension 8 every block size gives the
same answer as deep insertion, the root lattices being recovered by LLL
already.

Note that deep insertion is the better cost-quality point here, at 1.92x for
108 ms against BKZ-8's 2.06x for 337 ms. That is a recent inversion: before the
incremental Gram-Schmidt update deep insertion took 1008 ms and BKZ dominated
it outright. Use BKZ-8 when the extra quality is worth three times the time,
and raise `beta` only on evidence.

Entry points: **BKZreducedBasis** (`delta = 99/100`, no tour cap),
**BKZreducedBasisDelta** (explicit `delta` and tour cap, for an early abort),
and **IsBKZreduced**.


Slide reduction
---------------

`SlideReduction.h` implements the Gama-Nguyen slide reduction. BKZ enforces one
family of conditions on overlapping blocks and pays for it with a number of
tours that has no known polynomial bound. Slide reduction enforces **two**
families on blocks that do **not** overlap, and buys back the polynomial bound.
With `n = p k`, and indices 0-based:

* **primal**, for `i = 0..p-1`: the projected block on `[ik, ik+k-1]` is
  HKZ-reduced, which for a block of rank `k` is BKZ at block size `k`;
* **dual**, for `i = 0..p-2`: the projected block on `[ik+1, (i+1)k]` -- the
  same tiling shifted right by one -- has its **last** Gram-Schmidt norm
  `|b*_{(i+1)k}|` as large as possible.

The primal blocks tile the basis; the dual blocks are the same tiling shifted
by one, so each straddles exactly one primal boundary. That is the whole
design. A primal step rearranges vectors inside one block and so changes no
`D_{ik}`. A dual step on `[ik+1, (i+1)k]` raises `|b*_{(i+1)k}|` at fixed block
determinant, so it lowers the earlier norms in that block and strictly
decreases `D_{(i+1)k}`, touching no other boundary. Hence

```
Pi = prod_{i=1}^{p-1} D_{ik}
```

is a positive integer untouched by primal steps and strictly decreased by every
dual step: **the number of dual steps is polynomially bounded**, as in the LLL
analysis, where for BKZ no bound on the number of tours is known. The quality
bound is also slightly better than BKZ's at the same block size, the exponent
being `(n-k)/(2(k-1))` against `(n-1)/(2(k-1))`.

Two implementation points. The reduction interleaved after each step is a
**size** reduction (`IntegralSizeReduce` in `DeepLLL.h`) and not a full LLL:
size reduction changes no Gram-Schmidt norm, so it cannot disturb a block
condition just established, whereas LLL is free to swap across a boundary and
undo the work. And the dual step goes through the **reversed dual**: if the
block has Gram matrix `G` and `J` is the reversal, the reversed dual basis has
Gram matrix `J G^-1 J`, its first Gram-Schmidt norm is the reciprocal of the
block's last, and a transformation `U` of the reversed dual corresponds to
`J U^-T J` of the block. So maximising the last norm is finding a shortest
vector of `J adj(G) J` and putting it first -- the same enumeration the primal
step uses, and integral since the adjugate clears the denominators.

Entry points: **SlideReducedBasis** (explicit `k`, which must divide `n`),
**SlideReducedBasisAuto** (largest admissible `k` not exceeding a bound, via
**SlideBlockSize**), **SlideReducedBasisDelta**, and **IsSlideReduced**.

### Divide out the block content

A trap worth recording, because it cost a factor of more than ten thousand
before it was found. The projected block carries the factor `d_j`, the leading
minor of order `j`, which grows like a minor of the whole Gram matrix. Every
use of the block is invariant under a positive rescaling, so the factor looks
harmless. It is not: the dual step takes the adjugate of the block, raising the
factor to the power `k-1`, and the enumerator reduces internally through the
dual, raising it again. A `d_j` of twenty digits becomes entries of hundreds of
digits and the enumeration crawls. On one dimension-20 instance slide reduction
at `k = 5` had not finished after 300 seconds; dividing the content out of the
block brought it to 28 ms, with the Gram matrix entries themselves never
exceeding four bits throughout.

`BKZ_ProjectedBlockGram` therefore divides the content out, and BKZ's own loop
does the same, carrying the removed factor into the comparison so that the two
sides stay on one scale.

### Measured, dimension 20, 30 cases

| method | potential vs LLL | geo defect^2 | time | failures |
|---|---|---|---|---|
| LLL | 1.00x | 6650 | 60 ms | 3/30 |
| deep insertion | 1.92x | 6284 | 106 ms | 1/30 |
| BKZ-4 | 1.43x | 6227 | 146 ms | 2/30 |
| BKZ-8 | 2.06x | 6365 | 337 ms | 2/30 |
| BKZ-12 | 2.01x | 6320 | 655 ms | 1/30 |
| slide, k=4 | 1.69x | 6384 | 172 ms | 2/30 |
| slide, k=5 | 1.81x | 6320 | 212 ms | 1/30 |
| Seysen+LLL | 1.02x | **6170** | 192 ms | 1/30 |

Slide reduction at `k = 5` sits between BKZ-4 and BKZ-8 on the potential and
below both on time per unit of quality, and matches BKZ-12's failure count at a
third of its cost. Its real argument is not in this table, though: it is the
polynomial bound on the number of dual steps, which BKZ has no analogue of and
which is what one wants when the dimension grows past where a table like this
can be produced.

Reduction of a vector family
----------------------------

`VectFamilyReduction.h`, and the program **VectFamily_Reduction**, reduce a
family of vectors: a change of coordinates in the ambient space making the
coefficients of the family small. This is what is applied to an `EXT` matrix
before a dual description, and it is a different problem from reducing a
lattice basis.

The difference decides which reducer to use. What the consumer of the output
pays for is the size of every coefficient, and for the dual description the
relevant quantity is `sqr_estimate_facet_coefficients` of `norms.h`, the
Hadamard bound on the facet coefficients that will be produced. That depends
on all the vectors symmetrically, so a reduction aiming at one short vector,
which is what LLL does, is optimising the wrong thing here. Seysen's measure
is a much closer match, and deep insertion sometimes wins instead. Which of
them wins is not predictable from the input:

| instance | size | `direct` | best single | winner |
|---|---|---|---|---|
| ContactE8 | 240 x 9 | 1.46e8 | 1.00e8 | seysen |
| Perfect E7 | 63 x 28 | 2.43e24 | 1.41e24 | deep |
| ER35 | 35 x 8 | 6.0e4 | 6.0e4, L1 112 -> 105 | seysen_lll |
| CUT_7 | 64 x 22 | 9.07e22 | 9.07e22 | direct |
| CUT_K8 | 128 x 29 | 1.29e34 | 1.29e34, lower L1 | direct |
| 24cell | 24 x 5 | 81 | 81 | all tie |

So the useful method is **`best`**: run all the candidates and keep whichever
actually minimises the estimate, with the unreduced input among the candidates
so that the result is never worse than what was handed in. It costs three to
eighteen times a single reduction, which is tens to hundreds of milliseconds
on the instances above and negligible against the dual description that
follows.

Note that `dual` is consistently the worst of them on this measure, by one to
five orders of magnitude, which is worth knowing since it is one of the two
methods the program originally offered.

All thirteen single methods are available: `direct`, `dual`, `seysen`,
`seysen_best`, `seysen_lll`, `deep`, `deep5`, `deep10`, `bkz4`, `bkz8`,
`bkz12`, `slide4`, `slide8`, plus `best`. Running all thirteen costs 57 ms on
ContactE8 (240 x 9), 102 ms on CUT_7 (64 x 22) and 492 ms on CUT_K8
(128 x 29), which is nothing against the dual description that follows.

Adding BKZ and slide reduction did not change any winner on the corpus above:
`deep` still wins Perfect E7, at 1.41e24 against BKZ-12's 2.02e24 and BKZ-4's
2.06e24, and the block methods tie everywhere else. They are in the candidate
list because the winner varies by instance and there is no cost to having
more candidates, not because they have been observed to win here.

The underlying dispatch is `ReduceVectorFamilyGeneral`, returning the reduced
family, the change of coordinates, the method that won and its measured
quality. `ReduceVectorFamilyKernel` of `ClassicLLL.h` takes the Gram-matrix
reducer as a functor, so a new reducer is added by naming it in
`ReduceVectorFamilySingle` and nowhere else.
