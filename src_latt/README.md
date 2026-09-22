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
| deep insertion | 1.92x | 1008 ms | 1/30 |
| BKZ-4 | 1.43x | 146 ms | 2/30 |
| BKZ-8 | **2.06x** | 340 ms | 2/30 |
| BKZ-12 | 2.01x | 658 ms | 1/30 |

So BKZ-8 beats deep insertion on quality at a third of its cost, and BKZ-12
buys nothing over BKZ-8 at twice the price: on these instances the dial has
already saturated by `beta = 8`. At dimension 8 every block size gives the same
answer as deep insertion, the root lattices being recovered by LLL already.

Entry points: **BKZreducedBasis** (`delta = 99/100`, no tour cap),
**BKZreducedBasisDelta** (explicit `delta` and tour cap, for an early abort),
and **IsBKZreduced**.
