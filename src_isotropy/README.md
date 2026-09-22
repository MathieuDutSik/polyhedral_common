Indefinite form computations
============================

The functions in this subdirectory allow computing with indefinite forms.

The following can be done:
  * Compute the signature of a quadratic form.
  * Computing the Indefinite LLL reduction of a form.
  * Compute the classic LLL reduction of a form.
  * Compute the Seysen reduction of a form.
  * Compute the Schnorr-Euchner deep insertion reduction of a form.
  * Find isotropic vector of an indefinite form.
  * Test whether a quadratic form has an isotropic vector or not.
  * Find a vector of positive norm of a form

This is the basis of subsequent developments for indefinite forms
and lattice computations.

Algorithms
----------

A central method used is the Indefinite LLL for decomposing complicated forms.


Seysen reduction
----------------

`SeysenReduction.h` implements Seysen's reduction, which minimises

```
S(A) = sum_i A_ii (A^-1)_ii
```

over the bases of the lattice, `A` being the Gram matrix. The measure is at
least the dimension, with equality exactly for an orthogonal basis; it is
scale free; it depends on the basis only through `A`, so the orthogonal factor
never enters; and it satisfies `S(A) = S(A^-1)`, which is the point of the
method: a basis and its reciprocal are reduced by one and the same descent,
rather than by a primal pass and a dual pass.

The move is the transvection `b_i <- b_i + lambda b_j`, whose effect on the
reciprocal basis is `b^_j <- b^_j - lambda b^_i`, and the optimal integer
`lambda` is the rounded quotient

```
lambda = round( (A_jj A*_ij - A_ij A*_ii) / (2 A_jj A*_ii) ),   A* = A^-1.
```

Unlike the size reduction of LLL, this coefficient involves the geometry of
the reciprocal lattice, and the step is taken only if it lowers the measure.

Everything runs over the integers: the rational inverse is replaced by the
adjugate, `det(A)` cancels from `lambda`, and the descent decreases the
integer `Shat(A) = sum_i A_ii adj(A)_ii`, which is bounded below by
`n det(A)` with `det(A)` fixed. That gives a termination proof of the same
shape as the LLL one.

Entry points, all returning the same `LLLreduction` pair as the LLL reducers
and therefore interchangeable with them:
  * **SeysenReducedBasis** applies every improving pair as it is met. This is
    the one to use.
  * **SeysenReducedBasisBest** applies only the most profitable pair of each
    sweep. Slower, and in our measurements no better, but its trajectory does
    not depend on the order in which pairs are visited, which is what one
    wants when comparing implementations.
  * **SeysenLLLreducedBasis** alternates Seysen and LLL while the measure
    improves. Seysen keeps the basis balanced, LLL produces the short leading
    vector; they are complements.
  * **SeysenMeasure** and **SeysenIntegralPotential** expose the measure
    itself.

Validated by **TEST_SeysenReduction** in `src_latt`, which checks among other
things that the output really is a local minimum, by trying every transvection
with a small coefficient. That is a test of the closed form for `lambda`
rather than a restatement of it.

Reference:
  * Michael Seysen, Simultaneous reduction of a lattice basis and its
    reciprocal basis, Combinatorica 13 (1993) 363--376.


Deep insertion
--------------

`DeepLLL.h` implements the Schnorr-Euchner deep insertion. Where LLL enforces
the Lovasz condition at the single index `k-1`, a deep reduced basis satisfies

```
delta |b_i^*|^2 <= |pi_i(b_k)|^2    for every i < k,
```

`pi_i` being the projection orthogonal to the span of the first `i-1` vectors.
When it fails at some `i`, `b_k` is taken out and inserted at position `i`,
the vectors in between shifting up by one. The case `i = k-1` is the LLL swap,
so the condition is a strict strengthening of LLL's and a deep reduced basis
is in particular LLL reduced.

The test is free: running `i` upwards and maintaining

```
c_1 = |b_k|^2,   c_{i+1} = c_i - mu_{i,k}^2 |b_i^*|^2
```

gives `c_i = |pi_i(b_k)|^2` at each step, on data the size reduction has
already produced.

The arithmetic is **integral throughout**. Keeping the `mu` as rationals costs
dearly: numerators and denominators grow with the minors, and the data is
recomputed after every insertion. Instead we keep de Weger's integral form,
`d_i = D_i` and `lambda_{i,j} = d_j mu_{j,i}`. Three things become integer
operations: the Gram-Schmidt data itself, by a Bareiss-type recursion with
exact divisions; the projected norms, through `S_i = d_{i-1} |pi_i(b_k)|^2`
with `S_{i+1} = (S_i d_i - lambda_{k,i}^2) / d_{i-1}`, exact because `S_i` is
the Gram determinant of `(b_1, ..., b_{i-1}, b_k)`; and the deep test itself,
which with `delta = num/den` reads

```
|pi_i(b_k)|^2 < delta |b_i^*|^2    <=>    den * S_i < num * d_i,
```

a comparison of two integers. This is worth 4.6x: at dimension 20 the
benchmark run goes from 4.67 s to 1.01 s, with output bit-identical to the
rational implementation on all 75 cases of the benchmark.

Termination is **not** the LLL potential argument, and the difference is worth
knowing. Writing `D_j` for the leading principal minor of order `j`, an
insertion at `i` from `k` leaves `D_1..D_{i-1}` and `D_k..D_n` unchanged and
strictly decreases `D_i`, but says nothing about `D_j` for `i < j < k`, which
may rise. The product `prod_j D_j` therefore need not decrease, and only for
an adjacent swap, where those middle indices do not exist, does the classical
argument apply. What survives is that `(D_1, ..., D_{n-1})` decreases strictly
in the lexicographic order, which is well founded over the positive integers.
The bound it gives is exponential, and no polynomial bound on the number of
deep insertions is known: this is why the insertion depth is usually
restricted in practice.

Entry points, returning the same `LLLreduction` pair as the other reducers:
  * **DeepLLLreducedBasis** with unrestricted insertion depth.
  * **DeepLLLreducedBasisDepth** with a depth parameter `d`, trying only the
    positions `i < d` and `i >= k - d`. The tail positions include `i = k-1`,
    so the output is LLL reduced for every `d >= 1`.
  * **DeepLLLreducedBasisDepthDelta** additionally takes `delta` as a pair of
    integers; the others use 99/100, the value this package uses for LLL.
  * **IsDeepLLLreduced** tests a Gram matrix for the property, recomputing the
    integral Gram-Schmidt data from scratch.

Validated by **TEST_DeepLLL** in `src_latt`. Cost, on the benchmark there: at
dimension 8 deep insertion changes almost nothing, the root lattices being
already recovered by LLL; at dimension 20 it halves the LLL potential and cuts
the number of cases failing to reach the hidden presentation from 3 in 30 to
1, at about sixteen times LLL's time, down from seventy before the arithmetic
was made fraction free. What remains is the recomputation of the integral
Gram-Schmidt data after each insertion, which is where a further speedup would
have to start: it is recomputed from the insertion point down, and an
incremental update of the rows that actually move would avoid most of it.

Reference:
  * Claus-Peter Schnorr, M. Euchner, Lattice basis reduction: improved
    practical algorithms and solving subset sum problems, Mathematical
    Programming 66 (1994) 181--199.

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
| ER35 | 35 x 8 | 6.0e4 | 6.0e4, lower L1 | seysen_lll |
| CUT_7 | 64 x 22 | 9.07e22 | 9.07e22 | direct |
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

The underlying dispatch is `ReduceVectorFamilyGeneral`, returning the reduced
family, the change of coordinates, the method that won and its measured
quality. `ReduceVectorFamilyKernel` of `ClassicLLL.h` takes the Gram-matrix
reducer as a functor, so a new reducer is added by naming it in
`ReduceVectorFamilySingle` and nowhere else.
