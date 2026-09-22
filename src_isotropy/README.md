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
  * **IsDeepLLLreduced** tests a Gram matrix for the property, recomputing the
    Gram-Schmidt data from scratch.

Validated by **TEST_DeepLLL** in `src_latt`. Cost, on the benchmark there: at
dimension 8 deep insertion changes almost nothing, the root lattices being
already recovered by LLL; at dimension 20 it halves the LLL potential and cuts
the number of cases failing to reach the hidden presentation from 3 in 30 to
1, at about seventy times LLL's time. That cost is the exact rational
Gram-Schmidt recomputation after each insertion, not the search, and is where
a faster version would have to start.

Reference:
  * Claus-Peter Schnorr, M. Euchner, Lattice basis reduction: improved
    practical algorithms and solving subset sum problems, Mathematical
    Programming 66 (1994) 181--199.
