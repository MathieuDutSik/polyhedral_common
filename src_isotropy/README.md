Indefinite form computations
============================

The functions in this subdirectory allow computing with indefinite forms.

The following can be done:
  * Compute the signature of a quadratic form.
  * Computing the Indefinite LLL reduction of a form.
  * Compute the classic LLL reduction of a form.
  * Compute the Seysen reduction of a form.
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
