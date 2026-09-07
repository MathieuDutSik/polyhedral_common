# Gram matrices whose canonical form is very expensive

These are even positive definite Gram matrices of rank 14 and determinant 243
-- ordinary members of a genus one wants to enumerate, not contrived inputs --
on which `ComputeCanonicalForm` becomes unusable:

* `LATT_Canonicalize gmp slow_canonic_1.txt` does not terminate in 110 s;
* driven from `GENUS_Enumerate` the same computation aborted inside nauty with
  `malloc>E gtools: Cannot allocate memory`.

## The cause

`ExtractInvariantVectorFamilyZbasis` stops when the vectors GENERATE the
lattice over Z, not when they merely have full rank. On
`slow_canonic_1.txt` those two conditions are separated by one shell:

    norm <= 2 :    72 vectors, rank  9
    norm <= 4 :  1179 vectors, rank 14, index [L : <vectors>] = 2
    norm <= 6 : 12351 vectors, rank 14, index = 1

At norm 4 the family is already of full rank but generates a sublattice of
index 2, so the extraction takes the whole next shell. The family jumps from
1179 to 12351, and after the antipodal duplication that is 24702 -- the value
reported as `LSEC: nbRow=24702 n=14`, and the size of the graph handed to the
canonical labelling.

A single index-2 obstruction therefore multiplies the family by 10.5, which is
the difference between a graph that is canonicalised easily and one that
exhausts memory. `ExtractInvariantVectorFamilyFullRank` stops at 1179 (2358
after duplication) and is fine, which is why `ArithmeticEquivalence` remains
usable on these lattices while `ComputeCanonicalForm` does not: the equivalence
test is driven from the full-rank family.

Note that the full-rank family cannot simply be substituted for the canonical
form: a family generating a proper sublattice yields a canonical basis of that
sublattice rather than of the lattice. Completing it needs a canonically
determined choice, which is what makes this a design question rather than a
one-line change.
