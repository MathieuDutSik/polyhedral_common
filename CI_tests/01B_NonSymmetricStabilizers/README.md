This test checks the automorphism groups computed for non-symmetric
matrices, which go through the non-symmetric weight matrix (a graph on
2 n + 1 vertices) and, above 1000 vectors, through the heuristic scheme with
its transposed weights.

GRP_DirectMatrix_Stabilizer computes the permutations s with
M(s(i), s(j)) = M(i, j) of a square matrix M given directly:
- Paley tournaments for p = 7, 11, 19, 23: M(i, j) is the Legendre symbol of
  j - i. The group is {x -> a x + b, a a nonzero square}, of order
  p (p - 1) / 2. The symmetrized matrix of p = 7, the complete graph, is
  checked to give 7!.
- Circulants (j - i) mod n with n distinct values: the n translations.
- Matrices with all n^2 entries distinct, n = 16, 20, 22: the trivial group,
  with more distinct weights than an 8-bit weight index holds.
- The Petersen graph: order 120 (symmetric control).

GRP_LinPolytope_Automorphism_GramMat with a non-symmetric Gram matrix:
- The 1088 integer points of [-16,16]^2 minus the origin with the form
  [[1,1],[-1,1]]: the rotations by pi/2, order 4 (order 8 for the identity
  form, as control).
- PlusMinus_1100.ext, 550 random antipodal pairs with the same form: -Id
  only, order 2.
- Generic_20.ext, 20 generic vectors of dimension 4 with a non-symmetric
  form having 316 distinct values: the trivial group.

The group orders are known in advance, and each generator returned is
checked in GAP to preserve the matrix.
