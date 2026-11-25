Three drivers, each result being checked independently in GAP.

`FindPositiveVector.g`: LATT_FindPositiveVector on the indefinite forms of
LGramReflect, for the norms -1, 0 and 1, strict and non strict.

`TestReductions.g`: the programs of src_isotropy, on the indefinite forms of
LGramReflect (one in eight) and the root lattices of ListPositiveDefinite,
all scrambled by random unimodular matrices.
- LATT_IndefiniteLLL: B unimodular with Mred = B M B^T, and either an
  isotropic vector or the conditions of Simon's LLL with c = 7/8.
- LATT_IndefiniteReduction: B unimodular with Mred = B M B^T.
- LATT_LLLreduceBasis: the same and the LLL conditions with delta = 99/100.
- LATT_DetMinimization: Mred integral, P M P^T a multiple of Mred, and the
  determinant not increased.
- LATT_AutomorphyReduction: the same for M and for the primitive integral
  multiple of its inverse. The equivariant reduction can increase the
  determinant (for example -49 to -343), which is only printed.
- CheckPositiveSemiDefinite: definite, semidefinite and indefinite forms.
- LATT_SublatticeBasisReduction: the same lattice and no larger L1 norm.
- LATT_DiophantApprox: the guarantee of Lenstra, Lenstra, Lovasz
  (Proposition 1.39), q >= 1 and |q alpha_i - p_i| <= epsilon. The program
  works in double precision, which limits epsilon from below.

`TestShort.g`: the programs of src_short on the root lattices.
- SHORT_GetShortestVector: the minimum 2 and the shortest vectors of GAP.
- SHORT_AutomorphismGroup: the generators preserve the configuration of
  roots and generate a group of the order given in
  RootLatticeAutomorphismOrders.
- SHORT_ReduceVectorFamilyGAP: transformed copies of four configurations,
  rows permuted and signs changed, reduced to four classes.
- SHORT_SplitVectorFamily: the configurations of index i mod N.
