The Hecke matrices of the algebraic modular forms over Z
(GENUS_MPI_AlgebraicModularForms), from the classes of a genus, with and
without a level structure: E8 at levels 1 and 5, II_16, and a genus of 25
classes of DefLatDB at levels 1 and 11 (1365 states). Each result is
compared with the expected one and checked in GAP: rows summing to the
number of neighbors, self-adjoint for the mass inner product, the operators
at two primes commuting; one case also on two MPI processes, and the whole
again with the SANITY_CHECK macros. Over CM fields (TestHermitianAlgebraicModularForms.g, hermitian_cases.g):
the classes of a Hermitian genus by GENUS_HermitianEnumerate, then the
Hecke matrices by GENUS_MPI_AlgebraicModularForms with Ring = Hermitian, at
every prime above the given rational primes, with and without a level at an
inert or a split prime, over Q(sqrt -7), Q(sqrt -11), Z[sqrt -5] (the values
of GENUS_AlgebraicModularForms), two non-Galois CM quartics and a sextic CM
field over a non-Galois cubic; the masses are computed from certified
L-values (flint). At a split prime the operator is adjoint to the one at
the conjugate prime. Over the ring of integers of Q(sqrt -7),
Q(sqrt -11) and Q(sqrt -5) (GENUS_AlgebraicModularForms,
TestOrderAlgebraicModularForms.g): the neighbors at the ramified prime alone
close the genus of I_m at the mass formula, and the Hecke matrices at
ramified, split and inert primes commute, those at the primes equal to their
conjugate being self-adjoint.
