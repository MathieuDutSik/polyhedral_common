The Hecke matrices of the algebraic modular forms over Z
(GENUS_MPI_AlgebraicModularForms), from the classes of a genus, with and
without a level structure: E8 at levels 1 and 5, II_16, and a genus of 25
classes of DefLatDB at levels 1 and 11 (1365 states). Each result is
compared with the expected one and checked in GAP: rows summing to the
number of neighbors, self-adjoint for the mass inner product, the operators
at two primes commuting; one case also on two MPI processes, and the whole
again with the SANITY_CHECK macros.

Over CM fields (TestHermitianAlgebraicModularForms.g, hermitian_cases.g):
the classes of a Hermitian genus by GENUS_HermitianEnumerate, then the
Hecke matrices by GENUS_MPI_AlgebraicModularForms with Ring = Hermitian, at
every prime above the given rational primes, with and without a level, an
ideal of O_E: a prime (inert, split, ramified, the dyadic (1+i) of Z[i]) or
composite ((1+i)^4, (2+i)(2-i), 3(1+i), 4 over Q(sqrt -7)). Over Q(sqrt -7),
Q(sqrt -11), Z[sqrt -5] and Z[i] the values are those of the former
GENUS_AlgebraicModularForms (which the levels of su21_hecke.tex were
validated with against PARI); there are also two non-Galois CM quartics and
a sextic CM field over a non-Galois cubic. The masses are computed from
certified L-values (flint). At a split prime the operator is adjoint to the
one at the conjugate prime.
