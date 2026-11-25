Enumeration of a genus of lattices by Kneser p-neighbours, the class number
and the Smith-Minkowski-Siegel mass being checked against Hecke; and a
negative case with an unreachable mass, where the enumeration must report
itself incomplete rather than loop or claim completeness, since every result
of the program rests on that certificate.

The genera of Hermitian lattices over CM fields (GENUS_HermitianEnumerate,
TestHermitianEnumeration.g on hermitian_enumeration_cases.g): from their
genus symbols, the seed built and the genus enumerated by neighbours until
the mass, over abelian fields given by (n, H), where the mass is computed,
and over fields given by a polynomial, with the class numbers, the orders of
the unitary groups and the masses of Hecke.jl; for three CM fields that are
not abelian the program computes the mass itself (the zeta values of the
real subfield by Siegel's formula, the L-values of the relative character
numerically and certified with Arb).
