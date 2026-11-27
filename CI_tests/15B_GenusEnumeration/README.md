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
numerically and certified with Arb). The unimodular genera of
Schiemann's tables over imaginary quadratic fields (A. Schiemann, J. Symbolic
Comput. 26 (1998)), for 18 entries (field, rank) from Q(sqrt -3) to
Q(sqrt -163): the class numbers of the genera of an entry sum to his. Over
the fields where 2 ramifies the enumeration uses the neighbors at the dyadic
ramified prime first, and the test is run again with the SANITY_CHECK macros,
which check that every neighbor has the local symbols of its lattice.

The classes of indefinite genera of rank >= 3 (TestIndefiniteRepresentatives.g,
indefinite_cases.g): GENUS_IndefiniteRepresentatives gives one lattice per
improper spinor genus, on genera with 1, 2, 4 and 8 classes in ranks 3, 4,
5 and 10, again with the SANITY_CHECK macros. The numbers of classes are
those of SageMath's spinor kernel; the representatives of the genera with
4 and 8 classes were checked with Hecke.jl (in the genus, pairwise
non-isometric).
