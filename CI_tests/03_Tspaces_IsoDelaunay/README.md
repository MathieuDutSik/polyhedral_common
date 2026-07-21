This test code enumerates the iso-Delaunay / L-type domains via
`LATT_SerialLattice_IsoDelaunayDomain`. It is done in `TestEnumIsoDelaunayDomain.g`
and has two phases:

Phase 1: iso-Delaunay domains of stored T-spaces (Bravais / Coxeter).
The result had been published before, see page 160 of
Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, A generalization
of Voronoi's reduction theory and its application, preprint at arxiv:0601084,
Duke Mathematical Journal 142 (2008) 127--164.

Phase 2: L-type domains for T-spaces from algebraic theory (formerly the
separate CI entry `07B_L_domains_algebraic_rings`, merged here because it
relies on the same binary).
The test with `comm_choice:="Use_realimag"` corresponds to using the group `GL_n(Z[x])`.
If one uses `comm_choice:="Trivial"` then the group used is `GL_n(Z[x])` and the conjugation
operation.
Those cases are interesting but the number of cones generated is much higher than for
the stored T-spaces of phase 1.
