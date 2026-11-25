# CommonGramMat iso-Delaunay test (rigid lattices and their stars)

This test exercises the `CommonGramMat` support of the iso-Delaunay domain code
(`src_delaunay/IsoDelaunayDomains.h`). Given a fixed positive definite Gram
matrix `G`, the enumeration is restricted to the **star of `G`**: the
iso-Delaunay domains whose closure contains `G`, enumerated up to the subgroup
of arithmetic equivalences that fix `G`.

For each dimension the test performs the full pipeline:

1. **Enumerate** the top-dimensional iso-Delaunay domains of the classic
   `GL_d(Z)` T-space (`LATT_SerialLattice_IsoDelaunayDomain`), dumping every
   domain as a boost archive via the `PrefixIsoDelaunayDomains` option.
2. **Extract** the rigid lattices: for each domain, `LATT_AnalysisIsoDelaunay`
   with the `FileFullRankRays` QUERIES option emits the full-rank (rank `d`,
   positive definite) extreme rays. A lattice on such a ray is *rigid* — its
   iso-Delaunay domain is zero-dimensional.
3. **Deduplicate** the rays up to `GL_d(Z)`-equivalence with `LATT_Canonicalize`
   (grouping identical canonical Gram matrices). In dimension 5 this yields the
   7 rigid lattices.
4. **Star of each rigid lattice**: rerun the enumeration with `CommonGramMat`
   set to the (canonical, reduced) rigid form and count the domains of its star.
   A rigid lattice contained in an iso-Delaunay domain in two inequivalent ways
   shows up as a star count greater than 1.

The test checks the number of iso-Delaunay domains, the number of rigid
lattices, and the multiset of star counts against the reference values below.

The 5-dimensional iso-Delaunay classification (222 domains) is from:
Achill Schürmann, Mathieu Dutour Sikirić, Frank Vallentin, *A generalization of
Voronoi's reduction theory and its application*, Duke Mathematical Journal 142
(2008) 127--164, arXiv:math/0601084.

Note: this is a long-running test (the 5-dimensional enumeration alone takes
several minutes), appropriate for the monthly CI schedule.
