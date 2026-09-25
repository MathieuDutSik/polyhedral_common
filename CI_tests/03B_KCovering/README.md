This test exercises the order-k Delaunay tilings and the (L,k)-types of
`src_k_coverings`, in `TestKCovering.g`:

Phase 1 (Tiling): `LATT_SerialComputeKDelaunay` on classical lattices
(A2, Z^2, Z^3, A3, D4) for several k. The number of orbits of tiles, the
squared k-covering radius (OutFormat = GAP_Covering) and the volume identity
of the tiling (QUERIES/FileTilingVolume) are compared with stored values.
The values were checked by hand for A2 and Z^2 (kagome tiling for A2 with
k = 2).

Phase 2 (Enumeration): `LATT_SerialLattice_IsoKDelaunayDomain` on the space
of all forms (TypeTspace = Classic) in dimensions 2 and 3. The number of
(L,k)-types and the minimum over the types of the normalized k-covering
density (DATA/FileCoveringOptimum, the determinant maximization of
CoveringMaxdet.h) are compared with stored values. In the plane the optima
agree with Blundon's thinnest k-fold lattice coverings: 4 pi / sqrt(27) for
k = 2 and 25 pi / 18 (the square lattice) for k = 4.
