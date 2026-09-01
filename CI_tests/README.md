The CI tests
============

The CI tests are used to detect problems. They are run once per month
(each test has its own day). This avoids straining the GitHub usage
and it allows a more leisurely correction of the problem.

GAP interfacing
---------------

The polyhedral work started historically with the GAP code.
For reasons of performance, parallelism and expressivity,
things have moved to C++. But the GAP inheritance remains
for the CI tests and checking their correctness.

List of CI tests. Each directory is driven by the workflow named next
to it; `.github/workflows/ci_NN...` fires on day NN of the month.

* `01_RatIntAutomorphy` -- `ci_01_rat_int_automorphy`: Rational and integral automorphism groups of a configuration of vectors, the associated isomorphism and canonical form codes, and the decomposition of the rational group into double cosets of the integral one.
* `02A_FindPositiveVectors` -- `ci_02A_find_positive_vectors`: Finding positive vectors of indefinite forms.
* `02B_RealAlgebraicPolytope` -- `ci_02B_real_algebraic_polytope`: Automorphism group and dual description, by bb, cdd, lrs and normaliz, of the regular N-gons over the real algebraic field Q(sin(2*pi/N)).
* `03_Tspaces_IsoDelaunay` -- `ci_03_enum_isodelaunay`: Enumeration of isoDelaunay domains.
* `04A_LattAutIsomCan` -- `ci_04A_latt_aut_isom_can`: Lattice automorphism isomorphism and canonicalization.
* `05A_DirectVolumePolytope` -- `ci_05A_direct_volume`: Computation of the volume of polytopes.
* `05B_CddSkeletons` -- `ci_05B_cdd_skeletons`: CDD skeletons (DualDescriptionAdjacencies).
* `06A_GeometricallyUnique` -- `ci_06A_geometrically_unique`: Finding interior point in polytope defined by facets, which is invariant under transformation.
* `06B_CharacteristicVectorSet` -- `ci_06B_characteristic_vector_set`: Generation of characteristic vector sets.
* `07A_facelatticegen` -- `ci_07A_facelatticegen`: Computing the automorphism group by using the skeleton.
* `08_ConeIntersection` -- `ci_08A_cone_int`: This is checking different methods for computing the intersection of polyhedral cones.
* `09A_Isotropic` -- `ci_09A_isotropic`: Testing existence of isotropic in quadrtaic forms and computing them if proven to exist.
* `09B_PerfectContractingHomotopy` -- `ci_09B_contracting_homotopy`: For contracting homotopy computations.
* `11A_EquiStabDatabase` -- `ci_11A_equistabdatabase`: Memoization of equivalence and stabilizers.
* `11B_Robust_Covering` -- `ci_11B_robust_covering`: Robust covering density computations.
* `12A_Copositivity` -- `ci_12A_copos`: This is for checking the copositivity code.
* `12B_GenPolytopes` -- `ci_12B_gen_polytopes`: Generalized Polytopes.
* `13A_PeriodicDelaunay` -- `ci_13A_periodic_delaunay`: The periodic Delaunay and iso-Delaunay computations.
* `15_LinearProgramming` -- `ci_15_linear_programming`: Linear programming and LinearDetermineByInequalities.
* `16_EquivDualDesc` -- `ci_16_equiv_dual_desc`: Equivariant computation of dual description.
* `17_Laminations` -- `ci_17_two_laminations`: Computation of two laminations.
* `19_IndefiniteComp` -- `ci_19_indefinite_comp`: Computation of indefinite forms.
* `20_Reflective` -- `ci_20_reflective`: This is for using the edgewalk algorithm of Allcock for building the polyheral cone.
* `21A_SamplingFacets` -- `ci_21A_sampling_facets`: Sampling facets of polytopes.
* `21B_ShortRealizability` -- `ci_21B_shortest_realizability`: Sampling facets of polytopes.
* `22_CommonG_IsoDelaunay` -- `ci_22_commong_isodelaunay`: CommonGramMat iso-Delaunay (rigid lattices and stars).
* `23A_IntegralPoints` -- `ci_23A_integral_points`: Compute the integral points of some polytope.
* `23B_SimpleDualDesc` -- `ci_23B_simple_dual_desc`: This is for the code for computing the dual description of polyhedral cones.
* `24_GroupSimplification` -- `ci_24_group_simplification`: Group Simplification.
* `25_Redundancy` -- `ci_25_redundancy`: Redundancy tests.
* `26A_CTYPE_AdjSch` -- `ci_26A_enum_ctype`: Computation of the Ctypes in dimension 5 with the adjacency scheme.
* `26B_Matrix_group_finiteness` -- `ci_26B_finite_groups`: Test of group finiteness.
* `27A_OrbitShortestVector` -- `ci_27A_orbit_shortest`: Orbit of shortest vectors.
* `28A_WythoffH4` -- `ci_28A_wythoff_H4`: This is for computing the dual description of the facets of the orbit of `x W(H4)` for x a random vector (there are 4 orbits).
* `28B_LorentzianPerfStabEqui` -- `ci_28B_perfect_lorentzian`: Enumeration of perfect lorentzian cones, and stabilizer/equivalence of the resulting forms.
* `29_SerialDelaunayQueries` -- `ci_29_serial_delaunay_queries`: Serial Delaunay computations.
* `30_PerfectFormTspace` -- `ci_30_perfect_form_tspace`: Perfect Forms in T-spaces.

Directories under CI_tests that no workflow drives. They hold data or
past work rather than a runnable test, and are kept on purpose:
(only directories tracked by git are listed; an empty directory left in a
working tree is not part of the repository)

* `24B_BrandhortLattices`
* `DoubleCosets`
* `Rankin`

`DATA` holds inputs shared by several of the sections above.

More needs to be added:
* Computing the orbits of K-dim faces of polytopes (the G6 and G7 are good examples).
* Running of the various orbit splitting algorithm.
* Computing automorphism group of polytope with LinPolytope.
* Checking for pointedness of cones.
* Short vector enumeration.
