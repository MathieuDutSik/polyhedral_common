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
The "GitHub CI" times are those of the last successful scheduled run:
"build and run" is the test step (most workflows compile inside it), the
whole job adds the installation of the dependencies (nauty, GAP, ...).

* `01_RatIntAutomorphy` -- `ci_01A_rat_int_automorphy`: Rational and integral automorphism groups of a configuration of vectors, the associated isomorphism and canonical form codes, and the decomposition of the rational group into double cosets of the integral one.
  GitHub CI: build 8 min, run 2 min, whole job 12 min, 2026-11-01.
* `01B_NonSymmetricStabilizers` -- `ci_01B_direct_matrix_stabilizer`: The stabilizer of a matrix given directly, which need not be symmetric (GRP_DirectMatrix_Stabilizer: Paley tournaments, circulants with distinct values, matrices with all entries distinct, the Petersen graph), and the automorphism group of a configuration of vectors for a non-symmetric Gram matrix (GRP_LinPolytope_Automorphism_GramMat, including more than 1000 vectors, where the heuristic scheme runs). The group orders are known and each generator is checked in GAP.
* `02A_FindPositiveVectors` -- `ci_02A_find_positive_vectors`: Finding positive vectors of indefinite forms.
  GitHub CI: build and run 13 min, whole job 16 min, 2026-11-02.
* `02B_RealAlgebraicPolytope` -- `ci_02B_real_algebraic_polytope`: Automorphism group and dual description, by bb, cdd, lrs and normaliz, of the regular N-gons over the real algebraic field Q(2*sin(2*pi/N)); and the dual description of G553, a 7-dimensional cone on 150 generators over Q(2*sin(2*pi/5)) with a symmetry group, where the recursive dual description using the symmetry and the direct one ignoring it are checked against each other.
  GitHub CI: build and run 9 min, whole job 15 min, 2026-10-02.
* `03A_Tspaces_IsoDelaunay` -- `ci_03A_enum_isodelaunay`: Enumeration of isoDelaunay domains.
  GitHub CI: build and run 13 min, whole job 16 min, 2026-11-03.
* `03B_KCovering` -- `ci_03B_k_covering`: Order-k Delaunay tilings and (L,k)-types: the tilings of classical lattices (A2, Z^2, Z^3, A3, D4) with their k-covering radius and volume identity, and the enumeration of the (L,k)-types of the space of all forms in dimensions 2 and 3 with the minimal k-covering density, which in the plane agrees with Blundon's thinnest k-fold lattice coverings.
  GitHub CI: build and run 8 min, whole job 10 min, 2026-10-03.
* `04A_LattAutIsomCan` -- `ci_04A_latt_aut_isom_can`: Lattice automorphism isomorphism and canonicalization, on the classical lattices and on the well rounded lattices of dimension 10 realized from their shortest vectors.
  Local run: build 3 min, run 57 min (`AllTests.g`), 2026-11-08, machine under heavy load.
  GitHub CI: build and run 44 min, whole job 47 min, 2026-11-04.
* `05A_DirectVolumePolytope` -- `ci_05A_direct_volume`: Computation of the volume of polytopes.
  GitHub CI: build and run 20 s, whole job 3 min, 2026-11-05.
* `05B_CddSkeletons` -- `ci_05B_cdd_skeletons`: CDD skeletons (DualDescriptionAdjacencies).
  GitHub CI: build and run 23 s, whole job 3 min, 2026-10-05.
* `06A_GeometricallyUnique` -- `ci_06A_geometrically_unique`: Finding interior point in polytope defined by facets, which is invariant under transformation.
  Local run: about 3.3 h for the 3382 examples (38 min per 650), 2026-10-01.
  GitHub CI: build and run 1.9 h, whole job 2.0 h, 2026-11-06.
* `06B_CharacteristicVectorSet` -- `ci_06B_characteristic_vector_set`: Generation of characteristic vector sets.
  Local run: build 10 s, run 3 min (`ProcessEntries.g`), 2026-11-08, machine under heavy load.
  GitHub CI: build and run 4 min, whole job 8 min, 2026-10-06.
* `07A_facelatticegen` -- `ci_07A_facelatticegen`: Computing the automorphism group by using the skeleton.
  GitHub CI: build and run 1 min, whole job 5 min, 2026-11-07.
* `07B_StronglySemiEutactic` -- `ci_07B_strongly_semi_eutactic`: Determination of whether a positive definite quadratic form is strongly semi-eutactic, on the 34 strongly semi-eutactic perfect forms of dimension 9 and on lattices that are not strongly semi-eutactic.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-10-07.
* `08A_ConeIntersection` -- `ci_08A_cone_int`: This is checking different methods for computing the intersection of polyhedral cones.
  GitHub CI: build and run 23 s, whole job 3 min, 2026-11-08.
* `08B_AllZeroOneSolutions` -- `ci_08B_zero_one_solutions`: Enumeration of all the 0/1 solutions of a linear system A x = b.
  GitHub CI: build and run 17 min, whole job 20 min, 2026-10-08.
* `09A_Isotropic` -- `ci_09A_isotropic`: Testing existence of isotropic in quadrtaic forms and computing them if proven to exist.
  GitHub CI: build and run 32 min, whole job 35 min, 2026-11-09.
* `09B_PerfectPartialContractingHomotopy` -- `ci_09B_partial_contracting_homotopy`: For the partial contracting homotopy.
* `11A_EquiStabDatabase` -- `ci_11A_equistabdatabase`: Memoization of equivalence and stabilizers.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-11-11.
* `11B_Robust_Covering` -- `ci_11B_robust_covering`: Robust covering density computations.
* `12A_Copositivity` -- `ci_12A_copos`: This is for checking the copositivity code.
  GitHub CI: cancelled at the 2 h timeout on 2026-09-12 and 2026-11-12, no successful run in the retained history.
* `12B_GenPolytopes` -- `ci_12B_gen_polytopes`: Generalized Polytopes.
* `13A_PeriodicDelaunay` -- `ci_13A_periodic_delaunay`: The periodic Delaunay and iso-Delaunay computations.
  GitHub CI: build and run 7 min, whole job 10 min, 2026-11-13.
* `15A_LinearProgramming` -- `ci_15A_linear_programming`: Linear programming and LinearDetermineByInequalities, over the rational types from GAP driven data files.
  GitHub CI: build and run 54 s, whole job 4 min, 2026-11-15.
* `15B_GenusEnumeration` -- `ci_15B_genus_enumeration`: Enumeration of a genus of lattices by Kneser p-neighbours, the class number and the Smith-Minkowski-Siegel mass being checked against Hecke; and a negative case with an unreachable mass, where the enumeration must report itself incomplete rather than loop or claim completeness, since every result of the program rests on that certificate.
  GitHub CI: build and run 2 min, whole job 3 min, 2026-10-15.
* `16A_EquivDualDesc` -- `ci_16A_equiv_dual_desc`: Equivariant computation of dual description.
  GitHub CI: build and run 3 min, whole job 6 min, 2026-11-16.
* `16B_GramReductions` -- `ci_16B_gram_reductions`: The reductions of positive definite Gram matrices offered by LATT_Reduction (LLL, Seysen, deep insertion, BKZ, self-dual BKZ, slide, Minkowski and the best-of search), run on scrambled root lattices; the transformation is checked to be unimodular and to produce the returned form, and the condition each method claims is checked independently in GAP.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-10-16.
* `17A_Laminations` -- `ci_17A_two_laminations`: Computation of two laminations.
  GitHub CI: build and run 13 s, whole job 3 min, 2026-11-17.
* `17B_Erdahl` -- `ci_17B_erdahl`: Enumeration of the perfect Delaunay polyhedra with src_erdahl: the 4 of the full space in dimension 7 and the 3 centrally symmetric ones (functions of center e_1/2) in dimension 8, with the CUT_8 dual description heuristics of `16A_EquivDualDesc/CUT_K8`; and the perfect Delaunay polytopes only (method `polytopes`, from the Erdahl-Rybnikov polytope and its symmetrization) in the same two cases.
  GitHub CI: `TestErdahl_centered8` build and run 19 min, whole job 25 min; `TestErdahl_full7` build and run 23 min, whole job 29 min (jobs in parallel), 2026-10-17.
* `19A_IndefiniteComp` -- `ci_19A_indefinite_comp`: Computation of indefinite forms.
  GitHub CI: build and run 1.8 h, whole job 1.8 h, 2026-11-19.
* `19B_AlgebraicModularForms` -- `ci_19B_algebraic_modular_forms`: The Hecke matrices of the algebraic modular forms over Z (GENUS_MPI_AlgebraicModularForms), from the classes of a genus, with and without a level structure: E8 at levels 1 and 5, II_16, and a genus of 25 classes of DefLatDB at levels 1 and 11 (1365 states). Each result is compared with the expected one and checked in GAP: rows summing to the number of neighbors, self-adjoint for the mass inner product, the operators at two primes commuting; one case also on two MPI processes, and the whole again with the SANITY_CHECK macros.
* `20_Reflective` -- `ci_20_reflective`: This is for using the edgewalk algorithm of Allcock for building the polyheral cone.
  GitHub CI: build and run 45 min, whole job 48 min, 2026-11-20.
* `21A_SamplingFacets` -- `ci_21A_sampling_facets`: Sampling facets of polytopes.
  GitHub CI: build and run 2 min, whole job 5 min, 2026-11-21.
* `21B_GroupIndividualization` -- `ci_21B_group_individualization`: Automorphism group of point sets with a large block and a small group, through the orbit-stabilizer computation by individualization: the WythoffH4 polytope and orbits of B6 above the thresholds, and orbits of B4 and B5 with the lowered thresholds of `TEST_GROUP_THRESHOLD` and the sanity checks.
  GitHub CI: build and run 1 min, whole job 4 min, 2026-10-21.
* `22_CommonG_IsoDelaunay` -- `ci_22_commong_isodelaunay`: CommonGramMat iso-Delaunay (rigid lattices and stars).
  GitHub CI: build and run 20 min, whole job 27 min, 2026-11-22.
* `23A_IntegralPoints` -- `ci_23A_integral_points`: Compute the integral points of some polytope.
  GitHub CI: build and run 2 min, whole job 4 min, 2026-09-23.
* `23B_SimpleDualDesc` -- `ci_23B_simple_dual_desc`: This is for the code for computing the dual description of polyhedral cones.
* `24_GroupSimplification` -- `ci_24_group_simplification`: Group Simplification.
  GitHub CI: build and run 6 min, whole job 10 min, 2026-10-24.
* `25_Redundancy` -- `ci_25_redundancy`: Redundancy tests.
* `26A_CTYPE_AdjSch` -- `ci_26A_enum_ctype`: Computation of the Ctypes in dimension 5 with the adjacency scheme.
  GitHub CI: build and run 32 s, whole job 3 min, 2026-09-26.
* `26B_Matrix_group_finiteness` -- `ci_26B_finite_groups`: Test of group finiteness.
  GitHub CI: build and run 2 min, whole job 5 min, 2026-10-26.
* `27A_OrbitShortestVector` -- `ci_27A_orbit_shortest`: Orbit of shortest vectors.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-09-27.
* `27B_CoveringMaxdet` -- `ci_27B_covering_maxdet`: Covering density optimization by determinant maximization over the iso-Delaunay domains: the classical optima A_3^*, A_4^* and A_5^* over the lattice domains of dimensions 3, 4 and 5, the domains of the periodic point set Z^3 + {0, (1/3,1/3,1/3)}, and the random-walk record search on that same point set.
  GitHub CI: build and run 13 min, whole job 17 min, 2026-10-27.
* `28A_WythoffH4` -- `ci_28A_wythoff_H4`: This is for computing the dual description of the facets of the orbit of `x W(H4)` for x a random vector (there are 4 orbits).
  GitHub CI: build and run 7 min, whole job 9 min, 2026-09-28.
* `28B_LorentzianPerfStabEqui` -- `ci_28B_perfect_lorentzian`: Enumeration of perfect lorentzian cones, and stabilizer/equivalence of the resulting forms.
  GitHub CI: build and run 5 min, whole job 8 min, 2026-10-28.
* `29_SerialDelaunayQueries` -- `ci_29_serial_delaunay_queries`: Serial Delaunay computations.
  Local run: build 3 min, run 4 min (`TestSerialDelaunayQueries.g`) and 10 s (`TestPackingCoveringRecord.g`), 2026-11-08, machine under heavy load.
  GitHub CI: build and run 14 min, whole job 17 min, 2026-10-29.
* `30_PerfectFormTspace` -- `ci_30_perfect_form_tspace`: Perfect Forms in T-spaces.
  GitHub CI: build and run 7 min, whole job 9 min, 2026-10-30.

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
