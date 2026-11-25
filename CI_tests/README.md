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

* `01A_RatIntAutomorphy` -- `ci_01A_rat_int_automorphy`: Rational and integral automorphisms, isomorphisms and canonical forms of vector sets.
  GitHub CI: build 8 min, run 2 min, whole job 12 min, 2026-11-01.
* `01B_NonSymmetricStabilizers` -- `ci_01B_direct_matrix_stabilizer`: Stabilizers of non-symmetric matrices and Gram matrices.
  Local run: build 7 min, run 1 min (`TestNonSymmetricStabilizers.g`), 2026-11-24.
* `02A_FindPositiveVectors` -- `ci_02A_find_positive_vectors`: Positive vectors, indefinite LLL, reductions and shortest vector configurations.
  GitHub CI: build and run 13 min, whole job 16 min, 2026-11-02.
* `02B_RealAlgebraicPolytope` -- `ci_02B_real_algebraic_polytope`: Symmetries and dual descriptions over real algebraic fields.
  GitHub CI: build and run 9 min, whole job 15 min, 2026-10-02.
* `03A_Tspaces_IsoDelaunay` -- `ci_03A_enum_isodelaunay`: Enumeration of isoDelaunay domains.
  GitHub CI: build and run 13 min, whole job 16 min, 2026-11-03.
* `03B_KCovering` -- `ci_03B_k_covering`: Order-k Delaunay tilings and (L,k)-types.
  GitHub CI: build and run 8 min, whole job 10 min, 2026-10-03.
* `04A_LattAutIsomCan` -- `ci_04A_latt_aut_isom_can`: Lattice automorphisms, isomorphisms and canonical forms.
  Local run: build 3 min, run 57 min (`AllTests.g`), 2026-11-08, machine under heavy load.
  GitHub CI: build and run 44 min, whole job 47 min, 2026-11-04.
* `05A_DirectVolumePolytope` -- `ci_05A_direct_volume`: Computation of the volume of polytopes.
  GitHub CI: build and run 20 s, whole job 3 min, 2026-11-05.
* `05B_DoubleDescriptionSkeletons` -- `ci_05B_cdd_skeletons`: CDD skeletons (DualDescriptionAdjacencies).
  GitHub CI: build and run 23 s, whole job 3 min, 2026-10-05.
* `06A_GeometricallyUnique` -- `ci_06A_geometrically_unique`: Interior point of a polytope invariant under its transformations.
  Local run: build 10 s, run 2 min (`ProcessEntries.g`), 2026-11-25.
* `06B_CharacteristicVectorSet` -- `ci_06B_characteristic_vector_set`: Generation of characteristic vector sets.
  Local run: build 10 s, run 3 min (`ProcessEntries.g`), 2026-11-08, machine under heavy load.
  GitHub CI: build and run 4 min, whole job 8 min, 2026-10-06.
* `07A_facelatticegen` -- `ci_07A_facelatticegen`: Computing the automorphism group by using the skeleton.
  GitHub CI: build and run 1 min, whole job 5 min, 2026-11-07.
* `07B_StronglySemiEutactic` -- `ci_07B_strongly_semi_eutactic`: Whether a positive definite form is strongly semi-eutactic.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-10-07.
* `08A_ConeIntersection` -- `ci_08A_cone_int`: This is checking different methods for computing the intersection of polyhedral cones.
  GitHub CI: build and run 23 s, whole job 3 min, 2026-11-08.
* `08B_AllZeroOneSolutions` -- `ci_08B_zero_one_solutions`: Enumeration of all the 0/1 solutions of a linear system A x = b.
  GitHub CI: build and run 17 min, whole job 20 min, 2026-10-08.
* `09A_Isotropic` -- `ci_09A_isotropic`: Testing existence of isotropic in quadrtaic forms and computing them if proven to exist.
  GitHub CI: build and run 32 min, whole job 35 min, 2026-11-09.
* `09B_PerfectPartialContractingHomotopy` -- `ci_09B_partial_contracting_homotopy`: For the partial contracting homotopy.
  Local run: build 1 min, run 3 min (`TestPartialContractingHomotopy.g` and `TestHeckeOperators.g`), 2026-11-24.
* `11A_EquiStabDatabase` -- `ci_11A_equistabdatabase`: Memoization of equivalence and stabilizers.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-11-11.
* `11B_Robust_Covering` -- `ci_11B_robust_covering`: Robust covering density computations.
  Local run: fails after 10 s (`ComputeRobustCoveringDensity.g`), `Robust_ExactRobustCoveringDensity` stopping on the connectivity check of an orbit (`|lconn|=1 orbit_len=3`), 2026-11-24.
* `12A_Copositivity` -- `ci_12A_copos`: This is for checking the copositivity code.
  GitHub CI: cancelled at the 2 h timeout on 2026-09-12 and 2026-11-12, no successful run in the retained history.
  Local run: build 1 min, then `TestCompletePositivity.g` passes and `TestCopos.g` does 186 of its 330 cases in about 6 min together, but does not finish case 187 (`spn/7/cop7_14`), stopped after 2.75 h, 2026-11-24.
* `12B_GenPolytopes` -- `ci_12B_gen_polytopes`: Generalized Polytopes.
  Local run: build 1 min, run under 10 s (`TestGenPolytope.g`), 2026-11-24.
* `13A_PeriodicDelaunay` -- `ci_13A_periodic_delaunay`: The periodic Delaunay and iso-Delaunay computations.
  GitHub CI: build and run 7 min, whole job 10 min, 2026-11-13.
* `15A_LinearProgramming` -- `ci_15A_linear_programming`: Linear programming and LinearDetermineByInequalities.
  GitHub CI: build and run 54 s, whole job 4 min, 2026-11-15.
* `15B_GenusEnumeration` -- `ci_15B_genus_enumeration`: Genus enumeration by Kneser neighbours, checked by the mass.
  GitHub CI: build and run 2 min, whole job 3 min, 2026-10-15.
* `16A_EquivDualDesc` -- `ci_16A_equiv_dual_desc`: Equivariant computation of dual description.
  GitHub CI: build and run 3 min, whole job 6 min, 2026-11-16.
* `16B_GramReductions` -- `ci_16B_gram_reductions`: The reductions of positive definite Gram matrices of LATT_Reduction.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-10-16.
* `17A_Laminations` -- `ci_17A_two_laminations`: Computation of two laminations.
  GitHub CI: build and run 13 s, whole job 3 min, 2026-11-17.
* `17B_Erdahl` -- `ci_17B_erdahl`: Enumeration of the perfect Delaunay polyhedra in the Erdahl cone.
  GitHub CI: `TestErdahl_centered8` build and run 19 min, whole job 25 min; `TestErdahl_full7` build and run 23 min, whole job 29 min (jobs in parallel), 2026-10-17.
* `19A_IndefiniteComp` -- `ci_19A_indefinite_comp`: Stabilizers, equivalences and orbits of indefinite forms, on flint.
  GitHub CI: build and run 1.8 h, whole job 1.8 h, 2026-11-19.
* `19B_AlgebraicModularForms` -- `ci_19B_algebraic_modular_forms`: Hecke matrices of algebraic modular forms.
* `20_Reflective` -- `ci_20_reflective`: Allcock's edgewalk: reflective Lorentzian lattices and their isometries, on flint.
  GitHub CI: build and run 45 min, whole job 48 min, 2026-11-20.
* `21A_SamplingFacets` -- `ci_21A_sampling_facets`: Sampling facets of polytopes.
  GitHub CI: build and run 2 min, whole job 5 min, 2026-11-21.
* `21B_GroupIndividualization` -- `ci_21B_group_individualization`: Automorphism groups of point sets by individualization.
  GitHub CI: build and run 1 min, whole job 4 min, 2026-10-21.
* `22A_CommonG_IsoDelaunay` -- `ci_22_commong_isodelaunay`: CommonGramMat iso-Delaunay (rigid lattices and stars).
  GitHub CI: build and run 20 min, whole job 27 min, 2026-11-22.
* `23A_IntegralPoints` -- `ci_23A_integral_points`: Compute the integral points of some polytope.
  GitHub CI: build and run 2 min, whole job 4 min, 2026-09-23.
* `23B_SimpleDualDesc` -- `ci_23B_simple_dual_desc`: This is for the code for computing the dual description of polyhedral cones.
  Local run: build shared with `12B`, run 30 min (`TestSimpleDualDesc.g`), 2026-11-24.
* `24A_GroupSimplification` -- `ci_24_group_simplification`: Group Simplification.
  GitHub CI: build and run 6 min, whole job 10 min, 2026-10-24.
* `25A_Redundancy` -- `ci_25A_redundancy`: Redundancy tests.
  Local run: build 20 s, run 33 min (`TestCoherence.g`), peak memory 4.2 GB, 2026-11-24.
* `26A_CTYPE_AdjSch` -- `ci_26A_enum_ctype`: Computation of the Ctypes in dimension 5 with the adjacency scheme.
  GitHub CI: build and run 32 s, whole job 3 min, 2026-09-26.
* `26B_Matrix_group_finiteness` -- `ci_26B_finite_groups`: Test of group finiteness.
  GitHub CI: build and run 2 min, whole job 5 min, 2026-10-26.
* `27A_OrbitShortestVector` -- `ci_27A_orbit_shortest`: Orbit of shortest vectors.
  GitHub CI: build and run 1 min, whole job 3 min, 2026-09-27.
* `27B_CoveringMaxdet` -- `ci_27B_covering_maxdet`: Covering density optimization over iso-Delaunay domains.
  GitHub CI: build and run 13 min, whole job 17 min, 2026-10-27.
* `28A_WythoffH4` -- `ci_28A_wythoff_H4`: Dual description of the facets of the 4 orbits of x W(H4), x random.
  GitHub CI: build and run 7 min, whole job 9 min, 2026-09-28.
* `28B_LorentzianPerfStabEqui` -- `ci_28B_perfect_lorentzian`: Perfect Lorentzian cones and stabilizer/equivalence of the forms.
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
