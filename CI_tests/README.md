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

List of CI tests:
* `CTYPE_AdjSch`: Computation of the Ctypes in dimension 5 with the adjacency scheme.
* `Canonicalization`: This is for running the canonicalization of Gram
matrices.
* `ConeIntersection`: This is checking different methods for computing the
intersection of polyhedral cones.
* `Copositivity`: This is for checking the copositivity code.
* `DirectVolumePolytope`: computation of the volume of polytopes.
* `DoubleCosets`: Computation of double cosets (unfinished).
* `EnumLatticeDelaunays`: Enumeration of Delaunay polytopes of lattices.
* `EquiStabDatabase`: Memoization of equivalence and stabilizers.
* `EquivDualDesc`: Equivariant computation of dual description.
* `FindPositiveVectors`: Finding positive vectors of indefinite forms.
* `GeometricallyUnique`: Finding interior point in polytope defined by facets,
which is invariant under transformation.
* `GroupSkeleton`: Computing the automorphism group by using the skeleton.
* `Indefinite_Stab_Equi`: Finding automorphism group and equivalence of Lorentzian
matrices.
* `IntegralAutomorphy`: Computed the integral automorphism group of a
configuration of vectors.
* `IntegralPoints`: Compute the integral points of some polytope.
* `Isotropic`: Testing existence of isotropic in quadrtaic forms and computing them
if proven to exist.
* `RedundantCtype5`: We test 3 methods for reducing by redundancy.
* `Reflective`: This is for using the edgewalk algorithm of Allcock for
building the polyheral cone.
* `SimpleDualDesc`: This is for the code for computing the dual description
of polyhedral cones.
* `Tspace_IsoDelaunay`: Computing the IsoDelaunay spaces of T-spaces.
* `Tspace_StabEquiv`: Computing the stabilizer/equivalence of elements of T-spaces.
* `WythoffH4`: This is for computing the dual description of the facets of
the orbit of `x W(H4)` for x a random vector (there are 4 orbits).

More needs to be added:
* Computing the orbits of K-dim faces of polytopes (the G6 and G7 are good examples).
* Running of the various orbit splitting algorithm.
* Computing automorphism group of polytope with LinPolytope.
* Checking for pointedness of cones.
* Short vector enumeration.
