The CI tests
============

The CI tests are used to detect problems. They are run once per month
(each test has its own day). This avoids straining the GitHub usage
and it allows a more leisurely correction of the problem.

GAP interfacing
---------------

The polyhedral work started historically with the GAP code.
For reasons of performance, parallelism and expressivity,
things have moved to C++.. But the GAP inheritance remains
for the CI tests and checking their correctness.


List of CI tests:
* `Canonicalization`: This is for running the canonicalization of Gram
matrices.
* `ConeIntersection`: This is checking different methods for computing the
intersection of polyhedral cones.
* `Copositivity`: This is for checking the copositivity code.
* `Reflective`: This is for using the edgewalk algorithm of Allcock for
building the polyheral cone.
* `IntegralAutomorphy`: Computed the integral automorphism group of a
configuration of vectors.
* `SimpleDualDesc`: This is for the code for computing the dual description
of polyhedral cones.
* `WythoffH4`: This is for computing the dual description of the facets of
the orbit of `x W(H4)` for x a random vector (there are 4 orbits).
* `Redundant`: We test 3 methods for reducing by redundancy.

More needs to be added:
* Enumerating integral points in polytopes (the 8-dim perfect Delaunay polytopes
are good examples).
* Computing the orbits of K-dim faces of polytopes (the G6 and G7 are good examples).
* Running of the various orbit splitting algorithm.
* Computing automorphism group of polytope with LinPolytope.
* Computing automorphism group of skeletton of a polytope.
* Checking for pointedness of cones.
* The sampling of facets (using the Various CUT polytopes)
* Short vector enumeration.
