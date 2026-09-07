The dual description of polyhedral cones
========================================

The problem
-----------

A polyhedral cone in `R^d` admits two descriptions. It is the set of
nonnegative combinations of finitely many generators

    C = { sum_i lambda_i v_i : lambda_i >= 0 }

and it is the set of solutions of finitely many homogeneous linear
inequalities

    C = { x : <f_j, x> >= 0 for all j }

That the two notions coincide is the Minkowski-Weyl theorem. Passing
from one description to the other is the *dual description* problem,
also called representation conversion. The generators that cannot be
removed are the extreme rays, the inequalities that cannot be removed
are the facets, and the two families are exchanged when one replaces
the cone by its dual

    C* = { y : <y, x> >= 0 for all x in C }

so a single algorithm answers both questions. Polytopes are covered as
well: a polytope of dimension `d` is read as a cone in `R^{d+1}` by
homogenization, its vertices becoming extreme rays and its facets
staying facets.

The difficulty is one of size before it is one of algorithm. The number
of facets of a cone with `n` generators in dimension `d` can grow like
`n^{floor(d/2)}`, the bound of McMullen's upper bound theorem, and the
cones that arise below have facets in numbers that no machine will ever
list. The difficulty is not merely that of writing the output: counting
the vertices of the Voronoi cell of a lattice given by a basis is
already a #P-hard problem [5]. What makes the classical methods suffer
in practice is degeneracy, meaning extreme rays lying on many more
facets than the dimension requires, and degeneracy is exactly what the
cones coming from lattices and from combinatorial optimization have.

The classical methods answer the question directly: Fourier-Motzkin
elimination, the double description method of Motzkin, Raiffa, Thompson
and Thrall [A], the reverse search of Avis and Fukuda [B], and the
methods based on triangulations. Each is excellent in its own range,
and each stops well before the problems of interest here.


Symmetry
--------

What breaks the deadlock is that the cones one actually wants to
describe are highly symmetric. A finite group `G` of linear
transformations preserving the cone permutes its extreme rays and
permutes its facets, and the natural question becomes not the list of
facets but the list of *orbits* of facets. The gain is not a constant
factor. The contact polytope of the Leech lattice has

    1197362269604214277200 facets, falling into 232 orbits [6]

and the Birkhoff polytope of the Coxeter group `H4`, the convex hull of
its 14400 elements taken as matrices, has its facets in 1063 orbits [7].
In such a range the orbit list is the only description that can be
written down, and it is also the mathematically meaningful one: two
facets in the same orbit are the same inequality written in different
coordinates.

Working with orbits requires deciding when two faces are equivalent
under `G` and, better, attaching to each face a canonical representative
so that equality of representatives decides equivalence. It also
requires knowing `G` in the first place, which is a question of its own:
the symmetry groups of a polyhedron preserving its linear, projective or
combinatorial structure are studied in [3], where the computation of the
linear one is reduced to a graph automorphism problem.


The adjacency decomposition method
----------------------------------

The method that makes the symmetric problems tractable is the adjacency
decomposition, developed for combinatorial polytopes by Christof and
Reinelt [C] and treated together with the other symmetry-exploiting
conversion techniques in the survey [1]. It rests on a simple geometric
observation: a ridge, that is a face of codimension two, is contained in
exactly two facets. The facets of the cone therefore form a graph, in
which two facets are adjacent when they share a ridge, and one may
explore that graph instead of enumerating its vertices from scratch.

The exploration goes as follows. One facet is found, by linear
programming or by any direct method. That facet is itself a polyhedral
cone of one dimension less, and its own facets are precisely the ridges
lying on it; finding them is again a dual description problem, in a
smaller dimension and with the smaller group that stabilizes the facet.
Across each of those ridges there is exactly one other facet of the
original cone, obtained by a flip. The facets so produced are sorted
into orbits under `G`, the new orbits are put on a queue, and the
process repeats until no new orbit appears.

Two things make this work. The first is that the recursion is genuine:
the subproblems are dual descriptions of the same nature, so the method
calls itself, and the depth at which it stops calling itself and hands
over to a direct method is a matter of judgement rather than of
principle. The second is that the subproblems repeat. The facets of a
symmetric cone fall into few isomorphism types, and the same small cone
is met again and again in the course of the computation; recognizing it
by its canonical form and remembering the answer turns an exponential
amount of repeated work into a lookup.

The algorithms as they were assembled for a classification of real size
are described in [2].


Why the exploration is complete
-------------------------------

The correctness of the adjacency decomposition rests on the connectivity
of the ridge graph. This is Balinski's theorem [D]: the graph of a
`d`-dimensional polytope is `d`-connected. Applied to the dual polytope
it says that the graph whose vertices are the facets and whose edges are
the ridges is connected, indeed `d`-connected, so an exploration that
starts anywhere and never refuses an adjacency reaches every facet.

Balinski's theorem does more than guarantee correctness in the limit: it
gives a stopping criterion. Suppose the exploration has produced a
partial list of orbits, and let the extreme rays not covered by that
list be called undone. If the undone rays are too few, or span too small
a subspace, then removing them cannot disconnect the ridge graph, and
the `d`-connectivity forces the partial list to be complete. That
argument admits refinements by linear programming and by a rank
computation on the undone rays, and it can itself be applied recursively
to the faces. The technique, and its use to settle the facets of cut
polytopes over highly symmetric graphs, is described in [4].


Face lattices and sampling
--------------------------

Two neighbouring questions are answered by the same circle of ideas. The
full face lattice of a polytope, that is all the faces of all dimensions
with their inclusions, is obtained by iterating the passage from a face
to its own facets, again up to the group. And when the complete list of
facets is out of reach, one may still want a sample of them, for
instance to bound a linear functional over the cone or to produce
certificates; facets can be produced one at a time by linear programming
from a point outside the cone, which gives a usable, if incomplete,
description.


What has been reached this way
------------------------------

* The contact polytope of the Leech lattice, the convex hull of its
  196560 shortest vectors: its 1197362269604214277200 facets classified
  into 232 orbits [6].

* The Birkhoff polytopes of the Coxeter groups `F4` and `H4`, that is
  the convex hulls of their elements taken as matrices [7]. The answer
  is 2 orbits of facets for `F4`, which contradicts what had been
  published, and 1063 orbits for `H4`, which disproves a conjecture. The
  point is worth making: these computations are not only expensive, they
  are also the kind of thing about which the literature can be wrong
  until someone carries them out.

* The perfect lattices, by Voronoi's algorithm of 1908, which passes
  from one perfect form to its neighbours across the facets of its
  Voronoi domain and therefore needs a dual description at every step.
  The 10916 of dimension 8 [2], a classification that rests entirely on
  exploiting symmetry in the polyhedral computations, and then the
  2237251040 of dimension 9 [8]. The latter settles the lattice packing
  problem in that dimension: the laminated lattice `Lambda_9` is the
  densest, the Hermite constant `gamma_9` equals 2, and the possible
  kissing numbers are exactly `2 * {1, ..., 91, 99, 120, ..., 129, 136}`.
  An enumeration of that size is the strongest argument there is for the
  machinery described above.

* The facets of the cut polytopes over highly symmetric graphs with 15
  to 30 edges [4], where the Balinski-based criterion above is what
  allows the enumerations to be certified complete.

* The vertices of Voronoi cells of lattices in dimensions up to about
  12, for which the symmetric algorithm of [5] is what makes a #P-hard
  counting problem practical.


References
----------

The complete and current publication list, with links, is at

    https://mathieudutsik.github.io/Publications/index.html

### On the method

[1] D. Bremner, M. Dutour Sikirić, A. Schürmann, *Polyhedral
representation conversion up to symmetries*, CRM Proceedings **48**
(2009) 45–72, arXiv:math/0702239. A survey of the conversion techniques
that exploit symmetry: the decomposition methods reducing the problem to
lower-dimensional subproblems, an incremental method generalizing
Fourier-Motzkin elimination, and the use of pivots.

[2] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *Classification of
eight-dimensional perfect forms*, Electronic Research Announcements of
the American Mathematical Society **13** (2007) 21–32,
arXiv:math/0609388. Describes the algorithms that made the
classification possible.

[3] D. Bremner, M. Dutour Sikirić, D. V. Pasechnik, T. Rehn,
A. Schürmann, *Computing symmetry groups of polyhedra*, LMS Journal of
Computation and Mathematics **17-1** (2014) 565–581, arXiv:1210.0206.
The linear, projective and combinatorial symmetry groups, and the
reduction of the linear one to graph automorphism.

[4] M. Deza, M. Dutour Sikirić, *Enumeration of the facets of cut
polytopes over some highly symmetric graphs*, International Transactions
in Operational Research **23-5** (2016) 853–860, arXiv:1501.05407. The
source of the Balinski-based termination criterion used above.

[5] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *Complexity and
algorithms for computing Voronoi cells of lattices*, Mathematics of
Computation **78** (2009) 1713–1731, arXiv:0804.0036. The #P-hardness of
counting the vertices, and an algorithm suited to highly symmetric
lattices.

### Computations carried out with it

[6] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *The contact polytope
of the Leech lattice*, Discrete and Computational Geometry **44** (2010)
904–911, arXiv:0906.1427.

[7] M. Dutour Sikirić, *The Birkhoff polytope of the groups `F4` and
`H4`*, Proceedings of the 4th Croatian Combinatorial Days, edited by
T. Došlić, S. Majstorović and L. Podrug (2023) 21–26, arXiv:2212.08452.

[8] M. Dutour Sikirić, W. van Woerden, *The lattice packing problem in
dimension 9 by Voronoi's algorithm*, preprint arXiv:2508.20719, data at
Zenodo, record 15707640.

### Classical background

[A] T. S. Motzkin, H. Raiffa, G. L. Thompson, R. M. Thrall, *The double
description method*, in Contributions to the Theory of Games II, Annals
of Mathematics Studies **28**, Princeton University Press (1953) 51–73.

[B] D. Avis, K. Fukuda, *A pivoting algorithm for convex hulls and
vertex enumeration of arrangements and polyhedra*, Discrete and
Computational Geometry **8** (1992) 295–313.

[C] T. Christof, G. Reinelt, *Decomposition and parallelization
techniques for enumerating the facets of combinatorial polytopes*,
International Journal of Computational Geometry and Applications **11**
(2001) 423–437.

[D] M. L. Balinski, *On the graph structure of convex polyhedra in
n-space*, Pacific Journal of Mathematics **11** (1961) 431–434.
