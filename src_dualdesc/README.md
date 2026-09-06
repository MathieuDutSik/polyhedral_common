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
cones that arise in the classification problems below routinely have
generators in the thousands and facets in numbers that no machine will
ever list. Worse, the difficulty is not merely the output: no algorithm
is known whose running time is polynomial in the size of the input and
of the output taken together, and for general unbounded polyhedra the
associated decision problem is NP-hard. Degeneracy, meaning extreme
rays lying on many more facets than the dimension requires, is what
makes the classical methods suffer, and degeneracy is exactly what the
cones coming from lattices and from combinatorial optimization have.

The classical methods answer the question directly: Fourier-Motzkin
elimination, the double description method of Motzkin, Raiffa, Thompson
and Thrall, the reverse search of Avis and Fukuda, and the methods based
on triangulations. Each of them is excellent in its own range and each
of them stops well before the problems of interest here.


Symmetry
--------

What breaks the deadlock is that the cones one actually wants to
describe are highly symmetric. A finite group `G` of linear
transformations preserving the cone permutes its extreme rays and
permutes its facets, and the natural question is then not the list of
facets but the list of *orbits* of facets. That list can be smaller by
many orders of magnitude, and it is usually the mathematically
meaningful object: for a cut polytope or a hypermetric cone, two facets
in the same orbit are the same inequality written in different
coordinates.

Working with orbits requires being able to decide when two faces are
equivalent under `G`, and, better, to attach to each face a canonical
representative so that equality of representatives decides equivalence.
Both the computation of the symmetry group of a polyhedron and the
canonical form of the objects involved are questions in their own
right; they are treated in [1] and, for the positive definite matrices
that index many of these cones, in [2].


The adjacency decomposition method
----------------------------------

The method that makes the symmetric problems tractable is the adjacency
decomposition. It rests on a simple geometric observation: a ridge, that
is a face of codimension two, is contained in exactly two facets. The
facets of the cone therefore form a graph, in which two facets are
adjacent when they share a ridge, and one may hope to explore that graph
instead of enumerating its vertices from scratch.

The exploration goes as follows. One facet is found, by linear
programming or by any direct method. That facet is itself a polyhedral
cone of one dimension less, and its own facets are precisely the ridges
lying on it; finding them is again a dual description problem, in a
smaller dimension and with the smaller group that stabilizes the facet.
Across each of those ridges there is exactly one other facet of the
original cone, and it is obtained by a flip. The facets so produced are
sorted into orbits under `G`, the new orbits are put on a queue, and the
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

The method appears in [3] and its use for the extreme Delaunay polytopes
in [4]. The treatment of the symmetry it needs is that of [1].


Why the exploration is complete
-------------------------------

The correctness of the adjacency decomposition rests on the connectivity
of the ridge graph. This is Balinski's theorem [B]: the graph of a
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
to the faces. The technique, and its use to reach the facets of cut
polytopes over highly symmetric graphs, is described in [5].


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


Where this is used
------------------

The dual description of symmetric cones is not an end in itself. It is
the computational engine behind a number of classification results, and
the list gives a fair idea of the sizes involved.

* Lattices and quadratic forms. The classification of the eight
  dimensional perfect forms [6] rests on the description of the cones
  attached to the perfect forms; the generalization of Voronoi's
  reduction theory of [7] and the complete classification of the five
  dimensional Dirichlet-Voronoi polyhedra of translational lattices [8]
  are of the same nature, as is the classification of six dimensional
  iso-edge domains [9] and the general framework for them [10]. The
  Voronoi cells themselves, their complexity and their computation, are
  treated in [11], and the contact polytope of the Leech lattice, a
  polytope on 196560 vertices, in [12].

* Delaunay polytopes. The six dimensional Delaunay polytopes [13], the
  infinite series of extreme ones [14], the perfect Delaunay polytopes
  in low dimension [15], the seven dimensional perfect Delaunay polytopes
  [16], the rank of a Delaunay polytope [17] and the Delaunay polytopes
  derived from the Leech lattice [18] all lean on the same machinery.

* Metric and hypermetric geometry. The hypermetric cone on seven
  vertices [19] and on eight vertices [20], its decomposition into
  L-domains [21], the hypermetric cone and polytope on graphs [22], the
  cut polytopes over highly symmetric graphs [5] and their
  generalizations [23]. The general setting is that of [24].

* Parallelohedra and tilings. The Voronoi conjecture for combinatorially
  Voronoi parallelohedra [25], the sum of a parallelotope and a zonotope
  [26], zonotopes and parallelotopes [27], the periodic triangulations of
  `Z^n` [28] and the Voronoi polytopes for polyhedral norms [29].

* Further afield: the colouring of the Voronoi tessellation of lattices
  [30], the smoothness and singularities of the perfect form
  compactification of `A_g` [31], the Voronoi complexes in higher
  dimensions and the cohomology of `GL_N(Z)` [32], and rational
  copositive factorization [33].


References
----------

The papers below are those of M. Dutour Sikirić and coauthors; the
complete and current list, with the links, is at

    https://mathieudutsik.github.io/Publications/index.html

[B] M. L. Balinski, *On the graph structure of convex polyhedra in
n-space*, Pacific Journal of Mathematics **11** (1961) 431–434.

[1] D. Bremner, M. Dutour Sikirić, D. V. Pasechnik, T. Rehn,
A. Schürmann, *Computing symmetry groups of polyhedra*, LMS Journal of
Computation and Mathematics **17-1** (2014) 565–581.

[2] M. Dutour Sikirić, A. Haensch, J. Voight, W. van Woerden, *A
canonical form for positive definite matrices*, Proceedings of the
Fourteenth Algorithmic Number Theory Symposium (ANTS-XIV), Open Book
Series 4, Mathematical Sciences Publishers, 2020.

[3] D. Bremner, M. Dutour Sikirić, A. Schürmann, *Polyhedral
representation conversion up to symmetries*, CRM Proceedings **48**
(2009) 45–72.

[4] M. Dutour, *Adjacency method for extreme Delaunay polytopes*,
Voronoi's Impact on Modern Science, Book 3 (2005) 94–101.

[5] M. Deza, M. Dutour Sikirić, *Enumeration of the facets of cut
polytopes over some highly symmetric graphs*, International Transactions
in Operational Research **23-5** (2016) 853–860.

[6] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *Classification of
eight-dimensional perfect forms*, Electronic Research Announcements of
the American Mathematical Society **13** (2007) 21–32.

[7] A. Schürmann, M. Dutour Sikirić, F. Vallentin, *A generalization of
Voronoi's reduction theory and its application*, Duke Mathematical
Journal **142** (2008) 127–164.

[8] M. Dutour Sikirić, A. Garber, A. Schürmann, C. Waldmann, *The
complete classification of five-dimensional Dirichlet-Voronoi polyhedra
of translational lattices*, Acta Crystallographica A **72** (2016)
673–683.

[9] M. Dutour Sikirić, W. van Woerden, *Complete classification of
six-dimensional iso-edge domains*, Acta Crystallographica A **81** (2025)
9–15.

[10] M. Dutour Sikirić, M. Kummer, *Iso-edge domains*, Expositiones
Mathematicae **40-2** (2022) 302–314.

[11] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *Complexity and
algorithms for computing Voronoi cells of lattices*, Mathematics of
Computation **78** (2009) 1713–1731.

[12] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *The contact polytope
of the Leech lattice*, Discrete and Computational Geometry **44** (2010)
904–911.

[13] M. Dutour, *The six-dimensional Delaunay polytopes*, European
Journal of Combinatorics **25** (2004) 535–548.

[14] M. Dutour, *Infinite series of extreme Delaunay polytopes*, European
Journal of Combinatorics **26** (2005) 129–132.

[15] M. Dutour Sikirić, R. Erdahl, K. Rybnikov, *Perfect Delaunay
polytopes in low dimensions*, Integers **7** (2007) A39.

[16] M. Dutour Sikirić, *The seven dimensional perfect Delaunay polytopes
and Delaunay simplices*, Canadian Journal of Mathematics **69** (2017)
1143–1168.

[17] M. Dutour Sikirić, V. Grishukhin, *How to compute the rank of a
Delaunay polytope*, European Journal of Combinatorics **28** (2007)
762–773.

[18] M. Dutour Sikirić, K. Rybnikov, *Delaunay polytopes derived from the
Leech lattice*, Journal de Théorie des Nombres de Bordeaux **26-1** (2014)
85–101.

[19] M. Deza, M. Dutour Sikirić, *The hypermetric cone on seven
vertices*, Experimental Mathematics **12** (2004) 433–440.

[20] M. Deza, M. Dutour Sikirić, *The hypermetric cone on eight vertices
and some generalizations*, Journal of Symbolic Computation **88** (2018)
67–84.

[21] M. Dutour Sikirić, V. Grishukhin, *The decomposition of the
hypermetric cone into L-domains*, European Journal of Combinatorics **30**
(2009) 853–865.

[22] M. Dutour Sikirić, *The hypermetric cone and polytope on graphs*,
Chebyshevskii Sbornik **20-2** (2019) 160–168.

[23] M. Deza, M. Dutour Sikirić, *Generalized cut and metric polytopes of
graphs and simplicial complexes*, Optimization Letters **14** (2020)
273–289.

[24] M. Deza, M. Dutour Sikirić, E. Deza, *Generalizations of finite
metrics and cuts*, World Scientific, 2016.

[25] M. Dutour Sikirić, A. Garber, A. Magazinov, *On the Voronoi
conjecture for combinatorially Voronoi parallelohedra*, SIAM Journal on
Discrete Mathematics **34-4** (2020) 2481–2501.

[26] M. Dutour Sikirić, V. Grishukhin, A. Magazinov, *On the sum of a
parallelotope and a zonotope*, European Journal of Combinatorics **42**
(2014) 49–73.

[27] M. Dutour Sikirić, V. Grishukhin, *Zonotopes and parallelotopes*,
Southeast Asian Bulletin of Mathematics **41-2** (2017) 197–207.

[28] M. Dutour Sikirić, A. Garber, *Periodic triangulations of `Z^n`*,
Electronic Journal of Combinatorics **27** (2020) P2.36.

[29] M. Deza, M. Dutour Sikirić, *Voronoi polytopes for polyhedral norms
on lattices*, Discrete Applied Mathematics **197** (2015) 42–52.

[30] M. Dutour Sikirić, D. Madore, P. Moustrou, F. Vallentin, *Coloring
the Voronoi tessellation of lattices*, Journal of the London Mathematical
Society **104-2** (2021) 1135–1171.

[31] M. Dutour Sikirić, K. Hulek, A. Schürmann, *Smoothness and
singularities of the perfect form compactification of `A_g`*, Algebraic
Geometry **2-5** (2015) 642–653.

[32] M. Dutour Sikirić, P. Elbaz-Vincent, A. Kupers, J. Martinet,
*Voronoi complexes in higher dimensions, cohomology of `GL_N(Z)` for
`N >= 8` and the triviality of `K_8(Z)`*, Journal of the Institute of
Mathematics of Jussieu.

[33] M. Dutour Sikirić, A. Schürmann, F. Vallentin, *A simplex algorithm
for rational CP-factorization*, Mathematical Programming **187** (2020)
25–45.
