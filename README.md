polyhedral common
=================

This is the set of functionality for dealing with polytopes.

Overview
--------

Graphs are typically obtained as combinatorial data but for planar graphs and torus we are
typically interested in a graphical representation which helps reason about them.
This program allows to do exactly theat in a fairly efficient way. Features:

  * The graph can handle graphs having 1-gons and 2-gons 
  * For finite planar graph, it is possible to set the exterior face to be a face (standard case), an edge (in which case two half edges go to the external edge) or a vertex (in which case N directed edges go to the exterior vertex)
  * Iterative algorithms are used to find the coordinates.
  * The files are output to the SVG file format which thus can be used in Web pages or exported to .eps / .png / .pdf by using inkscape.

Finding coordinates
-------------------

Two main algorithms have been implemented:

  * The Primal Dual Circle Packing approach of Bojan Mohar.
  * The CaGe approach which handles better the small refinement.

References
----------

The algorithms used in that code has been described in following publications:

  * M. Deza, M. Dutour Sikirić, Enumeration of the facets of cut polytopes over some highly symmetric graphs, preprint at arxiv:1501.05407, to appear in International Transactions in Operational Research
  * M. Dutour Sikirić, A. Schürmann, F. Vallentin, The contact polytope of the Leech lattice, Discrete and Computational Geometry 44 (2010) 904--911, preprint at arxiv:0906.1427

Dependencies
------------

Following dependencies are needed for compiling the code:

  * Eigen: http://eigen.tuxfamily.org/
  * Boost: http://www.boost.org/
  * GNU MultiPrecision Library (GMP): https://gmplib.org/

Other dependencies are needed but are integrated into the repository:

  * bliss library: http://www.tcs.hut.fi/Software/bliss/
  * permlib: https://github.com/tremlin/PermLib

