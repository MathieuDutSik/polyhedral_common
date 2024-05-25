Polytope computations
=====================

The set of functions of this directory is directly related
to polytopes and used by other directories.

Geometric functions
-------------------

The following are provided:
  * **POLY_CreateAffineBasis** is for computingan affine basis of a polytope
  * **VectFamily_ColumnReduction** is for dropping columns and getting a full dimensional set of vectors.
  * **VectFamily_Reduction** is for reducing the size of the coefficients of a matrix by applying a LLL transformation

References for the affine basis:
  * Mathieu Dutour Sikirić, Viatcheslav Grishukhin, How to compute the rank of a Delaunay polytope, preprint at arxiv:0512193, European Journal of Combinatorics 28 (2007) 762--773


Linear Programming functions
----------------------------

The following functions are provided:
  * **POLY_GetFullRankFacetSet** Get a set of facets of a polytope that is full dimensional.
  * **POLY_LinearDetermineByInequalities** gets the linear space that is determined by a set of inequalities
  * **POLY_SolutionMatNonnegative** resolve equation **Ax = b** for **x >= 0** which is a linear program
  * **POLY_SolutionMatNonnegativeComplete** the same but provides a certificate if there is no solution.
  * **POLY_cdd_LinearProgramming** the simple call to **cdd** LP.
  * **POLY_IsPointedCone** for testing whether it is a pointed cone or not.

Another program that is a little bit more intricate is **TEST_GeometricallyUniquePoint** that computes a
geometrically unique integral point to a polyhedral cone defined by facet inequalities. But if you apply
an integral transformation to the defining inequalities then the obtained solution is changed accordingly.
In other words, there is no dependency on the chosen vertex and it can be used for isomorphism checks.

Other points:
  * For some problems, there are several direct and dual functions for addressing the problem. The choice is done via dimension.
  * When available, it is generally a good idea to use linear programming instead of dual description.
  * 

Dual description functions
--------------------------

The following functions are used for the dual description:
  * **POLY_SmallPolytope** for computing the dual description of polyhedral cones with n or n+1 facets
  * **POLY_sampling_facets** for sampling a set of facets of a polytope.
  * **POLY_dual_description** for computing the dual description of a polytope.
  * **POLY_dual_description_group** for computing the orbit representatives of a polytope (without using advanced techniques)
  * **POLY_cdd_skeletons** computing the dual description and the skeletons.

Other points:
  * For each dual description function, there is one function that returns a **vectface**, one function that returns a **MyMatrix<T>** and one function that returns a **std::vector<std::pair<Face,MyVector<T>>>**. The last function is taking the most memory.
  * The function provided are **lrs**, **lrs_ring**, **cdd**, **ppl**, **normaliz**.
  * Some of the function calls are external.
  * The method **pd_lrs** uses the primal dual technique by using lrs for the check.

Redundancy checks
-----------------

The following programs are provided:
  * **POLY_redundancy** for redundancy computation.
  * **POLY_redundancyClarksonCddlib** uses the Clarkson algorithm and double precision. Not 100% guaranteed to be correct.
  * **POLY_redundancyGroup** for using the group symmetries for the reduction of redundancies.

