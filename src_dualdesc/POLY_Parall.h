// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_PARALL_H_
#define SRC_DUALDESC_POLY_PARALL_H_

/* Let us try here to design a parallel enumeration structure.
   This is for the following enumerations:
   --- Facets of polytope
   --- Delaunay polytopes of lattice
   ... Perfect forms
   ... Perfect Lorentzian domains.
   What we want is some templatized code that can handle:
   --- Serial runs
   --- MPI runs
   --- C++11 threads

   Possible simple design:
   --- T: the types of the data being enumerated.
   The problem with that is that if we encode the faces of polytope this leads
   us to a std::vector<Face> which is suboptimal compared to vectface.
   --- VT (For "Vector T"): So we need a data type for the vector of face. This
   is the output of the function.
   --- DT (For "Domain T"): We need as well the domain of all T for processing
   the data.
   --- CT (For "Control T"): We need a control type for checking when things are
   finished or not.


   We do the following:
   --- In serial, we would have a single DT.
   --- In MPI, each process would have its own DT.
   --- In C++11 thread, we would have a std::vector<DT> and each process
   consider one

   It seems almost impossible to design a parallel code such that when doing the
   recursive calls, we have the full set of processors sill available. Unless
   maybe we have something like executors to use for the running of the code.

   In MPI for the computation of facets, if we use the recursive adjacency
   decomposition then we need to have a sending of data from all processors.


   There are many possible ways to make bad design decision, as happened to us
   before. What could be the right moves?
   --- Make the database bank a web service. This would cover the MPI and the
   different cone case. But we want this to work also in the standard case
   --- Offload a lot of code to the DatabaseOrbits.


   But we want a modularity in that respect.
   We need to avoid repetitions





 */

// clang-format off
#endif  // SRC_DUALDESC_POLY_PARALL_H_
// clang-format on
