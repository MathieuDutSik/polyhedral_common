#ifndef POLY_PARALL_H
#define POLY_PARALL_H


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
   --- VT: So we need a data type for the vector of face. This is the output
   of the function.
   --- DT: We need as well the domain of all T for processing the data.
   --- CT: We need a control type for storing the data.

   We do the following:
   --- In serial, we would have a single DT.
   --- In MPI, each process would have its own DT.
   --- In C++11 thread, we would have a std::vector<DT> and each process consider one

   


 */






#endif
