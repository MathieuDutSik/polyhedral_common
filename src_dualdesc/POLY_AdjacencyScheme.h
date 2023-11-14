// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_

/*
  We want here fully templatized code that allows to work with
  general code. Common features:
  ---Every entry has to be treated, no Balinski heuristic.
  ---Data can be written to file using boost.
  ---We need following functions:
     ---f_init for finding the initial.
     ---f_hash for getting the hash (of the canonicalized objects)
        The f_hash has to take into acount a seed, because we need
        to compute a hash for the partitioning and a hash for an HashMap.
     ---f_adj for getting the list of adjacent object (which should use
        the stabilizer computation)
     ---f_repr for testing the equivalence of the objects.
  ---We want to be able to restart the computation from saves.
  ---We want to be able to stop the communication in the moddle.
  ---No heuristic would be relevant

  Examples to be treatable:
  ---The Delaunay polytopes of a lattice.
     ---We can do canonicalization stuff if we have a spanning configuration
        of vectors, but no guarantee really.
     ---If we do not have a reasonable spanning family (which is the case
        for Niemeier) then this does not work. TODO: Have a canonicalization
        for
  ---The perfect forms of a T-space.
     ---There is likely no canonicalization for the T-space coming from number
        theory
  ---The perfect form of the hyperbolic cone.
     ---This might work similarly to the Delaunay. Right now we have
  ---The C-types of Z^n (but we do not have realistic examples)
     ---We have canonicalization stuff.

  Conclusion: The canonicalization cannot be assumed to exist in general.
  But we could work out with f_repr. And if by accident we have some canonicalization
  then we can have a corresponding hash function and then an equality test. So,
  there is no loss of generality in working that way.

  TODO:
  ---Resolve the scheme so that we can conclude when everything has been treated
  without any heuristic.
  ---When a node has everything processed, send a note to all others asking for their
  status. If all of them reply then things are clear.
  ---Have to learn about "distributed programming" from Leslie Lamport and others.
  Takada_M._Distributed_systems_for_fun_and_profit.epub (short one)
  Varela_C.A._Programming_distributed_computing_systems__a_foundational_approach.pdf
  (ok, teaches pi-calculus, join, etc.)



  
 */




/*
  input clear from preceding discussion.
 */

template<typename Tobj, typename Finit, typename Fadj, typename Fhash, typename Frepr>
bool compute_adjacency(boost::mpi::communicator &comm, std::string const& Prefix, int max_time_second,
                       Finit f_init, Fadj f_adj, Fhash f_hash, frepr f_repr) {

  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();

  std::vector<Tobj> V;
  std::vector<int> Vstatus;
  int n_obj = 0;
  while(true) {
    std::string FileI = Prefix + "_O_" + std::to_string(n_obj);
    if (IsExistingFile(FileI)) {
      break;
    }
    {
      std::ifstring isO(FileI);
      Tobj x;
      x << isO;
      V.emplace_back(std::move(x));
    }
    //
    {
      std::string FileS = Prefix + "_S_" + std::to_string(n_obj);
      std::ifstring isS(FileS);
      int status;
      status << isS;
      Vstatus.push_back(status);
    }
    //
    n_obj++;
  }
  //
  while(true) {
  }
}






// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
