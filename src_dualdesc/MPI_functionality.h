// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_MPI_FUNCTIONALITY_H_
#define SRC_DUALDESC_MPI_FUNCTIONALITY_H_


#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>



std::ofstream& get_standard_outstream(boost::mpi::communicator & comm) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string eFile = "err_" + std::to_string(i_rank) + "_" + std::to_string(n_proc);
  return std::ofstream(eFile);
}





// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
