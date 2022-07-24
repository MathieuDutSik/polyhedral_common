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

template<typename T>
std::vector<T> mpi_gather(boost::mpi::communicator & comm, T const& x, int const& i_proc) {
  int i_rank = comm.rank();
  std::vector<T> V;
  if (i_rank == i_proc) {
    boost::mpi::gather<T>(comm, x, V, i_proc);
  } else {
    boost::mpi::gather<T>(comm, x, i_proc);
  }
  return V;
}

template<typename T>
std::vector<T> mpi_allgather(boost::mpi::communicator & comm, T const& x) {
  std::vector<T> V;
  boost::mpi::all_gather<T>(comm, x, V);
  return V;
}




// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
