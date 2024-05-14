// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <iostream>

int main() {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  size_t irank = world.rank();
  if (irank == 0) {
    std::vector<int> V(100, 0);
    boost::mpi::request req = world.isend(1, 37, V);
    while (true) {
      boost::optional<boost::mpi::status> stat = req.test();
      if (stat) {
        if (stat->error() != 0) {
          std::cerr << "Rank 0: Non-zero error in stat\n";
          char error_string[10000];
          int length_of_error_string;
          MPI_Error_string(stat->error(), error_string,
                           &length_of_error_string);
          fprintf(stderr, "err: %s\n", error_string);
          return 1;
        }
        break;
      }
    }
  } else {
    while (true) {
      boost::optional<boost::mpi::status> prob = world.iprobe();
      if (prob) {
        if (prob->tag() == 37) {
          std::vector<int> V;
          world.recv(prob->source(), prob->tag(), V);
          std::cerr << "Rank 1: Correct receiving of the vector\n";
          break;
        }
      }
    }
  }
}
