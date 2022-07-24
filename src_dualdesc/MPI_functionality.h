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




struct request_status_list {
  size_t MaxNumberFlyingMessage;
  std::vector<boost::mpi::request> l_mpi_request;
  std::vector<int> l_mpi_status;
  request_status_list(size_t const& MaxNumberFlyingMessage) : MaxNumberFlyingMessage(MaxNumberFlyingMessage), l_mpi_request(MaxNumberFlyingMessage), l_mpi_status(MaxNumberFlyingMessage,0) {
  }
  boost::mpi::request& operator[](size_t const& pos) { l_mpi_status[pos]=1; return l_mpi_request[pos]; }
};


size_t GetFreeIndex(request_status_list & rsl) {
  size_t len = rsl.l_mpi_request.size();
  for (size_t i=0; i<len; i++) {
    if (rsl.l_mpi_status[i] == 1) {
      boost::optional<boost::mpi::status> stat = rsl.l_mpi_request[i].test();
      if (stat) {
        rsl.l_mpi_status[i] = 0;
      }
    }
  }
  for (size_t i = 0; i < len; i++)
    if (rsl.l_mpi_status[i] == 0)
      return i;
  if (rsl.MaxNumberFlyingMessage > 0)
    return std::numeri_limits<size_t>::max();
  rsl.l_mpi_request.push_back(boost::mpi::request());
  rsl.l_mpi_status.push_back(1);
  return len;
}


struct empty_message_management {
  boost::mpi::communicator & comm;
  request_status_list rsl;
  int expected_value; // Random, but identical on all process.
  int tag;
  empty_message_management(boost::mpi::communicator & comm, size_t const& MaxFly, int const& tag) : comm(comm), rsl(MaxFly), tag(tag) {
    int expected_value_pre = random();
    expected_value = boost::mpi::all_reduce(comm, expected_value_pre, mpi::minimum<int>());
  }
  void send_message(int dest) {
    size_t idx = rsl.GetFreeIndex();
    if (idx == std::numeric_limits<size_t>::max()) {
      std::cerr << "Failed to find an index\n";
      throw TerminalException{1};
    }
    rsl[idx] = comm.isend(dest, tag, expected_value);
  }
  void recv_message(int source) {
    int recv_message;
    comm.recv(source, tag, recv_message);
    if (recv_message != expected_value) {
      std::cerr << "The recv_message is incorrect\n";
      throw TerminalException{1};
    }
  }


};






// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
