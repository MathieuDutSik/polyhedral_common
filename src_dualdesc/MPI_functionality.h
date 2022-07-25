// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_MPI_FUNCTIONALITY_H_
#define SRC_DUALDESC_MPI_FUNCTIONALITY_H_


#include <string>
#include <vector>
#include "Timings.h"
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
  size_t clear_and_get_nb_undone() {
    size_t len = l_mpi_request.size();
    size_t nb_undone = 0;
    for (size_t i=0; i<len; i++) {
      if (l_mpi_status[i] == 1) {
        boost::optional<boost::mpi::status> stat = l_mpi_request[i].test();
        if (stat) {
          l_mpi_status[i] = 0;
        } else {
          nb_undone++;
        }
      }
    }
    return nb_undone;
  }
  size_t GetFreeIndex() {
    (void)clear_and_get_nb_undone();
    size_t len = l_mpi_request.size();
    for (size_t i = 0; i < len; i++)
      if (l_mpi_status[i] == 0)
        return i;
    if (MaxNumberFlyingMessage > 0)
      return std::numeric_limits<size_t>::max();
    l_mpi_request.push_back(boost::mpi::request());
    l_mpi_status.push_back(1);
    return len;
  }
  bool is_empty() {
    size_t n_undone = clear_and_get_nb_undone();
    return n_undone == 0;
  }
};



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




/* Managing exchanges at scale.
   We can T = int, T_vector = std::vector<T>
   or T = Face, T_vector = vectface
 */
template<typename T, typename T_vector>
struct buffered_T_exchanges {
  boost::mpi::communicator & comm;
  request_status_list rsl;
  int tag;
  T_vector empty;
  int n_proc;
  std::vector<T_vector> l_message;
  std::vector<T_vector> l_under_cons;
  buffered_T_exchanges(boost::mpi::communicator & comm, size_t const& MaxFly, int const& tag, T_vector const& empty) : comm(comm), rsl(MaxFly), tag(tag), empty(empty), n_proc(comm.size()), l_message(n_proc, empty), l_under_cons(MaxFly, empty) {
  }
  void insert_entry(size_t const& pos, T const& x) {
    l_message[pos].push_back(x);
  }
  T_vector recv_message(int source) {
    T_vector vf = empty;
    comm.recv(source, tag, vf);
    return vf;
  }
  void clear_one_entry() {
    int idx = rsl.GetFreeIndex();
    if (idx == -1)
      break;
    size_t max_siz = 0;
    int chosen_iproc = -1;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      size_t siz = ListFaceUnsent[i_proc].size();
      if (siz > max_siz) {
        max_siz = siz;
        chosen_iproc = i_proc;
      }
    }
    if (chosen_iproc == -1)
      break;
    l_under_cons[idx] = std::move(l_message[chosen_iproc]);
    l_message[chosen_iproc].clear();
    rsl[idx] = comm.isend(res, tag_new_facets, l_under_cons[idx]);
  }
  size_t get_unsent_size() const {
    size_t n_unsent = 0;
    for (auto & eEnt : l_message)
      n_unsent += eEnt.size();
    return n_unsent;
  }
  bool is_buffer_empty() {
    if (!rsl.is_empty())
      return false;
    return get_unsent_size() == 0;
  }
}



template<typename Tint>
struct database_balinski_info {
  boost::mpi::communicator & comm;
  int tag_request;
  int tag_info;
  int n_proc;
  int i_rank;
  request_status_list rsl;
  std::vector<UndoneOrbitInfo<Tint>> ListBalinski;
  // 0: not assigned
  // 1: Assigned but buffers are not free
  // 2: Assigned and buffers are free
  std::vector<int> ListStatus_Emptyness;
  Tint CritSiz;
  int expected_value;
  database_balinski_info(boost::mpi::communicator & comm, int tag_request, int tag_info, Tint const& CritSiz) : comm(comm),
    tag_request(tag_request), tag_info(tag_info), n_proc(comm.size()), i_rank(comm.rank()), rsl(n_proc),
    ListBalinski(n_proc), ListStatus_Emptyness(n_proc,0), CritSiz(CritSiz), expected_value(47) {
  }
  bool get_status() {
    for (int i_proc=0; i_proc<n_proc; i_proc++)
      if (ListStatus_Emptyness[i_proc] != 2)
        return false;
    UndoneOrbitInfo<Tint> uoi = CombineUndoneOrbitInfo(ListBalinski);
    return ComputeStatusUndone(uoi, CritSiz);
  }
  void reply_request(int source_dest, UndoneOrbitInfo<Tint> const& uoi, bool const& status) {
    int recv_message;
    comm.recv(source, tag, recv_message);
    if (recv_message != expected_value) {
      std::cerr << "The recv_message is incorrect\n";
      throw TerminalException{1};
    }
    //
    int idx = rsl.GetFreeIndex();
    if (idx == std::numeric_limits<size_t>::max())
      return;
    StatusUndoneOrbitInfo<Tint> suoi{status, uoi};
    rsl[idx] = comm.isend(res, tag_info, suoi);
  }
  void recv_message(int source) {
    StatusUndoneOrbitInfo<Tint> suoi;
    comm.recv(source, tag_info, recv_message);
    ListBalinski[source] = suoi.erec;
    ListStatus_Emptyness[source] = 1 + int(suoi.status);
  }
  void submit_request(int dest) {
    int idx = rsl.GetFreeIndex();
    if (idx == std::numeric_limits<size_t>::max())
      return;
    rsl[idx] = comm.isend(res, tag_request, expected_value);
  }
  template<typename F>
  void submit_uoi(UndoneOrbitInfo<Tint> const& uoi, F f) {
    // First checking natively
    ListBalinski[i_rank] = uoi;
    if (ComputeStatusUndone(uoi, CritSiz))
      return;
    ListStatus_Emptyness[i_rank] = 1 + int(f());
    // First checking for unassigned
    for (int i_proc=0; i_proc<n_proc; i_proc++)
      if (ListStatus_Emptyness[i_proc] == 0)
        return submit_request(i_proc);
    // First checking for non-empty buffer
    for (int i_proc=0; i_proc<n_proc; i_proc++)
      if (ListStatus_Emptyness[i_proc] == 1)
        return submit_request(i_proc);
    // Getting the highest value
    Tint max_val = 0;
    int chosen_idx = -1;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      Tint const& eval = ListBalinski[i_proc].nbUndone;
      if (eval > max_val) {
        max_val = eval;
        chosen_idx = i_proc;
      }
    }
    if (chosen_idx == -1)
      return;
    return submit_request(chosen_idx);
  }
};







// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
