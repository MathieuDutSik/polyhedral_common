// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_MPI_FUNCTIONALITY_H_
#define SRC_DUALDESC_MPI_FUNCTIONALITY_H_


#include <string>
#include <vector>
#include "Timings.h"
#include "Balinski_basic.h"
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>






template<typename T>
std::vector<T> my_mpi_gather(boost::mpi::communicator & comm, T const& x, int const& i_proc) {
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
std::vector<T> my_mpi_allgather(boost::mpi::communicator & comm, T const& x) {
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
    expected_value = boost::mpi::all_reduce(comm, expected_value_pre, boost::mpi::minimum<int>());
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
  int n_proc;
  std::vector<T_vector> l_message;
  std::vector<T_vector> l_under_cons;
  buffered_T_exchanges(boost::mpi::communicator & comm, size_t const& MaxFly, int const& tag) : comm(comm), rsl(MaxFly), tag(tag), n_proc(comm.size()), l_message(n_proc), l_under_cons(MaxFly) {
  }
  void insert_entry(size_t const& pos, T const& x) {
    l_message[pos].push_back(x);
  }
  T_vector recv_message(int source) {
    T_vector vf;
    comm.recv(source, tag, vf);
    return vf;
  }
  void clear_one_entry(std::ostream & os) {
    size_t idx = rsl.GetFreeIndex();
    os << "idx=" << idx << "\n";
    if (idx == std::numeric_limits<size_t>::max())
      return;
    size_t max_siz = 0;
    int chosen_iproc = -1;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      size_t siz = l_message[i_proc].size();
      if (siz > max_siz) {
        max_siz = siz;
        chosen_iproc = i_proc;
      }
    }
    os << "max_siz=" << max_siz << " chosen_iproc=" << chosen_iproc << "\n";
    if (chosen_iproc == -1)
      return;
    l_under_cons[idx] = std::move(l_message[chosen_iproc]);
    l_message[chosen_iproc].clear();
    rsl[idx] = comm.isend(chosen_iproc, tag, l_under_cons[idx]);
  }
  size_t get_unsent_size() const {
    size_t n_unsent = 0;
    for (auto & eEnt : l_message)
      n_unsent += eEnt.size();
    return n_unsent;
  }
  bool is_buffer_empty() {
    //    if (!rsl.is_empty())
    //      return false;
    return get_unsent_size() == 0;
  }
};



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
  // 1: Assigned and buffers are free
  std::vector<int> ListStatus_Emptyness;
  // 0: No message is in flight, free to send another one.
  // 1: A message is in flight, please wait before sending another one.
  std::vector<int> RequestNat;
  Tint CritSiz;
  int expected_value;
  std::vector<std::chrono::time_point<std::chrono::system_clock>> last_update;
  std::chrono::time_point<std::chrono::system_clock> last_database_update_time;
  Face ToBeAnswered;
  database_balinski_info(boost::mpi::communicator & comm, int tag_request, int tag_info, Tint const& CritSiz) :
    comm(comm), tag_request(tag_request), tag_info(tag_info), n_proc(comm.size()), i_rank(comm.rank()),
    rsl(n_proc), ListBalinski(n_proc), ListStatus_Emptyness(n_proc,0), RequestNat(n_proc,0), CritSiz(CritSiz),
    expected_value(47), last_update(n_proc, get_cpp_time(7, 1, 1974)),
    last_database_update_time(get_cpp_time(6, 1, 1974)), ToBeAnswered(n_proc) {
  }
  bool get_status(std::ostream & os) const {
    for (int i_proc=0; i_proc<n_proc; i_proc++)
      if (ListStatus_Emptyness[i_proc] == 0) {
        int val = ListStatus_Emptyness[i_proc];
        os << "Returning false at i_proc=" << i_proc << " val=" << val << "\n";
        return false;
      }
    UndoneOrbitInfo<Tint> uoi = CombineUndoneOrbitInfo(ListBalinski);
    os << "Merged uoi=" << uoi << "\n";
    return ComputeStatusUndone(uoi, CritSiz);
  }
  void read_request(int source) {
    int recv_message;
    comm.recv(source, tag_request, recv_message);
    if (recv_message != expected_value) {
      std::cerr << "The recv_message is incorrect\n";
      throw TerminalException{1};
    }
    ToBeAnswered[source] = 1;
  }
  void submit_info(int dest, UndoneOrbitInfo<Tint> const& uoi) {
    size_t idx = rsl.GetFreeIndex();
    if (idx == std::numeric_limits<size_t>::max())
      return;
    rsl[idx] = comm.isend(dest, tag_info, uoi);
    last_update[dest] = last_database_update_time;
    ToBeAnswered[dest] = 0;
  }
  void submit_request(int dest) {
    size_t idx = rsl.GetFreeIndex();
    if (idx == std::numeric_limits<size_t>::max())
      return;
    RequestNat[dest] = 1;
    rsl[idx] = comm.isend(dest, tag_request, expected_value);
  }
  void recv_info(int source) {
    UndoneOrbitInfo<Tint> uoi;
    comm.recv(source, tag_info, uoi);
    ListBalinski[source] = uoi;
    ListStatus_Emptyness[source] = 1;
    RequestNat[source] = 0;
  }
  void set_uoi_local(UndoneOrbitInfo<Tint> const& uoi_local) {
    ListBalinski[i_rank] = uoi_local;
    ListStatus_Emptyness[i_rank] = 1;
    last_database_update_time = std::chrono::system_clock::now();
  }
  void flush(std::ostream & os) {
    os << "Beginning of flush operation\n";
    UndoneOrbitInfo<Tint> const& uoi_local = ListBalinski[i_rank];
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      if (i_proc != i_rank) {
        bool test = last_database_update_time > last_update[i_proc];
        int dur_i = std::chrono::duration_cast<std::chrono::seconds>(last_database_update_time - last_update[i_proc]).count();
        os << "i_proc=" << i_proc << " ToBeAnswered=" << ToBeAnswered[i_proc] << " test=" << test << " dur_i=" << dur_i << "\n";
        if (ToBeAnswered[i_proc] == 1 && test) {
          os << "flush: Doing submit_info for i_proc=" << i_proc << "\n";
          submit_info(i_proc, uoi_local);
        }
      }
    }
  }
  void submit_request_uoi(std::ostream & os) {
    // First checking natively
    os << "Beginning of submit_uoi\n";
    // First checking for unassigned
    for (int i_proc=0; i_proc<n_proc; i_proc++)
      if (i_proc != i_rank && ListStatus_Emptyness[i_proc] == 0 && RequestNat[i_proc] == 0) {
        os << "submit_request_uoi (Case 0) at i_proc=" << i_proc << "\n";
        return submit_request(i_proc);
      }
    // Getting the highest value
    Tint max_val = 0;
    int chosen_idx = -1;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      if (i_proc != i_rank) {
        Tint const& eval = ListBalinski[i_proc].nbUndone;
        if (eval > max_val) {
          max_val = eval;
          chosen_idx = i_proc;
        }
      }
    }
    if (chosen_idx == -1)
      return;
    if (RequestNat[chosen_idx] == 0) {
      os << "submit_request_uoi (Case 1) at chosen_idx=" << chosen_idx << "\n";
      return submit_request(chosen_idx);
    }
  }
};







// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
