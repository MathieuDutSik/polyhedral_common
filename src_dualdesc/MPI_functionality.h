// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_MPI_FUNCTIONALITY_H_
#define SRC_DUALDESC_MPI_FUNCTIONALITY_H_

// clang-format off
#include "Balinski_basic.h"
#include "Boost_bitset_kernel.h"
#include "MAT_Matrix.h"
#include "MPI_basic.h"
#include "Namelist.h"
#include "Timings.h"
#include <limits>
#include <list>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  For the vector, simple syntactic sugar
 */

#ifdef DEBUG
#define DEBUG_MPI_FCT
#endif

template <typename T>
inline typename std::enable_if<!is_mymatrix<T>::value, std::vector<T>>::type
my_mpi_gather(boost::mpi::communicator &comm, T const &x, int const &i_proc) {
  int i_rank = comm.rank();
  std::vector<T> V;
  if (i_rank == i_proc) {
    boost::mpi::gather<T>(comm, x, V, i_proc);
  } else {
    boost::mpi::gather<T>(comm, x, i_proc);
  }
  return V;
}

template <typename T>
inline typename std::enable_if<!is_mymatrix<T>::value, std::vector<T>>::type
my_mpi_allgather(boost::mpi::communicator &comm, T const &x) {
  std::vector<T> V;
  boost::mpi::all_gather<T>(comm, x, V);
  return V;
}

/*
  For vectface, a little bit more advanced work
 */

vectface my_mpi_gather(boost::mpi::communicator &comm, vectface const &vf,
                       int i_proc) {
  int i_rank = comm.rank();
  size_t n_vert = vf.get_n();
  size_t n_face = vf.size();
  std::vector<uint8_t> const &V = vf.serial_get_std_vector_uint8_t();
  //
  std::vector<size_t> l_n_face = my_mpi_gather(comm, n_face, i_proc);
  std::vector<std::vector<uint8_t>> l_V = my_mpi_gather(comm, V, i_proc);
  if (i_rank == i_proc) {
    return vectface(n_vert, l_n_face, l_V);
  } else {
    return vectface(n_vert);
  }
}

vectface my_mpi_allgather(boost::mpi::communicator &comm, vectface const &vf) {
  size_t n_vert = vf.get_n();
  size_t n_face = vf.size();
  std::vector<uint8_t> const &V = vf.serial_get_std_vector_uint8_t();
  //
  std::vector<size_t> l_n_face = my_mpi_allgather(comm, n_face);
  std::vector<std::vector<uint8_t>> l_V = my_mpi_allgather(comm, V);
  return vectface(n_vert, l_n_face, l_V);
}

vectface merge_initial_samp(boost::mpi::communicator &comm, vectface const &vf,
                            std::string const &ansSamp,
                            [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: merge_initial_samp, |vf|=" << vf.size() << "\n";
#endif
  vectface vf_gather = my_mpi_allgather(comm, vf);
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: merge_initial_samp, |vf_gather|=" << vf_gather.size() << "\n";
#endif
  std::vector<std::string> ListStr = STRING_Split(ansSamp, ":");
  std::string ansOpt = ListStr[0];
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: ansSamp=" << ansSamp << " ansOpt=" << ansOpt << "\n";
#endif
  if (ansOpt == "lp_cdd_min") {
    vectface vf_ret = select_minimum_count(vf_gather);
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT: merge_initial_samp, |vf_ret|=" << vf_ret.size() << "\n";
#endif
    Face f_first = vf_ret[0];
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT: f_first.count=" << f_first.count() << "\n";
#endif
    return vf_ret;
  }
  return vf_gather;
}

/*
  For MyMatrix<T>, a little bit more advanced work. We merge the rows together
 */

template <typename T>
MyMatrix<T> MergeRows(int const &n_col, std::vector<int> const &l_n_rows,
                      std::vector<std::vector<T>> const &l_V) {
  int n_row = 0;
  for (auto &eVal : l_n_rows)
    n_row += eVal;
  MyMatrix<T> Mret(n_row, n_col);
  int pos = 0;
  for (size_t i = 0; i < l_n_rows.size(); i++) {
    int const &n_rows = l_n_rows[i];
    std::vector<T> const &V = l_V[i];
    size_t posV = 0;
    for (int i_row = 0; i_row < n_rows; i_row++) {
      for (int i_col = 0; i_col < n_col; i_col++) {
        Mret(pos, i_col) = V[posV];
        posV++;
      }
      pos++;
    }
  }
  return Mret;
}

template <typename T> std::vector<T> MatrixRowsAsVector(MyMatrix<T> const &M) {
  int n_rows = M.rows();
  int n_cols = M.cols();
  std::vector<T> V(n_rows * n_cols);
  size_t pos = 0;
  for (int i_row = 0; i_row < n_rows; i_row++) {
    for (int i_col = 0; i_col < n_cols; i_col++) {
      V[pos] = M(i_row, i_col);
      pos++;
    }
  }
  return V;
}

template <typename T>
MyMatrix<T> my_mpi_gather(boost::mpi::communicator &comm, MyMatrix<T> const &M,
                          int i_proc) {
  int i_rank = comm.rank();
  int n_rows = M.rows();
  int n_cols = M.cols();
  std::vector<T> V = MatrixRowsAsVector(M);
  //
  std::vector<int> l_n_rows = my_mpi_gather(comm, n_rows, i_proc);
  std::vector<std::vector<T>> l_V = my_mpi_gather(comm, V, i_proc);
  if (i_rank == i_proc) {
    return MergeRows(n_cols, l_n_rows, l_V);
  } else {
    return MyMatrix<T>(0, n_cols);
  }
}

template <typename T>
MyMatrix<T> my_mpi_allgather(boost::mpi::communicator &comm,
                             MyMatrix<T> const &M) {
  int n_rows = M.rows();
  int n_cols = M.cols();
  std::vector<T> V = MatrixRowsAsVector(M);
  std::vector<int> l_n_rows = my_mpi_allgather(comm, n_rows);
  std::vector<std::vector<T>> l_V = my_mpi_allgather(comm, V);
  return MergeRows(n_cols, l_n_rows, l_V);
}

/*
  For T, compute the sum of all the elements
 */

template <typename T>
T my_mpi_allreduce_sum(boost::mpi::communicator &comm, T const &x) {
  std::vector<T> V;
  boost::mpi::all_gather<T>(comm, x, V);
  T sum = V[0];
  for (size_t i = 1; i < V.size(); i++)
    sum += V[i];
  return sum;
}

template <typename T, typename Tgroup>
bool EvaluationConnectednessCriterion_KernelMPI_field(
    boost::mpi::communicator &comm, MyMatrix<T> const &FAC, Tgroup const &GRP,
    vectface const &vf_undone_loc, std::ostream &os) {
  vectface vf_undone_tot = my_mpi_allgather(comm, vf_undone_loc);
  MyMatrix<T> EXT_undone_loc = GetVertexSet_from_vectface(FAC, vf_undone_loc);
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 4.1\n";
#endif
  MyMatrix<T> EXT_undone_tot = my_mpi_allgather(comm, EXT_undone_loc);
  // Every processor is computing the adjacency stuff. Fairly inefficient but ok
  // for the time being.
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 5\n";
#endif
  size_t max_iter = 100;
  size_t n_iter = 0;
  auto f_recur = [&](const std::pair<size_t, Face> &pfr) -> bool {
    n_iter++;
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT:   f_recur n_iter=" << n_iter << "\n";
#endif
    if (n_iter == max_iter)
      return false;
    if (pfr.first > 1)
      return false;
    return true;
  };
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 6\n";
#endif
  return EvaluationConnectednessCriterion_Kernel(FAC, GRP, EXT_undone_tot,
                                                 vf_undone_tot, f_recur, os);
}

template <typename T, typename Tgroup>
inline typename std::enable_if<is_ring_field<T>::value, bool>::type
EvaluationConnectednessCriterion_KernelMPI(boost::mpi::communicator &comm,
                                           MyMatrix<T> const &FAC,
                                           Tgroup const &GRP,
                                           vectface const &vf_undone_loc,
                                           std::ostream &os) {
  return EvaluationConnectednessCriterion_KernelMPI_field(comm, FAC, GRP,
                                                          vf_undone_loc, os);
}

template <typename T, typename Tgroup>
inline typename std::enable_if<!is_ring_field<T>::value, bool>::type
EvaluationConnectednessCriterion_KernelMPI(boost::mpi::communicator &comm,
                                           MyMatrix<T> const &FAC,
                                           Tgroup const &GRP,
                                           vectface const &vf_undone_loc,
                                           std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> FACfield = UniversalMatrixConversion<Tfield, T>(FAC);
  return EvaluationConnectednessCriterion_KernelMPI_field(comm, FACfield, GRP,
                                                          vf_undone_loc, os);
}

template <typename TbasicBank>
bool EvaluationConnectednessCriterion_MPI(boost::mpi::communicator &comm,
                                          TbasicBank const &bb,
                                          std::ostream &os) {
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 1\n";
#endif
  using Tint = typename TbasicBank::Tint;
  // We need an heuristic to avoid building too large orbits.
  // A better system would have to balance out the cost of
  // doing that check with respect to the dual description itsef.
  //
  // We need to terminate the check if no orbit has been done.
  size_t nbOrbitDone_tot = 0, nbOrbitDone_loc = bb.foc.nbOrbitDone;
  all_reduce(comm, nbOrbitDone_loc, nbOrbitDone_tot,
             boost::mpi::maximum<size_t>());
  if (nbOrbitDone_tot == 0)
    return false;
  // In order for the check not to be too expensive, we limit ourselves to 1000
  Tint max_siz = 1000;
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 2 nbUndone="
     << bb.foc.nbUndone << "\n";
#endif
  Tint nbUndone_tot = my_mpi_allreduce_sum(comm, bb.foc.nbUndone);
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: nbUndone_tot=" << nbUndone_tot << "\n";
#endif
  if (nbUndone_tot > max_siz)
    return false;
#ifdef DEBUG_MPI_FCT
  os << "MPI_FCT: EvaluationConnectednessCriterion_MPI, step 3\n";
#endif
  vectface vf_undone_loc = ComputeSetUndone(bb);
  return EvaluationConnectednessCriterion_KernelMPI(comm, bb.EXT, bb.GRP,
                                                    vf_undone_loc, os);
}

/*
  The asynchronous exchanges require keeping track of when something has been
  sent and some exchange is done.
  We allow for:
  --- A number of MaxFly messages in fly which can be larger than n_proc.
  --- Two choices:
     --- n_proc == MaxFly : Only one message from A to B at the same time.
     --- MaxFly > n_proc : We can potentially have all the messages to node 0
  --- Methods:
     --- Access to the request
     --- clear_and_get_nb_done : clear the messages and return how much is left.
     --- GetFreeIndex : return an index for operation.
     --- is_empty : true if no pending messages.
 */
struct request_status_list {
  size_t n_proc;
  size_t MaxFly;
  std::vector<boost::mpi::request> l_mpi_request;
  std::vector<int> l_mpi_status;
  bool strict;
  request_status_list(size_t const &n_proc, size_t const &MaxFly)
      : n_proc(n_proc), MaxFly(MaxFly), l_mpi_request(MaxFly),
        l_mpi_status(MaxFly, 0) {
    if (n_proc == MaxFly) {
      strict = true;
    } else {
      strict = false;
    }
  }
  boost::mpi::request &operator[](size_t const &pos) {
    l_mpi_status[pos] = 1;
    return l_mpi_request[pos];
  }
  size_t clear_and_get_nb_undone() {
    size_t len = l_mpi_request.size();
    size_t nb_undone = 0;
    for (size_t i = 0; i < len; i++) {
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
  bool is_ok(size_t const &pos) {
    if (l_mpi_status[pos] == 0)
      return true;
    return false;
  }
  size_t GetFreeIndex(int dest) {
    (void)clear_and_get_nb_undone();
    if (strict) {
      if (is_ok(dest)) {
        return dest;
      } else {
        return std::numeric_limits<size_t>::max();
      }
    } else {
      size_t len = l_mpi_request.size();
      for (size_t i = 0; i < len; i++)
        if (l_mpi_status[i] == 0)
          return i;
      if (MaxFly > 0)
        return std::numeric_limits<size_t>::max();
      l_mpi_request.push_back(boost::mpi::request());
      l_mpi_status.push_back(1);
      return len;
    }
  }
  bool is_empty() {
    size_t n_undone = clear_and_get_nb_undone();
    return n_undone == 0;
  }
};

struct empty_message_management {
  boost::mpi::communicator &comm;
  request_status_list rsl;
  // expected_value is random, but identical on all process.
  int expected_value;
  int tag;
  empty_message_management(boost::mpi::communicator &comm, size_t const &MaxFly,
                           int const &tag)
      : comm(comm), rsl(comm.size(), MaxFly), tag(tag) {
    int expected_value_pre = random();
    expected_value = boost::mpi::all_reduce(comm, expected_value_pre,
                                            boost::mpi::minimum<int>());
  }
  void send_message(int dest) {
    size_t idx = rsl.GetFreeIndex(dest);
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

/* The operations is done via
   comm.isend(dest, tag, val);
   with val a reference. This means that data has to be preplaced before being
   sent.
   ---
   So, we have a number of arrays for sending data by mpi (the l_under_cons).
   And we have a number of pending messages to insert: (the l_messages)
   ---
   We have tzo kinds of template parameters
   * The T
   * The T_vector.
   Basically, the T_vector must have a .push_back function available to push
   entries. Eamplex:
   *  T = int, T_vector = std::vector<T>
   * T = Face, T_vector = vectface
 */
template <typename T, typename T_vector> struct buffered_T_exchanges {
  boost::mpi::communicator &comm;
  int n_proc;
  request_status_list rsl;
  int tag;
  std::vector<T_vector> l_message;
  std::vector<T_vector> l_under_cons;
  bool strict;
  buffered_T_exchanges(boost::mpi::communicator &comm, size_t const &MaxFly,
                       int const &tag)
      : comm(comm), n_proc(comm.size()), rsl(n_proc, MaxFly), tag(tag),
        l_message(n_proc), l_under_cons(MaxFly) {
    if (static_cast<size_t>(n_proc) == MaxFly) {
      strict = true;
    } else {
      strict = false;
    }
  }
  void insert_entry(size_t const &pos, T const &x) {
    l_message[pos].push_back(x);
  }
  T_vector recv_message(int source) {
    T_vector vf;
    comm.recv(source, tag, vf);
    return vf;
  }
  bool clear_one_entry([[maybe_unused]] std::ostream &os) {
    auto process = [&](int chosen_iproc, int idx) -> bool {
      if (chosen_iproc == -1)
        return false;
      l_under_cons[idx] = std::move(l_message[chosen_iproc]);
      l_message[chosen_iproc].clear();
      rsl[idx] = comm.isend(chosen_iproc, tag, l_under_cons[idx]);
      return true;
    };
    if (strict) {
#ifdef DEBUG_MPI_FCT
      size_t nb_undone = rsl.clear_and_get_nb_undone();
      os << "MPI_FCT: strict - clear_one_entry - nb_undone=" << nb_undone
         << "\n";
#else
      (void)rsl.clear_and_get_nb_undone();
#endif
      size_t max_siz = 0;
      int chosen_iproc = -1;
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        bool test = rsl.is_ok(i_proc);
#ifdef DEBUG_MPI_FCT
        os << "MPI_FCT: i_proc=" << i_proc << " test=" << test << "\n";
#endif
        if (test) {
          size_t siz = l_message[i_proc].size();
#ifdef DEBUG_MPI_FCT
          os << "MPI_FCT:  siz=" << siz << "\n";
#endif
          if (siz > max_siz) {
            max_siz = siz;
            chosen_iproc = i_proc;
          }
        }
      }
#ifdef DEBUG_MPI_FCT
      os << "MPI_FCT: strict: max_siz=" << max_siz
         << " chosen_iproc=" << chosen_iproc << "\n";
#endif
      return process(chosen_iproc, chosen_iproc);
    } else {
      size_t idx = rsl.GetFreeIndex(-1); // The dest is unused in that case
#ifdef DEBUG_MPI_FCT
      os << "MPI_FCT: idx=" << idx << "\n";
#endif
      if (idx == std::numeric_limits<size_t>::max())
        return false;
      size_t max_siz = 0;
      int chosen_iproc = -1;
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        size_t siz = l_message[i_proc].size();
        if (siz > max_siz) {
          max_siz = siz;
          chosen_iproc = i_proc;
        }
      }
#ifdef DEBUG_MPI_FCT
      os << "MPI_FCT: no_strict: max_siz=" << max_siz
         << " chosen_iproc=" << chosen_iproc << "\n";
#endif
      return process(chosen_iproc, idx);
    }
  }
  size_t get_unsent_size() const {
    size_t n_unsent = 0;
    for (auto &eEnt : l_message)
      n_unsent += eEnt.size();
    return n_unsent;
  }
  bool is_buffer_empty() {
    //    if (!rsl.is_empty())
    //      return false;
    return get_unsent_size() == 0;
  }
  bool is_completely_clear() {
    for (auto &eEnt : l_message) {
      if (eEnt.size() > 0) {
        return false;
      }
    }
    return rsl.is_empty();
  }
};

template <typename Tint> struct database_balinski_info {
  boost::mpi::communicator &comm;
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
  database_balinski_info(boost::mpi::communicator &comm, int tag_request,
                         int tag_info, Tint const &CritSiz)
      : comm(comm), tag_request(tag_request), tag_info(tag_info),
        n_proc(comm.size()), i_rank(comm.rank()), rsl(n_proc, n_proc),
        ListBalinski(n_proc), ListStatus_Emptyness(n_proc, 0),
        RequestNat(n_proc, 0), CritSiz(CritSiz), expected_value(47),
        last_update(n_proc, get_cpp_time(7, 1, 1974)),
        last_database_update_time(get_cpp_time(6, 1, 1974)),
        ToBeAnswered(n_proc) {}
  bool get_status([[maybe_unused]] std::ostream &os) const {
    for (int i_proc = 0; i_proc < n_proc; i_proc++)
      if (ListStatus_Emptyness[i_proc] == 0) {
#ifdef DEBUG_MPI_FCT
        int val = ListStatus_Emptyness[i_proc];
        os << "MPI_FCT: Returning false at i_proc=" << i_proc << " val=" << val
           << "\n";
#endif
        return false;
      }
    UndoneOrbitInfo<Tint> uoi = CombineUndoneOrbitInfo(ListBalinski);
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT: Merged uoi=" << uoi << "\n";
#endif
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
  void submit_info(int dest, UndoneOrbitInfo<Tint> const &uoi) {
    size_t idx = rsl.GetFreeIndex(dest);
    if (idx == std::numeric_limits<size_t>::max())
      return;
    rsl[idx] = comm.isend(dest, tag_info, uoi);
    last_update[dest] = last_database_update_time;
    ToBeAnswered[dest] = 0;
  }
  void submit_request(int dest) {
    size_t idx = rsl.GetFreeIndex(dest);
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
  void set_uoi_local(UndoneOrbitInfo<Tint> const &uoi_local) {
    ListBalinski[i_rank] = uoi_local;
    ListStatus_Emptyness[i_rank] = 1;
    last_database_update_time = std::chrono::system_clock::now();
  }
  void flush([[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT: Beginning of flush operation\n";
#endif
    UndoneOrbitInfo<Tint> const &uoi_local = ListBalinski[i_rank];
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      if (i_proc != i_rank) {
        bool test = last_database_update_time > last_update[i_proc];
#ifdef DEBUG_MPI_FCT
        int dur_i = std::chrono::duration_cast<std::chrono::seconds>(
                        last_database_update_time - last_update[i_proc])
                        .count();
        os << "MPI_FCT: i_proc=" << i_proc
           << " ToBeAnswered=" << ToBeAnswered[i_proc] << " test=" << test
           << " dur_i=" << dur_i << "\n";
#endif
        if (ToBeAnswered[i_proc] == 1 && test) {
#ifdef DEBUG_MPI_FCT
          os << "MPI_FCT: flush: Doing submit_info for i_proc=" << i_proc
             << "\n";
#endif
          submit_info(i_proc, uoi_local);
        }
      }
    }
  }
  void submit_request_uoi([[maybe_unused]] std::ostream &os) {
    // First checking natively
#ifdef DEBUG_MPI_FCT
    os << "MPI_FCT: Beginning of submit_uoi\n";
#endif
    // First checking for unassigned
    for (int i_proc = 0; i_proc < n_proc; i_proc++)
      if (i_proc != i_rank && ListStatus_Emptyness[i_proc] == 0 &&
          RequestNat[i_proc] == 0) {
#ifdef DEBUG_MPI_FCT
        os << "MPI_FCT: submit_request_uoi (Case 0) at i_proc=" << i_proc
           << "\n";
#endif
        return submit_request(i_proc);
      }
    // Getting the highest value
    Tint max_val = 0;
    int chosen_idx = -1;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      if (i_proc != i_rank) {
        Tint const &eval = ListBalinski[i_proc].nbUndone;
        if (eval > max_val) {
          max_val = eval;
          chosen_idx = i_proc;
        }
      }
    }
    if (chosen_idx == -1)
      return;
    if (RequestNat[chosen_idx] == 0) {
#ifdef DEBUG_MPI_FCT
      os << "MPI_FCT: submit_request_uoi (Case 1) at chosen_idx=" << chosen_idx
         << "\n";
#endif
      return submit_request(chosen_idx);
    }
  }
};

/*
  We want to handle a potentially unlimited number of pending requests.
  This forces having a std::list since a std::vector can deallocate and
  reallocate while being resized. This cannot happen for a std::list. Of course,
  there is never any deallocation as that would lead to runtime bugs
 */
struct unlimited_request {
  boost::mpi::communicator &comm;
  int n_proc;
  // First is false if not active.
  using Tpair = std::pair<bool, boost::mpi::request>;
  std::list<std::vector<Tpair>> tot_list;
  unlimited_request(boost::mpi::communicator &comm)
      : comm(comm), n_proc(comm.size()) {}
  size_t clear_and_get_nb_undone() {
    size_t n_undone = 0;
    for (auto &eList : tot_list) {
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        Tpair &ePair = eList[i_proc];
        if (ePair.first) {
          boost::optional<boost::mpi::status> stat = ePair.second.test();
          if (stat) {
            ePair.first = false;
          } else {
            n_undone++;
          }
        }
      }
    }
    return n_undone;
  }
  bool is_empty() { return clear_and_get_nb_undone() > 0; }
  boost::mpi::request &get_entry() {
    (void)clear_and_get_nb_undone();
    for (auto &eList : tot_list) {
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        Tpair &ePair = eList[i_proc];
        if (!ePair.first) {
          ePair.first = true;
          return ePair.second;
        }
      }
    }
    // That case should be rare indeed (and it is slow)
    std::vector<Tpair> V(n_proc);
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      V[i_proc].first = false;
    }
    tot_list.push_back(V);
    size_t siz = tot_list.size();
    auto iter = tot_list.begin();
    for (size_t u = 1; u < siz; u++) {
      iter++;
    }
    std::vector<Tpair> &eList = *iter;
    Tpair &ePair = eList[0];
    ePair.first = true;
    return ePair.second;
  }
};

bool ApplyStdUnitbuf(FullNamelist const &eFull) {
  SingleBlock BlockSYSTEM = eFull.get_block("SYSTEM");
  bool result = BlockSYSTEM.get_bool("ApplyStdUnitbuf");
  return result;
}

std::unique_ptr<std::ofstream>
get_mpi_log_stream(boost::mpi::communicator &comm, FullNamelist const &eFull) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string FileLog =
      "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::unique_ptr<std::ofstream> os = std::make_unique<std::ofstream>(FileLog);
  if (ApplyStdUnitbuf(eFull)) {
    *os << std::unitbuf;
    *os << "Apply UnitBuf\n";
  } else {
    *os << "Do not apply UnitBuf\n";
  }
  return os;
}

// clang-format off
#endif  // SRC_DUALDESC_MPI_FUNCTIONALITY_H_
// clang-format on
