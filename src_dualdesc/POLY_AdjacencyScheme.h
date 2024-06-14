// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_

// clang-format off
#include "boost_serialization.h"
#include "MPI_functionality.h"
#include "basic_datafile.h"
#include "Timings.h"
#include <unordered_map>
#include <utility>
#include <vector>
#include <set>
#include <functional>
// clang-format on

#ifdef DEBUG
#define DEBUG_ADJACENCY_SCHEME
#endif

#ifdef TIMINGS
#define TIMINGS_ADJACENCY_SCHEME
#endif

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
  ---The L-types of a T-space.
     ---Basic setting as for perfect forms.
  ---The perfect form of the hyperbolic cone.
     ---This might work similarly to the Delaunay. Right now we have
  ---The C-types of Z^n (but we do not have realistic examples)
     ---We have canonicalization stuff.

  Conclusion: The canonicalization cannot be assumed to exist in general.
  But we could work out with f_repr. And if by accident we have some
  canonicalization then we can have a corresponding hash function and then an
  equality test. So, there is no loss of generality in working that way.

  TODO:
  ---Resolve the scheme so that we can conclude when everything has been treated
  without any heuristic.
  ---When a node has everything processed, send a note to all others asking for
  their status. If all of them reply then things are clear.
  ---Have to learn about "distributed programming" from Leslie Lamport and
  others. Takada_M._Distributed_systems_for_fun_and_profit.epub (short one)
  Varela_C.A._Programming_distributed_computing_systems__a_foundational_approach.pdf
  (ok, teaches pi-calculus, join, etc.)
  Raynal_M._Distributed_Algorithms_for_Message-Passing_Systems.pdf
  (especially relevant to us)
  Fokkink_W._Distributed_Algorithms_-_An_Intuitive_Approach.pdf

  ---The Fokkink book provides a reasonable introduction to the question.


  Possible design:
  ---We can have two outcomes from a sending of data:
     ---For unlimited runtime: The termination is that all the entries in the
  vector<Tobj> have their adjacencies being computed.
     ---For specified runtime: The termination is that all the entries that have
  been sent have been inserted into the database.
  ---This is the definition of passive/active:
     ---It ensures that no message in flight is left.
     ---It ensures that no computation is finished.
  ---But that does not allow to determine that everything is finished.
  ---We need to keep a track of the number of operations done.
  ---So, we send a message to all the nodes and collect all the operation done.
  ---Then we call again and see if the nonce are the same.
  ---But the problem is that sending all the entries can introduce some
  deadlocks maybe.

 */

/*
  Possible extension to the scheme.
  ---Since in the scheme, we compute everything, if the first entry returned by
  the f_init is a big monster like the biggest entry in the list, then we are
  going to spend all the time on that one and the other nodes are not going to
  do anything.
  ---This is damaging since there are many well agreed use cases where we want
  to do a partial enumeration for this and that reason (mostly curiosity).
  ---So, there are some realistic scenario where we want to do the enumeration
  from the easiest to the hardest.

  How this can be implemented:
  ---We would simply need to have a function f_complexity that returns an size_t
  that encodes the complexity of the object considered.
  ---The task of having a good initial guess would be incumbent on the asker. A
  good f_init would be needed.
  ---(For example a good f_init could be obtained by an initial random walk).
  ---Therefore, we would have a std::map<size_t, ....> encoding the unused
  entries by complexity.
  ---Each node should have a value of first undone entry. We should also have
  something global.
  ---If node A and B have level X but B passes to X+1, in order to upgrade the
  global value, we need to know that there is no other node at X. So this
  forces having a std::vector<size_t> first_undone(m_proc)
  ---Any change of a node own level has to be followed by the emission to all
  nodes of the new level.
  ---
 */

template <typename T> void append_move(std::vector<T> &v1, std::vector<T> &v2) {
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()),
            std::make_move_iterator(v2.end()));
  v2.clear();
}

const size_t seed_partition = 10;
const size_t seed_hashmap = 20;

template <typename TadjI> struct entryAdjI {
  TadjI x;
  size_t hash_hashmap;
  int i_proc_orig;
  int i_orb_orig;
};

namespace boost::serialization {
template <class Archive, typename TadjI>
inline void serialize(Archive &ar, entryAdjI<TadjI> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("hash_hashmap", eRec.hash_hashmap);
  ar &make_nvp("i_proc_orig", eRec.i_proc_orig);
  ar &make_nvp("i_orb_orig", eRec.i_orb_orig);
}
} // namespace boost::serialization

template <typename TadjO> struct entryAdjO {
  TadjO x;
  int i_orb_orig;
};

namespace boost::serialization {
template <class Archive, typename TadjO>
inline void serialize(Archive &ar, entryAdjO<TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("i_orb_orig", eRec.i_orb_orig);
}
} // namespace boost::serialization

/*
  input clear from preceding discussion.
  The returned boolean is:
  true: if the enumeration finished. false otherwise.
  ---
  Types (cannot be set up via concept, but so be it)
  Tobj: The object type being created (like L-type domain)
  TadjI: The types returned by the spanning. Can be something like
    std::pair<Face, Tobj> with Face indicating the relevant Face.
  TadjO: The adjacency type after processing. Can be something like
    {Face, Trans, idx} with Face the corresponding Face, Trans the
    transformation realizing the equivalence and idx the equivalent
    object.
  All objects have to be transmitible via boost C++ mpi.
  ---
  Function types used:
  f_next() -> std::optional<std::pair<bool,Tobj>>: returns the entry
     to insert into the database together with its status. If a none
     if returned, then nothing else needs to be inserted.
  f_insert(Tobj) -> bool : insert the new object.
     If return true then early termination is triggered.
  f_adji_obj(TadjI) -> Tobj: should return the created object from the
     input adjacency.
  f_idx_obj(size_t) -> Tobj: should return the object from the database.
  f_save_status(int, bool) -> void : save the status in the database
  f_init() -> Tobj : get a starting element
  f_adj(i_orb) -> std::vector<TadjI> : get the adjacent object
     with i_orb the index of the orbit (used to assign for example the group)
  f_set_adj(int, std::vector<TadjO>) -> void : set the adjacencies to the ones
     computed.
  f_hash(size_t, Tobj) -> size_t : compute the hash from a specified seed.
  f_repr(Tobj, TadjI, int, int) -> std::optional<TadjO> : returns whether
    Tobj is equivalent to the spanned TadjI and find the equivalence if that
    is the case. Also take
  f_spann(TadjI, int, int) -> std::pair<Tobj, TabjO> : Generate from the
    equivalence the object to be inserted and the equivalence to be sent.
  ---
  The exchanges that we have:
  * Sending of entriesAdjI
  * Sending of entriesAdjO
  * Sending of nonce
  * termination
  ---
  The design for termination is the following:
  * The termination is done by checking the nonce of all the nodes.
  This is a blocking operation on the checker side, so only one
  can do it.
  * We want to avoid the busy loop of
  * When a node has nothing to do, if it is 0.

 */
template <typename Tobj, typename TadjI, typename TadjO, typename Fnext,
          typename Finsert, typename Fadji_obj, typename Fidx_obj,
          typename Fsave_status, typename Finit, typename Fadj,
          typename Fset_adj, typename Fhash, typename Frepr, typename Fspann>
bool compute_adjacency_mpi(boost::mpi::communicator &comm,
                           int const &max_time_second, Fnext f_next,
                           Finsert f_insert, Fadji_obj f_adji_obj,
                           Fidx_obj f_idx_obj, Fsave_status f_save_status,
                           Finit f_init, Fadj f_adj, Fset_adj f_set_adj,
                           Fhash f_hash, Frepr f_repr, Fspann f_spann,
                           [[maybe_unused]] std::ostream &os) {
  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  const int tag_initial = 35;
  const int tag_nonce_ask = 36;
  const int tag_nonce_reply = 37;
  const int tag_entriesadji_send = 38;
  const int tag_entriesadjo_send = 39;
  const int tag_termination = 40;
  const int tag_early_termination = 41;
  const int tag_no_operation = 42;
  unlimited_request ur(comm);
  //
  // The data sets
  //
  int n_obj = 0; // The number of objects generated
  // The map from the hash to the list of indices
  std::unordered_map<size_t, std::vector<size_t>> indices_by_hash;
  std::vector<size_t> undone; // The undone indices
  //
  // The entries of AdjI / AdjO
  //
  // The unsent entries by the processors.
  buffered_T_exchanges<entryAdjI<TadjI>, std::vector<entryAdjI<TadjI>>>
      buffer_entriesAdjI(comm, n_proc, tag_entriesadji_send);
  buffered_T_exchanges<entryAdjO<TadjO>, std::vector<entryAdjO<TadjO>>>
      buffer_entriesAdjO(comm, n_proc, tag_entriesadjo_send);
  std::vector<entryAdjI<TadjI>> unproc_entriesAdjI;
  // The mapping from the index to the list of adjacencices.
  std::unordered_map<int, std::pair<size_t, std::vector<TadjO>>> map_adjO;

  // The nonce is used so that a number is associated to a specific computation.
  // The function get_nonce returns 0 if there is something left to do and
  // nonzero if there is no pending computation. If the get_nonce returns
  // the same, then there was no computation going on.
  // The none is also used for the tracking of the adjacencies entries.
  size_t nonce = 1;
  // We have two ways of terminating a computation.
  // * One is by runtime elapsing. When this is reached, the processes need to
  //   finish the processing of the entries that are being treated.
  // * Another is by the early_termination. In that case, disorderly exit is
  // fine.
  //   We got something wrong, no point in continuing.
  // They behave in different ways
  bool early_termination = false;
  //
  // The lambda functions
  //
  auto send_early_termination = [&]() -> void {
    early_termination = true;
    int val_final = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      if (i_proc != i_rank) {
        ur.get_entry() = comm.isend(i_proc, tag_early_termination, val_final);
      }
    }
  };
  auto process_single_entryAdjI =
      [&](entryAdjI<TadjI> const &eI) -> std::pair<int, entryAdjO<TadjO>> {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Beginning of process_single_entryAdjI\n";
#endif
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    auto f_ret = [&](TadjO u) -> std::pair<int, entryAdjO<TadjO>> {
      entryAdjO<TadjO> eA{u, eI.i_orb_orig};
      return {eI.i_proc_orig, eA};
    };
    nonce++;
    std::vector<size_t> &vect = indices_by_hash[eI.hash_hashmap];
    for (auto &idx : vect) {
      Tobj x = f_idx_obj(idx);
#ifdef TIMINGS_ADJACENCY_SCHEME
      MicrosecondTime time_f_repr;
#endif
      std::optional<TadjO> opt = f_repr(x, eI.x, i_rank, idx);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |f_repr|=" << time_f_repr << "\n";
#endif
      if (opt) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: Conclude with an equivalence\n";
#endif
        return f_ret(*opt);
      }
    }
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |loop vect|=" << time << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Conclude with a new one\n";
#endif
    std::pair<Tobj, TadjO> pair = f_spann(eI.x, i_rank, n_obj);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |f_spann|=" << time << "\n";
#endif
    bool test = f_insert(pair.first);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |f_insert|=" << time << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_insert test=" << test << "\n";
#endif
    if (test) {
      send_early_termination();
    }
    undone.push_back(n_obj);
    vect.push_back(n_obj);
    f_save_status(n_obj, false);
    n_obj++;
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_single_entryAdjI, now n_obj=" << n_obj << "\n";
#endif
    return f_ret(pair.second);
  };
  auto insert_load = [&](Tobj const &x, bool const &is_treated) -> void {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    size_t hash_hashmap = f_hash(seed_hashmap, x);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |f_hash|=" << time << "\n";
#endif
    std::vector<size_t> &vect = indices_by_hash[hash_hashmap];
    vect.push_back(n_obj);
    if (!is_treated) {
      undone.push_back(n_obj);
    }
    n_obj++;
  };
  auto initial_insert = [&](Tobj const &x) -> void {
    nonce++;
    size_t hash_partition = f_hash(seed_partition, x);
    int i_proc_belong = static_cast<int>(hash_partition % size_t(n_proc));
    if (i_proc_belong != i_rank) {
      comm.send(i_proc_belong, tag_initial, x);
    } else {
      bool test = f_insert(x);
      if (test) {
        send_early_termination();
      }
      bool is_treated = false;
      f_save_status(n_obj, is_treated);
      insert_load(x, is_treated);
    }
  };
  auto process_entriesAdjO = [&](std::vector<entryAdjO<TadjO>> &v) -> void {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH:process_entriesAdjO : Before |v|=" << v.size()
       << " |map_adjO|=" << map_adjO.size() << " \n";
    for (auto &kv : map_adjO) {
      os << "ADJ_SCH:process_entriesAdjO : Before i_orb=" << kv.first
         << " siz=" << kv.second.first << " |l_adj|=" << kv.second.second.size()
         << "\n";
    }
#endif
    for (auto &eEnt : v) {
      map_adjO[eEnt.i_orb_orig].second.emplace_back(std::move(eEnt.x));
    }
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH:process_entriesAdjO : After |map_adjO|=" << map_adjO.size()
       << " \n";
    for (auto &kv : map_adjO) {
      os << "ADJ_SCH:process_entriesAdjO : After i_orb=" << kv.first
         << " siz=" << kv.second.first << " |l_adj|=" << kv.second.second.size()
         << "\n";
    }
#endif
    v.clear();
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |process_entriesAdjO|=" << time << "\n";
#endif
  };
  auto get_nonce = [&]() -> size_t {
    if (!buffer_entriesAdjI.is_completely_clear()) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because "
            "!buffer_entriesAdjI.is_completely_clear()\n";
#endif
      return 0;
    }
    if (!buffer_entriesAdjO.is_completely_clear()) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because "
            "!buffer_entriesAdjO.is_completely_clear()\n";
#endif
      return 0;
    }
    if (unproc_entriesAdjI.size() > 0) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because unproc_entriesAdjI.size() > 0\n";
#endif
      return 0;
    }
    if (map_adjO.size() > 0) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because map_adjO.size() > 0\n";
#endif
      return 0;
    }
    if (undone.size() > 0 && max_time_second == 0) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because undone.size() > 0 && "
            "max_time_second == 0\n";
#endif
      return 0;
    }
    if (undone.size() > 0 && max_time_second > 0 &&
        si(start) < max_time_second) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: max_time_second=" << max_time_second
         << " si(start)=" << si(start) << " |undone|=" << undone.size() << "\n";
      os << "ADJ_SCH: the_nonce=0, because max_time_second > 0 && si(start) < "
            "max_time_second\n";
#endif
      return 0;
    }
    return nonce;
  };
  auto process_mpi_status = [&](boost::mpi::status const &stat) -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_initial) {
      Tobj x;
      comm.recv(e_src, tag_initial, x);
      initial_insert(x);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 1|=" << time << "\n";
#endif
      return false;
    }
    if (e_tag == tag_entriesadji_send) {
      std::vector<entryAdjI<TadjI>> v;
      comm.recv(e_src, tag_entriesadji_send, v);
      append_move(unproc_entriesAdjI, v);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 2|=" << time << "\n";
#endif
      return false;
    }
    if (e_tag == tag_entriesadjo_send) {
      std::vector<entryAdjO<TadjO>> v;
      comm.recv(e_src, tag_entriesadjo_send, v);
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: Call process_entriesAdjO from process_mpi_status\n";
#endif
      process_entriesAdjO(v);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 3|=" << time << "\n";
#endif
      return false;
    }
    if (e_tag == tag_nonce_ask) {
      int val_recv;
      comm.recv(e_src, tag_nonce_ask, val_recv);
      size_t nonce = get_nonce();
      ur.get_entry() = comm.isend(e_src, tag_nonce_reply, nonce);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 4|=" << time << "\n";
#endif
      return false;
    }
    if (e_tag == tag_termination) {
      int val_recv;
      comm.recv(e_src, tag_termination, val_recv);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 5|=" << time << "\n";
#endif
      return true;
    }
    if (e_tag == tag_early_termination) {
      int val_recv;
      comm.recv(e_src, tag_termination, val_recv);
      early_termination = true;
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 6|=" << time << "\n";
#endif
      return true;
    }
    if (e_tag == tag_no_operation) {
      int val_recv;
      comm.recv(e_src, tag_no_operation, val_recv);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 7|=" << time << "\n";
#endif
      return false;
    }
    std::cerr << "The tag e_tag=" << e_tag << " is not matching\n";
    throw TerminalException{1};
  };
  auto get_undone_idx = [&]() -> size_t {
    size_t idx = undone[undone.size() - 1];
    undone.pop_back();
    return idx;
  };
  auto process_one_entry_obj = [&]() -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    if (max_time_second > 0 && si(start) > max_time_second) {
      return false;
    }
    if (undone.size() == 0) {
      return false;
    }
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: |undone|=" << undone.size() << " Before \n";
#endif
    size_t idx = get_undone_idx();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: |undone|=" << undone.size() << " After idx=" << idx << "\n";
#endif
    int idx_i = static_cast<int>(idx);
    f_save_status(idx, true);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |process_one_entry_obj f_save_status|=" << time << "\n";
#endif
    std::vector<TadjI> l_adj = f_adj(idx);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |process_one_entry_obj f_adj|=" << time << "\n";
#endif
    map_adjO[idx] = {l_adj.size(), {}};
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_one_entry_obj : idx=" << idx
       << " |l_adj|=" << l_adj.size() << "\n";
#endif
    nonce++;
    for (auto &x : l_adj) {
      Tobj x_obj = f_adji_obj(x);
      size_t hash_partition = f_hash(seed_partition, x_obj);
      size_t hash_hashmap = f_hash(seed_hashmap, x_obj);
      int i_proc_dest = static_cast<int>(hash_partition % size_t(n_proc));
      nonce++;
      entryAdjI<TadjI> e{x, hash_hashmap, i_rank, idx_i};
      buffer_entriesAdjI.insert_entry(i_proc_dest, e);
    }
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |process_one_entry_obj : dispatching|=" << time << "\n";
#endif
    return true;
  };
  auto flush_entriesAdjI = [&]() -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    std::vector<entryAdjI<TadjI>> &v = buffer_entriesAdjI.l_message[i_rank];
    append_move(unproc_entriesAdjI, v);
    bool test = buffer_entriesAdjI.clear_one_entry(os);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |flush_entriesAdjI|=" << time << "\n";
#endif
    return test;
  };
  auto flush_entriesAdjO = [&]() -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    std::vector<entryAdjO<TadjO>> &v = buffer_entriesAdjO.l_message[i_rank];
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Call process_entriesAdjO from flush_entriesAdjO\n";
#endif
    process_entriesAdjO(v);
    bool test = buffer_entriesAdjO.clear_one_entry(os);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |flush_entriesAdjO|=" << time << "\n";
#endif
    return test;
  };
  auto compute_entries_adjI = [&]() -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
    bool do_something = false;
    for (auto &eI : unproc_entriesAdjI) {
      std::pair<int, entryAdjO<TadjO>> pair = process_single_entryAdjI(eI);
      buffer_entriesAdjO.insert_entry(pair.first, pair.second);
      do_something = true;
    }
    unproc_entriesAdjI.clear();
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |compute_entries_adjI|=" << time << "\n";
#endif
    return do_something;
  };
  auto write_set_adj = [&]() -> bool {
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Begin of write_set_adj\n";
#endif
    std::vector<int> l_erase;
    bool do_something = false;
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: |map_adjO|=" << map_adjO.size() << "\n";
#endif
    for (auto &kv : map_adjO) {
      int i_orb = kv.first;
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: i_orb=" << i_orb << " kv.second.first=" << kv.second.first
         << " |kv.second.second|=" << kv.second.second.size() << "\n";
#endif
      if (kv.second.first == kv.second.second.size()) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: write_set_adj i_orb=" << i_orb
           << " |l_adj|=" << kv.second.second.size() << "\n";
#endif
        f_set_adj(i_orb, kv.second.second);
        l_erase.push_back(i_orb);
        do_something = true;
      }
    }
    for (auto &i_orb : l_erase) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: write_set_adj erasing i_orb=" << i_orb << "\n";
#endif
      map_adjO.erase(i_orb);
    }
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "ADJ_SCH: |write_set_adj|=" << time << "\n";
#endif
    return do_something;
  };
  // Returns true if something was actually done.
  auto f_clear_buffers = [&]() -> bool {
    // Transmitting the generated entriesAdjI
    bool test1 = flush_entriesAdjI();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffers : flush_entriesAdjI=" << test1 << "\n";
#endif
    if (test1) {
      return true;
    }
    // transmitting the generated entriesAdjO
    bool test2 = flush_entriesAdjO();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffers : flush_entriesAdjO=" << test2 << "\n";
#endif
    if (test2) {
      return true;
    }
    // compute the entryAdjO from entryAdjI
    bool test3 = compute_entries_adjI();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffers : compute_entries_adjI=" << test3 << "\n";
#endif
    if (test3) {
      return true;
    }
    // write down adjacencies if done
    bool test4 = write_set_adj();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffers : write_set_adj=" << test4 << "\n";
#endif
    if (test4) {
      return true;
    }
    // Nothing, then computing something new
    return process_one_entry_obj();
  };
#ifdef DEBUG_ADJACENCY_SCHEME
  auto print_status = [&](std::ostream &os_out) -> void {
    os_out << "ADJ_SCH: status n_obj=" << n_obj << "\n";
    os_out << "ADJ_SCH: status |map_adjO|=" << map_adjO.size() << "\n";
    for (auto &ent : map_adjO) {
      os_out << "ADJ_SCH: status index=" << ent.first << " value=("
             << ent.second.first << " | " << ent.second.second.size() << ")\n";
    }
    os_out << "ADJ_SCH: status |undone|=" << undone.size() << "\n";
    os_out << "ADJ_SCH: |buffer_entriesAdjI|="
           << buffer_entriesAdjI.get_unsent_size() << "\n";
    os_out << "ADJ_SCH: |buffer_entriesAdjO|="
           << buffer_entriesAdjO.get_unsent_size() << "\n";
  };
#endif
  int i_proc_termination = 0;
  auto f_consider_termination = [&]() -> void {
    int val_send = 42;
    ur.get_entry() = comm.isend(i_proc_termination, tag_no_operation, val_send);
  };
  auto f_terminate = [&]() -> bool {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: terminate, begin\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    if (i_rank > 0) {
      std::cerr << "The terminate cannot be called for a i_rank>0\n";
      std::cerr << "The reason is that the check for termination is a\n";
      std::cerr << "blocking check and so only one processor can do\n";
      std::cerr << "at a given time\n";
      throw TerminalException{1};
    }
#endif
    std::vector<size_t> l_nonce(n_proc - 1);
    for (int i_proc = 1; i_proc < n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_nonce_ask, val);
      size_t the_nonce;
      comm.recv(i_proc, tag_nonce_reply, the_nonce);
      if (the_nonce == 0) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: exiting at i_proc=" << i_proc << " for the_nonce=0\n";
#endif
        return false;
      }
      l_nonce[i_proc - 1] = the_nonce;
    }
    for (int i_proc = 1; i_proc < n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_nonce_ask, val);
      size_t the_nonce;
      comm.recv(i_proc, tag_nonce_reply, the_nonce);
      if (the_nonce != l_nonce[i_proc - 1]) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: exiting at i_proc=" << i_proc
           << " the_nonce=" << the_nonce << " l_nonce=" << l_nonce[i_proc - 1]
           << "\n";
#endif
        return false;
      }
    }
    // no operation done. Sending the termination notices
    for (int i_proc = 1; i_proc < n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_termination, val);
    }
    return true;
  };
  //
  // Loading the data
  //
  while (true) {
    std::optional<std::pair<bool, Tobj>> opt = f_next();
    if (opt) {
      std::pair<bool, Tobj> const &pair = *opt;
      Tobj const &x = pair.second;
      bool is_treated = pair.first;
      insert_load(x, is_treated);
    } else {
      break;
    }
  }
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: compute_adjacency_mpi, start : n_obj=" << n_obj
     << " n_undone=" << undone.size() << "\n";
#endif
  size_t n_orb_max = 0, n_orb_loc = n_obj;
  all_reduce(comm, n_orb_loc, n_orb_max, boost::mpi::maximum<size_t>());
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: beginning n_orb_max=" << n_orb_max
     << " n_orb_loc=" << n_orb_loc << "\n";
#endif
  if (n_orb_max == 0 && i_rank == 0) {
    Tobj x = f_init();
    initial_insert(x);
  }
  //
  // The infinite loop
  //
  while (true) {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Start while loop, early_termination=" << early_termination
       << "\n";
#endif
    if (early_termination) {
      break;
    }
    bool did_something;
    boost::optional<boost::mpi::status> opt = comm.iprobe();
    if (opt) {
      did_something = true;
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: prob is not empty\n";
#endif
      boost::mpi::status const &prob = *opt;
#ifdef TIMINGS_ADJACENCY_SCHEME
      MicrosecondTime time_process;
#endif
      bool test_process = process_mpi_status(prob);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 1|=" << time_process << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: process_mpi_status test_process=" << test_process << "\n";
#endif
      if (test_process) {
        break;
      }
    } else {
#ifdef TIMINGS_ADJACENCY_SCHEME
      MicrosecondTime time;
#endif
      did_something = f_clear_buffers();
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |f_clear_buffers|=" << time << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: f_clear_buffers did_something=" << did_something << "\n";
#endif
      if (!did_something) {
        if (i_rank != i_proc_termination) {
          f_consider_termination();
        } else {
#ifdef TIMINGS_ADJACENCY_SCHEME
          MicrosecondTime time_terminate;
#endif
          bool test_terminate = f_terminate();
#ifdef TIMINGS_ADJACENCY_SCHEME
          os << "ADJ_SCH: |f_terminate|=" << time_terminate << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
          os << "ADJ_SCH: test_terminate=" << test_terminate << "\n";
#endif
          if (test_terminate) {
            break;
          }
        }
      }
    }
    if (!did_something) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: Going to the blocking wait\n";
      os << "ADJ_SCH: time=" << timeanddate() << "\n";
      print_status(os);
#endif
      // Nothing was done, so we switch from iprobe to probe
#ifdef TIMINGS_ADJACENCY_SCHEME
      MicrosecondTime time_probe;
#endif
      boost::mpi::status prob = comm.probe();
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |probe|=" << time_probe << "\n";
#endif
      bool test = process_mpi_status(prob);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "ADJ_SCH: |process_mpi_status 2|=" << time_probe << "\n";
#endif
      if (test) {
        break;
      }
    }
  }
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: compute_adjacency_mpi, end : n_obj=" << n_obj
     << " n_undone=" << undone.size() << "\n";
#endif
  if (early_termination) {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: compute_adjacency_mpi, early_termination case\n";
#endif
    // We terminate early so by definition the enumeration is partial
    return false;
  } else {
    size_t n_undone_max = 0, n_undone_loc = undone.size();
    all_reduce(comm, n_undone_loc, n_undone_max, boost::mpi::maximum<size_t>());
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: compute_adjacency_mpi, n_undone_max=" << n_undone_max
       << " n_undone_loc=" << n_undone_loc << "\n";
#endif
    return n_undone_max == 0;
  }
}

template <typename Tobj, typename TadjI, typename TadjO, typename Fnext,
          typename Finsert, typename Fadji_obj, typename Fidx_obj,
          typename Fsave_status, typename Finit, typename Fadj,
          typename Fset_adj, typename Fhash, typename Frepr, typename Fspann>
bool compute_adjacency_serial(int const &max_time_second, Fnext f_next,
                              Finsert f_insert, Fadji_obj f_adji_obj,
                              Fidx_obj f_idx_obj, Fsave_status f_save_status,
                              Finit f_init, Fadj f_adj, Fset_adj f_set_adj,
                              Fhash f_hash, Frepr f_repr, Fspann f_spann,
                              [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: Beginning of compute_adjacency_serial\n";
#endif
  SingletonTime start;
  size_t n_obj = 0;
  std::unordered_map<size_t, std::vector<size_t>> indices_by_hash;
  std::vector<size_t> undone;
  bool early_termination = false;
  auto process_singleEntry_AdjI = [&](TadjI const &x_adjI) -> TadjO {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_singleEntry_AdjI before f_hash\n";
#endif
    size_t hash = f_hash(seed_hashmap, f_adji_obj(x_adjI));
    std::vector<size_t> &vect = indices_by_hash[hash];
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_singleEntry_AdjI hash=" << hash
       << " |vect|=" << vect.size() << "\n";
#endif
    for (auto &idx : vect) {
      Tobj y = f_idx_obj(idx);
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: process_singleEntry_AdjI Before f_repr idx=" << idx
         << "\n";
#endif
      std::optional<TadjO> opt = f_repr(y, x_adjI, idx);
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: process_singleEntry_AdjI After f_repr\n";
#endif
      if (opt) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: process_singleEntry_AdjI Find an equivalence\n";
#endif
        return *opt;
      }
    }
    std::pair<Tobj, TadjO> pair = f_spann(x_adjI, n_obj);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_singleEntry_AdjI after f_spann\n";
#endif
    vect.push_back(n_obj);
    bool test = f_insert(pair.first);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_singleEntry_AdjI after f_insert test=" << test
       << "\n";
#endif
    if (test) {
      early_termination = true;
    }
    bool is_treated = false;
    f_save_status(n_obj, is_treated);
    undone.push_back(n_obj);
    n_obj++;
    return pair.second;
  };
  auto insert_load = [&](Tobj const &x, bool const &is_treated) -> void {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    std::vector<size_t> &vect = indices_by_hash[hash_hashmap];
    vect.push_back(n_obj);
    if (!is_treated) {
      undone.push_back(n_obj);
    }
    n_obj++;
  };
  auto get_undone_idx = [&]() -> size_t {
    size_t idx = undone[undone.size() - 1];
    undone.pop_back();
    return idx;
  };
  auto treat_one_entry = [&]() -> void {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry beginning\n";
#endif
    size_t idx = get_undone_idx();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry idx=" << idx << "\n";
#endif
    f_save_status(idx, true);
    std::vector<TadjO> l_adj;
    for (auto &y : f_adj(idx)) {
      TadjO adj = process_singleEntry_AdjI(y);
      l_adj.push_back(adj);
    }
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry |l_adj|=" << l_adj.size() << "\n";
#endif
    f_set_adj(idx, l_adj);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry after f_set_adj\n";
#endif
  };
  while (true) {
    std::optional<std::pair<bool, Tobj>> opt = f_next();
    if (opt) {
      std::pair<bool, Tobj> const &pair = *opt;
      Tobj const &x = pair.second;
      bool is_treated = pair.first;
      insert_load(x, is_treated);
    } else {
      break;
    }
  }
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: after load loop n_obj=" << n_obj << "\n";
#endif
  if (n_obj == 0) {
    Tobj x = f_init();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: after f_init\n";
#endif
    bool is_treated = false;
    insert_load(x, is_treated);
    f_save_status(0, is_treated);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: before f_insert\n";
#endif
    bool test = f_insert(x);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: after f_insert\n";
#endif
    if (test) {
      early_termination = true;
    }
  }
  while (true) {
    if (early_termination) {
      return false;
    }
    if (undone.size() == 0) {
      return true;
    }
    if (max_time_second > 0 && si(start) > max_time_second) {
      return false;
    }
    treat_one_entry();
  }
}

template <typename T>
void PartialEnum_FullRead(std::string const &prefix, std::string const &suffix,
                          bool const &Saving, std::vector<T> &l_obj,
                          std::vector<uint8_t> &l_status, std::ostream &os) {
  std::string FileNb = prefix + "number_orbit" + suffix;
  std::string FileStatus = prefix + "orbit_status" + suffix;
  std::string FileDatabase = prefix + "database" + suffix;
  bool is_database_present = false;
  if (Saving) {
    is_database_present = IsExistingFile(FileNb);
  }
  if (is_database_present) {
    size_t n_orbit = FileNumber_Read(FileNb);
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "ADJ_SCH: reading database n_orbit=" << n_orbit << "\n";
#endif
    l_status = FileBool_FullRead(FileStatus, n_orbit, os);
    l_obj = FileData_FullRead<T>(FileDatabase, os);
#ifdef DEBUG_DELAUNAY_ENUMERATION
    os << "ADJ_SCH: reading database l_obj read\n";
#endif
    if (l_obj.size() != n_orbit) {
      std::cerr << "We have n_ent=" << l_obj.size() << " n_orbit=" << n_orbit
                << "\n";
      std::cerr << "But they should be matching\n";
      throw TerminalException{1};
    }
  }
}

template <typename T>
void PartialEnum_FullWrite(std::string const &prefix, std::string const &suffix,
                           bool const &Saving, std::vector<T> const &l_obj,
                           std::vector<uint8_t> const &l_status,
                           [[maybe_unused]] std::ostream &os) {
  std::string FileNb = prefix + "number_orbit" + suffix;
  std::string FileStatus = prefix + "orbit_status" + suffix;
  std::string FileDatabase = prefix + "database" + suffix;
  if (Saving) {
    size_t n_obj = l_obj.size();
    // The data
    FileData_FullWrite(FileDatabase, l_obj);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: writing database fdata written down\n";
#endif
    // The status
    FileBool_FullWrite(FileStatus, l_status);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: writing database FileStatus written down\n";
#endif
    // The number
    FileNumber_Write(FileNb, n_obj);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: writing database FileNumber written down\n";
#endif
  }
}

template <typename Tstor, typename Tout, typename F> struct NextIterator {
  std::vector<Tstor> &l_obj;
  std::vector<uint8_t> &l_status;
  F f;
  size_t pos_next;
  NextIterator(std::vector<Tstor> &_l_obj, std::vector<uint8_t> &_l_status,
               F _f)
      : l_obj(_l_obj), l_status(_l_status), f(_f), pos_next(0) {}
  std::optional<std::pair<bool, Tout>> f_next() {
    if (pos_next >= l_obj.size()) {
      return {};
    } else {
      bool is_treated = static_cast<bool>(l_status[pos_next]);
      Tout x = f(l_obj[pos_next]);
      std::pair<bool, Tout> pair{is_treated, x};
      pos_next++;
      return pair;
    }
  }
};

template <typename TadjO> struct AdjO_MPI {
  TadjO x;
  int iProc;
  int iOrb;
};

template <typename TadjO> struct AdjO_Serial {
  TadjO x;
  int iOrb;
};

template <typename TadjO>
void WriteEntryGAP(std::ostream &os_out, AdjO_MPI<TadjO> const &adj) {
  os_out << "rec(x:=";
  WriteEntryGAP(os_out, adj.x);
  os_out << ", iProc:=" << adj.iProc << ", iOrb:=" << adj.iOrb << ")";
}

template <typename TadjO>
void WriteEntryGAP(std::ostream &os_out, AdjO_Serial<TadjO> const &adj) {
  os_out << "rec(x:=";
  WriteEntryGAP(os_out, adj.x);
  os_out << ", iOrb:=" << adj.iOrb << ")";
}

namespace boost::serialization {
template <class Archive, typename TadjO>
inline void serialize(Archive &ar, AdjO_MPI<TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("iProc", eRec.iProc);
  ar &make_nvp("iOrb", eRec.iOrb);
}
template <class Archive, typename TadjO>
inline void serialize(Archive &ar, AdjO_Serial<TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("iOrb", eRec.iOrb);
}
} // namespace boost::serialization

template <typename Tobj, typename TadjO> struct DatabaseEntry_MPI {
  Tobj x;
  std::vector<AdjO_MPI<TadjO>> ListAdj;
};

template <typename Tobj, typename TadjO> struct DatabaseEntry_Serial {
  Tobj x;
  std::vector<AdjO_Serial<TadjO>> ListAdj;
};

template <typename Tobj, typename TadjO>
void WriteEntryGAP(std::ostream &os_out,
                   DatabaseEntry_MPI<Tobj, TadjO> const &dat_entry) {
  os_out << "rec(x:=";
  WriteEntryGAP(os_out, dat_entry.x);
  os_out << ", ListAdj:=[";
  bool IsFirst = true;
  for (auto &eAdj : dat_entry.ListAdj) {
    if (!IsFirst)
      os_out << ",";
    IsFirst = false;
    WriteEntryGAP(os_out, eAdj);
  }
}

template <typename Tobj, typename TadjO>
void WriteEntryGAP(std::ostream &os_out,
                   DatabaseEntry_Serial<Tobj, TadjO> const &dat_entry) {
  os_out << "rec(x:=";
  WriteEntryGAP(os_out, dat_entry.x);
  os_out << ", ListAdj:=[";
  bool IsFirst = true;
  for (auto &eAdj : dat_entry.ListAdj) {
    if (!IsFirst)
      os_out << ",";
    IsFirst = false;
    WriteEntryGAP(os_out, eAdj);
  }
}

namespace boost::serialization {
template <class Archive, typename Tobj, typename TadjO>
inline void serialize(Archive &ar, DatabaseEntry_MPI<Tobj, TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
template <class Archive, typename Tobj, typename TadjO>
inline void serialize(Archive &ar, DatabaseEntry_Serial<Tobj, TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
} // namespace boost::serialization

template <typename Tobj, typename TadjO>
std::vector<DatabaseEntry_Serial<Tobj, TadjO>>
my_mpi_gather(boost::mpi::communicator &comm,
              std::vector<DatabaseEntry_MPI<Tobj, TadjO>> const &blk,
              int const &i_proc_out) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  using T = typename std::vector<DatabaseEntry_MPI<Tobj, TadjO>>;
  std::vector<DatabaseEntry_Serial<Tobj, TadjO>> V;
  if (i_rank == i_proc_out) {
    std::vector<T> l_blk;
    boost::mpi::gather<T>(comm, blk, l_blk, i_proc_out);
    std::vector<int> l_sizes(n_proc), l_shift(n_proc);
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      l_sizes[i_proc] = l_blk[i_proc].size();
    }
    l_shift[0] = 0;
    for (int i_proc = 1; i_proc < n_proc; i_proc++) {
      l_shift[i_proc] = l_shift[i_proc - 1] + l_sizes[i_proc - 1];
    }
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      for (int u = 0; u < l_sizes[i_proc]; u++) {
        std::vector<AdjO_Serial<TadjO>> ListAdj;
        for (auto &eAdj : l_blk[i_proc][u].ListAdj) {
          int iOrb = eAdj.iOrb + l_shift[eAdj.iProc];
          AdjO_Serial<TadjO> adj{eAdj.x, iOrb};
          ListAdj.emplace_back(std::move(adj));
        }
        DatabaseEntry_Serial<Tobj, TadjO> entry{l_blk[i_proc][u].x,
                                                std::move(ListAdj)};
        V.emplace_back(std::move(entry));
      }
    }
  } else {
    boost::mpi::gather<T>(comm, blk, i_proc_out);
  }
  return V;
}

/*
  This code takes functions and then build new integrated functions for the MPI
  calls.
  It also manages the storage stuff.
 */
template <typename Tdata>
std::pair<
    bool,
    std::vector<DatabaseEntry_MPI<typename Tdata::Tobj, typename Tdata::TadjO>>>
EnumerateAndStore_MPI(boost::mpi::communicator &comm, Tdata &data,
                      std::string const &Prefix, bool const &Saving,
                      int const &max_runtime_second) {
  using Tobj = typename Tdata::Tobj;
  using TadjI = typename Tdata::TadjI;
  using TadjO = typename Tdata::TadjO;
  using TadjO_work = AdjO_MPI<TadjO>;
  std::ostream &os = data.get_os();
  auto f_init = [&]() -> Tobj { return data.f_init(); };
  auto f_hash = [&](size_t const &seed, Tobj const &x) -> size_t {
    return data.f_hash(seed, x);
  };
  auto f_repr = [&](Tobj const &x, TadjI const &y, int const &i_rank,
                    int const &i_orb) -> std::optional<TadjO_work> {
    std::optional<TadjO> opt = data.f_repr(x, y);
    if (opt) {
      TadjO_work ret{*opt, i_rank, i_orb};
      return ret;
    } else {
      return {};
    }
  };
  auto f_spann = [&](TadjI const &x, int i_rank,
                     int i_orb) -> std::pair<Tobj, TadjO_work> {
    std::pair<Tobj, TadjO> pair = data.f_spann(x);
    TadjO_work xo_work{pair.second, i_rank, i_orb};
    std::pair<Tobj, TadjO_work> pair_ret{pair.first, xo_work};
    return pair_ret;
  };
  std::vector<DatabaseEntry_MPI<Tobj, TadjO>> l_obj;
  std::vector<uint8_t> l_status;
  auto f_adj = [&](int const &i_orb) -> std::vector<TadjI> {
    Tobj &x = l_obj[i_orb].x;
    return data.f_adj(x);
  };
  auto f_set_adj = [&](int const &i_orb,
                       std::vector<TadjO_work> const &ListAdj) -> void {
    l_obj[i_orb].ListAdj = ListAdj;
  };
  auto f_adji_obj = [&](TadjI const &x) -> Tobj { return data.f_adji_obj(x); };
  auto f_idx_obj = [&](size_t const &idx) -> Tobj { return l_obj[idx].x; };
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string str_proc =
      "_nproc" + std::to_string(n_proc) + "_rank" + std::to_string(i_rank);
  PartialEnum_FullRead(Prefix, str_proc, Saving, l_obj, l_status, os);
  size_t pos_next = 0;
  auto f_next = [&]() -> std::optional<std::pair<bool, Tobj>> {
    if (pos_next >= l_obj.size()) {
      return {};
    } else {
      bool is_treated = static_cast<bool>(l_status[pos_next]);
      Tobj x = l_obj[pos_next].x;
      std::pair<bool, Tobj> pair{is_treated, x};
      pos_next++;
      return pair;
    }
  };
  auto f_insert = [&](Tobj const &x) -> bool {
    l_obj.push_back({x, {}});
    return false;
  };
  auto f_save_status = [&](size_t const &pos, bool const &val) -> void {
    uint8_t val_i = static_cast<uint8_t>(val);
    if (l_status.size() <= pos) {
      l_status.push_back(val_i);
    } else {
      l_status[pos] = val_i;
    }
  };
  bool test = compute_adjacency_mpi<
      Tobj, TadjI, TadjO_work, decltype(f_next), decltype(f_insert),
      decltype(f_adji_obj), decltype(f_idx_obj), decltype(f_save_status),
      decltype(f_init), decltype(f_adj), decltype(f_set_adj), decltype(f_hash),
      decltype(f_repr), decltype(f_spann)>(
      comm, max_runtime_second, f_next, f_insert, f_adji_obj, f_idx_obj,
      f_save_status, f_init, f_adj, f_set_adj, f_hash, f_repr, f_spann, os);
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: Termination test=" << test << "\n";
#endif
  PartialEnum_FullWrite(Prefix, str_proc, Saving, l_obj, l_status, os);
  return {test, std::move(l_obj)};
}

template <typename Tdata, typename Fincorrect>
std::optional<std::vector<
    DatabaseEntry_Serial<typename Tdata::Tobj, typename Tdata::TadjO>>>
EnumerateAndStore_Serial(Tdata &data, Fincorrect f_incorrect,
                         int const &max_runtime_second) {
  using Tobj = typename Tdata::Tobj;
  using TadjI = typename Tdata::TadjI;
  using TadjO = typename Tdata::TadjO;
  using TadjO_work = AdjO_Serial<TadjO>;
  std::ostream &os = data.get_os();
  auto f_init = [&]() -> Tobj { return data.f_init(); };
  auto f_hash = [&](size_t const &seed, Tobj const &x) -> size_t {
    return data.f_hash(seed, x);
  };
  auto f_repr = [&](Tobj const &x, TadjI const &y,
                    int const &i_orb) -> std::optional<TadjO_work> {
    std::optional<TadjO> opt = data.f_repr(x, y);
    if (opt) {
      TadjO_work ret{*opt, i_orb};
      return ret;
    } else {
      return {};
    }
  };
  auto f_spann = [&](TadjI const &x, int i_orb) -> std::pair<Tobj, TadjO_work> {
    std::pair<Tobj, TadjO> pair = data.f_spann(x);
    TadjO_work xo_work{pair.second, i_orb};
    std::pair<Tobj, TadjO_work> pair_ret{pair.first, xo_work};
    return pair_ret;
  };
  std::vector<DatabaseEntry_Serial<Tobj, TadjO>> l_obj;
  std::vector<uint8_t> l_status;
  auto f_adj = [&](int const &i_orb) -> std::vector<TadjI> {
    Tobj &x = l_obj[i_orb].x;
    return data.f_adj(x);
  };
  auto f_set_adj = [&](int const &i_orb,
                       std::vector<TadjO_work> const &ListAdj) -> void {
    l_obj[i_orb].ListAdj = ListAdj;
  };
  auto f_adji_obj = [&](TadjI const &x) -> Tobj { return data.f_adji_obj(x); };
  auto f_idx_obj = [&](size_t const &idx) -> Tobj { return l_obj[idx].x; };
  auto f_next = [&]() -> std::optional<std::pair<bool, Tobj>> { return {}; };
  auto f_insert = [&](Tobj const &x) -> bool {
    l_obj.push_back({x, {}});
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: EnumerateAndStore_Serial: |l_obj|=" << l_obj.size() << "\n";
#endif
    bool test = f_incorrect(x);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: EnumerateAndStore_Serial: test=" << test << "\n";
    size_t n_sum = 0;
    for (auto & status : l_status) {
      n_sum += static_cast<size_t>(status);
    }
    os << "ADJ_SCH: EnumerateAndStore_Serial: n_sum=" << n_sum << " test=" << test << "\n";
#endif
    return test;
  };
  auto f_save_status = [&](size_t const &pos, bool const &val) -> void {
    uint8_t val_i = static_cast<uint8_t>(val);
    if (l_status.size() <= pos) {
      l_status.push_back(val_i);
    } else {
      l_status[pos] = val_i;
    }
  };
  bool test = compute_adjacency_serial<
      Tobj, TadjI, TadjO_work, decltype(f_next), decltype(f_insert),
      decltype(f_adji_obj), decltype(f_idx_obj), decltype(f_save_status),
      decltype(f_init), decltype(f_adj), decltype(f_set_adj), decltype(f_hash),
      decltype(f_repr), decltype(f_spann)>(
      max_runtime_second, f_next, f_insert, f_adji_obj, f_idx_obj,
      f_save_status, f_init, f_adj, f_set_adj, f_hash, f_repr, f_spann, os);
  if (!test) {
    return {};
  }
  return l_obj;
}

template <typename Tobj, typename TadjO>
void WriteFamilyObjects(
    boost::mpi::communicator &comm, std::string const &OutFormat,
    std::string const &OutFile,
    std::vector<DatabaseEntry_MPI<Tobj, TadjO>> const &l_loc,
    [[maybe_unused]] std::ostream &os) {
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  int i_proc_out = 0;
  int i_rank = comm.rank();
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return;
  }
  if (OutFormat == "ObjectGAP") {
    std::vector<Tout> l_tot = my_mpi_gather(comm, l_loc, i_proc_out);
    if (i_rank == i_proc_out) {
      std::ofstream os_out(OutFile);
      os_out << "return [";
      size_t len = l_tot.size();
      for (size_t i = 0; i < len; i++) {
        if (i > 0)
          os_out << ",\n";
        WriteEntryGAP(os_out, l_tot[i].x);
      }
      os_out << "];\n";
    }
    return;
  }
  if (OutFormat == "ObjectAdjacencyGAP") {
    std::vector<Tout> l_tot = my_mpi_gather(comm, l_loc, i_proc_out);
    if (i_rank == i_proc_out) {
      std::ofstream os_out(OutFile);
      os_out << "return [";
      size_t len = l_tot.size();
      for (size_t i = 0; i < len; i++) {
        if (i > 0)
          os_out << ",\n";
        WriteEntryGAP(os_out, l_tot[i]);
      }
      os_out << "];\n";
    }
    return;
  }
  if (OutFormat == "AdjacencyGAP") {
    // It is actually a little suboptimal but until it is a problem,
    // let us leave it at that.
    std::vector<Tout> l_tot = my_mpi_gather(comm, l_loc, i_proc_out);
    if (i_rank == i_proc_out) {
      std::ofstream os_out(OutFile);
      std::set<int> set;
      size_t len = l_tot.size();
      os_out << "return [";
      for (size_t i = 0; i < len; i++) {
        if (i > 0) {
          os_out << ",\n";
        }
        set.clear();
        for (auto &eAdj : l_tot[i].ListAdj) {
          int pos = eAdj.iOrb + 1;
          set.insert(pos);
        }
        os_out << "[";
        bool IsFirst = true;
        for (auto &eAdj : set) {
          if (!IsFirst) {
            os_out << ",";
          }
          IsFirst = false;
          os_out << eAdj;
        }
        os_out << "]";
      }
      os_out << "];\n";
    }
    return;
  }
  if (OutFormat == "NumberGAP") {
    size_t nb_loc = l_loc.size();
    size_t nb_tot;
    all_reduce(comm, nb_loc, nb_tot, std::plus<size_t>());
    if (i_rank == i_proc_out) {
      std::ofstream os_out(OutFile);
      os_out << "return rec(nb:=" << nb_tot << ");\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
