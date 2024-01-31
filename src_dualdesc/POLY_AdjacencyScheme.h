// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_

#include "boost_serialization.h"
#include "MPI_functionality.h"

#ifdef DEBUG
#define DEBUG_ADJACENCY_SCHEME
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

template<typename T>
void append_move(std::vector<T> & v1, std::vector<T> & v2) {
  v1.insert(v1.end(),
            std::make_move_iterator(v2.begin()),
            std::make_move_iterator(v2.end()));
  v2.clear();
}

const size_t seed_partition = 10;
const size_t seed_hashmap = 20;


template<typename TadjI>
struct entryAdjI {
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
}

template<typename TadjO>
struct entryAdjO {
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
}

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
  f_obj(TadjI) -> Tobj: should return the created object from the
     input adjacency.
  f_save_status(int, bool) -> void : save the status in the database
  f_init() -> Tobj : get a starting element
  f_adj(Tobj, i_orb) -> std::vector<TadjI> : get the adjacent object
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
  --- Sending of entriesAdjI
  --- Sending of entriesAdjO
  --- Sending of nonce
  --- termination
 */
template <typename Tobj, typename TadjI, typename TadjO,
          typename Fnext, typename Finsert, typename Fobj,
          typename Fsave_status,
          typename Finit, typename Fadj, typename Fset_adj,
          typename Fhash,
          typename Frepr, typename Fspann>
bool compute_adjacency_mpi(boost::mpi::communicator &comm,
                           int const &max_time_second,
                           Fnext f_next, Finsert f_insert, Fobj f_obj,
                           Fsave_status f_save_status,
                           Finit f_init, Fadj f_adj, Fset_adj f_set_adj,
                           Fhash f_hash, Frepr f_repr, Fspann f_spann,
                           [[maybe_unused]] std::ostream & os) {
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
  unlimited_request ur(comm);
  //
  // The data sets
  //
  int n_obj = 0; // The number of objects generated
  // The map from the hash to the list of indices
  std::unordered_map<size_t, std::vector<size_t>> map;
  std::vector<Tobj> V; // The objects
  std::vector<size_t> undone; // The undone indices
  //
  // The entries of AdjI / AdjO
  //
  // The unsent entries by the processors.
  buffered_T_exchanges<entryAdjI<TadjI>, std::vector<entryAdjI<TadjI>>> buffer_entriesAdjI(comm, n_proc, tag_entriesadji_send);
  buffered_T_exchanges<entryAdjO<TadjO>, std::vector<entryAdjO<TadjO>>> buffer_entriesAdjO(comm, n_proc, tag_entriesadjo_send);
  std::vector<entryAdjI<TadjI>> unproc_entriesAdjI;
  // The mapping from the index to the list of adjacencices.
  std::unordered_map<int, std::pair<size_t, std::vector<TadjO>>> map_adjO;

  // The nonce is used so that a number is associated to a specific computation.
  // The function get_nonce returns 0 if there is something left to do and
  // nonzero if there is no pending computation. If the get_nonce returns
  // the same, then there was no computation going on.
  // The none is also used for the tracking of the adjacencies entries.
  size_t nonce = 1;
  bool early_termination = false;
  //
  // The lambda functions
  //
  auto send_early_termination=[&]() -> void {
    early_termination = true;
    int val_final = 0;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      if (i_proc != i_rank) {
        ur.get_entry() = comm.isend(i_proc, tag_early_termination, val_final);
      }
    }
  };
  auto process_single_entryAdjI = [&](entryAdjI<TadjI> const &eI) -> std::pair<int,entryAdjO<TadjO>> {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Beginning of process_single_entryAdjI\n";
#endif
    auto f_ret=[&](TadjO u) -> std::pair<int,entryAdjO<TadjO>> {
      entryAdjO<TadjO> eA{u, eI.i_orb_orig};
      return {eI.i_proc_orig, eA};
    };
    nonce++;
    std::vector<size_t> &vect = map[eI.hash_hashmap];
    for (auto &idx : vect) {
      Tobj &x = V[idx];
      std::optional<TadjO> opt = f_repr(x, eI.x, i_rank, idx);
      if (opt) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: Conclude with an equivalence\n";
#endif
        return f_ret(*opt);
      }
    }
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Conclude with a new one\n";
#endif
    std::pair<Tobj, TadjO> pair = f_spann(eI.x, i_rank, n_obj);
    bool test = f_insert(pair.first);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_insert test=" << test << "\n";
#endif
    if (test) {
      send_early_termination();
    }
    V.emplace_back(std::move(pair.first));
    undone.push_back(n_obj);
    vect.push_back(n_obj);
    f_save_status(n_obj, false);
    n_obj++;
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_single_entryAdjI, now n_obj=" << n_obj << "\n";
#endif
    return f_ret(pair.second);
  };
  auto insert_load=[&](Tobj const& x, bool const& is_treated) -> void {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    V.push_back(x);
    std::vector<size_t> &vect = map[hash_hashmap];
    vect.push_back(n_obj);
    if (!is_treated) {
      undone.push_back(n_obj);
    }
    n_obj++;
  };
  auto initial_insert=[&](Tobj const& x) -> void {
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
  auto process_entriesAdjO=[&](std::vector<entryAdjO<TadjO>> & v) -> void {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH:process_entriesAdjO : Before |v|=" << v.size() << " |map_adjO|=" << map_adjO.size() << " \n";
    for (auto & kv : map_adjO) {
      os << "ADJ_SCH:process_entriesAdjO : Before i_orb=" << kv.first << " siz=" << kv.second.first << " |l_adj|=" << kv.second.second.size() << "\n";
    }
#endif
    for (auto &eEnt : v) {
      map_adjO[eEnt.i_orb_orig].second.emplace_back(std::move(eEnt.x));
    }
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH:process_entriesAdjO : After |map_adjO|=" << map_adjO.size() << " \n";
    for (auto & kv : map_adjO) {
      os << "ADJ_SCH:process_entriesAdjO : After i_orb=" << kv.first << " siz=" << kv.second.first << " |l_adj|=" << kv.second.second.size() << "\n";
    }
#endif
    v.clear();
  };
  auto get_nonce = [&]() -> size_t {
    if (!buffer_entriesAdjI.is_completely_clear()) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because !buffer_entriesAdjI.is_completely_clear()\n";
#endif
      return 0;
    }
    if (!buffer_entriesAdjO.is_completely_clear()) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: the_nonce=0, because !buffer_entriesAdjO.is_completely_clear()\n";
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
      os << "ADJ_SCH: the_nonce=0, because undone.size() > 0 && max_time_second == 0\n";
#endif
      return 0;
    }
    if (undone.size() > 0 && max_time_second > 0 && si(start) < max_time_second) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: max_time_second=" << max_time_second << " si(start)=" << si(start) << " |undone|=" << undone.size() << "\n";
      os << "ADJ_SCH: the_nonce=0, because max_time_second > 0 && si(start) < max_time_second\n";
#endif
      return 0;
    }
    return nonce;
  };
  auto process_mpi_status = [&](boost::mpi::status const &stat) -> bool {
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_initial) {
      Tobj x;
      comm.recv(e_src, tag_initial, x);
      initial_insert(x);
      return false;
    }
    if (e_tag == tag_entriesadji_send) {
      std::vector<entryAdjI<TadjI>> v;
      comm.recv(e_src, tag_entriesadji_send, v);
      append_move(unproc_entriesAdjI, v);
      return false;
    }
    if (e_tag == tag_entriesadjo_send) {
      std::vector<entryAdjO<TadjO>> v;
      comm.recv(e_src, tag_entriesadjo_send, v);
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: Call process_entriesAdjO from process_mpi_status\n";
#endif
      process_entriesAdjO(v);
      return false;
    }
    if (e_tag == tag_nonce_ask) {
      int val_recv;
      comm.recv(e_src, tag_nonce_ask, val_recv);
      size_t nonce = get_nonce();
      ur.get_entry() = comm.isend(e_src, tag_nonce_reply, nonce);
      return false;
    }
    if (e_tag == tag_termination) {
      int val_recv;
      comm.recv(e_src, tag_termination, val_recv);
      return true;
    }
    if (e_tag == tag_early_termination) {
      int val_recv;
      comm.recv(e_src, tag_termination, val_recv);
      early_termination=true;
      return true;
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
    Tobj const& x = V[idx];
    std::vector<TadjI> l_adj = f_adj(x, idx);
    map_adjO[idx] = {l_adj.size(), {}};
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: process_one_entry_obj : idx=" << idx << " |l_adj|=" << l_adj.size() << "\n";
#endif
    nonce++;
    for (auto &x : l_adj) {
      Tobj x_obj = f_obj(x);
      size_t hash_partition = f_hash(seed_partition, x_obj);
      size_t hash_hashmap = f_hash(seed_hashmap, x_obj);
      int i_proc_dest = static_cast<int>(hash_partition % size_t(n_proc));
      nonce++;
      entryAdjI<TadjI> e{x, hash_hashmap, i_rank, idx_i};
      buffer_entriesAdjI.insert_entry(i_proc_dest, e);
    }
    return true;
  };
  auto flush_entriesAdjI=[&]() -> bool {
    std::vector<entryAdjI<TadjI>> & v = buffer_entriesAdjI.l_message[i_rank];
    append_move(unproc_entriesAdjI, v);
    return buffer_entriesAdjI.clear_one_entry(os);
  };
  auto flush_entriesAdjO=[&]() -> bool {
    std::vector<entryAdjO<TadjO>> & v = buffer_entriesAdjO.l_message[i_rank];
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Call process_entriesAdjO from flush_entriesAdjO\n";
#endif
    process_entriesAdjO(v);
    return buffer_entriesAdjO.clear_one_entry(os);
  };
  auto compute_entries_adjI=[&]() -> bool {
    bool do_something = false;
    for (auto & eI : unproc_entriesAdjI) {
      std::pair<int,entryAdjO<TadjO>> pair = process_single_entryAdjI(eI);
      buffer_entriesAdjO.insert_entry(pair.first, pair.second);
      do_something = true;
    }
    unproc_entriesAdjI.clear();
    return do_something;
  };
  auto write_set_adj=[&]() -> bool {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Begin of write_set_adj\n";
#endif
    std::vector<int> l_erase;
    bool do_something = false;
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: |map_adjO|=" << map_adjO.size() << "\n";
#endif
    for (auto & kv : map_adjO) {
      int i_orb = kv.first;
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: i_orb=" << i_orb << " kv.second.first=" << kv.second.first << " |kv.second.second|=" << kv.second.second.size() << "\n";
#endif
      if (kv.second.first == kv.second.second.size()) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: write_set_adj i_orb=" << i_orb << " |l_adj|=" << kv.second.second.size() << "\n";
#endif
        f_set_adj(i_orb, kv.second.second);
        l_erase.push_back(i_orb);
        do_something = true;
      }
    }
    for (auto & i_orb : l_erase) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: write_set_adj erasing i_orb=" << i_orb << "\n";
#endif
      map_adjO.erase(i_orb);
    }
    return do_something;
  };
  auto f_clear_buffers=[&]() -> bool {
    // Transmitting the generated entriesAdjI
    bool test1 = flush_entriesAdjI();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffer : flush_entriesAdjI=" << test1 << "\n";
#endif
    if (test1) {
      return true;
    }
    // transmitting the generated entriesAdjO
    bool test2 = flush_entriesAdjO();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffer : flush_entriesAdjO=" << test2 << "\n";
#endif
    if (test2) {
      return true;
    }
    // compute the entryAdjO from entryAdjI
    bool test3 = compute_entries_adjI();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffer : compute_entries_adjI=" << test3 << "\n";
#endif
    if (test3) {
      return true;
    }
    // write down adjacencies if done
    bool test4 = write_set_adj();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: f_clear_buffer : write_set_adj=" << test4 << "\n";
#endif
    if (test4) {
      return true;
    }
    // Nothing, then computing something new
    return process_one_entry_obj();
  };
  auto terminate = [&]() -> bool {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: terminate, begin\n";
#endif
    if (i_rank > 0) {
      return false;
    }
    std::vector<size_t> l_nonce(n_proc-1);
    for (int i_proc=1; i_proc<n_proc; i_proc++) {
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
      l_nonce[i_proc-1] = the_nonce;
    }
    for (int i_proc=1; i_proc<n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_nonce_ask, val);
      size_t the_nonce;
      comm.recv(i_proc, tag_nonce_reply, the_nonce);
      if (the_nonce != l_nonce[i_proc-1]) {
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: exiting at i_proc=" << i_proc << " the_nonce=" << the_nonce << " l_nonce=" << l_nonce[i_proc-1] << "\n";
#endif
        return false;
      }
    }
    // no operation done. Sending the termination notices
    for (int i_proc=1; i_proc<n_proc; i_proc++) {
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
      std::pair<bool, Tobj> const& pair = *opt;
      Tobj const& x = pair.second;
      bool is_treated = pair.first;
      insert_load(x, is_treated);
    } else {
      break;
    }
  }
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: compute_adjacency_mpi, start : n_obj=" << n_obj << " n_undone=" << undone.size() << "\n";
#endif
  size_t n_orb_max = 0, n_orb_loc = V.size();
  all_reduce(comm, n_orb_loc, n_orb_max, boost::mpi::maximum<size_t>());
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: beginning n_orb_max=" << n_orb_max << " n_orb_loc=" << n_orb_loc << "\n";
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
    os << "ADJ_SCH: Start while loop, early_termination=" << early_termination << "\n";
#endif
    if (early_termination) {
      break;
    }
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: prob is not empty\n";
#endif
      bool test = process_mpi_status(*prob);
      if (test) {
        break;
      }
    } else {
      bool test = f_clear_buffers();
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: f_clear_buffer test=" << test << "\n";
#endif
      if (!test) {
        bool test_terminate=terminate();
#ifdef DEBUG_ADJACENCY_SCHEME
        os << "ADJ_SCH: test_terminate=" << test_terminate << "\n";
#endif
        if (test_terminate) {
          break;
        }
      }
    }
  }
#ifdef DEBUG_ADJACENCY_SCHEME
  os << "ADJ_SCH: compute_adjacency_mpi, end : n_obj=" << n_obj << " n_undone=" << undone.size() << "\n";
#endif
  if (early_termination) {
    // We terminate early so by definition the enumeration is partial
    return false;
  } else {
    size_t n_undone_max = 0, n_undone_loc = undone.size();
    all_reduce(comm, n_undone_loc, n_undone_max, boost::mpi::maximum<size_t>());
    return n_undone_max > 0;
  }
}

template <typename Tobj, typename TadjI, typename TadjO,
          typename Fnext, typename Finsert, typename Fobj,
          typename Fsave_status,
          typename Finit, typename Fadj, typename Fset_adj,
          typename Fhash,
          typename Frepr, typename Fspann>
bool compute_adjacency_serial(int const &max_time_second,
                              Fnext f_next, Finsert f_insert, Fobj f_obj,
                              Fsave_status f_save_status,
                              Finit f_init, Fadj f_adj, Fset_adj f_set_adj,
                              Fhash f_hash, Frepr f_repr, Fspann f_spann,
                              [[maybe_unused]] std::ostream& os) {
  SingletonTime start;
  size_t n_obj = 0;
  std::vector<Tobj> V;
  std::unordered_map<size_t, std::vector<size_t>> map;
  std::vector<size_t> undone;
  bool early_termination = false;
  auto process_singleEntry_AdjI = [&](TadjI const &x_adjI) -> TadjO {
    size_t hash = f_hash(seed_hashmap, f_obj(x_adjI));
    std::vector<size_t> &vect = map[hash];
    for (auto &idx : vect) {
      Tobj &y = V[idx];
      std::optional<TadjO> opt = f_repr(y, x_adjI, idx);
      if (opt) {
        return *opt;
      }
    }
    std::pair<Tobj, TadjO> pair = f_spann(x_adjI, n_obj);
    V.push_back(pair.first);
    vect.push_back(n_obj);
    bool test = f_insert(pair.first);
    if (test) {
      early_termination = true;
    }
    bool is_treated = false;
    f_save_status(n_obj, is_treated);
    undone.push_back(n_obj);
    n_obj++;
    return pair.second;
  };
  auto insert_load=[&](Tobj const& x, bool const& is_treated) -> void {
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    V.push_back(x);
    std::vector<size_t> &vect = map[hash_hashmap];
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
  auto treat_one_entry=[&]() -> void {
    size_t idx = get_undone_idx();
    f_save_status(idx, true);
    Tobj const &x = V[idx];
    std::vector<TadjO> l_adj;
    for (auto &y : f_adj(x, idx)) {
      TadjO adj = process_singleEntry_AdjI(y);
      l_adj.push_back(adj);
    }
    f_set_adj(idx, l_adj);
  };
  while (true) {
    std::optional<std::pair<bool, Tobj>> opt = f_next();
    if (opt) {
      std::pair<bool, Tobj> const& pair = *opt;
      Tobj const& x = pair.second;
      bool is_treated = pair.first;
      insert_load(x, is_treated);
    } else {
      break;
    }
  }
  if (n_obj == 0) {
    Tobj x = f_init();
    bool is_treated = false;
    insert_load(x, is_treated);
    f_save_status(0, is_treated);
    bool test = f_insert(x);
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

// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
