// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_


// We clear the mpi requests and return the total length.
// Total length should be 0 in order to consider exiting.
size_t clear_mpi_request(std::vector<boost::mpi::request> & rsl) {
  size_t len = rsl.size();
  std::vector<boost::mpi::request> new_rsl;
  for (size_t u=0; u<len; u++) {
    boost::optional<boost::mpi::status> stat = rsl[u].test();
    if (!stat) {
      new_rsl.emplace_back(std::move(rsl[u]));
    }
  }
  rsl = std::move(new_rsl);
  return rsl.size();
}

template<typename T>
void append_move(std::vector<T> & v1, std::vector<T> & v2) {
  v1.insert(v1.end(),
            std::make_move_iterator(v2.begin()),
            std::make_move_iterator(v2.end()));
  v2.clear();
}


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

const size_t seed_partition = 10;
const size_t seed_hashmap = 20;

/*
  input clear from preceding discussion.
  The returned boolean is:
  true: if the enumeration finished. false otherwise.
  ---
  Types (cannot be set up via concept, but so be it)
  Tobj: The object type being created (like L-type domain)
  TadjI: The types returned by the spanning. Can be something like
    std::pair<Face, Tobj> with Face indicating the relevant Face.
    Doing the .obj should return the object.
  TadjO: The adjacency type after processing. Can be something like
    {Face, Trans, idx} with Face the corresponding Face, Trans the
    transformation realizing the equivalence and idx the equivalent
    object.
  All objects have to be transmitible via boost C++ mpi.
  ---
  Function types used:
  f_exist(int) -> bool : whether there is an entry at this level
  f_insert(Tobj) -> void : insert the new object.
  f_load(int) -> Tobj : get one object from the database
  f_save_status(int, bool) -> void : save the status in the database
  f_load_status(int) -> bool : get the status in the database
  f_init() -> Tobj : get a starting element
  f_adj(Tobj) -> std::vector<TadjI> : get the adjacent object
  f_set_adj(int, std::vector<TadjO>) -> void : set the adjacencies to the
  f_hash(size_t, Tobj) -> size_t : compute the hash from a specified seed.
  f_repr(Tobj, TadjI, int, int) -> std::optional<TadjO> : returns whether
    Tobj is equivalent to the spanned TadjI and find the equivalence if that
    is the case. Also take
  f_spann(TadjI, int, int) -> std::pair<Tobj, TabjO> : Generate from the
    equivalence the object to be inserted and the equivalence to be sent.
 */
template <typename Tobj, typename TadjI, typename TadjO,
          typename Fexist, typename Finsert, typename Fload,
          typename Fsave_status, typename Fload_status,
          typename Finit, typename Fadj, typename Fset_adj,
          typename Fhash,
          typename Frepr, typename Fspann>
bool compute_adjacency_mpi(boost::mpi::communicator &comm,
                           std::ostream & os,
                           int const &max_time_second,
                           Fexist f_exist, Finsert f_insert, Fload f_load,
                           Fsave_status f_save_status, Fload_status f_load_status,
                           Finit f_init, Fadj f_adj, Fset_adj f_set_adj,
                           Fhash f_hash, Frepr f_repr, Fspann f_spann) {
  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  const int tag_new_object = 34;
  const int tag_indicate_processed = 35;
  const int tag_nonce_ask = 36;
  const int tag_nonce_reply = 37;
  const int tag_entriesadji_send = 38;
  const int tag_entriesadjo_send = 38;
  const int tag_termination = 39;
  std::vector<boost::mpi::request> rsl_comp, rsl_admin;
  //
  // The combined data types
  //
  using Ttrack = std::pair<size_t,int>;
  struct entryObj {
    Tobj x; // The object 
    size_t hash_hashmap;
    bool is_treated;
    int i_proc;
  };
  struct entryAdjI {
    TadjI x;
    size_t hash_hashmap;
    Ttrack track;
    int i_proc_orig;
    int i_proc_dest;
  };
  struct entryAdjO {
    TadjO x;
  };
  //
  // The data sets
  //
  int n_obj = 0; // The number of objects generated
  // The map from the hash to the list of indices
  std::unordered_map<size_t, std::vector<size_t>> map;
  std::vector<entryObj> V; // The objects
  std::vector<size_t> undone; // The undone indices
  //
  // The tracking information
  //
  // The unsent entries by the processors.
  std::vector<std::vector<entryAdjI>> unsent_entriesAdjI;
  std::vector<std::vector<entryAdjO>> unsent_entriesAdjO;
  std::vector<entryAdjI> unproc_entriesAdjI;
  std::vector<entryAdjO> unproc_entriesAdjO;
  // The mapping from the index to the list of adjacencices.
  std::unordered_map<int, std::pair<size_t, std::vector<TadjO>>> map_adjO;


  
  // The nonce is used so that a number is associated to a specific computation.
  // The function get_nonce returns 0 if there is something left to do and
  // nonzero if there is no pending computation. If the get_nonce returns
  // the same, then there was no computation going on.
  // The none is also used for the tracking of the adjacencies entries.
  size_t nonce = 1;
  bool is_past_time = false;
  bool has_pending_work = true;
  //
  // The lambda functions
  //
  auto insert_local = [&](entryAdjI const &eI) -> void {
    std::vector<size_t> &vect = map[eI.hash_hashmap];
    for (auto &idx : vect) {
      entry &f = V[idx];
      std::optional<TadjO> opt = f_repr(f.x, eI.x, i_rank, idx);
      if (opt) {
        if (eI.track.second != i_proc) {
          l_ack.push_back(e.track);
        }
      }
    }
    V.push_back(e);
    vect.push_back(n_obj);
    if (e.track.second != i_proc) {
      l_ack.push_back(e.track);
    }
    n_obj++;
    nonce++;
  };
  auto insert_local_and_save = [&](entryI const &eI) -> void {
    f_insert(n_obj, eI.x);
    f_save_status(n_obj, eI.is_treated);
    insert_local(eI);
  };
  auto insert_entry = [&](entryI const &eI) -> void {
    if (eI.track.second == i_proc) {
      insert_local(eI);
    } else {
      l_ack_to_send.push_back(e.track);
      rsl_comp.push_back(comm.isend(res, tag_new_object, e));
    }
  };
  auto get_nonce = [&]() -> size_t {
    if (s_ack_waiting.size() > 0) {
      return 0;
    }
    (void)clear_mpi_request(rsl_admin);
    size_t left_oper = clear_mpi_request(rsl_comp);
    if (left_oper > 0) {
      return 0;
    }
    if (undone.size() == 0) {
      return nonce;
    }
    if (max_time_second > 0 && si(start) > max_time_second) {
      return nonce;
    }
    return 0;
  };
  auto process_received_nonce = [&](track const &t) -> void {
    std::unordered_set<Ttrack>::iterator iter = s_ack_waiting.find(t);
    if (t == s_ack_waiting.end()) {
      std::cerr << "We have a consistency error\n";
      throw TerminalException{1};
    }
    s_ack_waiting.erase(iter);
  };
  auto process_mpi_status = [&](boost::mpi::status const &stat) -> void {
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_new_object) {
      entryI e;
      comm.recv(e_src, tag_new_object, e);
      insert_local(e);
      return false;
    }
    if (e_tag == tag_entriesadji_send) {
      std::vector<entryAdjI>> v;
      comm.recv(e_src, tag_entriesadji_send, v);
      append_move(unproc_entriesAdjI, v);
      return false;
    }
    if (e_tag == tag_nonce_ask) {
      int val;
      comm.recv(e_src, tag_nonce_ask, val);
      size_t nonce = get_nonce();
      rsl_admin.push_back(comm.isend(e_src, tag_nonce_reply, nonce));
      return false;
    }
    if (e_tag == tag_termination) {
      int val;
      comm.recv(e_src, tag_termination, val);
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
  auto process_one_entry_obj = [&]() -> void {
    size_t idx = get_undone_idx();
    entryObj &e = V[idx];
    std::vector<TadjI> l_adj = f_adj(e.x);
    nonce++;
    for (auto &x : l_adj) {
      size_t hash_partition = f_hash(seed_partition, x.obj);
      size_t hash_hashmap = f_hash(seed_hashmap, x.obj);
      int res = static_cast<int>(hash_partition % size_t(n_proc));
      Track track{nonce, res};
      nonce++;
      bool is_treated = false;
      entryAdjI e{x, hash_hashmap, is_treated, track};
      unsent_entriesAdjI[res].push_back(entryAdjI);
    }
  };
  auto flush_entriesAdjI=[&]() -> bool {
    bool do_something = false;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      std::vector<entryAdjI> & v = unsent_entriesAdjI[i_proc];
      if (v.size() > 0) {
        do_something = true;
        if (i_proc == i_rank) {
          append_move(unproc_entriesAdjI, v);
        } else {
          // This does not work because v is passed as reference and
          // so we cannot clear it. We need to transfer to a vector
          rsl_comp.push_back(comm.isend(i_proc, tag_entriesadji_send, v));
        }
      }
    }
    return do_something;
  };
  auto flush_entriesAdjO=[&]() -> bool {
    bool do_something = false;
    for (int i_proc=0; i_proc<n_proc; i_proc++) {
      std::vector<entryAdjO> & v = unsent_entriesAdjO[i_proc];
      if (v.size() > 0) {
        do_something = true;
        if (i_proc == i_rank) {
          for (auto &eEnt : v) {
            
          }
          v.clear();
        } else {
          // This does not work because v is passed as reference and
          // so we cannot clear it. We need to transfer to a vector
          rsl_comp.push_back(comm.isend(i_proc, tag_entriesadji_send, v));
        }
      }
    }
    return do_something;
  };
  auto write_set_adj=[&]() -> bool {
    
  };
  auto terminate = [&]() -> bool {
    std::vector<size_t> l_nonce(n_proc-1);
    for (int i_proc=1; i_proc<n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_nonce_ask, val);
      size_t the_nonce;
      comm.recv(e_src, tag_nonce_reply, the_nonce);
      if (the_nonce == 0) {
        return false;
      }
      l_nonce[i_proc-1] = the_nonce;
    }
    for (int i_proc=1; i_proc<n_proc; i_proc++) {
      int val = 0;
      comm.send(i_proc, tag_nonce_ask, val);
      size_t the_nonce;
      comm.recv(e_src, tag_nonce_reply, the_nonce);
      if (the_nonce != l_nonce[i_proc-1]) {
        return false;
      }
    }
    return true;
  };
  //
  // Loading the data
  //
  while (true) {
    if (!f_exist(n_obj)) {
      break;
    }
    Tobj x = f_load(n_obj);
    size_t hash_hashmap = f_hash(seed_hashmap, x);
    bool is_treated = f_load_status(n_obj);
    Ttrack track{nonce, i_rank};
    entry e{x, hash_hashmap, is_treated, track};
    f_insert_local_and_save(e);
  }
  //
  // The infinite loop
  //
  while (true) {
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      os << "prob is not empty\n";
      bool test = process_mpi_status(*prob);
      if (test) {
        break;
      }
    } else {
      if (l_ack_to_send.size() > 0) {
        for (auto &track : l_ack_to_send) {
          rsl_comp.push_back(comm.isend(e_src, tag_indicate_processed, track));
        }
      } else {
        if (undone.size() > 0) {
          process_one_entry_obj();
        } else {
          if (terminate()) {
            return is_past_time;
          }
        }
      }
    }
  }
}

template <typename Tobj, typename Finit, typename Fadj, typename Fhash,
          typename Frepr>
bool compute_adjacency_serial(std::ostream& os,
                              int const &max_time_second, Fexist f_exist,
                              Finsert f_insert, Fload f_load,
                              Fsave_status f_save_status,
                              Fload_status f_load_status, Finit f_init,
                              Fadj f_adj, Fhash f_hash, frepr f_repr) {
  SingletonTime start;
  size_t n_obj = 0;
  std::vector<Tobj> V;
  std::unordered_map<size_t, std::vector<size_t>> map;
  std::vector<int> Vstatus;
  std::vector<size_t> undone;
  auto f_insert = [&](Tobj const &x, bool is_treated,
                      bool const &save_if_new) -> void {
    size_t hash = f_hash(seed_hashmap, x);
    std::vector<size_t> &vect = map[hash];
    for (auto &idx : vect) {
      Tobj &y = V[idx];
      if (f_repr(x, y)) {
        return;
      }
    }
    V.push_back(x);
    vect.push_back(n_obj);
    if (save_if_new) {
      f_save(n_obj, x);
      f_save_status(n_obj, is_new);
    }
    if (is_treated) {
      undone.push_back(n_obj);
    }
    n_obj++;
  };
  auto get_undone_idx = [&]() -> size_t {
    size_t idx = undone[undone.size() - 1];
    undone.pop_back();
    return idx;
  };
  while (true) {
    if (!f_exist(n_obj)) {
      break;
    }
    Tobj x = f_load(n_obj);
    bool is_treated = f_load_status(n_obj);
    bool save_if_new = false;
    f_insert(x, is_treated, save_if_new);
  }
  if (n_obj == 0) {
    Tobj x = f_init();
    f_insert(x, 0, true);
  }
  while (true) {
    if (unordered_set.size() == 0) {
      return true;
    }
    if (max_time_second > 0 && si(start) > max_time_second) {
      return false;
    }
    size_t idx = get_undone_idx();
    Tobj const &x = V[idx];
    for (auto &y : f_adj(x)) {
      bool is_treated = false;
      bool save_if_new = false;
      f_insert(y, is_treated, save_if_new);
    }
  }
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
