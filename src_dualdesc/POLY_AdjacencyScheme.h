// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_

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
  But we could work out with f_repr. And if by accident we have some canonicalization
  then we can have a corresponding hash function and then an equality test. So,
  there is no loss of generality in working that way.

  TODO:
  ---Resolve the scheme so that we can conclude when everything has been treated
  without any heuristic.
  ---When a node has everything processed, send a note to all others asking for their
  status. If all of them reply then things are clear.
  ---Have to learn about "distributed programming" from Leslie Lamport and others.
  Takada_M._Distributed_systems_for_fun_and_profit.epub (short one)
  Varela_C.A._Programming_distributed_computing_systems__a_foundational_approach.pdf
  (ok, teaches pi-calculus, join, etc.)
  Raynal_M._Distributed_Algorithms_for_Message-Passing_Systems.pdf
  (especially relevant to us)
  Fokkink_W._Distributed_Algorithms_-_An_Intuitive_Approach.pdf

  ---The Fokkink book provides a reasonable introduction to the question.


  Possible design:
  ---We can have two outcomes from a sending of data:
     ---For unlimited runtime: The termination is that all the entries in the vector<Tobj>
     have their adjacencies being computed.
     ---For specified runtime: The termination is that all the entries that have been sent
     have been inserted into the database.
  ---This is the definition of passive/active:
     ---It ensures that no message in flight is left.
     ---It ensures that no computation is finished.
  ---But that does not allow to determine that everything is finished.
  ---We need to keep a track of the number of operations done.
  ---So, we send a message to all the nodes and collect all the operation done.
  ---Then we call again and see if the nonce are the same.
  ---But the problem is that sending all the entries can introduce some deadlocks maybe.

 */

const size_t seed_partition = 10;
const size_t seed_hashmap = 20;





/*
  input clear from preceding discussion.
  The returned boolean is:
  true: if the enumeration finished
 */
template<typename Tobj, typename Finit, typename Fadj, typename Fhash, typename Frepr>
bool compute_adjacency_mpi(boost::mpi::communicator &comm,
                           int const& max_time_second,
                           Fexist f_exist,
                           Fsave f_save,
                           Fload f_load,
                           Fsave_status f_save_status,
                           Fload_status f_load_status,
                           Finit f_init, Fadj f_adj, Fhash f_hash, frepr f_repr) {
  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  const int tag_new_object = 34;
  const int tag_indicate_processed = 35;
  const int tag_query_n_oper_ask = 36;
  const int tag_query_n_oper_reply = 37;
  const int tag_termination = 38;
  //
  // The data sets
  //
  int n_obj = 0;
  std::unordered_map<size_t, std::vector<size_t>> map;
  using Ttrack = std::pair<size_t, int>;
  struct entry {
    Tobj x;
    size_t hash_hashmap;
    bool is_treated;
    Ttrack track;
  };
  std::vector<entry> V;
  std::vector<size_t> undone;
  //
  // The tracking information
  //
  std::vector<Ttrack> l_ack_to_send;
  std::unordered_set<Ttrack> s_ack_waiting;
  int nonce = 1; // This is so that never ever two objects get generated with the same id.
  bool is_past_time = false;
  //
  // The lambda functions
  //
  auto f_insert_local=[&](entry const& e) -> void {
    std::vector<size_t>& vect = map[hash];
    for (auto & idx : vect) {
      entry & f = V[idx];
      if (f_repr(e.x, f.x)) {
        if (e.track.second != i_proc) {
          l_ack.push_back(pair);
        }
      }
    }
    V.push_back(e);
    vect.push_back(n_obj);
    if (pair.second != i_proc) {
      l_ack.push_back(pair);
    }
    n_obj++;
    nonce++;
  };
  auto f_insert_local_and_save=[&](entry const& e) -> void {
    f_save(n_obj, e.x);
    f_save_status(n_obj, e.is_treated);
    f_insert_local(e);
  };
  auto f_insert=[&](entry const& e) -> void {
    if (res == i_proc) {
      f_insert_local(e);
    } else {
      l_ack_to_send.push_back(e.track);
      // We drop the return status here
      comm.isend(res, tag_new_object, e);
    }
  };
  auto get_nonce = [&]() -> size_t {
    if (s_ack_waiting.size() > 0) {
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
  auto process_received_nonce=[&](track const& t) -> void {
    std::unordered_set<Ttrack>::iterator iter =  s_ack_waiting.find(t);
    if (t == s_ack_waiting.end()) {
      std::cerr << "We have a consistency error\n";
      throw TerminalException{1};
    }
    s_ack_waiting.erase(iter);
  }
  auto process_mpi_status = [&](boost::mpi::status const &stat) -> void {
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_new_object) {
      entry e;
      comm.recv(e_src, tag_new_object, e);
      f_insert_local(e);
      return false;
    }
    if (e_tag == tag_indicate_processed) {
      track t;
      comm.recv(e_src, tag_indicate_processed, t);
      process_received_nonce(t);
      return false;
    }
    if (e_tag == tag_query_n_oper_ask) {
      int val;
      comm.recv(e_src, tag_query_n_oper_ask, val);
      size_t nonce = get_nonce();
      comm.isend(e_src, tag_query_n_oper_reply, nonce);
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
  auto get_undone_idx=[&]() -> size_t {
    size_t idx = undone[undone.size() - 1];
    undone.pop_back();
    return idx;
  };
  auto process_one_entry=[&]() -> void {
    size_t idx = get_undone_idx();
    entry & e = V[idx];
    std::vector<Tobj> l_adj = f_adj(e.x);
    for (auto & x : l_adj) {
      size_t hash_partition = f_hash(seed_partition, x);
      size_t hash_hashmap = f_hash(seed_hashmap, x);
      int res = static_cast<int>(hash_partition % size_t(n_proc));
      Track track{nonce, res};
      nonce++;
      bool is_treated = false;
      entry e{x, hash_hashmap, is_treated, track};
      f_insert(e);
    }
  };
  //
  // Loading the data
  //
  while(true) {
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
  while(true) {
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      os << "prob is not empty\n";
      bool test = process_mpi_status(*prob);
      if (test) {
        break;
      }
    } else {
      if (l_ack_to_send.size() > 0) {
        for (auto & track : l_ack_to_send) {
          comm.isend(e_src, tag_indicate_processed, track);
        }
      } else {
        if (undone.size() > 0) {
          process_one_entry();
        }
      }
    }
  }
}

template<typename Tobj, typename Finit, typename Fadj, typename Fhash, typename Frepr>
bool compute_adjacency_serial(int const& max_time_second,
                              Fexist f_exist,
                              Fsave f_save,
                              Fload f_load,
                              Fsave_status f_save_status,
                              Fload_status f_load_status,
                              Finit f_init, Fadj f_adj, Fhash f_hash, frepr f_repr) {
  SingletonTime start;
  size_t n_obj = 0;
  std::vector<Tobj> V;
  std::unordered_map<size_t, std::vector<size_t>> map;
  std::vector<int> Vstatus;
  std::vector<size_t> undone;
  auto f_insert=[&](Tobj const& x, bool is_treated, bool const& save_if_new) -> void {
    size_t hash = f_hash(seed_hashmap, x);
    std::vector<size_t>& vect = map[hash];
    for (auto & idx : vect) {
      Tobj & y = V[idx];
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
  auto get_undone_idx=[&]() -> size_t {
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
    Tobj const& x = V[idx];
    for (auto & y : f_adj(x)) {
      bool is_treated = false;
      bool save_if_new = false;
      f_insert(y, is_treated, save_if_new);
    }
  }
}


// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
