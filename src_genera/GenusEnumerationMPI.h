// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATIONMPI_H_
#define SRC_GENERA_GENUSENUMERATIONMPI_H_

// clang-format off
#include "GenusEnumerationInput.h"
#include "MPI_functionality.h"
#include <boost/mpi.hpp>
#include <cstdint>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_GENUS_ENUMERATION_MPI
#endif

#ifdef DEBUG
# define DEBUG_GENUS_ENUMERATION_MPI
#endif

/*
  The enumeration of a genus over MPI (GENUS_MPI_Enumerate).

  The classes are distributed by their key: the canonical Gram matrix
  (recognition canonical) or the BV invariant (recognition invariant). A
  class is owned by the node hash(key) mod n_proc, which keeps its Gram
  matrix, its automorphism group and its Hecke bookkeeping, and expands it.
  The neighbors built by an expansion get their key on the node that builds
  them, the expensive canonical form being thus spread with the expansions,
  and are sent to their owner, with their Hecke contribution (the orbit
  size over |Aut| of the expanded class), so that the deficits stay exact.
  The owner recognizes them, computes the automorphism group of the new
  ones, and expands them in turn. Every node chooses the next of its own
  classes to expand by the order of the enumeration (deficit, mixed, aut),
  on its own classes.

  Rank 0 coordinates.
  * The mass. Every node reports its local mass, the sum of 1/|Aut| over
    the classes it owns, when it changes. A class is owned by exactly one
    node and a report is a sum over classes found, so the sum of the
    latest reports never exceeds the true total, which never exceeds the
    mass of the genus: when it reaches the mass, the enumeration is
    complete, and rank 0 asks every node to stop. Stale reports only delay
    that.
  * The idleness. When rank 0 has nothing to do, it asks every node for its
    state: busy or not, and the number of entries sent and received. Two
    successive rounds with every node idle, as many entries received as
    sent, and the same counts, mean that nothing is in flight or to do
    (the four-counter method). Then, if the mass is open, iterate brings in
    the next prime and out_prime ends incomplete.
  * The end. After the stop the nodes drain the entries still in flight,
    which rank 0 detects by the same rounds, then the classes are gathered
    on rank 0.

  The strategies out_prime and iterate are supported, and the recognitions
  canonical and invariant (isomorphism needs every class of a bucket on one
  node, which the hash of the key does not give).
 */

// A neighbor sent to the node owning its class.
template <typename T, typename Tmass> struct GenusMPIEntry {
  // the canonical Gram matrix (canonical) or the neighbor (invariant)
  MyMatrix<T> Gram;
  // the index of the prime of the expansion, -1 for a seed
  int i_prime;
  // the orbit size over |Aut| of the expanded class
  Tmass contrib;
  // the BV invariant (invariant)
  T bv_max_norm;
  uint64_t bv_n_vect;
  uint64_t bv_hash;
};

// The state of a node, as reported to the coordinator.
struct GenusMPIStatus {
  int busy;
  uint64_t n_sent;
  uint64_t n_recv;
};

namespace boost::serialization {
template <class Archive, typename T, typename Tmass>
inline void serialize(Archive &ar, GenusMPIEntry<T, Tmass> &e,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("Gram", e.Gram);
  ar &make_nvp("i_prime", e.i_prime);
  ar &make_nvp("contrib", e.contrib);
  ar &make_nvp("bv_max_norm", e.bv_max_norm);
  ar &make_nvp("bv_n_vect", e.bv_n_vect);
  ar &make_nvp("bv_hash", e.bv_hash);
}
template <class Archive>
inline void serialize(Archive &ar, GenusMPIStatus &s,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("busy", s.busy);
  ar &make_nvp("n_sent", s.n_sent);
  ar &make_nvp("n_recv", s.n_recv);
}
} // namespace boost::serialization

// The classes owned by a node at the end, for the gathering on rank 0.
template <typename T, typename TintGroup> struct GenusMPILocalClasses {
  std::vector<MyMatrix<T>> ListGram;
  std::vector<TintGroup> ListAutOrder;
  uint64_t n_expansion;
  uint64_t n_neighbor;
};

namespace boost::serialization {
template <class Archive, typename T, typename TintGroup>
inline void serialize(Archive &ar, GenusMPILocalClasses<T, TintGroup> &c,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("ListGram", c.ListGram);
  ar &make_nvp("ListAutOrder", c.ListAutOrder);
  ar &make_nvp("n_expansion", c.n_expansion);
  ar &make_nvp("n_neighbor", c.n_neighbor);
}
} // namespace boost::serialization

// A hash of a matrix that is the same on every node.
template <typename T> size_t GenusMPIHashMatrix(MyMatrix<T> const &M) {
  std::ostringstream str;
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      str << M(i, j) << " ";
    }
  }
  return std::hash<std::string>{}(str.str());
}

/*
  The distributed enumeration of the even genus of the seeds ListSeed, of
  mass TotalMass, from the classes of which the result keeps those owned by
  this node; complete tells whether the mass closed, on every node.
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<GenusEnumerationResult<T, Tgroup>, bool>
GenusEnumerationMPI(boost::mpi::communicator &comm,
                    std::vector<MyMatrix<T>> const &ListSeed,
                    typename GenusEnumerationResult<T, Tgroup>::Tmass const
                        &TotalMass,
                    int prime, GenusSpec const &spec, int max_runtime_second,
                    std::ostream &os) {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  using TintGroup = typename Tgroup::Tint;
  using Tentry = GenusMPIEntry<T, Tmass>;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  const int tag_entry = 201;
  const int tag_mass = 202;
  const int tag_ask = 203;
  const int tag_reply = 204;
  const int tag_prime = 205;
  const int tag_stop = 206;
  const int tag_final = 207;
  GenusMethods const &methods = spec.methods;
  bool by_canonical = methods.recognition == GenusRecognition::Canonical;
  if (!by_canonical && methods.recognition != GenusRecognition::Invariant) {
    std::cerr << "GENUS_MPI: the recognition must be canonical or invariant\n";
    throw TerminalException{1};
  }
  if (spec.strategy != NeighborStrategy::OutPrime &&
      spec.strategy != NeighborStrategy::Iterate) {
    std::cerr << "GENUS_MPI: the strategy must be out_prime or iterate\n";
    throw TerminalException{1};
  }
  int n = ListSeed[0].rows();
  GenusEnumerationState<T, Tint, Tgroup> state(n, methods, TotalMass,
                                               spec.hecke, spec.order, os);
  state.distributed = true;
  auto &result = state.result;
  buffered_T_exchanges<Tentry, std::vector<Tentry>> buffer(comm, n_proc,
                                                           tag_entry);
  uint64_t n_sent = 0, n_recv = 0;
  bool final = false;
  // the primes in play, the classes expanded at each
  std::vector<int> &ListPrimeUsed = result.ListPrimeUsed;
  std::vector<std::vector<bool>> is_expanded;
  auto add_prime = [&](int p) -> void {
    ListPrimeUsed.push_back(p);
    result.ListNewClassPerPrime.push_back(0);
    is_expanded.push_back(std::vector<bool>());
  };
  add_prime(spec.strategy == NeighborStrategy::OutPrime ? prime : 2);
  //
  // The keys and the owners
  //
  auto make_entry = [&](MyMatrix<T> const &G, int i_prime,
                        Tmass const &contrib) -> Tentry {
    Tentry e;
    e.i_prime = i_prime;
    e.contrib = contrib;
    e.bv_max_norm = T(0);
    e.bv_n_vect = 0;
    e.bv_hash = 0;
    if (by_canonical) {
      e.Gram =
          GenusCanonicalGram<T, Tint, Tgroup>(G, methods.canonical, os);
    } else {
      e.Gram = G;
      LatticeBVInvariant<T> bv = LatticeBV<T, Tint>(G, methods.bv, os);
      e.bv_max_norm = bv.max_norm;
      e.bv_n_vect = bv.n_vect;
      e.bv_hash = bv.hash;
    }
    return e;
  };
  auto owner = [&](Tentry const &e) -> int {
    size_t h = by_canonical ? GenusMPIHashMatrix(e.Gram)
                            : static_cast<size_t>(e.bv_hash);
    return static_cast<int>(h % static_cast<size_t>(n_proc));
  };
  // The insertion of an entry owned here; the Hecke contribution of a
  // remote expansion is added to the class (the local expansions add theirs
  // themselves, in expand_at).
  auto insert_entry = [&](Tentry const &e,
                          bool add_contrib) -> std::pair<bool, size_t> {
    std::pair<bool, size_t> pair =
        by_canonical
            ? state.insert_canonical(e.Gram)
            : state.insert_with_bv(
                  e.Gram, LatticeBVInvariant<T>{e.bv_max_norm,
                                                static_cast<size_t>(e.bv_n_vect),
                                                e.bv_hash});
    if (add_contrib && state.hecke && e.i_prime >= 0) {
      size_t i = static_cast<size_t>(e.i_prime);
      if (state.hecke_in.size() <= i) {
        state.hecke_in.resize(i + 1);
      }
      state.hecke_in[i].resize(state.ListGramWork.size(), Tmass(0));
      state.hecke_in[i][pair.second] += e.contrib;
    }
    // the new classes of the local expansions are counted by expand_at
    if (add_contrib && pair.first && e.i_prime >= 0) {
      result.ListNewClassPerPrime[static_cast<size_t>(e.i_prime)]++;
    }
    return pair;
  };
  //
  // The coordinator (rank 0)
  //
  std::vector<Tmass> ListMass(n_proc, Tmass(0));
  Tmass last_reported(0);
  bool coord_stop_sent = false;
  bool coord_complete = false;
  bool round_in_progress = false;
  int n_reply = 0;
  std::vector<GenusMPIStatus> ListStatus(n_proc);
  std::optional<std::pair<uint64_t, uint64_t>> prev_idle_counts;
  auto local_busy = [&]() -> bool {
    if (!buffer.is_completely_clear()) {
      return true;
    }
    if (state.stop_requested) {
      return false;
    }
    for (size_t i = 0; i < ListPrimeUsed.size(); i++) {
      for (size_t pos = 0; pos < state.ListGramWork.size(); pos++) {
        if (pos >= is_expanded[i].size() || !is_expanded[i][pos]) {
          return true;
        }
      }
    }
    return false;
  };
  auto local_status = [&]() -> GenusMPIStatus {
    return {local_busy() ? 1 : 0, n_sent, n_recv};
  };
  // The control messages are small and blocking: no buffer has to outlive
  // the call.
  auto send_all = [&](int tag, int val) -> void {
    for (int i_proc = 1; i_proc < n_proc; i_proc++) {
      comm.send(i_proc, tag, val);
    }
  };
  auto coordinator_stop = [&](bool complete) -> void {
    if (coord_stop_sent) {
      return;
    }
    coord_stop_sent = true;
    coord_complete = complete;
#ifdef DEBUG_GENUS_ENUMERATION_MPI
    os << "GENUS_MPI: the coordinator stops, complete=" << complete << "\n";
#endif
    send_all(tag_stop, complete ? 1 : 0);
    state.stop_requested = true;
    prev_idle_counts.reset();
  };
  auto coordinator_check_mass = [&]() -> void {
    if (coord_stop_sent) {
      return;
    }
    Tmass total(0);
    for (auto const &m : ListMass) {
      total += m;
    }
    if (total == TotalMass) {
      coordinator_stop(true);
    }
  };
  //
  // The messages
  //
  auto poll = [&]() -> void {
    while (buffer.clear_one_entry(os)) {
    }
    while (true) {
      boost::optional<boost::mpi::status> st =
          comm.iprobe(boost::mpi::any_source, boost::mpi::any_tag);
      if (!st) {
        break;
      }
      int src = st->source();
      int tag = st->tag();
      if (tag == tag_entry) {
        std::vector<Tentry> v = buffer.recv_message(src);
        for (auto const &e : v) {
          n_recv++;
          if (!state.stop_requested) {
            (void)insert_entry(e, true);
          }
        }
      } else if (tag == tag_mass) {
        Tmass m;
        comm.recv(src, tag_mass, m);
        ListMass[src] = m;
      } else if (tag == tag_ask) {
        int val;
        comm.recv(src, tag_ask, val);
        GenusMPIStatus status = local_status();
        comm.send(0, tag_reply, status);
      } else if (tag == tag_reply) {
        GenusMPIStatus status;
        comm.recv(src, tag_reply, status);
        ListStatus[src] = status;
        n_reply++;
      } else if (tag == tag_prime) {
        int p;
        comm.recv(src, tag_prime, p);
        add_prime(p);
      } else if (tag == tag_stop) {
        int val;
        comm.recv(src, tag_stop, val);
        state.stop_requested = true;
      } else if (tag == tag_final) {
        int val;
        comm.recv(src, tag_final, val);
        final = true;
      }
#ifdef SANITY_CHECK_GENUS_ENUMERATION_MPI
      else {
        std::cerr << "GENUS_MPI: unexpected tag " << tag << "\n";
        throw TerminalException{1};
      }
#endif
    }
    // the mass report
    if (result.accumulated_mass != last_reported) {
      last_reported = result.accumulated_mass;
      if (i_rank == 0) {
        ListMass[0] = last_reported;
      } else {
        comm.send(0, tag_mass, last_reported);
      }
    }
    if (i_rank == 0) {
      coordinator_check_mass();
      if (max_runtime_second > 0 && state.is_finished()) {
        coordinator_stop(false);
      }
    }
  };
  // The neighbors of the expansions: inserted here or sent to their owner.
  state.f_route = [&](MyMatrix<T> const &G, size_t pos, size_t i_prime,
                      size_t orbit_size)
      -> std::optional<std::pair<bool, size_t>> {
    Tmass order_mass =
        UniversalScalarConversion<Tmass, TintGroup>(result.ListAutOrder[pos]);
    Tentry e = make_entry(G, static_cast<int>(i_prime),
                          Tmass(orbit_size) / order_mass);
    int dest = owner(e);
    if (dest == i_rank) {
      return insert_entry(e, false);
    }
    buffer.insert_entry(dest, e);
    n_sent++;
    poll();
    return {};
  };
  // An idle step of the coordinator: an idleness round, then its outcome.
  auto coordinator_idle_step = [&]() -> void {
    if (!round_in_progress) {
      round_in_progress = true;
      n_reply = 0;
      send_all(tag_ask, 0);
      return;
    }
    if (n_reply < n_proc - 1) {
      return;
    }
    round_in_progress = false;
    ListStatus[0] = local_status();
    bool all_idle = true;
    uint64_t sum_sent = 0, sum_recv = 0;
    for (auto const &status : ListStatus) {
      if (status.busy) {
        all_idle = false;
      }
      sum_sent += status.n_sent;
      sum_recv += status.n_recv;
    }
    if (!all_idle || sum_sent != sum_recv) {
      // the nodes are working: wait before the next round
      prev_idle_counts.reset();
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      return;
    }
    std::pair<uint64_t, uint64_t> counts{sum_sent, sum_recv};
    if (!prev_idle_counts || *prev_idle_counts != counts) {
      prev_idle_counts = counts;
      return;
    }
    prev_idle_counts.reset();
    // nothing in flight, nothing to do
    if (coord_stop_sent) {
#ifdef DEBUG_GENUS_ENUMERATION_MPI
      os << "GENUS_MPI: drained, the end\n";
#endif
      send_all(tag_final, 0);
      final = true;
      return;
    }
    if (spec.strategy == NeighborStrategy::Iterate) {
      int p_last =
          *std::max_element(ListPrimeUsed.begin(), ListPrimeUsed.end());
      int p_new = NextPrimeAfter(p_last);
#ifdef DEBUG_GENUS_ENUMERATION_MPI
      os << "GENUS_MPI: every class expanded, the prime " << p_new
         << " comes in\n";
#endif
      send_all(tag_prime, p_new);
      add_prime(p_new);
      return;
    }
    coordinator_stop(false);
  };
  //
  // The seeds, inserted by their owners
  //
  for (auto const &eSeed : ListSeed) {
    Tentry e = make_entry(eSeed, -1, Tmass(0));
    if (owner(e) == i_rank) {
      (void)insert_entry(e, false);
    }
  }
  //
  // The main loop
  //
  while (!final) {
    poll();
    if (final) {
      break;
    }
    if (!state.stop_requested) {
      // the first prime in play with an unexpanded class of this node
      size_t n_class = state.ListGramWork.size();
      std::optional<std::pair<size_t, std::vector<size_t>>> work;
      for (size_t i = 0; i < ListPrimeUsed.size() && !work; i++) {
        is_expanded[i].resize(n_class, false);
        std::vector<size_t> unexpanded;
        for (size_t pos = 0; pos < n_class; pos++) {
          if (!is_expanded[i][pos]) {
            unexpanded.push_back(pos);
          }
        }
        if (!unexpanded.empty()) {
          work = std::make_pair(i, unexpanded);
        }
      }
      if (work) {
        size_t i = work->first;
        int p = ListPrimeUsed[i];
        size_t pos = state.pick(work->second, p, i);
        is_expanded[i][pos] = true;
        state.expand_at(pos, p, i);
        continue;
      }
    }
    // nothing to do here
    if (i_rank == 0) {
      coordinator_idle_step();
      if (round_in_progress && n_reply < n_proc - 1) {
        (void)comm.probe(boost::mpi::any_source, boost::mpi::any_tag);
      }
    } else if (buffer.is_completely_clear()) {
      (void)comm.probe(boost::mpi::any_source, boost::mpi::any_tag);
    }
  }
  // The entries are drained by now (the last rounds found as many received as
  // sent), and the control messages are blocking.
  bool complete = coord_complete;
  boost::mpi::broadcast(comm, complete, 0);
  result.complete = complete;
  return {result, complete};
}

/*
  The classes of all the nodes, gathered on rank 0 into one result (the
  counters summed); on the other ranks the result is empty.
 */
template <typename T, typename Tgroup>
GenusEnumerationResult<T, Tgroup>
GenusMPIGather(boost::mpi::communicator &comm,
               GenusEnumerationResult<T, Tgroup> const &local, bool complete) {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  using TintGroup = typename Tgroup::Tint;
  GenusMPILocalClasses<T, TintGroup> mine{local.ListGram, local.ListAutOrder,
                                          local.n_expansion,
                                          local.n_neighbor_full};
  std::vector<GenusMPILocalClasses<T, TintGroup>> all;
  boost::mpi::gather(comm, mine, all, 0);
  GenusEnumerationResult<T, Tgroup> result = local;
  result.ListGram.clear();
  result.ListAutOrder.clear();
  result.accumulated_mass = Tmass(0);
  result.n_expansion = 0;
  result.n_neighbor_full = 0;
  result.complete = complete;
  for (auto const &c : all) {
    for (size_t i = 0; i < c.ListGram.size(); i++) {
      result.ListGram.push_back(c.ListGram[i]);
      result.ListAutOrder.push_back(c.ListAutOrder[i]);
      result.accumulated_mass +=
          Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(
                         c.ListAutOrder[i]);
    }
    result.n_expansion += c.n_expansion;
    result.n_neighbor_full += c.n_neighbor;
  }
  return result;
}

// clang-format off
#endif  // SRC_GENERA_GENUSENUMERATIONMPI_H_
// clang-format on
