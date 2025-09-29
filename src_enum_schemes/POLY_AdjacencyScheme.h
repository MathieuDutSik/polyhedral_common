// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
#define SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_

// clang-format off
#include "boost_serialization.h"
#include "basic_datafile.h"
#include "Timings.h"
#include <unordered_map>
#include <utility>
#include <vector>
#include <set>
#include <functional>
#include <string>
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
        for Niemeier) then this does not work.
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
#ifdef TIMINGS_ADJACENCY_SCHEME
      MicrosecondTime time_f_repr;
#endif
      std::optional<TadjO> opt = f_repr(y, x_adjI, idx);
#ifdef TIMINGS_ADJACENCY_SCHEME
      os << "|ADJ_SCH: f_repr|=" << time_f_repr << "\n";
#endif
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
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time_f_spann;
#endif
    std::pair<Tobj, TadjO> pair = f_spann(x_adjI, n_obj);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "|ADJ_SCH: f_spann|=" << time_f_spann << "\n";
#endif
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
    n_obj += 1;
    return pair.second;
  };
  auto insert_load = [&](Tobj const &x, bool const &is_treated) -> void {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: Before f_hash\n";
#endif
    size_t hash_hashmap = f_hash(seed_hashmap, x);
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: After f_hash\n";
#endif
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
#ifdef TIMINGS_ADJACENCY_SCHEME
    MicrosecondTime time;
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry beginning\n";
#endif
    size_t idx = get_undone_idx();
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry idx=" << idx << "\n";
#endif
    bool is_treated = true;
    f_save_status(idx, is_treated);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "|ADJ_SCH: get_undone_idx / f_save_status|=" << time << "\n";
#endif
    std::vector<TadjI> l_adj_i = f_adj(idx);
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "|ADJ_SCH: f_adj|=" << time << "\n";
#endif
    std::vector<TadjO> l_adj_o;
    for (auto &adj_i : l_adj_i) {
      TadjO adj_o = process_singleEntry_AdjI(adj_i);
      l_adj_o.push_back(adj_o);
    }
#ifdef TIMINGS_ADJACENCY_SCHEME
    os << "|ADJ_SCH: l_adj_o|=" << time << "\n";
#endif
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: treat_one_entry |l_adj_o|=" << l_adj_o.size() << "\n";
#endif
    f_set_adj(idx, l_adj_o);
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
    os << "ADJ_SCH: after f_insert test=" << test << "\n";
#endif
    if (test) {
      early_termination = true;
    }
  }
  while (true) {
#ifdef DEBUG_ADJACENCY_SCHEME
    os << "ADJ_SCH: early_termination=" << early_termination << " n_obj=" << n_obj << " |undone|=" << undone.size() << "\n";
#endif
    if (early_termination) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: returning false due to early_termination = true\n";
#endif
      return false;
    }
    if (undone.size() == 0) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: returning true due to |undone| = 0\n";
#endif
      return true;
    }
    if (max_time_second > 0 && si(start) > max_time_second) {
#ifdef DEBUG_ADJACENCY_SCHEME
      os << "ADJ_SCH: returning false due to si(start) > max_time_second > 0\n";
#endif
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

template <typename TadjO> struct AdjO_Serial {
  TadjO x;
  int iOrb;
};

template <typename TadjO>
void WriteEntryGAP(std::ostream &os_out, AdjO_Serial<TadjO> const &adj) {
  os_out << "rec(x:=";
  WriteEntryGAP(os_out, adj.x);
  os_out << ", iOrb:=" << adj.iOrb << ")";
}

template <typename TadjO>
void WriteEntryPYTHON(std::ostream &os_out, AdjO_Serial<TadjO> const &adj) {
  os_out << "{\"x\":";
  WriteEntryPYTHON(os_out, adj.x);
  os_out << ", \"iOrb\":" << adj.iOrb << "}";
}

template <typename Tobj, typename TadjO> struct DatabaseEntry_Serial {
  Tobj x;
  std::vector<AdjO_Serial<TadjO>> ListAdj;
};

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
  os_out << "])";
}

template <typename Tobj, typename TadjO>
void WriteEntryPYTHON(std::ostream &os_out,
                   DatabaseEntry_Serial<Tobj, TadjO> const &dat_entry) {
  os_out << "{\"x\":";
  WriteEntryPYTHON(os_out, dat_entry.x);
  os_out << ", \"ListAdj\":[";
  bool IsFirst = true;
  for (auto &eAdj : dat_entry.ListAdj) {
    if (!IsFirst)
      os_out << ",";
    IsFirst = false;
    WriteEntryPYTHON(os_out, eAdj);
  }
  os_out << "]}";
}

namespace boost::serialization {
template <class Archive, typename TadjO>
inline void serialize(Archive &ar, AdjO_Serial<TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("iOrb", eRec.iOrb);
}
template <class Archive, typename Tobj, typename TadjO>
inline void serialize(Archive &ar, DatabaseEntry_Serial<Tobj, TadjO> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("x", eRec.x);
  ar &make_nvp("ListAdj", eRec.ListAdj);
}
} // namespace boost::serialization

template <typename Tdata, typename Fincorrect>
std::vector<DatabaseEntry_Serial<typename Tdata::Tobj, typename Tdata::TadjO>>
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
    for (auto &status : l_status) {
      n_sum += static_cast<size_t>(status);
    }
    os << "ADJ_SCH: EnumerateAndStore_Serial: |l_status|=" << l_status.size()
       << " n_sum=" << n_sum << " test=" << test << "\n";
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
  (void)compute_adjacency_serial<
      Tobj, TadjI, TadjO_work, decltype(f_next), decltype(f_insert),
      decltype(f_adji_obj), decltype(f_idx_obj), decltype(f_save_status),
      decltype(f_init), decltype(f_adj), decltype(f_set_adj), decltype(f_hash),
      decltype(f_repr), decltype(f_spann)>(
      max_runtime_second, f_next, f_insert, f_adji_obj, f_idx_obj,
      f_save_status, f_init, f_adj, f_set_adj, f_hash, f_repr, f_spann, os);
  return l_obj;
}

template <typename TadjO>
std::vector<int> UniqueAdjacencies(std::vector<AdjO_Serial<TadjO>> const& ListAdj, int const& offset) {
  std::set<int> set;
  for (auto &eAdj : ListAdj) {
    int new_val = eAdj.iOrb + offset;
    set.insert(new_val);
  }
  std::vector<int> new_v;
  for (auto & new_val: set) {
    new_v.push_back(new_val);
  }
  return new_v;
}

template <typename Tdata, typename Tobj, typename TadjO>
bool WriteFamilyObjects(
    Tdata const& data,
    std::string const &OutFormat,
    std::ostream& os_out,
    std::vector<DatabaseEntry_Serial<Tobj, TadjO>> const &l_tot,
    std::ostream &os) {
  if (OutFormat == "nothing") {
    std::cerr << "No output\n";
    return false;
  }
  if (OutFormat == "NumberGAP") {
    size_t len = l_tot.size();
    os_out << "return rec(nb:=" << len << ");\n";
    return false;
  }
  if (OutFormat == "ObjectGAP") {
    os_out << "return [";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",\n";
      WriteEntryGAP(os_out, l_tot[i].x);
    }
    os_out << "];\n";
    return false;
  }
  if (OutFormat == "DetailedObjectGAP") {
    os_out << "return rec(ListEntry:=[\n";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",\n";
      WriteDetailedEntryGAP(os_out, data, l_tot[i].x, os);
    }
    os_out << "], n_obj:=" << len << ");\n";
    return false;
  }
  if (OutFormat == "ObjectPYTHON") {
    os_out << "[";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",";
      WriteEntryPYTHON(os_out, l_tot[i].x);
    }
    os_out << "]";
    return false;
  }
  if (OutFormat == "ObjectFullAdjacencyGAP") {
    os_out << "return [";
    // INCORRECT CODE
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",\n";
      WriteEntryGAP(os_out, l_tot[i]);
    }
    os_out << "];\n";
    return false;
  }
  if (OutFormat == "ObjectFullAdjacencyPYTHON") {
    os_out << "[";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",";
      WriteEntryPYTHON(os_out, l_tot[i]);
    }
    os_out << "]";
    return false;
  }
  if (OutFormat == "ObjectReducedAdjacencyGAP") {
    os_out << "return [";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",\n";
      os_out << "rec(obj:=";
      WriteEntryGAP(os_out, l_tot[i]);
      os_out << ", LAdj:=";
      std::vector<int> v = UniqueAdjacencies(l_tot[i].ListAdj, 1);
      os_out << StringStdVectorGAP(v) << ")";
    }
    os_out << "];\n";
    return false;
  }
  if (OutFormat == "ObjectReducedAdjacencyPYTHON") {
    os_out << "[";
    size_t len = l_tot.size();
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os_out << ",";
      os_out << "{\"obj\":";
      WriteEntryPYTHON(os_out, l_tot[i]);
      os_out << ", \"LAdj\":";
      std::vector<int> v = UniqueAdjacencies(l_tot[i].ListAdj, 0);
      os_out << StringStdVectorPYTHON(v) << "}";
    }
    os_out << "]";
    return false;
  }
  if (OutFormat == "AdjacencyGAP") {
    size_t len = l_tot.size();
    os_out << "return [";
    for (size_t i = 0; i < len; i++) {
      if (i > 0) {
        os_out << ",\n";
      }
      std::vector<int> v = UniqueAdjacencies(l_tot[i].ListAdj, 1);
      os_out << StringStdVectorGAP(v);
    }
    os_out << "];\n";
    return false;
  }
  if (OutFormat == "AdjacencyPYTHON") {
    size_t len = l_tot.size();
    os_out << "[";
    for (size_t i = 0; i < len; i++) {
      if (i > 0) {
        os_out << ",";
      }
      std::vector<int> v = UniqueAdjacencies(l_tot[i].ListAdj, 0);
      os_out << StringStdVectorPYTHON(v);
    }
    os_out << "]";
    return false;
  }
  return true;
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_ADJACENCYSCHEME_H_
// clang-format on
