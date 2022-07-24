// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
#define SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_

#include "POLY_RecursiveDualDesc.h"
#include "MPI_functionality.h"
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct message_facet {
  size_t e_hash;
  vectface vf; // List of vectface by the DatabaseBank
};

struct message_query {
  size_t e_hash;
  uint8_t query;
  // 0: for TerminationInfo
  // 1: For destroying the databank and sending back the vectface
};

template <typename Tint> struct HashUndoneOrbitInfo {
  size_t e_hash;
  UndoneOrbitInfo<Tint> erec;
};

namespace boost::serialization {

template <class Archive>
inline void serialize(Archive &ar, message_facet &mesg,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("hash", mesg.e_hash);
  ar &make_nvp("vf", mesg.vf);
}

template <class Archive>
inline void serialize(Archive &ar, message_query &mesg,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("hash", mesg.e_hash);
  ar &make_nvp("query", mesg.query);
}

template <class Archive, typename Tint>
inline void serialize(Archive &ar, HashUndoneOrbitInfo<Tint> &mesg,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("hash", mesg.e_hash);
  ar &make_nvp("nborbitdone", mesg.erec.nbOrbitDone);
  ar &make_nvp("nbundone", mesg.erec.nbUndone);
  ar &make_nvp("setundone", mesg.erec.eSetUndone);
}

} // namespace boost::serialization





/*
template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface MPI_DUALDESC_AdjacencyDecomposition_General(
    Tbank &TheBank, MyMatrix<T> const &EXT, Tgroup const &GRP,
    boost::mpi::communicator &comm,
    PolyHeuristicSerial<typename Tgroup::Tint> const &AllArr,
    std::string const &ePrefix) {
  //  using Tgr = GraphListAdj;
  using Tint = typename Tgroup::Tint;
  int irank = comm.rank();
  int size = comm.size();
  // New Facets (processed asynchronously)
  const int tag_new_facets =
      36; // New facets to be added, the most common request
  // Message query being received
  const int tag_message_query = 37;
  // undone information for Balinski termination
  const int tag_nbundone_balinski = 38;
  // undone information for Balinski termination
  const int tag_terminate_send_vf = 38;
  // undone information for Balinski termination
  const int tag_setup_databank = 39;

  using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
  //  using DataFacet = typename TbasicBank::DataFacet;
  struct BigRecordEntry {
    std::vector<size_t> subset_index_proc;
    DatabaseOrbits<TbasicBank> databank;
    bool did_something;
    std::vector<UndoneOrbitInfo<Tint>> list_undoneinfo;
    int initiating_proc;
  };

  struct SetupDatabank {
    MyMatrix<T> EXT;
    Tgroup GRP;
    size_t e_hash;
  };

  // We can have DatabseBank created on disjoint processes.
  // The order will not be the same between processors.
  std::vector<BigRecordEntry> ListRPL;
  TbasicBank bb(EXT, GRP);
  ListRPL.emplace_back({
      get_subset_index_rev(EXT.rows()),
      DatabaseOrbits<TbasicBank>(bb, ePrefix, AllArr.Saving, std::cerr), false,
      std::vector<UndoneOrbitInfo<Tint>>(size,
                                         get_default_undoneinfo(EXT.rows())),
      -1 // no initiating for the main one
  });
  std::unordered_map<size_t, uint8_t>
      map_databank; // mapping from the hash of database orbit to the index in
                    // ListRPL
  std::vector<size_t> map_databank_rev; // mapping from index of ListRPL to hash
  map_databank[0] = 0;
  map_databank_rev.push_back(0);
  uint8_t selected_pos = 0;
  size_t selected_hash = 0;
  auto get_pos = [&](auto &x) -> uint8_t { return map_databank[x.e_hash]; };
  auto set_selected_pos = [&]() -> void {
    int min_dim = std::numeric_limits<int>::max();
    for (size_t i = 0; i < ListRPL.size(); i++) {
      int dim = ListRPL[i].databank.EXT.cols();
      if (dim < min_dim) {
        selected_pos = i;
        min_dim = dim;
      }
    }
    selected_hash = map_databank_rev[selected_pos];
  };
  auto remove_databank = [&](const size_t &e_hash, const uint8_t &pos) -> void {
    // Correct ListRPL
    ListRPL.erase(pos);
    // Correct map_databank
    map_databank.erase(e_hash);
    for (auto &kv : map_databank) {
      if (kv.second > pos)
        map_databank[kv.first] = kv.second - 1;
    }
    // Correct map_databank
    for (size_t i = pos + 1; i < map_databank_rev.size(); i++)
      map_databank_rev[i - 1] = map_databank_rev[i];
    map_databank_rev.pop_back();
    // Recomputed working_pos
    set_selected_pos();
  };
  auto insert_databank = [&](int init_proc,
                             const SetupDatabank &setup_db) -> void {
    // Update ListRPL
    TbasicBank bb(setup_db.EXT,
                  setup_db.GRP); // Bugged as we are passing a reference to a
                                 // temporary object.
    ListRPL.emplace_back(
        {get_subset_index_rev(setup_db.EXT.rows()),
         DatabaseOrbits<TbasicBank>(bb, ePrefix, AllArr.Saving, std::cerr),
         false,
         std::vector<UndoneOrbitInfo<Tint>>(
             size, get_default_undoneinfo(setup_db.EXT.rows())),
         init_proc});
    uint8_t pos = ListRPL.size();
    // map_databank
    map_databank[setup_db.e_hash] = pos;
    // map_databank_rev
    map_databank_rev.push_back(setup_db.e_hash);
    // Recomputed working_pos
    set_selected_pos();
  };
  // The Buffers in output and receive
  std::vector<message_facet> ListEntries_OUT(size);
  std::vector<message_facet> ListEntries_IN;
  // The queries
  std::vector<std::pair<int, message_query>> ListMesgQuery;
  // The entry nbundone
  std::vector<std::pair<int, HashUndoneOrbitInfo<Tint>>> ListMesgUndone;
  // The partial dual desc (at the end of computation)
  std::vector<std::pair<int, message_facet>> ListMesgPartDualDesc;
  // The bank setup
  std::vector<std::pair<int, SetupDatabank>> ListSetupDatabank;
  // The MPI related stuff
  std::vector<boost::mpi::request> ListRequest;
  std::vector<int> RequestStatus;
  auto GetFreeIndex = [&]() -> size_t {
    size_t len = RequestStatus.size();
    for (size_t i = 0; i < len; i++) {
      if (RequestStatus[i] == 1) {
        boost::optional<boost::mpi::status> stat = ListRequest[i].test();
        if (stat) { // that request has ended. Let's read it.
          RequestStatus[i] = 0;
        }
      }
    }
    for (size_t i = 0; i < len; i++)
      if (RequestStatus[i] == 0)
        return i;
    ListRequest.push_back(boost::mpi::request());
    RequestStatus.push_back(1);
    return len;
  };
  // The infinite loop to do the enumeration
  while (true) {
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      // We only do data reception processing is done afterwards
      if (prob->tag() == tag_new_facets) {
        message_facet e_mesg;
        comm.recv(prob->source(), prob->tag(), e_mesg);
        ListEntries_IN.push_back(e_mesg);
      }
      if (prob->tag() == tag_message_query) {
        message_query e_mesg;
        comm.recv(prob->source(), prob->tag(), e_mesg);
        ListMesgQuery.push_back({prob->source(), e_mesg});
      }
      if (prob->tag() == tag_nbundone_balinski) {
        HashUndoneOrbitInfo<Tint> e_mesg;
        comm.recv(prob->source(), prob->tag(), e_mesg);
        ListMesgUndone.push_back({prob->source(), e_mesg});
      }
      if (prob->tag() == tag_terminate_send_vf) {
        message_facet e_mesg;
        comm.recv(prob->source(), prob->tag(), e_mesg);
        ListMesgPartDualDesc.push_back({prob->source(), std::move(e_mesg)});
      }
      if (prob->tag() == tag_setup_databank) {
        SetupDatabank e_mesg;
        comm.recv(prob->source(), prob->tag(), e_mesg);
        ListSetupDatabank.push_back({prob->source(), std::move(e_mesg)});
      }
    } else {
      // First clearing the buffers of facets
      for (auto &eEntry_IN : ListEntries_IN) {
        uint8_t pos = get_pos(eEntry_IN);
        for (auto &eFace : eEntry_IN.vf)
          ListRPL[pos].databank.FuncInsert(eFace);
      }
      ListEntries_IN.clear();
      // Processing the message queries
      for (auto &eMesgQuery : ListMesgQuery) {
        uint8_t pos = get_pos(eMesgQuery.second);
        if (eMesgQuery.second.query == 0) {
          UndoneOrbitInfo<Tint> recundone =
              ListRPL[pos].databank.GetTerminationInfo();
          HashUndoneOrbitInfo<Tint> hashrecundone{eMesgQuery.second.e_hash,
                                                  recundone};
          size_t u = GetFreeIndex();
          ListRequest[u] = comm.isend(0, tag_nbundone_balinski, hashrecundone);
          RequestStatus[u] = 1;
        }
        if (eMesgQuery.second.query == 1) {
          message_facet e_mesg{eMesgQuery.second.e_hash,
                               ListRPL[pos].databank.FuncListOrbitIncidence()};
          size_t u = GetFreeIndex();
          ListRequest[u] = comm.isend(0, tag_terminate_send_vf, e_mesg);
          RequestStatus[u] = 1;
          remove_databank(eMesgQuery.second.e_hash, pos);
        }
      }
      ListMesgQuery.clear();
      // Update information regarding the undone status of orbits
      for (auto &e_ent : ListMesgUndone) {
        int jrank = e_ent.first;
        uint8_t pos = get_pos(e_ent.second);
        ListRPL[pos].list_undoneinfo[jrank] = e_ent.second.erec;
      }
      ListMesgUndone.clear();
      // Create new databanks as required
      for (auto &e_ent : ListSetupDatabank) {
        int jrank = e_ent.first;
        insert_databank(jrank, e_ent.second);
      }
      ListSetupDatabank.clear();
      // Now treating the selected_block
      if (ListRPL[selected_pos].did_something) {
        ListRPL[selected_pos].did_something = false;
        if (irank == 0) {
          ListRPL[selected_pos].list_undoneinfo[0] =
              ListRPL[selected_pos].databank.GetTerminationInfo();
          if (MonotonicCheckStatusUndone(
                  ListRPL[selected_pos].list_undoneinfo[0],
                  ListRPL[selected_pos].databank.CritSiz)) {
            UndoneOrbitInfo<Tint> undoneinfo =
                CombineUndoneOrbitInfo(ListRPL[selected_pos].list_undoneinfo);
            if (ComputeStatusUndone(undoneinfo)) {
              if (selected_pos ==
                  0) { // This is the main databank. Work is finished
                break;
              } else { // Not the main one, more work needed. Relocation to the
                       // initiating processor
              }
            } else {
              message_query mq{selected_hash, 0};
              for (int jrank = 1; jrank < size; jrank++) {
                size_t u = GetFreeIndex();
                ListRequest[u] = comm.isend(0, tag_message_query, mq);
                RequestStatus[u] = 1;
              }
            }
          }
        }
      } else {
        ListRPL[selected_pos].did_something = true;
      }
    }
  }
  return ListRPL[0].FuncListOrbitIncidence();
}
*/




vectface mpi_gather(vectface const& vf, boost::mpi::communicator & comm, int i_proc) {
  int i_rank = comm.rank();
  size_t n_vert = vf.get_n();
  size_t n_face = vf.size();
  std::vector<uint8_t> const& V = vf.serial_get_std_vector_uint8_t();
  //
  std::vector<size_t> l_n_face = mpi_gather(n_face, comm, i_proc);
  std::vector<std::vector<uint8_t>> l_V = mpi_gather(V, comm, i_proc);
  if (i_rank == i_proc) {
    return vectface(n_vert, l_n_face, l_V);
  } else {
    return vectface(n_vert);
  }
}



vectface mpi_allgather(vectface const& vf, boost::mpi::communicator & comm) {
  size_t n_vert = vf.get_n();
  size_t n_face = vf.size();
  std::vector<uint8_t> const& V = vf.serial_get_std_vector_uint8_t();
  //
  std::vector<size_t> l_n_face = mpi_allgather(comm, n_face, i_proc);
  std::vector<std::vector<uint8_t>> l_V = mpi_allgather(V, comm, i_proc);
  return vectface(n_vert, l_n_face, l_V);
}




/*
  This is the code for the MPI parallelization.
  Since the first attempt (see above) was maybe too complicated.
  --
  Designs:
  ---At the beginning, read the existing databases.
     ---Compute the total number of existing entries with a mpi.allreduce.
     ---If zero, then all the process do sampling and insert into databases.
     ---If not zero, read it and start.
  ---A needed signal of balinski termination has to be computed from time to time.
  ---A termination check (from runtime or Ctrl-C) has to be handled.
  ---Computation of existing with a call to the serial code. Should be modelled
  on the C-type code.


 */
template <typename Tbank, typename T, typename Tgroup, typename Tidx_value>
vectface MPI_Kernel_DUALDESC_AdjacencyDecomposition(
    Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> const &AllArr,
    boost::mpi::communicator &comm,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> const &TheMap) {
  using DataFacet = typename TbasicBank::DataFacet;
  //  using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
  SingletonTime start;
  int irank = comm.rank();
  int n_proc = comm.size();
  std::ostream& os = get_standard_outstream(comm);
  DatabaseOrbits<TbasicBank> RPL(bb, ePrefix, AllArr.Saving, os);
  int n_vert = bb.nbRow;
  //
  // The types of exchanges
  //
  // New facets to be added, the most common request
  const int tag_new_facets = 36;
  // tag_balinski
  const int tag_balinski_info = 38;
  // undone information for Balinski termination
  const int tag_termination = 38;
  // undone information for Balinski termination
  const int tag_requestinfo_balinski = 39;
  const int tag_balinski_info = 39;
  //
  // Reading the input
  //
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");

  std::vector<std::optional<UndoneOrbitInfo<Tint>>> ListBalinski(n_proc);
  int MaxFly = 4 * n_proc;


  //
  // The parallel MPI classes
  //
  empty_message_management emm_termin(comm, 0, tag_termination);
  empty_message_management emm_balinski(comm, 0, tag_requestinfo_balinski);
  buffered_T_exchanges<Face,vectface> bte_facet(comm, MaxFly, tag_new_facets, vectface(n_vert));


  std::vector<int> StatusNeighbors(n_proc, 0);
  auto fInsertUnsent = [&](Face const &face) -> void {
    int res = IntegerDiscriminantInvariant(face) % n_proc;
    if (res == irank) {
      RPL.FuncInsert(face);
    } else {
      bte_facet.insert_entry(res, face);
    }
  };
  //
  // Initial invocation of the synchronization code
  //
  size_t n_orb_tot = 0, n_orb_loc = RPL.FuncNumberOrbit();
  all_reduce(comm, n_orb_loc, n_orb_tot, mpi::minimum<size_t>());
  if (n_orb_tot == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    for (auto &face : RPL.ComputeInitialSet(ansSamp))
      fInsertUnsent(face);
  }
  //
  // The infinite loop
  //
  while (true) {
    bool SendTerminationCriterion = true;
    bool StopComputation = MaxRunTimeSecond > 0 && si(start) > MaxRunTimeSecond;
    bool NeedComputeBalinski = false;
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      if (prob->tag() == tag_new_facets) {
        StatusNeighbors[prob->source()] = 0;
        vectface l_recv_face = bte_facet.recv_message(prob->source())
        for (auto & face : l_recv_face)
          RPL.FuncInsert(face);
      }
      if (prob->tag() == tag_termination) {
        StatusNeighbors[prob->source()] = 1;
        emm_termin.recv_message(prob->source());
      }
      if (prob->tag() == tag_requestinfo_balinski) {
        int recv_message;
        emm_balinski.recv_message(prob->source());
        UndoneOrbitInfo<Tint> uoi = RPL.GetTerminationInfo();
      }
      if (prob->tag() == tag_balinski_info) {
        UndoneOrbitInfo<Tint> uoi
        comm.recv(prob->source(), prob->tag(), recv_message);
        ListBalinski[prob->source()] = uoi;
      }
    } else {
      if (!StopComputation) {
        DataFacet df = RPL.FuncGetMinimalUndoneOrbit();
        size_t SelectedOrbit = df.SelectedOrbit;
        std::string NewPrefix =
          ePrefix + "ADM" + std::to_string(SelectedOrbit) + "_";
        vectface TheOutput =
          DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(TheBank, df.FF.EXT_face, df.Stab, AllArr, NewPrefix);
        for (auto &eOrbB : TheOutput) {
          Face eFlip = df.flip(eOrbB);
          fInsertUnsent(eFlip);
        }
      }
    }
    os << "Before ClearUnsentAsPossible\n";
    bte_facet.clear_one_entry();
    os << "After ClearUnsentAsPossible\n";
    //
    // Determine Balinski stuff
    //
    if (NeedComputeBalinski) {
      auto get_value=[&]() -> bool {
        std::vector<UndoneOrbitInfo<Tint>> ListPart;
        for (auto & eComp : ListBalinski) {
          if (opt) {
            ListPart.push_back(*opt);
          } else {
            return false;
          }
        }
        UndoneOrbitInfo<Tint> uoi = CombineUndoneOrbitInfo(ListPart);
        return ComputeStatusUndone(uoi, CritSiz);
      };
      SendTerminationCriterion = get_value();
    }
    //
    // Sending termination criterion
    //
    if (SendTerminationCriterion) {
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        if (i_proc == irank) {
          StatusNeighbors[irank] = 1;
        } else {
          emm_termin.send_message(i_proc);
        }
      }
    }
    //
    // Termination criterion
    //
    int nb_finished = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++)
      nb_finished += StatusNeighbors[i_pes];
    if (nb_finished == n_proc)
      break;
  }
  return RPL.FuncListOrbitIncidence();
}










// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
// clang-format on
