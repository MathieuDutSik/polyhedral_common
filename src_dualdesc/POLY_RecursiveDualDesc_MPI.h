// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
#define SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_

#include "POLY_RecursiveDualDesc.h"
#include "MPI_functionality.h"
#include "Balinski_basic.h"
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
template <typename Tbank, typename TbasicBank, typename T, typename Tgroup, typename Tidx_value>
vectface MPI_Kernel_DUALDESC_AdjacencyDecomposition(
    boost::mpi::communicator &comm,
    Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> const &TheMap, std::ostream& os) {
  using DataFacet = typename TbasicBank::DataFacet;
  using Tint = typename TbasicBank::Tint;
  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string lPrefix = ePrefix + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  DatabaseOrbits<TbasicBank> RPL(bb, lPrefix, AllArr.Saving, AllArr.AdvancedTerminationCriterion, os);
  Tint CritSiz = RPL.CritSiz;
  bool LocalBalinskiStatus = false;
  UndoneOrbitInfo<Tint> uoi_local;
  os << "DirectFacetOrbitComputation, step 1\n";
  int n_vert = bb.nbRow;
  int n_vert_div8 = (n_vert + 7) / 8;
  std::vector<uint8_t> V_hash(n_vert_div8,0);
  auto get_hash=[&](Face const& x) -> size_t {
    for (int i_vert=0; i_vert<n_vert; i_vert++) {
      setbit(V_hash, i_vert, x[i_vert]);
    }
    uint32_t seed = 0x1b873560;
    return robin_hood_hash_bytes(V_hash.data(), n_vert_div8, seed);
  };
  //
  // The types of exchanges
  //
  // New facets to be added, the most common request
  const int tag_new_facets = 36;
  // undone information for Balinski termination
  const int tag_termination = 38;
  // undone information for Balinski termination
  const int tag_balinski_request = 39;
  const int tag_balinski_info = 40;
  //
  // Reading the input
  //
  size_t MaxBuffered = 10000 * n_proc;
  int MaxFly = 4 * n_proc;


  //
  // The parallel MPI classes
  //
  os << "DirectFacetOrbitComputation, step 2\n";
  empty_message_management emm_termin(comm, 0, tag_termination);
  os << "DirectFacetOrbitComputation, step 3\n";
  buffered_T_exchanges<Face,vectface> bte_facet(comm, MaxFly, tag_new_facets);
  os << "DirectFacetOrbitComputation, step 4\n";
  database_balinski_info dbi(comm, tag_balinski_request, tag_balinski_info, CritSiz);
  os << "DirectFacetOrbitComputation, step 5\n";

  std::vector<int> StatusNeighbors(n_proc, 0);
  auto fInsertUnsent = [&](Face const &face) -> void {
    int res = int(get_hash(face) % size_t(n_proc));
    os << "|face|=" << face.count() << " res=" << res << "\n";
    if (res == i_rank) {
      RPL.FuncInsert(face);
    } else {
      bte_facet.insert_entry(res, face);
    }
  };
  //
  // Initial invocation of the synchronization code
  //
  os << "Compute initital\n";
  size_t n_orb_tot = 0, n_orb_loc = RPL.FuncNumberOrbit();
  all_reduce(comm, n_orb_loc, n_orb_tot, boost::mpi::maximum<size_t>());
  os << "n_orb_loc=" << n_orb_loc << " n_orb_tot=" << n_orb_tot << "\n";
  if (n_orb_tot == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    os << "ansSamp=" << ansSamp << "\n";
    for (auto &face : RPL.ComputeInitialSet(ansSamp, os)) {
      fInsertUnsent(face);
    }
  }
  os << "DirectFacetOrbitComputation, step 6\n";
  //
  // The infinite loop
  //
  auto balinski_flush=[&]() -> void {
    os << "Beginning of balinski_flush\n";
    uoi_local = RPL.GetTerminationInfo();
    LocalBalinskiStatus = ComputeStatusUndone(uoi_local, CritSiz);
    os << "We have uoi_local=" << uoi_local << " LocalBalinskiStatus=" << LocalBalinskiStatus << "\n";
    if (LocalBalinskiStatus) {
      dbi.set_uoi_local(uoi_local);
      os << "last_database_update has been set\n";
    }
  };
  auto process_mpi_status=[&](boost::mpi::status const& stat) -> void {
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_new_facets) {
      os << "RECV of tag_new_facets from " << e_src << "\n";
      StatusNeighbors[e_src] = 0;
      vectface l_recv_face = bte_facet.recv_message(e_src);
      os << "|l_recv_face|=" << l_recv_face.size() << "\n";
      for (auto & face : l_recv_face)
        RPL.FuncInsert(face);
      balinski_flush();
    }
    if (e_tag == tag_termination) {
      os << "RECV of tag_termination from " << e_src << "\n";
      StatusNeighbors[e_src] = 1;
      emm_termin.recv_message(e_src);
    }
    if (e_tag == tag_balinski_request) {
      os << "RECV of tag_balinski_request from " << e_src << "\n";
      dbi.read_request(e_src);
      dbi.flush(os);
      os << "dbi.flush of uoi_local=" << uoi_local << "\n";
    }
    if (e_tag == tag_balinski_info) {
      os << "RECV of tag_balinski_info from " << e_src << "\n";
      dbi.recv_info(e_src);
    }
    os << "Exiting process_mpi_status\n";
  };
  auto process_database=[&]() -> void {
    os << "process_database, begin\n";
    DataFacet df = RPL.FuncGetMinimalUndoneOrbit();
    os << "process_database, we have df\n";
    size_t SelectedOrbit = df.SelectedOrbit;
    std::string NewPrefix =
      ePrefix + "PROC" + std::to_string(i_rank) + "_ADM" + std::to_string(SelectedOrbit) + "_";
    try {
      os << "Before call to DUALDESC_AdjacencyDecomposition\n";
      vectface TheOutput =
        DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(TheBank, df.FF.EXT_face, df.Stab, AllArr, NewPrefix, os);
      os << "We have TheOutput, |TheOutput|=" << TheOutput.size() << "\n";
      for (auto &eOrbB : TheOutput) {
        Face eFlip = df.flip(eOrbB);
        fInsertUnsent(eFlip);
      }
      RPL.FuncPutOrbitAsDone(SelectedOrbit);
    } catch (RuntimeException const &e) {
      os << "The computation of DUALDESC_AdjacencyDecomposition has ended by runtime exhaustion\n";
    }
    os << "process_database, EXIT\n";
    balinski_flush();
  };
  while (true) {
    os << "DirectFacetOrbitComputation, inf loop, start\n";
    bool HasReachedBalinskiConclusion = false;
    bool MaxRuntimeReached = AllArr.max_runtime > 0 && si(start) > AllArr.max_runtime;
    bool SomethingToDo = !MaxRuntimeReached && !RPL.IsFinished();
    os << "DirectFacetOrbitComputation, HasRechedBalinskiConclusion=" << HasReachedBalinskiConclusion << " MaxRuntimeReached=" << MaxRuntimeReached << " SomethingToDo=" << SomethingToDo << "\n";
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      os << "prob is not empty\n";
      process_mpi_status(*prob);
    } else {
      if (SomethingToDo) {
        os << "Case something to do\n";
        // we have to clear our buffers sometimes while running
        // otherwise we will only do it in the very end
        // which could lead to extremely large buffers
        // and possibly treating orbits with unnecessary large incidence
        if(bte_facet.get_unsent_size() >= MaxBuffered) {
            os << "Calling clear_one_entry after reaching MaxBuffered\n";
            bte_facet.clear_one_entry(os);
            StatusNeighbors[i_rank] = 0;
        }

        process_database();
      } else {
        if (!bte_facet.is_buffer_empty()) {
          os << "Calling clear_one_entry\n";
          bte_facet.clear_one_entry(os);
          StatusNeighbors[i_rank] = 0;
        } else {
          os << "Nothing to do, so we do a blocking wait (This avoids busy wait)\n";
          boost::mpi::status stat = comm.probe();
          process_mpi_status(stat);
        }
      }
    }
    //
    bool test_buffer = bte_facet.is_buffer_empty();
    os << "LocalBalinskiStatus=" << LocalBalinskiStatus << " test_buffer=" << test_buffer << "\n";
    if (LocalBalinskiStatus && test_buffer) {
      os << "Trying our chance at terminating the run\n";
      dbi.flush(os);
      os << "dbi.flush has been called\n";
      if (!HasReachedBalinskiConclusion) {
        dbi.submit_request_uoi(os);
        os << "dbi.submit_request_uoi has been called\n";
      }
    }
    //
    // Determine Balinski stuff
    //
    HasReachedBalinskiConclusion = dbi.get_status(os);
    os << "Currently HasReachedBalinskiConclusion=" << HasReachedBalinskiConclusion << "\n";
    //
    // Sending termination criterion
    //
    if ((MaxRuntimeReached || HasReachedBalinskiConclusion) && StatusNeighbors[i_rank]==0) {
      os << "Sending messages for terminating the run\n";
      for (int i_proc = 0; i_proc < n_proc; i_proc++) {
        if (i_proc == i_rank) {
          StatusNeighbors[i_rank] = 1;
        } else {
          os << "Before send_message at i_proc=" << i_proc << "\n";
          emm_termin.send_message(i_proc);
          os << "After send_message\n";
        }
      }
    }
    //
    // Termination criterion
    //
    int nb_finished = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      os << "i_proc=" << i_proc << " status=" << StatusNeighbors[i_proc] << "\n";
      nb_finished += StatusNeighbors[i_proc];
    }
    os << "nb_finished=" << nb_finished << "\n";
    if (nb_finished == n_proc)
      break;
    os << "End of the while loop, continuing\n";
  }
  os << "We just exited the infinite loop\n";
  return RPL.FuncListOrbitIncidence();
}



template<typename T>
void Reset_Directories(boost::mpi::communicator & comm, PolyHeuristicSerial<T> & AllArr) {
  int n_proc = comm.size();
  int i_rank = comm.rank();
  std::string postfix = "_nproc" + std::to_string(n_proc) + "_rank" + std::to_string(i_rank);
  auto update_string=[&](std::string & str_ref) -> void {
    size_t len = str_ref.size();
    std::string part1 = str_ref.substr(0,len-1);
    std::string part2 = str_ref.substr(len-1,1);
    if (part2 != "/") {
      std::cerr << "str_ref=" << str_ref << "\n";
      std::cerr << "Last character should be a /\n";
      throw TerminalException{1};
    }
    str_ref = part1 + postfix + part2;
    CreateDirectory(str_ref);
  };
  update_string(AllArr.BANK_Prefix);
  update_string(AllArr.DD_Prefix);
}




template <typename T, typename Tgroup, typename Tidx_value>
void MPI_MainFunctionDualDesc(boost::mpi::communicator & comm, FullNamelist const &eFull) {
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tkey = MyMatrix<T>;
  using Tval = PairStore<Tgroup>;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  //
  std::string FileLog = "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::cerr << "We have moved. See the log in FileLog=" << FileLog << "\n";
  std::ofstream os(FileLog);
  os << std::unitbuf;
  //  std::ostream& os = std::cerr;
  os << "Initial writing of the log\n";
  os.flush();
  //
  MyMatrix<T> EXT = Get_EXT_DualDesc<T,Tidx>(eFull, os);
  Tgroup GRP = Get_GRP_DualDesc<Tgroup>(eFull, os);
  PolyHeuristicSerial<Tint> AllArr = Read_AllStandardHeuristicSerial<Tint>(eFull, os);
  Reset_Directories(comm, AllArr);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  //
  using Tbank = DataBankClient<Tkey, Tval>;
  Tbank TheBank(AllArr.port);

  using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
  TbasicBank bb(EXTred, GRP);
  std::map<std::string, Tint> TheMap = ComputeInitialMap<Tint>(EXTred, GRP);
  vectface vf = MPI_Kernel_DUALDESC_AdjacencyDecomposition<Tbank, TbasicBank, T, Tgroup, Tidx_value>(comm, TheBank, bb, AllArr, AllArr.DD_Prefix, TheMap, os);
  os << "We have vf\n";
  int i_proc_ret = 0;
  vectface vf_tot = my_mpi_gather(comm, vf, i_proc_ret);
  os << "We have vf_tot\n";
  if (comm.rank() == i_proc_ret)
    OutputFacets(vf_tot, AllArr.OUTfile, AllArr.OutFormat);
  os << "We have done our output\n";
}









// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
// clang-format on

