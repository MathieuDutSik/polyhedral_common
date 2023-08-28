// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
#define SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_

// clang-format off
#include "POLY_RecursiveDualDesc.h"
#include "Databank_mpi.h"
#include "MPI_functionality.h"
#include "Balinski_basic.h"
#include <chrono>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#include <map>
#include <atomic>
// clang-format on

#ifdef TIMINGS
# define TIMINGS_RECURSIVE_DUAL_DESC_MPI
#endif

struct message_facet {
  size_t e_hash;
  // vf: List of vectface by the DatabaseBank
  vectface vf;
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

// clang-format off
}  // namespace boost::serialization
// clang-format on

size_t mpi_get_hash_kernel(Face const &x, size_t const &n_vert,
                           int const &n_vert_div8,
                           std::vector<uint8_t> &V_hash) {
  for (size_t i_vert = 0; i_vert < n_vert; i_vert++)
    setbit_vector(V_hash, i_vert, x[i_vert]);
  uint32_t seed = 0x1b873560;
  return robin_hood_hash_bytes(V_hash.data(), n_vert_div8, seed);
}

size_t mpi_get_hash(Face const &x, size_t const &n_vert) {
  int n_vert_div8 = (n_vert + 7) / 8;
  std::vector<uint8_t> V_hash(n_vert_div8, 0);
  return mpi_get_hash_kernel(x, n_vert, n_vert_div8, V_hash);
}

// When we upgrade the canonicalization scheme, we need to recompute the hashes
// and redistribute the data.
vectface mpi_shuffle(boost::mpi::communicator &comm, vectface &&vf,
                     size_t const &n) {
  int n_proc = comm.size();
  size_t n_face = vf.size();
  size_t siz_f = vf.get_n(); // the full size of the entries
  int n_div8 = (n + 7) / 8;
  std::vector<uint8_t> V_hash(n_div8, 0);
  std::vector<std::vector<uint8_t>> list_in_vect(n_proc);
  std::vector<size_t> list_in_cnt(n_proc, 0);
  for (size_t i_face = 0; i_face < n_face; i_face++) {
    Face f = vf[i_face];
    size_t hash = mpi_get_hash_kernel(f, n, n_div8, V_hash);
    int res = static_cast<int>(hash % size_t(n_proc));
    size_t pos = list_in_cnt[res];
    size_t needed_len = (siz_f * (pos + 1) + 7) / 8;
    size_t curr_len = list_in_vect[res].size();
    for (size_t u = curr_len; u < needed_len; u++) {
      list_in_vect[res].push_back(0);
    }
    size_t pos_vect = pos * siz_f;
    for (size_t u = 0; u < siz_f; u++) {
      bool val = f[u];
      setbit_vector(list_in_vect[res], pos_vect, val);
      pos_vect++;
    }
    list_in_cnt[res] = pos + 1;
  }
  std::vector<std::vector<uint8_t>> list_out_vect(n_proc);
  std::vector<size_t> list_out_cnt(n_proc);
  boost::mpi::all_to_all(comm, list_in_cnt, list_out_cnt);
  boost::mpi::all_to_all(comm, list_in_vect, list_out_vect);
  vectface vf_ret(siz_f);
  Face f(siz_f);
  for (int i_proc = 0; i_proc < n_proc; i_proc++) {
    size_t cnt_out = list_out_cnt[i_proc];
    size_t pos = 0;
    std::vector<uint8_t> const &V = list_out_vect[i_proc];
    for (size_t i_face = 0; i_face < cnt_out; i_face++) {
      for (size_t u = 0; u < siz_f; u++) {
        bool val = getbit_vector(V, pos);
        f[u] = val;
        pos++;
      }
      vf_ret.push_back(f);
    }
  }
  return vf_ret;
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
  ---A needed signal of balinski termination has to be computed from time to
  time.
  ---A termination check (from runtime or Ctrl-C) has to be handled.
  ---Computation of existing with a call to the serial code. Should be modelled
  on the C-type code.
 */
template <typename Tbank, typename T, typename Tgroup, typename Tidx_value,
          typename TbasicBank, typename Finsert, typename Fcomm>
void DUALDESC_AdjacencyDecomposition_and_insert_commthread(
    Tbank &TheBank, TbasicBank &bb, DataFacet<T, Tgroup> & df,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr, Finsert f_insert,
    Fcomm f_comm, std::string const &ePrefix, std::ostream &os) {
  using Tint = typename Tgroup::Tint;
  CheckTermination<Tgroup>(AllArr);
  std::map<std::string, Tint> TheMap =
      ComputeInitialMap<Tint>(df.FF.EXT_face, df.Stab, AllArr);
  std::string ansSplit = HeuristicEvaluation(TheMap, AllArr.Splitting);
  std::string ansCommThread = HeuristicEvaluation(TheMap, AllArr.CommThread);
  bool launch_comm_thread = (ansCommThread == "yes");
  std::thread comm_thread;
  std::atomic_bool done = false;
  auto start_comm_thread = [&]() -> void {
    if (launch_comm_thread) {
      os << "Start Thread" << std::endl;
      comm_thread = std::thread(f_comm, std::ref(done));
    }
  };
  auto stop_comm_thread = [&]() -> void {
    if (launch_comm_thread) {
      done = true;
      MicrosecondTime time_join;
      os << "Join thread" << std::endl;
      comm_thread.join();
      os << "|join|=" << time_join << "\n";
    }
  };
  if (ansSplit != "split") {
    auto EXT = df.FF.EXT_face;
    auto Stab = df.Stab;

    start_comm_thread();
    std::string ansProg = AllArr.DualDescriptionProgram.get_eval(TheMap);
    vectface TheOutput = DirectFacetOrbitComputation(EXT, Stab, ansProg, os);
    AllArr.DualDescriptionProgram.pop(os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
    MicrosecondTime time_full;
    os << "|outputsize|=" << TheOutput.size() << "\n";
#endif
    stop_comm_thread();

    for (auto &eOrb : TheOutput) {
      Face eFlipPre = df.FlipFace(eOrb);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      MicrosecondTime time;
#endif
      Face eFlip = bb.operation_face(eFlipPre);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|operation_face1|=" << time << "\n";
#endif
      f_insert(eFlip);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|insert1|=" << time << "\n";
#endif
    }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
    os << "|outputtime|=" << time_full << "\n";
#endif
  } else {
    start_comm_thread();
    try {
      vectface TheOutput =
          DUALDESC_AdjacencyDecomposition<Tbank, T, Tgroup, Tidx_value>(
              TheBank, df.FF.EXT_face, df.Stab, TheMap, AllArr, ePrefix, os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      MicrosecondTime time_full;
      os << "|outputsize|=" << TheOutput.size() << "\n";
#endif
      stop_comm_thread();

      for (auto &eOrb : TheOutput) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
        MicrosecondTime time;
#endif
        Face eFlipPre = df.FlipFace(eOrb);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
        os << "|FlipFace2|=" << time << "\n";
#endif
        Face eFlip = bb.operation_face(eFlipPre);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
        os << "|operation_face2|=" << time << "\n";
#endif
        f_insert(eFlip);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
        os << "|insert2|=" << time << "\n";
#endif
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|outputtime|=" << time_full << "\n";
#endif
    } catch (RuntimeException const &e) {
      os << "RuntimeException, join comm thread\n";
      stop_comm_thread();
      throw; // rethrow
    }
  }
}

template <typename Tbank, typename TbasicBank, typename T, typename Tgroup,
          typename Tidx_value>
vectface MPI_Kernel_DUALDESC_AdjacencyDecomposition(
    boost::mpi::communicator &comm, Tbank &TheBank, TbasicBank &bb,
    PolyHeuristicSerial<typename Tgroup::Tint> &AllArr,
    std::string const &ePrefix,
    std::map<std::string, typename Tgroup::Tint> const &TheMap,
    std::ostream &os) {
  using Tint = typename TbasicBank::Tint;
  SingletonTime start;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  std::string lPrefix = ePrefix + "database";
  DatabaseOrbits<TbasicBank> RPL(bb, lPrefix, AllArr.Saving,
                                 AllArr.AdvancedTerminationCriterion, os);
  auto set_up = [&]() -> void {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
    MicrosecondTime time;
#endif
    std::string ansChoiceCanonic =
        HeuristicEvaluation(TheMap, AllArr.ChoiceCanonicalization);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
    os << "|HeuristicEvaluation|=" << time
       << " ansChoiceCanonic=" << ansChoiceCanonic << "\n";
#endif
    int action = RPL.determine_action_database(ansChoiceCanonic);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
    os << "|determine_action_database|=" << time << " action=" << action
       << "\n";
#endif
    auto f_recompute = [&](int const &method) -> void {
      size_t n_orbit = RPL.preload_nb_orbit();
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|n_orbit|=" << time << "\n";
#endif
      vectface vfo = RPL.ReadDatabase(n_orbit);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|ReadDatabase|=" << time << "\n";
#endif
      RPL.set_method(method);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|set_method|=" << time << "\n";
#endif
      vectface_update_method(vfo, bb, os);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|method update|=" << time << "\n";
#endif
      vectface vfb = mpi_shuffle(comm, std::move(vfo), bb.nbRow);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|mpi_shuffle|=" << time << "\n";
#endif
      RPL.DirectAppendDatabase(std::move(vfb));
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|DirectAppendDatabase|=" << time << "\n";
#endif
    };
    if (action == DATABASE_ACTION__SIMPLE_LOAD) {
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "Before RPL.LoadDatabase()\n";
#endif
      return RPL.LoadDatabase();
    }
    if (action == DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE) {
      int method = bb.convert_string_method(ansChoiceCanonic);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "Before f_recompute, method=" << method
         << " ansChoiceCanonic=" << ansChoiceCanonic << "\n";
#endif
      return f_recompute(method);
    }
    if (action == DATABASE_ACTION__GUESS) {
      vectface vf = RPL.get_runtime_testcase();
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|get_runtime_testcase|=" << time << "\n";
#endif
      int method = RPL.bb.evaluate_method_mpi(comm, vf);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC_MPI
      os << "|evaluate_method_serial|=" << time << "\n";
#endif
      if (method == bb.the_method) {
        return RPL.LoadDatabase();
      } else {
        return f_recompute(method);
      }
    }
    std::cerr << "Failed to find a matching entry for action=" << action
              << "\n";
    throw TerminalException{1};
  };
  set_up();
  bool HasReachedRuntimeException = false;
  int n_vert = bb.nbRow;
  int n_vert_div8 = (n_vert + 7) / 8;
  bool use_f_insert_pair = bb.use_f_insert_pair();
  std::vector<uint8_t> V_hash(n_vert_div8, 0);
  auto get_hash = [&](Face const &x) -> size_t {
    return mpi_get_hash_kernel(x, n_vert, n_vert_div8, V_hash);
  };
  if (AllArr.max_runtime < 0) {
    std::cerr << "The MPI version requires a strictly positive runtime\n";
    throw TerminalException{1};
  }
  //
  // The types of exchanges
  //
  // New facets to be added, the most common request
  const int tag_new_facets = 36;
  // undone information for Balinski termination
  const int tag_termination = 38;
  //
  // Reading the input
  //
  // Grow max buffer like sqrt(n_proc) instead of n_proc to
  // 1. Get O(sqrt(n_proc)) memory usage per thread
  // 2. Send first orbits earlier, to prevent empty threads
  //    to do nothing for many hours at the start.
  size_t MaxBuffered = 10000 * size_t(std::sqrt(n_proc));
  int MaxFly;
  if (AllArr.SimpleExchangeScheme) {
    MaxFly = n_proc;
  } else {
    MaxFly = 4 * n_proc;
  }
  //
  // The parallel MPI classes
  //
  empty_message_management emm_termin(comm, 0, tag_termination);
  os << "emm_termin has been created\n";
  buffered_T_exchanges<Face, vectface> bte_facet(comm, MaxFly, tag_new_facets);
  os << "bte_facet has been created\n";

  std::vector<int> StatusNeighbors(n_proc, 0);
  auto FuncInsertGeneral = [&](Face const &face) -> void {
    if (use_f_insert_pair)
      RPL.FuncInsertPair(face);
    else
      RPL.FuncInsert(face);
  };
  auto f_insert = [&](Face const &face) -> void {
    int res = static_cast<int>(get_hash(face) % size_t(n_proc));
    if (res == i_rank) {
      FuncInsertGeneral(face);
    } else {
      bte_facet.insert_entry(res, face);
    }
  };
  //
  // Initial invocation of the synchronization code
  //
  os << "Compute initital\n";
  size_t n_orb_max = 0, n_orb_loc = RPL.FuncNumberOrbit();
  all_reduce(comm, n_orb_loc, n_orb_max, boost::mpi::maximum<size_t>());
  os << "n_orb_loc=" << n_orb_loc << " n_orb_max=" << n_orb_max << "\n";
  if (n_orb_max == 0) {
    std::string ansSamp = HeuristicEvaluation(TheMap, AllArr.InitialFacetSet);
    os << "ansSamp=" << ansSamp << "\n";
    vectface vf_init = DirectComputationInitialFacetSet(bb.EXT, ansSamp, os);
    vectface vf_init_merge = merge_initial_samp(comm, vf_init, ansSamp, os);
    for (auto &face : vf_init_merge) {
      Face face_can = bb.operation_face(face);
      f_insert(face_can);
    }
  }
  os << "DirectFacetOrbitComputation, step 6\n";
  //
  // The infinite loop
  //
  auto process_mpi_status = [&](boost::mpi::status const &stat) -> void {
    int e_tag = stat.tag();
    int e_src = stat.source();
    if (e_tag == tag_new_facets) {
      os << "RECV of tag_new_facets from " << e_src << "\n";
      StatusNeighbors[e_src] = 0;
      vectface l_recv_face = bte_facet.recv_message(e_src);
      os << "|l_recv_face|=" << l_recv_face.size() << "\n";
      for (auto &face : l_recv_face)
        FuncInsertGeneral(face);
    }
    if (e_tag == tag_termination) {
      os << "RECV of tag_termination from " << e_src << "\n";
      StatusNeighbors[e_src] = 1;
      emm_termin.recv_message(e_src);
    }
    os << "Exiting process_mpi_status\n";
  };
  auto get_maxruntimereached = [&]() -> bool {
    if (HasReachedRuntimeException)
      return true;
    return si(start) > AllArr.max_runtime;
  };
  auto send_termination_notice = [&]() -> void {
    // This function should be passed only one time.
    // It says to the receiving nodes: "You wil never ever received anything
    // more from me"
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
  };
  auto get_nb_finished = [&]() -> int {
    int nb_finished = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      os << "get_nb_finished : i_proc=" << i_proc
         << " status=" << StatusNeighbors[i_proc] << "\n";
      nb_finished += StatusNeighbors[i_proc];
    }
    os << "nb_finished=" << nb_finished << " n_proc=" << n_proc << "\n";
    return nb_finished;
  };
  auto get_nb_finished_oth = [&]() -> int {
    int nb_finished_oth = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++) {
      if (i_proc != i_rank) {
        os << "get_nb_finished_oth : i_proc=" << i_proc
           << " status=" << StatusNeighbors[i_proc] << "\n";
        nb_finished_oth += StatusNeighbors[i_proc];
      }
    }
    os << "nb_finished_oth=" << nb_finished_oth << " n_proc=" << n_proc << "\n";
    return nb_finished_oth;
  };
  auto short_wait = [&]() -> void {
    int n_milliseconds = 10;
    std::this_thread::sleep_for(std::chrono::milliseconds(n_milliseconds));
  };

  auto wait = [&]() -> void {
    int n_milliseconds = 1000;
    std::this_thread::sleep_for(std::chrono::milliseconds(n_milliseconds));
  };
  auto f_comm = [&](std::atomic_bool &done) -> void {
    // start with no work
    size_t cnt = 0;
    while (!done && cnt < 10000) {
      cnt++;
      short_wait();
    }
    if (!done)
      os << "Start Probing" << std::endl;
    // Start probing
    while (!done) {
      boost::optional<boost::mpi::status> prob = comm.iprobe();
      if (prob) {
        os << "process_mpi_status" << std::endl;
        process_mpi_status(*prob);
      } else {
        wait();
      }
    }
  };
  auto process_database = [&]() -> void {
    DataFacet<T, Tgroup> df = RPL.FuncGetMinimalUndoneOrbit();
    size_t SelectedOrbit = df.SelectedOrbit;
    std::string NewPrefix = ePrefix + "PROC" + std::to_string(i_rank) + "_ADM" +
                            std::to_string(SelectedOrbit) + "_";
    try {
      DUALDESC_AdjacencyDecomposition_and_insert_commthread<
          Tbank, T, Tgroup, Tidx_value, TbasicBank, decltype(f_insert),
          decltype(f_comm)>(TheBank, bb, df, AllArr, f_insert, f_comm,
                            NewPrefix, os);

      RPL.FuncPutOrbitAsDone(SelectedOrbit);
    } catch (RuntimeException const &e) {
      HasReachedRuntimeException = true;
      os << "The computation of DUALDESC_AdjacencyDecomposition has ended by "
            "runtime exhaustion\n";
    }
  };
  bool HasSendTermination = false;
  while (true) {
    os << "DirectFacetOrbitComputation, inf loop, start\n";
    bool MaxRuntimeReached = get_maxruntimereached();
    if (MaxRuntimeReached && !HasSendTermination) {
      if (bte_facet.get_unsent_size() == 0) {
        HasSendTermination = true;
        send_termination_notice();
      }
    }
    bool SomethingToDo = !MaxRuntimeReached && !RPL.IsFinished();
    os << "DirectFacetOrbitComputation, MaxRuntimeReached=" << MaxRuntimeReached
       << " SomethingToDo=" << SomethingToDo << "\n";
    os << "get_unsent_size()=" << bte_facet.get_unsent_size()
       << " MaxBuffered=" << MaxBuffered << "\n";
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      os << "prob is not empty\n";
      process_mpi_status(*prob);
    } else {
      os << "prob is empty\n";
      if (SomethingToDo) {
        os << "Case something to do\n";
        // we have to clear our buffers sometimes while running
        // otherwise we will only do it in the very end
        // which could lead to extremely large buffers
        // and possibly treating orbits with unnecessary large incidence
        if (bte_facet.get_unsent_size() >= MaxBuffered) {
          os << "Calling clear_one_entry after reaching MaxBuffered\n";
          if (!bte_facet.clear_one_entry(os))
            wait();
        }
        process_database();
      } else {
        bool test = bte_facet.is_buffer_empty();
        os << "Case nothing to do test=" << test << "\n";
        if (test) {
          int nb_finished_oth = get_nb_finished_oth();
          os << "Nothing to do, entering the busy loop status="
             << get_maxruntimereached()
             << " get_nb_finished_oth()=" << nb_finished_oth << "\n";
          // do while, such that we also sleep after maxruntimereached
          // to prevent log spam when other threads are still working
          do {
            boost::optional<boost::mpi::status> prob = comm.iprobe();
            if (prob) {
              process_mpi_status(*prob);
              break;
            }
            wait();
          } while (!get_maxruntimereached());
        } else {
          bool test = bte_facet.clear_one_entry(os);
          os << "Calling clear_one_entry test=" << test << "\n";
          if (!test)
            wait();
        }
      }
    }
    //
    // Termination criterion
    //
    if (get_nb_finished() == n_proc)
      break;
    os << "End of the while loop. From start=" << si(start) << "\n";
  }
  os << "We just exited the infinite loop\n";
  bool test_termination = EvaluationConnectednessCriterion_MPI(comm, bb, os);
  os << "We have test_termination=" << test_termination << "\n";
  if (test_termination) {
    os << "Correct termination, returning the database\n";
    FaceOrbitsizeTableContainer<Tint> fotc = RPL.GetListFaceOrbitsize();
    return fotc.GetListFaces();
  } else {
    os << "RuntimeException, terminating the computation\n";
    throw RuntimeException{1};
  }
}

void update_path_using_nproc_iproc(std::string &str_ref, int const &n_proc,
                                   int const &i_rank) {
  std::string postfix =
      "_nproc" + std::to_string(n_proc) + "_rank" + std::to_string(i_rank);
  size_t len = str_ref.size();
  std::string part1 = str_ref.substr(0, len - 1);
  std::string part2 = str_ref.substr(len - 1, 1);
  if (part2 != "/") {
    std::cerr << "str_ref=" << str_ref << "\n";
    std::cerr << "Last character should be a /\n";
    throw TerminalException{1};
  }
  str_ref = part1 + postfix + part2;
}

template <typename T>
void Reset_Directories(boost::mpi::communicator &comm,
                       PolyHeuristicSerial<T> &AllArr) {
  int n_proc = comm.size();
  int i_rank = comm.rank();
  auto update_string = [&](std::string &str_ref) -> void {
    if (i_rank < n_proc - 1 ||
        AllArr.bank_parallelization_method != "bank_mpi") {
      update_path_using_nproc_iproc(str_ref, n_proc, i_rank);
      CreateDirectory(str_ref);
    }
  };
  if (AllArr.bank_parallelization_method == "serial")
    update_string(AllArr.BANK_Prefix);
  else
    CreateDirectory(AllArr.BANK_Prefix);
  update_string(AllArr.DD_Prefix);
}

template <typename T, typename Tgroup, typename Tidx_value>
void MPI_MainFunctionDualDesc(boost::mpi::communicator &comm,
                              FullNamelist const &eFull) {
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tkey = MyMatrix<T>;
  using Tval = TripleStore<Tgroup>;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  int pos_generator = 0;
  //
  std::string FileLog =
      "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::cerr << "We have moved. See the log in FileLog=" << FileLog << "\n";
  std::ofstream os(FileLog);
  if (ApplyStdUnitbuf(eFull))
    os << std::unitbuf;
  os << "Initial writing of the log\n";
  os.flush();
  //
  MyMatrix<T> EXT = Get_EXT_DualDesc<T, Tidx>(eFull, os);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  Tgroup GRP = Get_GRP_DualDesc<Tgroup>(eFull, os);
  PolyHeuristicSerial<Tint> AllArr =
      Read_AllStandardHeuristicSerial<T, Tint>(eFull, EXTred, os);
  srand(time(NULL) + 12345 * i_rank);
  Reset_Directories(comm, AllArr);
  size_t n_rows = EXTred.rows();
  if (AllArr.bank_parallelization_method == "bank_mpi" && n_proc < 2) {
    std::cerr << "For the bank_mpi we need at least 2 nodes. n_proc=" << n_proc
              << "\n";
    throw TerminalException{1};
  }
  int proc_bank = n_proc - 1;
  int i_proc_ret = 0;
  auto msg_term_bank = [&]() -> void {
    if (AllArr.bank_parallelization_method == "bank_mpi" &&
        i_rank == i_proc_ret) {
      os << "sending bank_mpi termination signal\n";
      comm.send(proc_bank, tag_mpi_bank_end, val_mpi_bank_end);
    }
  };
  if (AllArr.bank_parallelization_method == "bank_mpi" && i_rank == n_proc - 1)
    pos_generator = 1;
  boost::mpi::communicator comm_work = comm.split(pos_generator);
  //
  using TbasicBank = DatabaseCanonic<T, Tint, Tgroup>;
  TbasicBank bb(EXTred, GRP, os);
  std::map<std::string, Tint> TheMap =
      ComputeInitialMap<Tint>(EXTred, GRP, AllArr);
  //
  auto get_vectface = [&]() -> vectface {
    if (AllArr.bank_parallelization_method == "serial") {
      using Tbank = DataBank<Tkey, Tval>;
      Tbank TheBank(AllArr.BANK_IsSaving, AllArr.BANK_Prefix, os);
      return MPI_Kernel_DUALDESC_AdjacencyDecomposition<Tbank, TbasicBank, T,
                                                        Tgroup, Tidx_value>(
          comm, TheBank, bb, AllArr, AllArr.DD_Prefix, TheMap, os);
    }
    if (AllArr.bank_parallelization_method == "bank_asio") {
      using Tbank = DataBankAsioClient<Tkey, Tval>;
      Tbank TheBank(AllArr.port);
      return MPI_Kernel_DUALDESC_AdjacencyDecomposition<Tbank, TbasicBank, T,
                                                        Tgroup, Tidx_value>(
          comm, TheBank, bb, AllArr, AllArr.DD_Prefix, TheMap, os);
    }
    if (AllArr.bank_parallelization_method == "bank_mpi") {
      using Tbank = DataBankMpiClient<Tkey, Tval>;
      Tbank TheBank(comm);
      if (i_rank < proc_bank) {
        return MPI_Kernel_DUALDESC_AdjacencyDecomposition<Tbank, TbasicBank, T,
                                                          Tgroup, Tidx_value>(
            comm_work, TheBank, bb, AllArr, AllArr.DD_Prefix, TheMap, os);
      } else {
        DataBankMpiServer<Tkey, Tval>(comm, AllArr.BANK_IsSaving,
                                      AllArr.BANK_Prefix, os);
        return vectface(n_rows);
      }
    }
    std::cerr << "No match for bank_parallelization_method\n";
    std::cerr << "AllArr.bank_parallelization_method="
              << AllArr.bank_parallelization_method << "\n";
    std::cerr << "Allowed methods are serial, bank_asio, bank_mpi\n";
    throw TerminalException{1};
  };

  try {
    vectface vf = get_vectface();
    if(AllArr.bank_parallelization_method == "bank_mpi" and i_rank == proc_bank) {
      os << "Closed DataBankMpiServer" << std::endl;
      return;
    }
    msg_term_bank();
    // output
    os << "We have vf |vf|=" << vf.size() << " / " << vf.get_n() << "\n";
    vectface vf_tot = my_mpi_gather(comm_work, vf, i_proc_ret);
    os << "We have vf_tot |vf_tot|=" << vf_tot.size() << " / " << vf_tot.get_n()
       << " i_proc_ret=" << i_proc_ret << "\n";
    if (i_rank == i_proc_ret)
      OutputFacets(EXT, GRP, vf_tot, AllArr.OUTfile, AllArr.OutFormat);
    os << "We have done our output\n";
  } catch (RuntimeException const &e) {
    msg_term_bank();
    throw RuntimeException{1};
  }
}

// clang-format off
#endif  // SRC_DUALDESC_POLY_RECURSIVEDUALDESC_MPI_H_
// clang-format on
