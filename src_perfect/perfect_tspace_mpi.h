// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECT_TSPACE_MPI_H_
#define SRC_PERFECT_PERFECT_TSPACE_MPI_H_

// clang-format off
#include "perfect_tspace.h"
#include "POLY_MPI_AdjacencyScheme.h"
#include "Tspace_Namelist.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_TSPACE_MPI
#endif

template <typename T, typename Tint, typename Tgroup>
void ComputePerfectTspace_mpi(boost::mpi::communicator &comm,
                             FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
#ifdef DEBUG_PERFECT_TSPACE_MPI
  os << "PERFTSPACEMPI: ComputePerfectTspace_mpi, beginning\n";
#endif
  SingleBlock const& BlockDATA = eFull.get_block("DATA");
  SingleBlock const& BlockSTORAGE = eFull.get_block("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.get_bool("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.get_int("max_runtime_second");
#ifdef DEBUG_PERFECT_TSPACE_MPI
  os << "PERFTSPACEMPI: max_runtime_second=" << max_runtime_second << "\n";
#endif
  SingleBlock const& BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
#ifdef DEBUG_PERFECT_TSPACE_MPI
  os << "PERFTSPACEMPI: LinSpa dimension=" << LinSpa.n << "\n";
#endif
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  std::string OutFile = BlockDATA.get_string("OutFile");
#ifdef DEBUG_PERFECT_TSPACE_MPI
  os << "PERFTSPACEMPI: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
#endif
  //
  int n = LinSpa.n;
  int dimEXT = n * (n + 1) / 2;
  using TintGroup = typename Tgroup::Tint;
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  //
  DataPerfectTspace<T, Tint, Tgroup> data{n, LinSpa, std::move(rddo)};
  using Tdata = DataPerfectTspaceFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  //
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
#ifdef DEBUG_PERFECT_TSPACE_MPI
  os << "PERFTSPACEMPI: We now have max_runtime_second=" << max_runtime_second << "\n";
  os << "PERFTSPACEMPI: We now have IsFinished=" << pair.first << "\n";
  os << "PERFTSPACEMPI: We now have |ListPerfect|=" << pair.second.size() << "\n";
#endif
  //
  if (pair.first) {
#ifdef DEBUG_PERFECT_TSPACE_MPI
    os << "PERFTSPACEMPI: Doing some output\n";
#endif
    auto f_print=[&](std::ostream& os_out) -> void {
      bool result = WriteFamilyObjects_MPI<DataPerfectTspace<T, Tint, Tgroup>, Tobj, TadjO>(
          comm, data, OutFormat, os_out, pair.second, os);
      if (result) {
        std::cerr << "PERFTSPACEMPI: Failed to find a matching entry for OutFormat=" 
                  << OutFormat << "\n";
        throw TerminalException{1};
      }
    };
    print_stderr_stdout_file(OutFile, f_print);
  } else {
#ifdef DEBUG_PERFECT_TSPACE_MPI
    os << "PERFTSPACEMPI: No output being done\n";
#endif
  }
}

template <typename Tdata, typename Tobj, typename TadjO>
bool WriteFamilyObjects_MPI([[maybe_unused]] boost::mpi::communicator const &comm, 
                          [[maybe_unused]] Tdata const &data,
                          std::string const &OutFormat,
                          std::ostream &os_out,
                          std::vector<DatabaseEntry_MPI<Tobj, TadjO>> const &l_obj,
                          [[maybe_unused]] std::ostream &os) {
  if (OutFormat == "GAP") {
    os_out << "[";
    bool IsFirst = true;
    for (auto &ent : l_obj) {
      if (!IsFirst)
        os_out << ",";
      IsFirst = false;
      WriteEntryGAP(os_out, ent);
    }
    os_out << "]";
    return false;
  }
  if (OutFormat == "PYTHON") {
    os_out << "[";
    bool IsFirst = true;
    for (auto &ent : l_obj) {
      if (!IsFirst)
        os_out << ",";
      IsFirst = false;
      WriteEntryPYTHON(os_out, ent);
    }
    os_out << "]";
    return false;
  }
  if (OutFormat == "nothing") {
    return false;
  }
  return true;
}


// clang-format off
#endif  // SRC_PERFECT_PERFECT_TSPACE_MPI_H_
// clang-format on