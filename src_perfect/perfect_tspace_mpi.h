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
bool WriteFamilyObjects_MPI(boost::mpi::communicator const &comm, 
                          Tdata const &data,
                          std::string const &OutFormat,
                          std::ostream &os_out,
                          std::vector<DatabaseEntry_MPI<Tobj, TadjO>> const &l_obj,
                          std::ostream &os) {
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

FullNamelist NAMELIST_GetStandard_ENUMERATE_PERFECT_TSPACE() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  SingleBlock BlockDATA;
  BlockDATA.setListIntValues(ListIntValues1);
  BlockDATA.setListBoolValues(ListBoolValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // STORAGE
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, std::string> ListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  SingleBlock BlockSTORAGE;
  BlockSTORAGE.setListBoolValues(ListBoolValues2);
  BlockSTORAGE.setListStringValues(ListStringValues2);
  ListBlock["STORAGE"] = BlockSTORAGE;
  // TSPACE
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, std::string> ListStringValues3;
  ListStringValues3["TypeTspace"] = "File";
  ListStringValues3["FileLinSpa"] = "unset.linspa";
  ListStringValues3["SuperMatMethod"] = "NotNeeded";
  ListStringValues3["ListComm"] = "Trivial";
  ListStringValues3["PtGroupMethod"] = "Trivial";
  ListStringValues3["FileListSubspaces"] = "unset";
  ListIntValues3["RealImagDim"] = 0;
  ListIntValues3["RealImagSum"] = 0;
  ListIntValues3["RealImagProd"] = 0;
  SingleBlock BlockTSPACE;
  BlockTSPACE.setListIntValues(ListIntValues3);
  BlockTSPACE.setListStringValues(ListStringValues3);
  ListBlock["TSPACE"] = BlockTSPACE;
  // Merging all data
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_PERFECT_PERFECT_TSPACE_MPI_H_
// clang-format on