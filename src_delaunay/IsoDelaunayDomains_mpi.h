// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "POLY_MPI_AdjacencyScheme.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_ISO_DELAUNAY_DOMAINS_MPI
#endif

#ifdef DISABLE_DEBUG_ISO_DELAUNAY_DOMAINS_MPI
#undef DEBUG_ISO_DELAUNAY_DOMAINS_MPI
#endif

template <typename T, typename Tint, typename Tgroup>
void ComputeLatticeIsoDelaunayDomains_MPI(boost::mpi::communicator &comm,
                                          FullNamelist const &eFull) {
  using TintGroup = typename Tgroup::Tint;
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  int dimEXT = LinSpa.n + 1;
  //
  bool STORAGE_Saving = BlockSYSTEM.get_bool("Saving");
  std::string STORAGE_Prefix = BlockSYSTEM.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
#ifdef DEBUG_ISO_DELAUNAY_DOMAINS_MPI
  os << "ISODELMPI: max_runtime_second=" << max_runtime_second << "\n";
#endif
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
#ifdef DEBUG_ISO_DELAUNAY_DOMAINS_MPI
  os << "ISODELMPI: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
#endif
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  DataIsoDelaunayDomains<T, Tint, Tgroup> data =
      get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, os);
  //
  using Tdata = DataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
  if (pair.first) {
    // Need to reread after the move
    DataIsoDelaunayDomains<T, Tint, Tgroup> data =
        get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, os);
    auto f_print = [&](std::ostream &os_out) -> void {
      bool result =
          WriteFamilyObjects_MPI<DataIsoDelaunayDomains<T, Tint, Tgroup>, Tobj,
                                 TadjO>(comm, data, OutFormat, os_out,
                                        pair.second, os);
      if (result) {
        std::cerr << "ISODELMPI: Failed to find a matching entry for OutFormat="
                  << OutFormat << "\n";
        throw TerminalException{1};
      }
    };
    print_stderr_stdout_file(OutFile, f_print);
  }
}

// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_
// clang-format on
