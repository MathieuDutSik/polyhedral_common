// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_
#define SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_

// clang-format off
#include "IsoDelaunayDomains.h"
#include "POLY_MPI_AdjacencyScheme.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void ComputeLatticeIsoDelaunayDomains_MPI(boost::mpi::communicator &comm,
                                          FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  //
  bool STORAGE_Saving = BlockDATA.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockDATA.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  os << "ISODELMPI: max_runtime_second=" << max_runtime_second << "\n";
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  os << "ISODELMPI: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";

  DataIsoDelaunayDomains<T, Tint, Tgroup> data = get_data_isodelaunay_domains<T,Tint,Tgroup>(eFull, os);

  using Tdata = DataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
  if (pair.first) {
    std::ofstream os_out(OutFile);
    bool result = WriteFamilyObjects_MPI<Tobj, TadjO>(comm, OutFormat, os_out, pair.second, os);
    if (result) {
      std::cerr	<< "ISODELMPI: Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    }
  }
}

// clang-format off
#endif  // SRC_LATT_ISODELAUNAYDOMAINS_MPI_H_
// clang-format on
