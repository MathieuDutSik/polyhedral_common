// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_
#define SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_

// clang-format off
#include "lorentzian_perfect.h"
#include "POLY_MPI_AdjacencyScheme.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void ComputePerfectLorentzian(boost::mpi::communicator &comm,
                              FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockSTORAGE = eFull.ListBlock.at("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" << max_runtime_second << "\n";
  std::string LorMatFile = BlockDATA.ListStringValues.at("LorMatFile");
  MyMatrix<T> LorMat = ReadMatrixFile<T>(LorMatFile);
  //
  std::string TheOption_str = BlockDATA.ListStringValues.at("TheOption");
  auto get_option = [&]() -> int {
    if (TheOption_str == "isotropic") {
      return LORENTZIAN_PERFECT_OPTION_ISOTROP;
    }
    if (TheOption_str == "total") {
      return LORENTZIAN_PERFECT_OPTION_TOTAL;
    }
    std::cerr << "Failed to find a matching entry for TheOption\n";
    throw TerminalException{1};
  };
  int TheOption = get_option();
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  //
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  std::string FileDualDesc =
      BlockDATA.ListStringValues.at("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  //
  DataPerfectLorentzian<T, Tint, Tgroup> data{n, LorMat, TheOption,
                                              std::move(rddo)};
  using Tdata = DataPerfectLorentzianFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  //
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "LORPERF: We now have IsFinished=" << pair.first << "\n";
  os << "LORPERF: We now have |ListPerfect|=" << pair.second.size() << "\n";
#endif
  //
  if (pair.first) {
    WriteFamilyObjects<Tobj, TadjO>(comm, OutFormat, OutFile, pair.second, os);
  }
}

// clang-format off
#endif  // SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_
// clang-format on
