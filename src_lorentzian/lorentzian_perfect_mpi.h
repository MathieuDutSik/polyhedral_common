// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_
#define SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_

// clang-format off
#include "lorentzian_perfect.h"
#include "Isotropic.h"
#include "POLY_MPI_AdjacencyScheme.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_LORENTZIAN_PERFECT_MPI
#endif

template <typename T, typename Tint, typename Tgroup>
void ComputePerfectLorentzian(boost::mpi::communicator &comm,
                              FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
  os << "LORPERFMPI: ComputePerfectLorentzian, beginning\n";
#endif
  SingleBlock const& BlockDATA = eFull.get_block("DATA");
  SingleBlock const& BlockSTORAGE = eFull.get_block("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.get_bool("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.get_int("max_runtime_second");
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
  os << "LORPERFMPI: max_runtime_second=" << max_runtime_second << "\n";
#endif
  std::string LorMatFile = BlockDATA.get_string("LorMatFile");
  MyMatrix<T> LorMat = ReadMatrixFile<T>(LorMatFile);
  check_correctness_lorentzian_perfect(LorMat);
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
  os << "LORPERFMPI: Pass the lorentzian correctness check\n";
#endif
  //
  std::string TheOption_str = BlockDATA.get_string("Option");
  auto get_option = [&]() -> int {
    if (TheOption_str == "isotropic") {
      return LORENTZIAN_PERFECT_OPTION_ISOTROP;
    }
    if (TheOption_str == "total") {
      return LORENTZIAN_PERFECT_OPTION_TOTAL;
    }
    std::cerr << "Failed to find a matching entry for TheOption_str=" << TheOption_str << " allowed: total, isotropic\n";
    throw TerminalException{1};
  };
  int TheOption = get_option();
  if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
    bool test = is_isotropic(LorMat, os);
    if (!test) {
      std::cerr << "LORPERFMPI: We have a request with isotropic\n";
      std::cerr << "LORPERFMPI: However, the matrix is not isotropic\n";
      throw TerminalException{1};
    }
  }
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  std::string OutFile = BlockDATA.get_string("OutFile");
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
  os << "LORPERFMPI: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
#endif
  //
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
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
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
  os << "LORPERFMPI: We now have max_runtime_second=" << max_runtime_second << "\n";
  os << "LORPERFMPI: We now have IsFinished=" << pair.first << "\n";
  os << "LORPERFMPI: We now have |ListPerfect|=" << pair.second.size() << "\n";
#endif
  //
  if (pair.first) {
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
    os << "LORPERFMPI: Doing some output\n";
#endif
    auto f_print=[&](std::ostream& os_out) -> void {
      bool result = WriteFamilyObjects_MPI<DataPerfectLorentzian<T, Tint, Tgroup>, Tobj, TadjO>(comm, data, OutFormat, os_out, pair.second, os);
      if (result) {
        std::cerr	<< "LORPERFMPI: Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
        throw TerminalException{1};
      }
    };
    print_stderr_stdout_file(OutFile, f_print);
  } else {
#ifdef DEBUG_LORENTZIAN_PERFECT_MPI
    os << "LORPERFMPI: No output being done\n";
#endif
  }
}

// clang-format off
#endif  // SRC_LORENTZIAN_LORENTZIAN_PERFECT_MPI_H_
// clang-format on
