// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_LORENTZIAN_PERFECT_SERIAL_H_
#define SRC_LORENTZIAN_LORENTZIAN_PERFECT_SERIAL_H_

// clang-format off
#include "lorentzian_perfect.h"
#include "Isotropic.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_LORENTZIAN_PERFECT_SERIAL
#endif

/*
  The serial counterpart of ComputePerfectLorentzian_mpi
  (lorentzian_perfect_mpi.h): it enumerates the perfect Lorentzian
  domains through the serial adjacency scheme rather than the MPI one, so
  it needs no MPI toolchain. The namelist and the output formats are the
  same, so a namelist written for the MPI program runs here unchanged.
 */
template <typename T, typename Tint, typename Tgroup>
void ComputePerfectLorentzian_serial(FullNamelist const &eFull,
                                     std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
  os << "LORPERFSER: ComputePerfectLorentzian_serial, beginning\n";
#endif
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  //
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
  os << "LORPERFSER: max_runtime_second=" << max_runtime_second << "\n";
#endif
  std::string LorMatFile = BlockDATA.get_string("LorMatFile");
  MyMatrix<T> LorMat = ReadMatrixFile<T>(LorMatFile);
  check_correctness_lorentzian_perfect(LorMat, os);
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
  os << "LORPERFSER: Pass the lorentzian correctness check\n";
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
    std::cerr << "Failed to find a matching entry for TheOption_str="
              << TheOption_str << " allowed: total, isotropic\n";
    throw TerminalException{1};
  };
  int TheOption = get_option();
  if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
    bool test = is_isotropic(LorMat, os);
    if (!test) {
      std::cerr << "LORPERFSER: We have a request with isotropic\n";
      std::cerr << "LORPERFSER: However, the matrix is not isotropic\n";
      throw TerminalException{1};
    }
  }
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
  os << "LORPERFSER: OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
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
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  //
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::optional<std::vector<Tout>> opt_l_obj =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
  os << "LORPERFSER: We now have IsFinished=" << opt_l_obj.has_value() << "\n";
#endif
  //
  if (opt_l_obj) {
    std::vector<Tout> const &l_obj = *opt_l_obj;
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
    os << "LORPERFSER: |ListPerfect|=" << l_obj.size() << ", doing output\n";
#endif
    auto f_print = [&](std::ostream &os_out) -> void {
      bool result = WriteFamilyObjects<DataPerfectLorentzian<T, Tint, Tgroup>,
                                       Tobj, TadjO>(data_func.data, OutFormat,
                                                    os_out, l_obj, os);
      if (result) {
        std::cerr << "LORPERFSER: Failed to find a matching entry for "
                  << "OutFormat=" << OutFormat << "\n";
        throw TerminalException{1};
      }
    };
    FILE_PrintStderrStdoutFile(OutFile, f_print);
  } else {
#ifdef DEBUG_LORENTZIAN_PERFECT_SERIAL
    os << "LORPERFSER: No output being done\n";
#endif
  }
}

// clang-format off
#endif  // SRC_LORENTZIAN_LORENTZIAN_PERFECT_SERIAL_H_
// clang-format on
