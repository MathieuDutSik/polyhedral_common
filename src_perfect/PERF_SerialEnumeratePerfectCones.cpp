// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "perfect_tspace.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process_A(FullNamelist const &eFull) {
  using TintGroup = typename Tgroup::Tint;
  SingleBlock const& BlockDATA = eFull.get_block("DATA");
  SingleBlock const& BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const& BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
  int n = LinSpa.n;
  int dimEXT = n * (n + 1) / 2;
  //
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  std::string STORAGE_Prefix = BlockSYSTEM.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT, std::cerr);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, std::cerr);
  //
  bool keep_generators = false;
  bool reduce_gram_matrix = false;
  DataPerfectTspace<T, Tint, Tgroup> data{LinSpa,
                                          OnlineHierarchicalMatrixReduction<Tint>(n, std::cerr),
                                          keep_generators,
                                          reduce_gram_matrix,
                                          std::move(rddo)};
  using Tdata = DataPerfectTspaceFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::vector<Tout> l_tot = EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
      data_func, f_incorrect, max_runtime_second);
  auto f_print=[&](std::ostream& os_out) -> void {
    bool result = WriteFamilyObjects(data, OutFormat, os_out, l_tot, std::cerr);
    if (result) {
      std::cerr << "PERFSERIAL: Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    }
  };
  print_stderr_stdout_file(OutFile, f_print);
}

template <typename T, typename Tint>
void process_B(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  return process_A<T, Tint, Tgroup>(eFull);
}

template <typename T>
void process_C(FullNamelist const &eFull) {
  std::string arithmetic_Tint =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_Tint");
  if (arithmetic_Tint == "gmp_integer") {
    using Tint = mpz_class;
    return process_B<T, Tint>(eFull);
  }
  std::cerr << "PERF_SerialEnumeratePerfectCones: Failed to find a matching type for "
               "arithmetic_Tint="
            << arithmetic_Tint << "\n";
  std::cerr << "Available types: gmp_integer\n";
  throw TerminalException{1};
}

void process_D(FullNamelist const &eFull) {
  std::string arithmetic_T =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_T");
  if (arithmetic_T == "gmp_rational") {
    using T = mpq_class;
    return process_C<T>(eFull);
  }
  std::cerr << "PERF_SerialEnumeratePerfectCones: Failed to find a matching type for "
               "arithmetic_T="
            << arithmetic_T << "\n";
  std::cerr << "Available types: gmp_rational\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT_TSPACE();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PERF_SerialEnumeratePerfectCones [file.nml]\n";
      std::cerr << "\n";
      std::cerr << "        --- or with a custom template ---\n";
      std::cerr << "\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_D(eFull);
    std::cerr << "Normal termination of PERF_SerialEnumeratePerfectCones\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERF_SerialEnumeratePerfectCones\n";
    exit(e.eVal);
  } catch (std::runtime_error const &e) {
    std::cerr << "Runtime error in PERF_SerialEnumeratePerfectCones\n";
    std::cerr << "error: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
  runtime(time);
}
