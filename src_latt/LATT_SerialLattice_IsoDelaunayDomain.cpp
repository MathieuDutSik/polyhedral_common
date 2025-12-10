// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint> void process_A(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa =
      ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
  int dimEXT = LinSpa.n + 1;
  //
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  std::string STORAGE_Prefix = BlockSYSTEM.get_string("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         std::cerr);

  DataIsoDelaunayDomains<T, Tint, Tgroup> data =
      get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, std::cerr);
  //
  using Tdata = DataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::vector<Tout> l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  std::ofstream os_out(OutFile);
  bool result = WriteFamilyObjects(data, OutFormat, os_out, l_tot, std::cerr);
  if (result) {
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
}

template <typename T> void process_B(FullNamelist const &eFull) {
  std::string arithmetic_Tint =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_Tint");
  if (arithmetic_Tint == "gmp_integer") {
    using Tint = mpz_class;
    return process_A<T, Tint>(eFull);
  }
  std::cerr << "LATT_SerialLattice_IsoDelaunayDomain B: Failed to find a "
               "matching type for arithmetic_Tint="
            << arithmetic_Tint << "\n";
  std::cerr << "Available types: gmp_integer\n";
  throw TerminalException{1};
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic_T =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic_T");
  if (arithmetic_T == "gmp_rational") {
    using T = mpq_class;
    return process_B<T>(eFull);
  }
  std::cerr << "LATT_SerialLattice_IsoDelaunayDomain A: Failed to find a "
               "matching type for arithmetic_T="
            << arithmetic_T << "\n";
  std::cerr << "Available types: gmp_rational\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull =
        NAMELIST_GetStandard_COMPUTE_LATTICE_IsoDelaunayDomains();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialLattice_IsoDelaunayDomain [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    //
    std::cerr << "Normal termination of LATT_SerialLattice_IsoDelaunayDomain\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialLattice_IsoDelaunayDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}
