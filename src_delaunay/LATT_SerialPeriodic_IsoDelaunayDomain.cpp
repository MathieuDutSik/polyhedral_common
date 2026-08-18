// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
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
  //
  // The periodic point set, with the validation of the input. Those are
  // checks of user input, not internal invariants, so they are performed
  // unconditionally.
  //
  std::string FileCosets = BlockDATA.get_string("FileCosets");
  MyMatrix<T> Cosets = ReadMatrixFile<T>(FileCosets);
  if (Cosets.cols() != LinSpa.n) {
    std::cerr << "LATT_SerialPeriodic_IsoDelaunayDomain: the cosets have "
              << Cosets.cols() << " columns while the T-space is of dimension "
              << LinSpa.n << "\n";
    throw TerminalException{1};
  }
  std::optional<PeriodicPointSet<Tint>> opt_pps =
      PeriodicPointSetFromRational_Opt<Tint, T>(Cosets);
  if (!opt_pps) {
    std::cerr << "LATT_SerialPeriodic_IsoDelaunayDomain: the cosets form a "
                 "group modulo Z^n, so the point set is a lattice: it has to "
                 "be treated by LATT_SerialLattice_IsoDelaunayDomain over "
                 "that lattice instead\n";
    throw TerminalException{1};
  }
  PeriodicPointSet<Tint> pps = *opt_pps;
  for (auto &eGen : LinSpa.PtStabGens) {
    if (!PeriodicAffineExtension<Tint, T>(pps, eGen)) {
      std::cerr << "LATT_SerialPeriodic_IsoDelaunayDomain: a generator of the "
                   "pointwise stabilizer of the T-space admits no affine "
                   "extension preserving the point set. The point stabilizer "
                   "has to stabilize the cosets for the enumeration to be "
                   "correct. The offending generator is\n";
      WriteMatrix(std::cerr, eGen);
      throw TerminalException{1};
    }
  }
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         std::cerr);
  std::optional<MyMatrix<T>> CommonGramMat;
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, std::cerr);
  LinSpaceMatrix<Tint> LinSpaRing = LINSPA_GetRingVersion(LinSpa);
  std::optional<MyMatrix<Tint>> CommonGramMatRing =
      GetCommonGramMatRing<Tint, T>(CommonGramMat);
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{
      std::move(LinSpa), std::move(LinSpaRing), std::move(rddo), CommonGramMat,
      CommonGramMatRing};

  using Tdata = PeriodicDataIsoDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data), std::move(pps)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  std::optional<std::vector<Tout>> opt_l_tot =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  std::vector<Tout> l_tot =
      unfold_opt(opt_l_tot, "EnumerateAndStore_Serial (periodic iso-Delaunay)");

  std::ofstream os_out(OutFile);
  bool result =
      WriteFamilyObjects(data_func.data, OutFormat, os_out, l_tot, std::cerr);
  if (result) {
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  }
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_A<T, Tint>(eFull);
  }
  std::cerr << "LATT_SerialPeriodic_IsoDelaunayDomain: Failed to find a "
               "matching type for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull =
        NAMELIST_GetStandard_COMPUTE_PERIODIC_IsoDelaunayDomains();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialPeriodic_IsoDelaunayDomain [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    //
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    //
    std::cerr
        << "Normal termination of LATT_SerialPeriodic_IsoDelaunayDomain\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialPeriodic_IsoDelaunayDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}
