// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "lt_complex.h"
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
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         std::cerr);

  bool compute_full_dimensional =
      BlockDATA.get_bool("compute_full_dimensional");

  DataIsoDelaunayDomains<T, Tint, Tgroup> data =
      get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, std::cerr);

  // Always start from the iso-Delaunay enumeration of the top-dimensional
  // cells; the `compute_full_dimensional` flag only decides what we then
  // do with that database.
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

  if (!compute_full_dimensional) {
    // Backward-compatible path: just emit the top-dimensional iso-Delaunay
    // domains from l_tot.
    std::ofstream os_out(OutFile);
    bool result =
        WriteFamilyObjects(data_func.data, OutFormat, os_out, l_tot, std::cerr);
    if (result) {
      std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    }
    return;
  }
  // compute_full_dimensional = T: enumerate every cell of the L-type
  // complex (top-dim + lower-dim) and keep only the cells whose relative
  // interior contains positive-definite forms.
  LtypeComplexOptions opts{compute_full_dimensional};
  FullLtypeComplexEnumeration<T> flce =
      get_full_ltype_complex_enumeration_from_ltot<T, Tint, Tgroup>(
          l_tot, data_func.data.LinSpa, opts, std::cerr);
  std::ofstream os_out(OutFile);
  if (OutFormat == "GAP") {
    WriteFullLtypeComplexEnumerationGAP(os_out, flce);
  } else if (OutFormat == "NumberGAP") {
    size_t total = 0;
    os_out << "return rec(per_idx:=[";
    int n_level = flce.levels.size();
    for (int i_level = 0; i_level < n_level; i_level++) {
      if (i_level > 0)
        os_out << ", ";
      size_t cnt = flce.levels[i_level].l_cells.size();
      total += cnt;
      os_out << "rec(idx:=" << flce.levels[i_level].idx << ", count:=" << cnt
             << ")";
    }
    os_out << "], nb:=" << total << ");\n";
  } else {
    std::cerr << "LATT_SerialLattice_IsoDelaunayDomain: Unsupported OutFormat="
              << OutFormat << " for compute_full_dimensional=T\n";
    std::cerr << "Supported: GAP, NumberGAP\n";
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
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_A<T, Tint>(eFull);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_A<T, Tint>(eFull);
  }
  if (arithmetic == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_A<T, Tint>(eFull);
  }
  std::cerr << "LATT_SerialLattice_IsoDelaunayDomain: Failed to find a "
               "matching type for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp, gmp_boost, multi_boost\n";
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
