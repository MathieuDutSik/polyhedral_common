// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "IsoKDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_oarchive.hpp>
// clang-format on

/*
  Enumeration of the top-dimensional iso-k-Delaunay domains ((L,k)-types)
  of a T-space, up to the equivalence of the space, with the optional
  minimization of the k-covering density over each domain.
 */

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
  DataIsoKDelaunayDomains<T, Tint, Tgroup> data =
      get_data_iso_k_delaunay_domains<T, Tint, Tgroup>(eFull, AllArr,
                                                       std::cerr);
  int k = data.k;
  using Tdata = DataIsoKDelaunayDomainsFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
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
      unfold_opt(opt_l_tot, "EnumerateAndStore_Serial (iso-k-Delaunay)");
  //
  std::string PrefixIsoKDel =
      BlockDATA.get_string("PrefixIsoKDelaunayDomains");
  if (PrefixIsoKDel != "unset") {
    for (size_t i = 0; i < l_tot.size(); i++) {
      std::string FileName = PrefixIsoKDel + std::to_string(i);
      std::ofstream ofs(FileName);
      boost::archive::text_oarchive oa(ofs);
      oa << l_tot[i].x.dom;
    }
  }
  //
  std::string FileCoveringOptimum = BlockDATA.get_string("FileCoveringOptimum");
  if (FileCoveringOptimum != "null") {
    std::ofstream os_out(FileCoveringOptimum);
    os_out << "return [";
    for (size_t i = 0; i < l_tot.size(); i++) {
      if (i > 0) {
        os_out << ",\n";
      }
      covering_maxdet::MaxdetResult<double> res =
          OptimizeKCovering<T, Tint, Tgroup>(l_tot[i].x.dom, LinSpa, std::cerr);
      WriteKCoveringOptimumRecordGAP(os_out, k, res);
    }
    os_out << "];\n";
  }
  //
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
#ifdef ENABLE_BOOST_TYPES
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
#endif
#ifdef ENABLE_FLINT_SUPPORT
  if (arithmetic == "flint") {
    using T = fmpq_class;
    using Tint = fmpz_class;
    return process_A<T, Tint>(eFull);
  }
#endif
  std::cerr << "LATT_SerialLattice_IsoKDelaunayDomain: Failed to find a "
               "matching type for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp";
#ifdef ENABLE_BOOST_TYPES
  std::cerr << ", gmp_boost, multi_boost";
#endif
#ifdef ENABLE_FLINT_SUPPORT
  std::cerr << ", flint\n";
#else
  std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull =
        NAMELIST_GetStandard_COMPUTE_LATTICE_IsoKDelaunayDomains();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialLattice_IsoKDelaunayDomain [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_SerialLattice_IsoKDelaunayDomain\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialLattice_IsoKDelaunayDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}
