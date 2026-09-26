// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "IsoKDelaunayDomains.h"
#include "Tspace_Generation.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_oarchive.hpp>
// clang-format on

/*
  The order-k Delaunay tiling of a lattice given by its Gram matrix: the
  tiles (P_-, P_0) up to the affine isometries of the lattice, the
  k-covering radius and density (OutFormat = GAP_Covering), and the queries
  on the tiling:
  * FileTilingVolume: the volume identity of the tiling, a check of the
    enumeration.
  * FileIsoKDelaunayDomain: the (L,k)-domain of the form in the space of
    all symmetric matrices, as a boost text-archive.
  * FileCoveringOptimum: the minimization of the k-covering density over
    that domain.
 */

FullNamelist NAMELIST_GetStandard_SERIAL_COMPUTE_K_DELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["GRAMfile"] = "unset.gram";
    ListStringValues["SVRfile"] = "unset.svr";
    ListStringValues["choice_initial"] = "direct";
    ListStringValues["FileDualDescription"] = "unset";
    ListStringValues["CacheFile"] = "none";
    std::map<std::string, int> ListIntValues;
    ListIntValues["k"] = 2;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // QUERIES
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["FileTilingVolume"] = "null";
    ListStringValues["FileIsoKDelaunayDomain"] = "null";
    ListStringValues["FileCoveringOptimum"] = "null";
    SingleBlock BlockQUERIES;
    BlockQUERIES.setListStringValues(ListStringValues);
    ListBlock["QUERIES"] = BlockQUERIES;
  }
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint, typename Tgroup>
void process_A(FullNamelist const &eFull, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockQUERIES = eFull.get_block("QUERIES");
  //
  std::string GRAMfile = BlockDATA.get_string("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  int n = GramMat.rows();
  int dimEXT = n + 1;
  int k = BlockDATA.get_int("k");
  if (k < 1) {
    std::cerr << "LATT_SerialComputeKDelaunay: k=" << k
              << " should be at least 1\n";
    throw TerminalException{1};
  }
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  DataLattice<T, Tint, Tgroup> data =
      get_data_lattice<T, Tint, Tgroup>(eFull, AllArr, os);
  //
  std::string CacheFile = BlockDATA.get_string("CacheFile");
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  KDelaunayTesselation<Tint, Tgroup> DT =
      get_k_delaunay_tessellation_serial<T, Tint, Tgroup>(
          data, k, CacheFile, max_runtime_second, os);
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  auto f = [&](std::ostream &os_out) -> void {
    WriteKDelaunayTesselation(OutFormat, os_out, GramMat, DT);
  };
  FILE_PrintStderrStdoutFile(OutFile, f);
  //
  std::string FileTilingVolume = BlockQUERIES.get_string("FileTilingVolume");
  if (FileTilingVolume != "null") {
    KTilingVolumeCheck<T, Tint, Tgroup> check =
        CheckKDelaunayTesselationVolume<T, Tint, Tgroup>(GramMat, DT, os);
    std::ofstream os_out(FileTilingVolume);
    WriteKTilingVolumeCheckGAP(os_out, check);
  }
  //
  std::string FileIsoKDelaunayDomain =
      BlockQUERIES.get_string("FileIsoKDelaunayDomain");
  std::string FileCoveringOptimum =
      BlockQUERIES.get_string("FileCoveringOptimum");
  if (FileIsoKDelaunayDomain != "null" || FileCoveringOptimum != "null") {
    // The domain of the form in the space of all symmetric matrices. It is
    // full dimensional only when no tile is co-spherical by accident, i.e.
    // when every P_0 is a simplex.
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
    std::vector<std::vector<Tint>> ListGramRing =
        GetListGramRing(LinSpa.ListLineMat);
    if (IsKDelaunayTesselationInducingEqualities(DT, ListGramRing, os)) {
      std::cerr << "LATT_SerialComputeKDelaunay: a tile of the tiling has "
                   "more than n+1 points on its sphere, so the form is on a "
                   "wall and its (L,k)-domain is not full dimensional\n";
      throw TerminalException{1};
    }
    IsoKDelaunayDomain<T, Tint, Tgroup> dom =
        BuildIsoKDelaunayDomain<T, Tint, Tgroup>(DT, GramMat, LinSpa,
                                                 ListGramRing, os);
    if (FileIsoKDelaunayDomain != "null") {
      std::ofstream ofs(FileIsoKDelaunayDomain);
      boost::archive::text_oarchive oa(ofs);
      oa << dom;
    }
    if (FileCoveringOptimum != "null") {
      covering_maxdet::MaxdetResult<double> res =
          OptimizeKCovering<T, Tint, Tgroup>(dom, LinSpa, os);
      std::ofstream os_out(FileCoveringOptimum);
      os_out << "return ";
      WriteKCoveringOptimumRecordGAP(os_out, k, res);
      os_out << ";\n";
    }
  }
}

template <typename T, typename Tint> void process_B(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return process_A<T, Tint, Tgroup>(eFull, std::cerr);
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic = GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(eFull);
  }
#ifdef ENABLE_BOOST_TYPES
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_B<T, Tint>(eFull);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_B<T, Tint>(eFull);
  }
#endif
#ifdef ENABLE_FLINT_SUPPORT
  if (arithmetic == "flint") {
    using T = fmpq_class;
    using Tint = fmpz_class;
    return process_B<T, Tint>(eFull);
  }
#endif
  std::cerr << "LATT_SerialComputeKDelaunay: Failed to find a matching type "
               "for arithmetic="
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
    FullNamelist eFull = NAMELIST_GetStandard_SERIAL_COMPUTE_K_DELAUNAY();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialComputeKDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_SerialComputeKDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialComputeKDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
