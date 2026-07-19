// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_iarchive.hpp>
// clang-format on

/*
  Reads an iso-Delaunay domain (boost text_oarchive of IsoDelaunayDomain<T,
  Tint, Tgroup> as written by LATT_SerialComputeDelaunay's FileIsoDelaunayDomain
  QUERIES option or by LATT_SerialLattice_IsoDelaunayDomain's per-domain dump)
  and reports:
    * |ListIneq|   — full set of defining inequalities
                     (ComputeDefiningIneqIsoDelaunayDomain)
    * |ListIrred|  — irredundant inequalities (get_non_redundant_indices)
    * For each extreme ray r of the cone, the rank of
        sum_u r_u * LinSpa.ListMat[u]
      and a tally of how many rays fall in each rank bucket.

  Driven by a namelist (DATA: arithmetic, FileIsoDelaunay; SYSTEM: OutFile;
  QUERIES: FileFullRankRays). The QUERIES entry FileFullRankRays, when not
  "null", writes the full-rank (rank n, positive definite) extreme rays — the
  "rigid" forms on the domain's boundary — as a GAP list of Gram matrices.
  Collecting these over every iso-Delaunay domain and reducing up to arithmetic
  equivalence gives the rigid lattices of the T-space (7 in the classic
  5-dimensional case).
 */

template <typename T, typename Tint>
void process(std::string const &FileIso, std::string const &FileFullRankRays,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ostream &os = std::cerr;
  //
  IsoDelaunayDomain<T, Tint, Tgroup> x;
  {
    std::ifstream ifs(FileIso);
    boost::archive::text_iarchive ia(ifs);
    ia >> x;
  }
  int n = x.GramMat.rows();
  int sym_dim = (n * (n + 1)) / 2;
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
  os << "ANA: n=" << n << " sym_dim=" << sym_dim
     << " |DT|=" << x.DT.l_dels.size() << "\n";
  // Same sequence as IsoDelaunayDomains.h::CountNonFullRankRays:
  // (1) build the full inequality set,
  // (2) eliminate redundancy with get_non_redundant_indices,
  // (3) one dual description on the irredundant FACred.
  HumanTime t1;
  std::vector<FullAdjInfo<T>> ListIneq =
      ComputeListIneqFromTesselationIneq<T, Tgroup>(x.DT);
  int n_ineq = ListIneq.size();
  os << "ANA: ComputeDefiningIneqIsoDelaunayDomain done n_ineq=" << n_ineq
     << " |elapsed|=" << t1 << "\n";
  //
  HumanTime t2;
  MyMatrix<T> FAC = GetFACineq(ListIneq);
  std::vector<int> ListIrred = get_non_redundant_indices(FAC, os);
  int n_irred = ListIrred.size();
  MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
  os << "ANA: get_non_redundant_indices done n_irred=" << n_irred
     << " |elapsed|=" << t2 << "\n";
  //
  HumanTime t3;
  MyMatrix<T> EXT = DirectDualDescription_mat(FACred, os);
  int n_ray = EXT.rows();
  os << "ANA: DirectDualDescription_mat done n_ray=" << n_ray
     << " |elapsed|=" << t3 << "\n";
  //
  std::map<int, int> rank_tally;
  // The full-rank (rank n, positive definite) extreme rays: the "rigid" forms
  // living on this domain's boundary. Kept in primitive (fraction-free) form so
  // that a ray and its positive multiples collapse to a single matrix.
  std::vector<MyMatrix<T>> full_rank_rays;
  int dimSpace = LinSpa.ListMat.size();
  for (int i_row = 0; i_row < n_ray; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      RayMat += EXT(i_row, u) * LinSpa.ListMat[u];
    }
    int rk = RankMat(RayMat);
    rank_tally[rk]++;
    if (rk == n && IsPositiveDefinite(RayMat, os)) {
      full_rank_rays.push_back(RemoveFractionMatrix(RayMat));
    }
  }
  //
  os_out << "n=" << n << " sym_dim=" << sym_dim << "\n";
  os_out << "n_ineq=" << n_ineq << "\n";
  os_out << "n_irred=" << n_irred << "\n";
  os_out << "n_ray=" << n_ray << "\n";
  os_out << "rank_tally=[";
  bool first = true;
  for (auto &kv : rank_tally) {
    if (!first) {
      os_out << ", ";
    }
    first = false;
    os_out << "(rank=" << kv.first << ", count=" << kv.second << ")";
  }
  os_out << "]\n";
  // Optional QUERIES output: the full-rank extreme rays (rigid forms) of this
  // domain, as a GAP list of Gram matrices. Skipped when the entry is "null".
  if (FileFullRankRays != "null") {
    std::string FileName = FileFullRankRays;
    WriteListMatrixFileGAP(FileName, full_rank_rays);
  }
}

FullNamelist NAMELIST_GetStandard_ANALYSIS_ISODELAUNAY() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["FileIsoDelaunay"] = "unset";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // QUERIES
  {
    std::map<std::string, std::string> ListStringValues;
    // The full-rank (rigid) extreme rays of the domain, as a GAP list of Gram
    // matrices. "null" (the default) skips the computation.
    ListStringValues["FileFullRankRays"] = "null";
    SingleBlock BlockQUERIES;
    BlockQUERIES.setListStringValues(ListStringValues);
    ListBlock["QUERIES"] = BlockQUERIES;
  }
  return FullNamelist(ListBlock);
}

void process_C(FullNamelist const &eFull) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockQUERIES = eFull.get_block("QUERIES");
  std::string arith = BlockDATA.get_string("arithmetic");
  std::string FileIso = BlockDATA.get_string("FileIsoDelaunay");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  std::string FileFullRankRays = BlockQUERIES.get_string("FileFullRankRays");
  auto f = [&](std::ostream &os_out) -> void {
    if (arith == "gmp") {
      using T = mpq_class;
      using Tint = mpz_class;
      return process<T, Tint>(FileIso, FileFullRankRays, os_out);
    }
    if (arith == "gmp_boost") {
      using T = boost::multiprecision::mpq_rational;
      using Tint = boost::multiprecision::mpz_int;
      return process<T, Tint>(FileIso, FileFullRankRays, os_out);
    }
    if (arith == "multi_boost") {
      using T = boost::multiprecision::cpp_rational;
      using Tint = boost::multiprecision::cpp_int;
      return process<T, Tint>(FileIso, FileFullRankRays, os_out);
    }
    if (arith == "safe") {
      using T = Rational<SafeInt64>;
      using Tint = SafeInt64;
      return process<T, Tint>(FileIso, FileFullRankRays, os_out);
    }
    std::cerr << "LATT_AnalysisIsoDelaunay: Failed to find a matching type for "
                 "arithmetic="
              << arith << "\n";
    std::cerr << "Available types: gmp, gmp_boost, multi_boost, safe\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ANALYSIS_ISODELAUNAY();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_AnalysisIsoDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_AnalysisIsoDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_AnalysisIsoDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
