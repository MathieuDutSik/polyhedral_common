// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "PeriodicStructures.h"
#include "CoveringMaxdet.h"
#include "QuantizerLtypeExport.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_iarchive.hpp>
// clang-format on

/*
  Reads an iso-Delaunay domain (boost text_oarchive of IsoDelaunayDomain<T,
  Tint, Tgroup> as written by LATT_SerialComputeDelaunay's FileIsoDelaunayDomain
  QUERIES option or by LATT_SerialLattice_IsoDelaunayDomain's per-domain dump)
  and always reports the cheap combinatorial data:
    * |ListIneq|   — full set of defining inequalities
                     (ComputeDefiningIneqIsoDelaunayDomain)
    * |ListIrred|  — irredundant inequalities (get_non_redundant_indices)

  Everything that needs the extreme rays of the cone (the dual description) or
  is otherwise expensive is a QUERIES option, computed only when its entry is
  set to something other than "null":

    * FileFacets            — writes the irredundant defining inequalities of
                              the L-type cone (FACred, i.e. its facets) as a
                              matrix (CPP format, one facet inequality per row).
                              Cheap: it needs only the redundancy elimination,
                              not the dual description.
    * FileFullRankRays      — writes the full-rank (rank n, positive definite)
                              extreme rays — the "rigid" forms on the domain's
                              boundary — as a GAP list of Gram matrices.
                              Collecting these over every iso-Delaunay domain and
                              reducing up to arithmetic equivalence gives the
                              rigid lattices of the T-space (7 in the classic
                              5-dimensional case).
    * FileRankTally         — writes, as a GAP list of [rank, count] pairs, the
                              tally of how many extreme rays r fall in each rank
                              bucket of sum_u r_u * LinSpa.ListMat[u].
    * FileInnerGramMat      — writes a Gram matrix in the relative interior of
                              the L-type cone (get_interior_gram_matrix_lp).
    * FileQuantizerLtypeJSON — writes the (simplex-Delaunay) L-type as JSON for
                              the Julia downstream (formerly the standalone
                              LATT_ExportQuantizerLtype program).
    * FileCoveringOptimum   — minimizes the sphere covering density over the
                              domain by the determinant maximization ("MAXDET")
                              formulation of the L-type covering problem (see
                              CoveringMaxdet.h) and writes the optimum as a GAP
                              record. The solve is numerical, so the covering
                              density it reports is an estimate rather than a
                              certificate. Works verbatim on the iso-Delaunay
                              domains of periodic point sets, which are stored
                              in the same IsoDelaunayDomain objects.

  FileFullRankRays and FileRankTally share the extreme-ray enumeration, so the
  dual description is computed once if either is requested.

  Driven by a namelist (DATA: arithmetic, FileIsoDelaunay, FileLinSpace;
  SYSTEM: OutFile; QUERIES: the entries above). The T-space has to be the one
  the iso-Delaunay domain was computed with: by default the canonical GL_n(Z)
  one (ComputeCanonicalSpace), or the one read from FileLinSpace.
 */

template <typename T, typename Tint>
void process(FullNamelist const &eFull, std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ostream &os = std::cerr;
  //
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockQUERIES = eFull.get_block("QUERIES");
  std::string FileIso = BlockDATA.get_string("FileIsoDelaunay");
  std::string FileFacets = BlockQUERIES.get_string("FileFacets");
  std::string FileFullRankRays = BlockQUERIES.get_string("FileFullRankRays");
  std::string FileRankTally = BlockQUERIES.get_string("FileRankTally");
  std::string FileInnerGramMat = BlockQUERIES.get_string("FileInnerGramMat");
  std::string FileQuantizerLtypeJSON =
      BlockQUERIES.get_string("FileQuantizerLtypeJSON");
  std::string FileCoveringOptimum =
      BlockQUERIES.get_string("FileCoveringOptimum");
  std::string FileLinSpace = BlockDATA.get_string("FileLinSpace");
  std::string FileCosets = BlockDATA.get_string("FileCosets");
  //
  IsoDelaunayDomain<T, Tint, Tgroup> x;
  {
    std::ifstream ifs(FileIso);
    boost::archive::text_iarchive ia(ifs);
    ia >> x;
  }
  int n = x.GramMat.rows();
  int sym_dim = (n * (n + 1)) / 2;
  // The T-space the domain was computed with. The defining inequalities of
  // the domain are expressed in its coordinates, so a mismatch would make
  // every cone computation below meaningless; it is checked right after.
  LinSpaceMatrix<T> LinSpa = FileLinSpace == "null"
                                 ? ComputeCanonicalSpace<T>(n)
                                 : ReadLinSpaceFile<T>(FileLinSpace, os);
  os << "ANA: n=" << n << " sym_dim=" << sym_dim
     << " |DT|=" << x.DT.l_dels.size() << "\n";
  // Always-computed combinatorial data:
  // (1) build the full inequality set,
  // (2) eliminate redundancy with get_non_redundant_indices.
  HumanTime t1;
  std::vector<FullAdjInfo<Tint>> ListIneq =
      ComputeListIneqFromTesselationIneq(x.DT);
  int n_ineq = ListIneq.size();
  os << "ANA: ComputeDefiningIneqIsoDelaunayDomain done n_ineq=" << n_ineq
     << " |elapsed|=" << t1 << "\n";
  int dimSpaceIneq = ListIneq[0].eIneq.size();
  if (dimSpaceIneq != static_cast<int>(LinSpa.ListMat.size())) {
    std::cerr << "LATT_AnalysisIsoDelaunay: the defining inequalities live in "
                 "dimension "
              << dimSpaceIneq << " while the T-space has dimension "
              << LinSpa.ListMat.size()
              << ". Set FileLinSpace to the T-space the iso-Delaunay domain "
                 "was computed with\n";
    throw TerminalException{1};
  }
  //
  HumanTime t2;
  MyMatrix<T> FAC = UniversalMatrixConversion<T, Tint>(GetFACineq(ListIneq));
  std::vector<int> ListIrred = get_non_redundant_indices(FAC, os);
  int n_irred = ListIrred.size();
  MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
  os << "ANA: get_non_redundant_indices done n_irred=" << n_irred
     << " |elapsed|=" << t2 << "\n";
  //
  os_out << "n=" << n << " sym_dim=" << sym_dim << "\n";
  os_out << "n_ineq=" << n_ineq << "\n";
  os_out << "n_irred=" << n_irred << "\n";
  //
  // QUERIES option: the facets of the L-type cone, i.e. the irredundant
  // defining inequalities FACred (one facet inequality per row). Cheap: it
  // reuses the redundancy elimination already done above.
  if (FileFacets != "null") {
    WriteMatrixFile(FileFacets, FACred);
  }
  //
  // QUERIES options over the extreme rays (rank tally and/or full-rank rays).
  // Both need the extreme rays, so the (potentially expensive) dual description
  // on FACred is computed once, only when at least one of them is requested.
  bool need_rank_tally = FileRankTally != "null";
  bool need_full_rank_rays = FileFullRankRays != "null";
  if (need_rank_tally || need_full_rank_rays) {
    HumanTime t3;
    MyMatrix<T> EXT = DirectDualDescription_mat(FACred, os);
    int n_ray = EXT.rows();
    os << "ANA: DirectDualDescription_mat done n_ray=" << n_ray
       << " |elapsed|=" << t3 << "\n";
    os_out << "n_ray=" << n_ray << "\n";
    std::map<int, int> rank_tally;
    // The full-rank (rank n, positive definite) extreme rays: the "rigid" forms
    // living on this domain's boundary. Kept in primitive (fraction-free) form
    // so that a ray and its positive multiples collapse to a single matrix.
    std::vector<MyMatrix<T>> full_rank_rays;
    int dimSpace = LinSpa.ListMat.size();
    for (int i_row = 0; i_row < n_ray; i_row++) {
      MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
      for (int u = 0; u < dimSpace; u++) {
        RayMat += EXT(i_row, u) * LinSpa.ListMat[u];
      }
      int rk = RankMat(RayMat);
      if (need_rank_tally) {
        rank_tally[rk]++;
      }
      if (need_full_rank_rays && rk == n && IsPositiveDefinite(RayMat, os)) {
        full_rank_rays.push_back(RemoveFractionMatrix(RayMat));
      }
    }
    if (need_rank_tally) {
      // GAP list of [rank, count] pairs.
      std::ofstream os_rt(FileRankTally);
      os_rt << "return [";
      bool first = true;
      for (auto &kv : rank_tally) {
        if (!first) {
          os_rt << ", ";
        }
        first = false;
        os_rt << "[" << kv.first << ", " << kv.second << "]";
      }
      os_rt << "];\n";
    }
    if (need_full_rank_rays) {
      std::string FileName = FileFullRankRays;
      WriteListMatrixFileGAP(FileName, full_rank_rays);
    }
  }
  //
  // QUERIES option: the interior (inner) Gram matrix of the iso-Delaunay domain,
  // i.e. a Gram matrix in the relative interior of the L-type cone.
  if (FileInnerGramMat != "null") {
    MyMatrix<T> InnerGram = get_interior_gram_matrix_lp(LinSpa, FAC, os);
    WriteMatrixFile(FileInnerGramMat, InnerGram);
  }
  //
  // QUERIES option: JSON export of the (simplex-Delaunay) L-type for the Julia
  // downstream (formerly the standalone LATT_ExportQuantizerLtype program).
  if (FileQuantizerLtypeJSON != "null") {
    quantizer_export::WriteQuantizerLtypeJSON<T, Tint, Tgroup>(
        x, LinSpa, FileQuantizerLtypeJSON, os);
  }
  //
  // QUERIES option: minimization of the sphere covering density over the
  // domain by the determinant maximization formulation (CoveringMaxdet.h).
  // The problem is built exactly and solved numerically, so the reported
  // density is an estimate. FACred rather than FAC is used: the redundant
  // inequalities would only add barrier terms without changing the feasible
  // set, at the cost of speed and of conditioning.
  if (FileCoveringOptimum != "null") {
    using Tfloat = double;
    HumanTime t4;
    // The covering density formula differs between a lattice and a periodic
    // point set only by the constant factor m / N^n of CoveringMaxdet.h,
    // recovered here from the very coset file the periodic enumeration was
    // run with.
    T point_density(1);
    if (FileCosets != "null") {
      MyMatrix<T> Cosets = ReadMatrixFile<T>(FileCosets);
      if (Cosets.cols() != n) {
        std::cerr << "LATT_AnalysisIsoDelaunay: the cosets have "
                  << Cosets.cols() << " columns while the domain is of "
                  << "dimension " << n << "\n";
        throw TerminalException{1};
      }
      PeriodicPointSet<Tint> pps =
          PeriodicPointSetFromRational<Tint, T>(Cosets);
      T N_T = UniversalScalarConversion<T, Tint>(pps.N);
      T denom(1);
      for (int i = 0; i < n; i++) {
        denom *= N_T;
      }
      point_density = T(pps.cosets_num.rows()) / denom;
      os << "ANA: periodic point set, n_coset=" << pps.cosets_num.rows()
         << " N=" << pps.N << " point_density=" << point_density << "\n";
    }
    covering_maxdet::CoveringData<T> cd = covering_maxdet::BuildCoveringData(
        x, LinSpa, FACred, point_density, os);
    covering_maxdet::CoveringData<Tfloat> cd_f =
        covering_maxdet::ConvertCoveringData<Tfloat, T>(cd);
    MyVector<Tfloat> x_start =
        covering_maxdet::GetStartingPoint<T, Tfloat>(cd, cd_f, os);
    covering_maxdet::MaxdetOptions<Tfloat> opts =
        covering_maxdet::GetDefaultMaxdetOptions<Tfloat>();
    covering_maxdet::MaxdetResult<Tfloat> res =
        covering_maxdet::SolveCoveringMaxdet(cd_f, x_start, opts, os);
    os << "ANA: covering optimization done success=" << res.success
       << " message=" << res.message << " |elapsed|=" << t4 << "\n";
    covering_maxdet::WriteCoveringOptimumGAP(FileCoveringOptimum, res);
    std::streamsize prec = os_out.precision();
    os_out << std::setprecision(17);
    os_out << "covering_success=" << res.success << "\n";
    os_out << "covering_density=" << res.cov_density << "\n";
    os_out << "covering_radius_sq=" << res.cov_radius_sq << "\n";
    os_out << "covering_gap_bound=" << res.gap_bound << "\n";
    os_out << std::setprecision(prec);
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
    // The T-space the iso-Delaunay domain was computed with, in the format
    // of WriteLinSpaceFile. "null" (the default) uses the canonical GL_n(Z)
    // T-space, which is the right one for the enumerations run over the full
    // space of quadratic forms.
    ListStringValues["FileLinSpace"] = "null";
    // For a domain coming from LATT_SerialPeriodic_IsoDelaunayDomain, the
    // very coset file that enumeration was run with. It is what tells the
    // covering density apart from the lattice one. "null" (the default)
    // treats the domain as a lattice one.
    ListStringValues["FileCosets"] = "null";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // QUERIES
  {
    std::map<std::string, std::string> ListStringValues;
    // The facets (irredundant defining inequalities FACred) of the L-type cone,
    // as a matrix in CPP format. "null" (the default) skips the writing.
    ListStringValues["FileFacets"] = "null";
    // The full-rank (rigid) extreme rays of the domain, as a GAP list of Gram
    // matrices. "null" (the default) skips the computation.
    ListStringValues["FileFullRankRays"] = "null";
    // The tally of extreme-ray ranks, as a GAP list of [rank, count] pairs.
    // "null" (the default) skips the computation.
    ListStringValues["FileRankTally"] = "null";
    // A Gram matrix in the relative interior of the L-type cone. "null" (the
    // default) skips the computation.
    ListStringValues["FileInnerGramMat"] = "null";
    // JSON export of the (simplex-Delaunay) L-type for the Julia downstream.
    // "null" (the default) skips the computation.
    ListStringValues["FileQuantizerLtypeJSON"] = "null";
    // The covering density optimum over the domain, as a GAP record.
    // "null" (the default) skips the computation.
    ListStringValues["FileCoveringOptimum"] = "null";
    SingleBlock BlockQUERIES;
    BlockQUERIES.setListStringValues(ListStringValues);
    ListBlock["QUERIES"] = BlockQUERIES;
  }
  return FullNamelist(ListBlock);
}

void process_C(FullNamelist const &eFull) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  std::string arith = BlockDATA.get_string("arithmetic");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  auto f = [&](std::ostream &os_out) -> void {
    if (arith == "gmp") {
      using T = mpq_class;
      using Tint = mpz_class;
      return process<T, Tint>(eFull, os_out);
    }
#ifdef ENABLE_ALL_NUMERICAL_TYPES
    if (arith == "gmp_boost") {
      using T = boost::multiprecision::mpq_rational;
      using Tint = boost::multiprecision::mpz_int;
      return process<T, Tint>(eFull, os_out);
    }
    if (arith == "multi_boost") {
      using T = boost::multiprecision::cpp_rational;
      using Tint = boost::multiprecision::cpp_int;
      return process<T, Tint>(eFull, os_out);
    }
#endif
    std::cerr << "LATT_AnalysisIsoDelaunay: Failed to find a matching type for "
                 "arithmetic="
              << arith << "\n";
    std::cerr << "Available types: gmp, gmp_boost, multi_boost\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f);
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
