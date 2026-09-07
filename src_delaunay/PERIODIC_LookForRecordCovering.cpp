// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "CoveringRecordSearch.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  Random-walk search for a periodic point set whose covering density beats
  the best known lattice covering of its dimension.

  Given a periodic point set Z^n + {c_1, ..., c_m} through its coset file,
  the parameter space is the cone of positive definite forms, subdivided
  into iso-Delaunay domains. Over each of them the covering density is
  minimized by the determinant maximization of CoveringMaxdet.h, and the
  covering density of the point set is the minimum over every domain. The
  domains being far too numerous to enumerate in the dimensions of interest,
  the search walks the adjacency graph instead, descending on the per-domain
  optimum and jumping randomly out of the local minima; see
  CoveringRecordSearch.h.

  Not finding a record is the expected outcome: for 3 <= n <= 5 the best
  covering is conjectured to be the lattice one, A_n^*. The program says so
  and terminates normally, reporting the best density it saw.

  Driven by a namelist:
    SYSTEM   max_runtime_second is the search budget, Prefix the directory
             where the Gram matrix of every improvement is written, OutFile
             the GAP record of the outcome.
    DATA     arithmetic, FileDualDescription, FileCosets.
    SEARCH   RecordToBeat and n_walk_steps.
    TSPACE   the T-space, which has to be the one the point set lives in.
 */

template <typename T, typename Tint>
void process_A(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ostream &os = std::cerr;

  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockSEARCH = eFull.get_block("SEARCH");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  int n = LinSpa.n;
  int dimEXT = n + 1;
  //
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  std::string Prefix = BlockSYSTEM.get_string("Prefix");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  int n_walk_steps = BlockSEARCH.get_int("n_walk_steps");
  std::string RecordToBeat = BlockSEARCH.get_string("RecordToBeat");
  //
  // The record. "auto" takes the least dense known lattice covering of the
  // dimension, which is what a periodic covering has to beat.
  //
  double record;
  if (RecordToBeat == "auto") {
    std::optional<double> opt =
        covering_record::GetKnownLatticeCoveringRecord(n);
    if (!opt) {
      std::cerr << "PERIODIC_LookForRecordCovering: no known lattice covering "
                   "record is recorded for dimension "
                << n
                << ", so RecordToBeat has to be set to an explicit value\n";
      throw TerminalException{1};
    }
    record = *opt;
  } else {
    record = ParseScalar<double>(RecordToBeat);
  }
  os << "PERIODIC_RECORD: n=" << n << " record=" << record << "\n";
  //
  // The periodic point set, with the same validation as the enumeration:
  // cosets that form a group modulo Z^n describe a lattice, and the
  // pointwise stabilizer of the T-space has to stabilize the point set.
  //
  std::string FileCosets = BlockDATA.get_string("FileCosets");
  MyMatrix<T> Cosets = ReadMatrixFile<T>(FileCosets);
  if (Cosets.cols() != n) {
    std::cerr << "PERIODIC_LookForRecordCovering: the cosets have "
              << Cosets.cols() << " columns while the T-space is of dimension "
              << n << "\n";
    throw TerminalException{1};
  }
  std::optional<PeriodicPointSet<Tint>> opt_pps =
      PeriodicPointSetFromRational_Opt<Tint, T>(Cosets);
  if (!opt_pps) {
    std::cerr << "PERIODIC_LookForRecordCovering: the cosets form a group "
                 "modulo Z^n, so the point set is a lattice and cannot beat "
                 "the lattice record by construction\n";
    throw TerminalException{1};
  }
  PeriodicPointSet<Tint> pps = *opt_pps;
  for (auto &eGen : LinSpa.PtStabGens) {
    if (!PeriodicAffineExtension<Tint, T>(pps, eGen)) {
      std::cerr << "PERIODIC_LookForRecordCovering: a generator of the "
                   "pointwise stabilizer of the T-space admits no affine "
                   "extension preserving the point set. The point stabilizer "
                   "has to stabilize the cosets for the walk to be correct. "
                   "The offending generator is\n";
      WriteMatrix(std::cerr, eGen);
      throw TerminalException{1};
    }
  }
  // The covering density of a periodic point set is the lattice formula
  // multiplied by m / N^n, see CoveringMaxdet.h.
  T N_T = UniversalScalarConversion<T, Tint>(pps.N);
  T denom(1);
  for (int i = 0; i < n; i++) {
    denom *= N_T;
  }
  T point_density = T(pps.cosets_num.rows()) / denom;
  os << "PERIODIC_RECORD: n_coset=" << pps.cosets_num.rows()
     << " N=" << pps.N << " point_density=" << point_density << "\n";
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  std::optional<MyMatrix<T>> CommonGramMat;
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  LinSpaceMatrix<Tint> LinSpaRing = LINSPA_GetRingVersion(LinSpa);
  std::optional<MyMatrix<Tint>> CommonGramMatRing =
      GetCommonGramMatRing<Tint, T>(CommonGramMat);
  LinSpaceMatrix<T> LinSpaCopy = LinSpa;
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{
      std::move(LinSpa), std::move(LinSpaRing), std::move(rddo), CommonGramMat,
      CommonGramMatRing};
  //
  // The stabilizer of a domain for the point set rather than for the
  // lattice: this is the single place where the walk differs from the
  // lattice one, the wall and flip kernels below it being shared.
  //
  using Tdom = IsoDelaunayDomain<T, Tint, Tgroup>;
  auto f_stab_gens =
      [&](Tdom const &x) -> std::vector<MyMatrix<Tint>> {
    Result_ComputeStabilizer_SHV<Tint, Tgroup> result =
        LINSPA_ComputeStabilizer_SHV_Periodic<Tint, Tint, Tgroup>(
            data.LinSpaRing, pps, x.GramMat, x.SHV, data.CommonGramMatRing, os);
    return result.get_list_matrix(x.SHV, x.GramMat, data.LinSpaRing, os);
  };
  //
  if (Prefix != "/irrelevant/") {
    FILE_CreateDirectory(Prefix);
  } else {
    Prefix = "";
  }
  covering_record::RecordSearchOptions opts{record, n_walk_steps,
                                            max_runtime_second, Prefix};
  IsoDelaunayDomain<T, Tint, Tgroup> start =
      GetInitialPeriodicIsoDelaunayDomain(data, pps);
  covering_record::RecordSearchResult<Tint> res =
      covering_record::LookForRecordCovering<T, Tint, Tgroup,
                                             decltype(f_stab_gens)>(
          data, start, LinSpaCopy, point_density, f_stab_gens, opts, os);
  //
  covering_record::WriteRecordSearchGAP(OutFile, res, record);
  os << "PERIODIC_RECORD: " << res.message << "\n";
  if (res.found_record) {
    os << "PERIODIC_RECORD: FOUND a periodic covering of density "
       << res.best_density << " beating the record " << record << "\n";
  } else if (res.has_best) {
    os << "PERIODIC_RECORD: no record found, the best density seen is "
       << res.best_density << " against the record " << record << "\n";
  } else {
    os << "PERIODIC_RECORD: no domain could be optimized at all\n";
  }
  os << "PERIODIC_RECORD: n_iter=" << res.n_iter << " n_walk=" << res.n_walk
     << " n_domain_evaluated=" << res.n_domain_evaluated
     << " n_domain_failed=" << res.n_domain_failed
     << " runtime_second=" << res.runtime_second << "\n";
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic = GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_A<T, Tint>(eFull);
  }
  std::cerr << "PERIODIC_LookForRecordCovering: Failed to find a matching "
               "type for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp\n";
  throw TerminalException{1};
}

FullNamelist NAMELIST_GetStandard_PERIODIC_RECORD_COVERING() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["FileDualDescription"] = "unset";
    // The file with the rational coset matrix of the periodic point set:
    // one coset per row, the zero coset included.
    ListStringValues["FileCosets"] = "unset";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // SEARCH
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    // The covering density to beat. "auto" (the default) takes the least
    // dense known lattice covering of the dimension: A_n^* up to dimension
    // 5, then the L^c_n of Schuermann and Vallentin.
    ListStringValues["RecordToBeat"] = "auto";
    // The number of random adjacency jumps taken to leave a local minimum.
    ListIntValues["n_walk_steps"] = 20;
    SingleBlock BlockSEARCH;
    BlockSEARCH.setListStringValues(ListStringValues);
    BlockSEARCH.setListIntValues(ListIntValues);
    ListBlock["SEARCH"] = BlockSEARCH;
  }
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  return FullNamelist(ListBlock);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_PERIODIC_RECORD_COVERING();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PERIODIC_LookForRecordCovering [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    unsigned seed = get_random_seed();
    std::cerr << "seed=" << seed << "\n";
    srand(seed);
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    //
    std::cerr << "Normal termination of PERIODIC_LookForRecordCovering\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PERIODIC_LookForRecordCovering\n";
    exit(e.eVal);
  }
  runtime(time);
}
