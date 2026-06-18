// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizerLtypeExport.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

FullNamelist NAMELIST_GetStandard_EXPORT_QUANTIZER_LTYPE() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["OutFile"] = "ltype.json";
    SingleBlock BlockSYSTEM;
    BlockSYSTEM.setListStringValues(ListStringValues);
    ListBlock["SYSTEM"] = BlockSYSTEM;
  }
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["FileDualDescription"] = "unset";
    ListStringValues["FileGramMat"] = "unset";
    ListStringValues["CVPmethod"] = "SVexact";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  // TSPACE
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  return FullNamelist(ListBlock);
}

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

  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  std::string FileGramMat = BlockDATA.get_string("FileGramMat");
  std::cerr << "QuantExport: OutFile=" << OutFile << "\n";

  int dimEXT = LinSpa.n + 1;
  std::cerr << "QuantExport: reading dual-desc heuristics, dimEXT=" << dimEXT
            << "\n";
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         std::cerr);

  // For the quantizer use case we want to allow the full T-space S^n (which
  // GetInitialIsoDelaunayDomain rejects via is_integrally_saturated_matrix_space),
  // and we need a sample Q_0 interior to the L-type cone of interest (e.g.
  // V^n's principal type, where all Delaunay polytopes are simplices).
  MyMatrix<T> GramMat;
  if (FileGramMat != "unset") {
    GramMat = ReadMatrixFile<T>(FileGramMat);
    std::cerr << "QuantExport: read sample Q_0 from " << FileGramMat << "\n";
  } else {
    GramMat = LinSpa.SuperMat;
    std::cerr << "QuantExport: using LinSpa.SuperMat as sample Q_0 "
                 "(may not be in the desired L-type cone)\n";
  }

  DataLattice<T, Tint, Tgroup> data_lattice =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);
  auto f_incorrect =
      [&]([[maybe_unused]] Delaunay_Obj<T, Tgroup> const &x) -> bool {
    return false;
  };
  std::cerr << "QuantExport: enumerating Delaunay polytopes for this Gram matrix\n";
  std::optional<DelaunayTesselation<T, Tgroup>> opt_DT =
      EnumerationDelaunayPolytopes<T, Tint, Tgroup, decltype(f_incorrect)>(
          data_lattice, f_incorrect, /*max_runtime_second=*/0);
  if (!opt_DT) {
    std::cerr << "QuantExport: Delaunay enumeration returned nothing\n";
    throw TerminalException{1};
  }
  DelaunayTesselation<T, Tgroup> DT = std::move(*opt_DT);
  std::cerr << "QuantExport: have |DT.l_dels|=" << DT.l_dels.size() << "\n";

  // Build a minimal IsoDelaunayDomain. The SHV_T field is not used by the
  // exporter; provide an empty matrix to satisfy the struct.
  MyMatrix<T> SHV_T(0, LinSpa.n);
  IsoDelaunayDomain<T, Tint, Tgroup> IDD{std::move(DT), GramMat,
                                         std::move(SHV_T)};

  quantizer_export::WriteQuantizerLtypeJSON<T, Tint, Tgroup>(
      IDD, LinSpa, OutFile, std::cerr);
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
  std::cerr << "LATT_ExportQuantizerLtype: Unknown arithmetic=" << arithmetic
            << "\n";
  std::cerr << "Available types: gmp, gmp_boost\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_EXPORT_QUANTIZER_LTYPE();
    if (argc != 2) {
      std::cerr << "Usage: LATT_ExportQuantizerLtype [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_ExportQuantizerLtype\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ExportQuantizerLtype\n";
    exit(e.eVal);
  }
  runtime(time);
}
