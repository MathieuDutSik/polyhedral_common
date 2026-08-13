// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  Enumeration of the Delaunay cells of a periodic point set Z^n + {c_i}
  at a fixed Gram matrix, the periodic counterpart of the plain Delaunay
  enumeration of LATT_SerialComputeDelaunay.
 */

FullNamelist NAMELIST_GetStandard_COMPUTE_PeriodicDelaunay() {
  std::map<std::string, SingleBlock> ListBlock;
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  {
    std::map<std::string, std::string> ListStringValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["GRAMfile"] = "unset";
    // The file with the rational coset matrix of the periodic point set:
    // one coset per row, the zero coset included.
    ListStringValues["FileCosets"] = "unset";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint> void process_A(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  std::ostream &os = std::cerr;
  //
  std::string GRAMfile = BlockDATA.get_string("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  int n = GramMat.rows();
  // The input validation, performed unconditionally.
  if (!IsSymmetricMatrix(GramMat) || !IsPositiveDefinite(GramMat, os)) {
    std::cerr << "LATT_SerialPeriodicDelaunay: the Gram matrix should be "
                 "symmetric positive definite\n";
    throw TerminalException{1};
  }
  std::string FileCosets = BlockDATA.get_string("FileCosets");
  MyMatrix<T> Cosets = ReadMatrixFile<T>(FileCosets);
  if (Cosets.cols() != n) {
    std::cerr << "LATT_SerialPeriodicDelaunay: the cosets have "
              << Cosets.cols() << " columns while the Gram matrix has order "
              << n << "\n";
    throw TerminalException{1};
  }
  std::optional<PeriodicPointSet<Tint>> opt_pps =
      PeriodicPointSetFromRational_Opt<Tint, T>(Cosets);
  if (!opt_pps) {
    std::cerr << "LATT_SerialPeriodicDelaunay: the cosets form a group "
                 "modulo Z^n, so the point set is a lattice: it has to be "
                 "treated by LATT_SerialComputeDelaunay over that lattice\n";
    throw TerminalException{1};
  }
  PeriodicPointSet<Tint> pps = *opt_pps;
  //
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  MyMatrix<Tint> Graver = GetGraverBasis<T, Tint>(GramMat);
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(n + 1, os);
  PeriodicDataDelaunay<T, Tint, Tgroup> data =
      GetPeriodicDataDelaunay<T, Tint, Tgroup>(GramMat, pps, SHV_T, Graver,
                                               AllArr, os);
  PeriodicDataDelaunayFunc<T, Tint, Tgroup> data_func{data};
  using Tobj = PeriodicDelaunay_Obj<Tint, Tgroup>;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  auto opt = EnumerateAndStore_Serial(data_func, f_incorrect,
                                      max_runtime_second);
  if (!opt) {
    std::cerr << "LATT_SerialPeriodicDelaunay: the enumeration did not "
                 "terminate\n";
    throw TerminalException{1};
  }
  auto const &l_ent = *opt;
  //
  std::ofstream os_out(OutFile);
  if (OutFormat == "NumberGAP") {
    os_out << "return rec(nb:=" << l_ent.size() << ");\n";
    return;
  }
  if (OutFormat == "SummaryGAP") {
    os_out << "return rec(nb:=" << l_ent.size() << ", ListRec:=[";
    bool IsFirst = true;
    for (auto &eEnt : l_ent) {
      if (!IsFirst)
        os_out << ",\n";
      IsFirst = false;
      os_out << "rec(nVert:=" << eEnt.x.EXT.rows()
             << ", ordStab:=" << eEnt.x.GRP.size()
             << ", nAdj:=" << eEnt.ListAdj.size() << ")";
    }
    os_out << "]);\n";
    return;
  }
  std::cerr << "LATT_SerialPeriodicDelaunay: Unsupported OutFormat="
            << OutFormat << "\n";
  std::cerr << "Supported: NumberGAP, SummaryGAP\n";
  throw TerminalException{1};
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_A<T, Tint>(eFull);
  }
  std::cerr << "LATT_SerialPeriodicDelaunay: Failed to find a matching type "
               "for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_COMPUTE_PeriodicDelaunay();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialPeriodicDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_SerialPeriodicDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialPeriodicDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
