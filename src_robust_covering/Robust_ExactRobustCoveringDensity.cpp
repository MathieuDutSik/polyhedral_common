// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "enum_robust_covering.h"
#include "pyvista_json.h"
// clang-format on

template <typename T, typename Tint>
void process_B(std::string const &MatFile, std::string const &OutFormat,
               std::string const &OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  //
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  DataLattice<T, Tint, Tgroup> eData =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);

  std::vector<PVoronoi<T, Tint>> l_ppoly = compute_all_p_polytopes(eData);
  T sqr_dist =
      square_robust_covering_radius_from_ppoly(l_ppoly, GramMat, std::cerr);

  T det = DeterminantMat(GramMat);
  int dim = GramMat.rows();
  ResultCov<T> rc = ComputeCoveringDensityFromDimDetCov(dim, det, sqr_dist);
  auto f_print = [&](std::ostream &osf) -> void {
    if (OutFormat == "GAP") {
      osf << "return ";
      osf << to_stringGAP(rc);
      osf << ";\n";
      return;
    }
    if (OutFormat == "GAP_extend") {
      // The covering result together with the P-polytopes. Each PVoronoi record
      // carries its generalized polytope in its "gp" field.
      osf << "return rec(cov:=";
      osf << to_stringGAP(rc);
      osf << ", l_ppoly:=[";
      for (size_t i = 0; i < l_ppoly.size(); i++) {
        if (i > 0) {
          osf << ",\n";
        }
        WriteEntryGAP(osf, l_ppoly[i]);
      }
      osf << "]);\n";
      return;
    }
    if (OutFormat == "PyVista_json") {
      if (dim != 3) {
        std::cerr << "The PyVista_json output is only for dimension 3\n";
        throw TerminalException{1};
      }
      write_pyvista_json(osf, l_ppoly, GramMat, std::cerr);
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: GAP, GAP_extend, PyVista_json\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}

void process_A(std::string const &arithmetic, std::string const &MatFile,
               std::string const &OutFormat, std::string const &OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(MatFile, OutFormat, OutFile);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_B<T, Tint>(MatFile, OutFormat, OutFile);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_B<T, Tint>(MatFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_ExactRobustCoveringDensity [arith] [MatFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_ExactRobustCoveringDensity [arith] [MatFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp, gmp_boost, multi_boost\n";
      std::cerr << "  gmp         : T = mpq_class, Tint = mpz_class\n";
      std::cerr << "  gmp_boost   : T = boost::multiprecision::mpq_rational, "
                   "Tint = boost::multiprecision::mpz_int\n";
      std::cerr << "  multi_boost : T = boost::multiprecision::cpp_rational, "
                   "Tint = boost::multiprecision::cpp_int\n";
      std::cerr << "OutFormat: GAP, GAP_extend, PyVista_json\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string MatFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    process_A(arithmetic, MatFile, OutFormat, OutFile);
    std::cerr << "Normal termination of Robust_ExactRobustCoveringDensity\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_ExactRobustCoveringDensity\n";
    exit(e.eVal);
  }
  runtime(time);
}
