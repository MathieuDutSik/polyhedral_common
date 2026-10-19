// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "enum_robust_covering.h"
#include "pyvista_json.h"
#include "Arithmetic_dispatch_types.h"
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

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string MatFile;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process_B<T, Tint>(MatFile, OutFormat, OutFile);
        });
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_ExactRobustCoveringDensity [arith] [MatFile] "
                << "[OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_ExactRobustCoveringDensity [arith] [MatFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "OutFormat: GAP, GAP_extend, PyVista_json\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    MatFile = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    dispatch.run(arithmetic);
    std::cerr << "Normal termination of Robust_ExactRobustCoveringDensity\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_ExactRobustCoveringDensity\n";
    exit(e.eVal);
  }
  runtime(time);
}
