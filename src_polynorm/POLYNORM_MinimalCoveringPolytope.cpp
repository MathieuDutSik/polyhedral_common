// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_MinimalCovering.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeMinimalCovering(std::string const &FileEXT,
                            std::string const &OutFormat, size_t max_step,
                            std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  PolyNormMinimalCovering<T> res =
      PolyNorm_MinimalCoveringPolytope<T, Tint, Tgroup>(EXT, max_step,
                                                        std::cerr);
  if (!res.is_covering) {
    std::cerr << "The polytope is not covering, its covering radius is mu="
              << res.mu_input << " > 1\n";
    throw TerminalException{1};
  }
  if (OutFormat == "text") {
    os_out << "Covering radius of the input mu=" << res.mu_input << "\n";
    os_out << "Reductions after the shrinking by mu: " << res.n_remove
           << " vertices removed, " << res.n_push << " facets pushed, "
           << res.n_cut << " vertices cut\n";
    os_out << "Search: " << res.n_node << " nodes, " << res.n_lp
           << " linear programs\n";
    if (res.is_minimal) {
      os_out << "The polytope is minimal covering\n";
    } else {
      os_out << "The limit of " << max_step << " reductions was reached, the "
             << "polytope is covering but not proved minimal\n";
    }
    os_out << "Vertices (" << res.EXT.rows() << "):\n";
    WriteMatrix(os_out, res.EXT);
    os_out << "Facets b + a.x >= 0 (" << res.FAC.rows() << "):\n";
    WriteMatrix(os_out, res.FAC);
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(mu_input:=" << res.mu_input
           << ", is_minimal:=" << GAP_logical(res.is_minimal)
           << ", n_remove:=" << res.n_remove << ", n_push:=" << res.n_push
           << ", n_cut:=" << res.n_cut << ", EXT:=";
    WriteMatrixGAP(os_out, res.EXT);
    os_out << ", FAC:=";
    WriteMatrixGAP(os_out, res.FAC);
    os_out << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileEXT;
    std::string OutFormat = "text";
    std::string OutFile = "stderr";
    size_t max_step = 1000;
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      ComputeMinimalCovering<T>(FileEXT, OutFormat, max_step, *os_ptr);
    });
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_MinimalCoveringPolytope [arith] [FileEXT] "
                   "[OutFormat] [OutFile] [MaxStep]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_MinimalCoveringPolytope [arith] [FileEXT] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_MinimalCoveringPolytope [arith] [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "FileEXT   : the vertices of the polytope P, one per row\n";
      std::cerr << "            in homogeneous coordinates (1, v)\n";
      std::cerr << "OutFormat : text (default) or GAP\n";
      std::cerr << "OutFile   : the output file, or stderr / stdout\n";
      std::cerr << "MaxStep   : the maximal number of reductions (default 1000)\n";
      std::cerr << "\n";
      std::cerr << "P must be covering, P + Z^n = R^n. Returns a minimal\n";
      std::cerr << "covering polytope Q contained in P: Q + Z^n = R^n and\n";
      std::cerr << "no convex body strictly contained in Q is covering\n";
      return -1;
    }
    std::string arith = argv[1];
    FileEXT = argv[2];
    if (argc >= 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (argc == 6) {
      max_step = ParseScalar<size_t>(argv[5]);
    }
    auto prt = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_MinimalCoveringPolytope\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_MinimalCoveringPolytope\n";
    exit(e.eVal);
  }
  runtime(time);
}
