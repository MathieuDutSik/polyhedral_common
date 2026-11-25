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

// The result of either method, for the output.
template <typename T> struct MinimalCoveringOutput {
  T mu_input;
  bool is_covering;
  bool is_minimal;
  std::string description;
  MyMatrix<T> EXT;
  MyMatrix<T> FAC;
};

template <typename T>
void ComputeMinimalCovering(std::string const &FileEXT,
                            std::string const &OutFormat,
                            std::string const &method,
                            std::vector<int> const &ListParam,
                            std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  MinimalCoveringOutput<T> out;
  if (method == "grid") {
    int N = ListParam.size() >= 1 ? ListParam[0] : 1;
    int n_round = ListParam.size() >= 2 ? ListParam[1] : 3;
    if (N < 1 || n_round < 1) {
      std::cerr << "The grid method needs N >= 1 and n_round >= 1\n";
      throw TerminalException{1};
    }
    PolyNormMinimalCoveringGrid<T> res =
        PolyNorm_MinimalCoveringPolytopeGrid<T, Tint, Tgroup>(EXT, N, n_round,
                                                              std::cerr);
    std::ostringstream os;
    os << "Grid method: " << res.n_round << " rounds, last grid (1/"
       << res.N_last << ") Z^n with " << res.n_point << " candidate points, "
       << res.n_drop << " points dropped, " << res.n_covering_test
       << " covering tests, " << res.n_node << " nodes and " << res.n_lp
       << " linear programs for the vertex cuts\n";
    if (res.is_tile) {
      os << "The polytope has volume 1, it is a lattice tile\n";
    }
    if (!res.is_minimal) {
      os << "Some vertex can be cut: the polytope is covering but not "
         << "minimal, a finer grid or more rounds are needed\n";
    }
    out = {res.mu_input, res.is_covering, res.is_minimal,
           os.str(),     res.EXT,         res.FAC};
  } else if (method == "reduction") {
    size_t max_step = ListParam.size() >= 1 ? ListParam[0] : 1000;
    PolyNormMinimalCovering<T> res =
        PolyNorm_MinimalCoveringPolytope<T, Tint, Tgroup>(EXT, max_step,
                                                          std::cerr);
    std::ostringstream os;
    os << "Reductions after the shrinking by mu: " << res.n_remove
       << " vertices removed, " << res.n_push << " facets pushed, " << res.n_cut
       << " vertices cut\n";
    os << "Search: " << res.n_node << " nodes, " << res.n_lp
       << " linear programs\n";
    if (!res.is_minimal) {
      os << "The limit of " << max_step << " reductions was reached, the "
         << "polytope is covering but not proved minimal\n";
    }
    out = {res.mu_input, res.is_covering, res.is_minimal,
           os.str(),     res.EXT,         res.FAC};
  } else {
    std::cerr << "Failed to find a matching entry for method=" << method
              << ", allowed: grid and reduction\n";
    throw TerminalException{1};
  }
  if (!out.is_covering) {
    std::cerr << "The polytope is not covering, its covering radius is mu="
              << out.mu_input << " > 1\n";
    throw TerminalException{1};
  }
  T volume = PolyNorm_MinCovVolume(out.EXT, std::cerr);
  if (OutFormat == "text") {
    os_out << "Covering radius of the input mu=" << out.mu_input << "\n";
    os_out << out.description;
    if (out.is_minimal) {
      os_out << "The polytope is minimal covering\n";
    }
    os_out << "Volume=" << volume << "\n";
    os_out << "Vertices (" << out.EXT.rows() << "):\n";
    WriteMatrix(os_out, out.EXT);
    os_out << "Facets b + a.x >= 0 (" << out.FAC.rows() << "):\n";
    WriteMatrix(os_out, out.FAC);
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(mu_input:=" << out.mu_input
           << ", is_minimal:=" << GAP_logical(out.is_minimal)
           << ", volume:=" << volume << ", EXT:=";
    WriteMatrixGAP(os_out, out.EXT);
    os_out << ", FAC:=";
    WriteMatrixGAP(os_out, out.FAC);
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
    std::string method = "grid";
    std::vector<int> ListParam;
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      ComputeMinimalCovering<T>(FileEXT, OutFormat, method, ListParam, *os_ptr);
    });
    if (argc != 3 && argc != 5 && argc < 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_MinimalCoveringPolytope [arith] [FileEXT] "
                   "[OutFormat] [OutFile] [method] [param ...]\n";
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
      std::cerr << "method    : grid (default) or reduction\n";
      std::cerr << "            grid [N] [n_round]: the candidate points in\n";
      std::cerr << "            (1/D) Z^n, then (1/2D) Z^n, (1/3D) Z^n, ...\n";
      std::cerr << "            for at most n_round rounds, with D = N times\n";
      std::cerr << "            the lcm of the denominators of the vertices\n";
      std::cerr << "            of mu(P) P (defaults N=1, n_round=3)\n";
      std::cerr << "            reduction [MaxStep]: vertex removals, facet\n";
      std::cerr << "            pushes and vertex cuts, at most MaxStep of\n";
      std::cerr << "            them (default 1000)\n";
      std::cerr << "\n";
      std::cerr << "P must be covering, P + Z^n = R^n. Returns a covering\n";
      std::cerr << "polytope Q contained in P and says whether it is proved\n";
      std::cerr << "minimal: no convex body strictly contained in Q is\n";
      std::cerr << "covering\n";
      return -1;
    }
    std::string arith = argv[1];
    FileEXT = argv[2];
    if (argc >= 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (argc >= 6) {
      method = argv[5];
    }
    for (int i = 6; i < argc; i++) {
      ListParam.push_back(ParseScalar<int>(argv[i]));
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
