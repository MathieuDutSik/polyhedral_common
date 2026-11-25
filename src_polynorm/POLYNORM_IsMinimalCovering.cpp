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
void ComputeIsMinimalCovering(std::string const &FileEXT,
                              std::string const &OutFormat,
                              std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT_input = ReadMatrixFile<T>(FileEXT);
  int n = EXT_input.cols() - 1;
  for (int i = 0; i < EXT_input.rows(); i++) {
    if (EXT_input(i, 0) != 1) {
      std::cerr << "The vertices should be in homogeneous coordinates with "
                << "first entry 1, row i=" << i << " has " << EXT_input(i, 0)
                << "\n";
      throw TerminalException{1};
    }
  }
  if (EXT_input.rows() == 0 || RankMat(EXT_input) != n + 1) {
    std::cerr << "The polytope should be full dimensional\n";
    throw TerminalException{1};
  }
  // Repeated rows and rows that are not vertices are dropped.
  MyMatrix<T> EXT = SelectRow(
      EXT_input, get_distinct_extreme_rays_indices(EXT_input, std::cerr));
  MyMatrix<T> FAC = PolyNorm_MinCovFacets(EXT, std::cerr);
  PolyNormIsMinimalCovering<T> res =
      PolyNorm_IsMinimalCovering<T, Tint, Tgroup>(EXT, FAC, std::cerr);
  auto get_vertex = [&](int iv) -> MyVector<T> {
    MyVector<T> v(n);
    for (int k = 0; k < n; k++) {
      v(k) = EXT(iv, k + 1);
    }
    return v;
  };
  if (OutFormat == "text") {
    os_out << "Dimension n=" << n << " with " << EXT.rows() << " vertices and "
           << FAC.rows() << " facets\n";
    os_out << "is_covering=" << res.is_covering << "\n";
    os_out << "is_minimal=" << res.is_minimal << "\n";
    if (res.is_covering) {
      os_out << "Vertices that can be cut keeping the covering ("
             << res.ListCuttable.size() << "):\n";
      for (auto &iv : res.ListCuttable) {
        os_out << "  " << StringVector(get_vertex(iv)) << "\n";
      }
      os_out << "Local tests: " << res.n_lp << " linear programs\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(is_covering:=" << GAP_logical(res.is_covering)
           << ", is_minimal:=" << GAP_logical(res.is_minimal)
           << ", ListCuttable:=[";
    bool IsFirst = true;
    for (auto &iv : res.ListCuttable) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      WriteVectorGAP(os_out, get_vertex(iv));
    }
    os_out << "]);\n";
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
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      ComputeIsMinimalCovering<T>(FileEXT, OutFormat, *os_ptr);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_IsMinimalCovering [arith] [FileEXT] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_IsMinimalCovering [arith] [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "FileEXT   : the vertices of the polytope P, one per row\n";
      std::cerr << "            in homogeneous coordinates (1, v)\n";
      std::cerr << "OutFormat : text (default) or GAP\n";
      std::cerr << "OutFile   : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Decides whether P is minimal covering: P + Z^n = R^n\n";
      std::cerr << "and no convex body strictly contained in P is covering.\n";
      std::cerr << "When P is covering, lists the vertices of P that can be\n";
      std::cerr << "cut keeping the covering\n";
      return -1;
    }
    std::string arith = argv[1];
    FileEXT = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_IsMinimalCovering\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_IsMinimalCovering\n";
    exit(e.eVal);
  }
  runtime(time);
}
