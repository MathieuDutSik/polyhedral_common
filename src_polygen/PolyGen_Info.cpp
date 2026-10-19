// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "generalized_polytopes.h"
#include "Arithmetic_dispatch.h"
// clang-format on

// The structural information on a generalized polytope given as a list of
// vertex sets: whether the pieces have disjoint interiors, the connected
// components (pieces sharing a facet are adjacent), the number of pieces
// after merging the convex unions, and the boundary hyperplanes.
template <typename T>
void process(std::string const &GenPolyFile, std::string const &OutFormat,
             std::string const &OutFile) {
  std::vector<MyMatrix<T>> list_ext = ReadListMatrixFile<T>(GenPolyFile);
  int dim = list_ext[0].cols();
  GeneralizedPolytope<T> gp = list_ext_to_generalizedpolytope(dim, list_ext);
  size_t n_polytope = gp.size();
  T volume = volume_gp(gp, std::cerr);
  bool is_disjoint = check_pairwise_intersection(gp, std::cerr);
  std::vector<GeneralizedPolytope<T>> l_conn =
      connected_components_decomposition(gp, std::cerr);
  std::vector<size_t> l_conn_size;
  for (auto &gp_conn : l_conn) {
    l_conn_size.push_back(gp_conn.size());
  }
  std::sort(l_conn_size.begin(), l_conn_size.end());
  BoundaryGeneralizedPolytope<T> bnd =
      find_generalized_polytope_boundary(gp, std::cerr);
  size_t n_boundary_hyperplane = bnd.size();
  GeneralizedPolytope<T> gp_simp = gp;
  size_t n_simplified = simplify_generalized_polytopes(gp_simp, std::cerr);
  T volume_simplified = volume_gp(gp_simp, std::cerr);

  auto f_print = [&](std::ostream &os_out) -> void {
    if (OutFormat == "GAP") {
      os_out << "return rec(n_polytope:=" << n_polytope
             << ", volume:=" << volume << ", is_disjoint:="
             << GAP_logical(is_disjoint) << ", n_component:=" << l_conn.size()
             << ", component_sizes:=" << StdVectorToStringGAP(l_conn_size)
             << ", n_boundary_hyperplane:=" << n_boundary_hyperplane
             << ", n_simplified:=" << n_simplified
             << ", volume_simplified:=" << volume_simplified << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string PolyFile;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    ArithmeticDispatch dispatch;
    dispatch.add("mpq_class", [&]() -> void {
      process<mpq_class>(PolyFile, OutFormat, OutFile);
    });
#ifdef ENABLE_BOOST_TYPES
    dispatch.add("mpq_rational", [&]() -> void {
      process<boost::multiprecision::mpq_rational>(PolyFile, OutFormat, OutFile);
    });
    dispatch.add("cpp_rational", [&]() -> void {
      process<boost::multiprecision::cpp_rational>(PolyFile, OutFormat, OutFile);
    });
#endif
#ifdef ENABLE_FLINT_SUPPORT
    dispatch.add("fmpq_class", [&]() -> void {
      process<fmpq_class>(PolyFile, OutFormat, OutFile);
    });
#endif
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "PolyGen_Info arith [ListEXT] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "PolyGen_Info arith [ListEXT]\n";
      std::cerr << "\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: " << dispatch.choices() << "\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arith = argv[1];
    PolyFile = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    dispatch.run(arith);
    std::cerr << "Normal termination PolyGen_Info\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PolyGen_Info\n";
    exit(e.eVal);
  }
  runtime(time);
}
