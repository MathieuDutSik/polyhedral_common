// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "erdahl_polytope_enumeration.h"
#include "Namelist.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The perfect Delaunay polyhedra of Z^n relative to a space W of degree 2
  functions, up to the affine integral transformations preserving W.
  See the namelist below for the options.
 */

FullNamelist NAMELIST_GetStandard_ENUMERATE_PERFECT() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    // The dimension.
    ListIntValues["n"] = -1;
    // "full": all the functions a + b.x + C[x].
    // "centered": the functions a + C[x - c], whose zero sets are the
    // Delaunay polyhedra symmetric around c, read from FileCenter in the
    // format of ReadVectorFile (n followed by the n entries, rational entries
    // written as p/q).
    ListStringValues["space"] = "full";
    ListStringValues["FileCenter"] = "unset";
    // "recursive": all the perfect Delaunay polyhedra, by the recursive
    // adjacency decomposition through the degenerate polyhedra.
    // "polytopes": only the perfect Delaunay polytopes (or the perfect
    // polyhedra P + L with dim L <= MaxDimL), by flips between them from a
    // starting perfect polytope; complete only if the graph restricted to
    // them is connected.
    ListStringValues["method"] = "recursive";
    // For the method "polytopes": "auto" starts from the Erdahl-Rybnikov
    // polytope P_ER(n) (full space) or its symmetrization from dimension
    // n-1 (centered space), see erdahl_initial_polytope; otherwise a file
    // in the format of ReadMatrixFile with the function of a Delaunay
    // polytope of the space, made perfect by moves if it is not already.
    ListStringValues["FileStart"] = "auto";
    // For the method "polytopes": the perfect polyhedra P + L with
    // dim L <= MaxDimL are kept, the polytopes only for 0.
    ListIntValues["MaxDimL"] = 0;
    // The namelist of the heuristics of the dual descriptions, e.g.
    // CI_tests/16A_EquivDualDesc/CUT_K8/input.nml, or "unset" for the
    // default ones.
    ListStringValues["FileDualDesc"] = "unset";
    // The output, in GAP format: a file, "stderr" or "stdout".
    ListStringValues["OutFile"] = "stderr";
    SingleBlock BlockDATA;
    BlockDATA.setListIntValues(ListIntValues);
    BlockDATA.setListStringValues(ListStringValues);
    ListBlock["DATA"] = BlockDATA;
  }
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  int n = BlockDATA.get_int("n");
  std::string space = BlockDATA.get_string("space");
  std::string FileCenter = BlockDATA.get_string("FileCenter");
  std::string method = BlockDATA.get_string("method");
  std::string FileStart = BlockDATA.get_string("FileStart");
  int MaxDimL = BlockDATA.get_int("MaxDimL");
  std::string FileDualDesc = BlockDATA.get_string("FileDualDesc");
  std::string OutFile = BlockDATA.get_string("OutFile");
  std::ostream &os = std::cerr;
  if (n < 1) {
    std::cerr << "ERDAHL_EnumeratePerfect: n should be at least 1\n";
    throw TerminalException{1};
  }
  ErdahlFunctionSpace<T> W = [&]() -> ErdahlFunctionSpace<T> {
    if (space == "full") {
      return erdahl_full_space<T>(n);
    }
    if (space == "centered") {
      MyVector<T> c = ReadVectorFile<T>(FileCenter);
      if (c.size() != n) {
        std::cerr << "ERDAHL_EnumeratePerfect: the center should be of "
                     "length n\n";
        throw TerminalException{1};
      }
      return erdahl_centered_space(c);
    }
    std::cerr << "ERDAHL_EnumeratePerfect: space should be full or centered, "
                 "not "
              << space << "\n";
    throw TerminalException{1};
  }();
  std::vector<DelaunayPolyhedron<T, Tint>> l_perf;
  std::string summary;
  if (method == "recursive") {
    l_perf = erdahl_enumerate_perfect<T, Tint, Tgroup>(W, FileDualDesc, os);
  } else if (method == "polytopes") {
    std::optional<DelaunayPolyhedron<T, Tint>> opt_D0 =
        [&]() -> std::optional<DelaunayPolyhedron<T, Tint>> {
      if (FileStart == "auto") {
        return erdahl_initial_polytope<T, Tint>(W, os);
      }
      MyMatrix<T> F = ReadMatrixFile<T>(FileStart);
      if (F.rows() != n + 1 || F.cols() != n + 1) {
        std::cerr << "ERDAHL_EnumeratePerfect: the function of FileStart "
                     "should be of size n+1\n";
        throw TerminalException{1};
      }
      if (!erdahl_is_in_space(W, F)) {
        std::cerr << "ERDAHL_EnumeratePerfect: the function of FileStart is "
                     "not in the space\n";
        throw TerminalException{1};
      }
      return erdahl_polyhedron_from_function<T, Tint>(F, os);
    }();
    if (!opt_D0) {
      // No perfect polytope in that space and dimension.
      summary = " no perfect polytope";
    } else {
      DelaunayPolyhedron<T, Tint> Dstart =
          erdahl_perfect_polytope_from<T, Tint>(W, *opt_D0, os);
      ErdahlPerfectPolytopes<T, Tint> res =
          erdahl_enumerate_perfect_polytopes<T, Tint, Tgroup>(
              W, Dstart, MaxDimL, FileDualDesc, os);
      l_perf = res.l_perfect;
      summary = " n_flip=" + std::to_string(res.n_flip) +
                " n_degenerate=" + std::to_string(res.n_degenerate);
    }
  } else {
    std::cerr << "ERDAHL_EnumeratePerfect: method should be recursive or "
                 "polytopes, not "
              << method << "\n";
    throw TerminalException{1};
  }
  auto f = [&](std::ostream &os_out) -> void {
    os_out << "return [";
    for (size_t i = 0; i < l_perf.size(); i++) {
      if (i > 0) {
        os_out << ",\n";
      }
      DelaunayPolyhedron<T, Tint> const &D = l_perf[i];
      os_out << "rec(EXT:=" << StringMatrixGAP(D.EXT)
             << ", L:=" << StringMatrixGAP(D.L)
             << ", F:=" << StringMatrixGAP(D.F) << ")";
    }
    os_out << "];\n";
  };
  if (OutFile == "stderr") {
    f(std::cerr);
  } else {
    if (OutFile == "stdout") {
      f(std::cout);
    } else {
      std::ofstream os_out(OutFile);
      f(os_out);
    }
  }
  std::cerr << "ERDAHL_EnumeratePerfect: n=" << n << " space=" << space
            << " method=" << method << " |l_perf|=" << l_perf.size() << summary
            << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_PERFECT();
    if (argc != 2) {
      std::cerr << "ERDAHL_EnumeratePerfect [file.nml]\n";
      std::cerr << "with file.nml of the form\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    process<T, Tint, Tgroup>(eFull);
    std::cerr << "Normal termination of ERDAHL_EnumeratePerfect\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ERDAHL_EnumeratePerfect\n";
    exit(e.eVal);
  }
  runtime(time);
}
