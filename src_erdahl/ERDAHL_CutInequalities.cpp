// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "erdahl_cuts.h"
#include "Namelist.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The facets of the cone CUT_E(n) of the Erdahl cuts phi(phi - 1), see
  erdahl_cuts.h, from the centrally symmetric perfect Delaunay polyhedra
  of Z^{n+1} of center e_1/2. Each perfect function of the full space of
  dimension n is also separated from CUT_E(n) by the cutting plane method,
  which checks the facets independently.
 */

FullNamelist NAMELIST_GetStandard_CUT_INEQUALITIES() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    // The dimension n of the cone CUT_E(n).
    ListIntValues["n"] = -1;
    // The enumeration of the centrally symmetric perfect polyhedra of
    // dimension n+1: "recursive" (complete) or "polytopes" (the perfect
    // polyhedra with dim L <= MaxDimL connected to a starting polytope),
    // see ERDAHL_EnumeratePerfect.
    ListStringValues["method"] = "recursive";
    ListIntValues["MaxDimL"] = 0;
    // The perfect functions of the full space of dimension n separated from
    // CUT_E(n) are those of the method "polytopes" with dim L <= this.
    // -1 for none.
    ListIntValues["TestMaxDimL"] = 1;
    // The heuristics of the dual descriptions, "unset" for the default.
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
  std::string method = BlockDATA.get_string("method");
  int MaxDimL = BlockDATA.get_int("MaxDimL");
  int TestMaxDimL = BlockDATA.get_int("TestMaxDimL");
  std::string FileDualDesc = BlockDATA.get_string("FileDualDesc");
  std::string OutFile = BlockDATA.get_string("OutFile");
  std::ostream &os = std::cerr;
  if (n < 1) {
    std::cerr << "ERDAHL_CutInequalities: n should be at least 1\n";
    throw TerminalException{1};
  }
  // The centrally symmetric perfect polyhedra of Z^{n+1}.
  MyVector<T> c = ZeroVector<T>(n + 1);
  c(0) = T(1) / T(2);
  ErdahlFunctionSpace<T> Wc = erdahl_centered_space(c);
  std::vector<DelaunayPolyhedron<T, Tint>> l_cs;
  if (method == "recursive") {
    l_cs = erdahl_enumerate_perfect<T, Tint, Tgroup>(Wc, FileDualDesc, os);
  } else if (method == "polytopes") {
    std::optional<DelaunayPolyhedron<T, Tint>> opt =
        erdahl_initial_polytope<T, Tint>(Wc, os);
    if (opt) {
      l_cs = erdahl_enumerate_perfect_polytopes<T, Tint, Tgroup>(
                 Wc, *opt, MaxDimL, FileDualDesc, os)
                 .l_perfect;
    }
  } else {
    std::cerr << "ERDAHL_CutInequalities: method should be recursive or "
                 "polytopes\n";
    throw TerminalException{1};
  }
  std::ostringstream out;
  out << "return rec(n:=" << n << ",\n";
  out << "facets:=[";
  size_t n_facet_orbit = 0;
  size_t n_lifted = 0;
  bool first = true;
  for (auto &D : l_cs) {
    if (D.L.rows() > 0) {
      // D = D' x Z^d: its facets are lifted from CUT_E(n - d).
      n_lifted++;
      std::cerr << "ERDAHL_CutInequalities: |EXT|=" << D.EXT.rows()
                << " dim L=" << D.L.rows() << ": lifted from CUT_E("
                << n - D.L.rows() << ")\n";
      continue;
    }
    std::vector<ErdahlCutFacet<T, Tint>> l_facet =
        erdahl_cut_facets_of_polytope<T, Tint, Tgroup>(Wc, D, os);
    for (auto &fac : l_facet) {
      n_facet_orbit++;
      int dim = ((n + 1) * (n + 2)) / 2;
      std::cerr << "ERDAHL_CutInequalities: facet from |EXT|=" << D.EXT.rows()
                << " orbit of " << fac.orbit_size << " antipodal pairs"
                << " n_tight=" << fac.n_tight << " rank_tight="
                << fac.rank_tight << " (facet if " << dim - 1 << ")\n";
      if (!first) {
        out << ",\n";
      }
      first = false;
      out << "rec(M:=" << StringMatrixGAP(fac.M)
          << ", n_vert:=" << D.EXT.rows()
          << ", orbit_size:=" << fac.orbit_size
          << ", n_tight:=" << fac.n_tight
          << ", rank_tight:=" << fac.rank_tight
          << ", EXT:=" << StringMatrixGAP(fac.D.EXT) << ")";
    }
  }
  out << "],\n";
  // The separation of the perfect functions of the full space.
  out << "separations:=[";
  if (TestMaxDimL >= 0) {
    ErdahlFunctionSpace<T> W = erdahl_full_space<T>(n);
    std::vector<MyMatrix<T>> l_F;
    {
      // The strip, which is a cut: the minimum is 0.
      MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
      F(1, 1) = 1;
      F(0, 1) = T(-1) / T(2);
      F(1, 0) = T(-1) / T(2);
      l_F.push_back(F);
    }
    std::optional<DelaunayPolyhedron<T, Tint>> opt =
        erdahl_initial_polytope<T, Tint>(W, os);
    if (opt) {
      for (auto &D : erdahl_enumerate_perfect_polytopes<T, Tint, Tgroup>(
                         W, *opt, TestMaxDimL, FileDualDesc, os)
                         .l_perfect) {
        l_F.push_back(D.F);
      }
    }
    for (size_t i = 0; i < l_F.size(); i++) {
      MyMatrix<T> const &F = l_F[i];
      auto res_F = erdahl_zero_set<T, Tint>(F, os);
      int n_vert = res_F.is_ok() ? res_F.get_ok().EXT.rows() : -1;
      int dimL = res_F.is_ok() ? res_F.get_ok().L.rows() : -1;
      ErdahlCutSeparation<T, Tint> sep = erdahl_cut_separation<T, Tint>(F, os);
      std::cerr << "ERDAHL_CutInequalities: separation of the function of "
                << "|EXT|=" << n_vert << " dim L=" << dimL
                << ": value=" << sep.value << " tight |EXT|="
                << sep.tight.EXT.rows() << " dim L=" << sep.tight.L.rows()
                << " n_iter=" << sep.n_iter
                << " n_constraint=" << sep.n_constraint << "\n";
      if (i > 0) {
        out << ",\n";
      }
      out << "rec(F:=" << StringMatrixGAP(F) << ", n_vert:=" << n_vert
          << ", dimL:=" << dimL << ", value:=" << sep.value
          << ", M:=" << StringMatrixGAP(RemoveFractionMatrix(sep.M))
          << ", tight_EXT:=" << StringMatrixGAP(sep.tight.EXT)
          << ", tight_L:=" << StringMatrixGAP(sep.tight.L) << ")";
    }
  }
  out << "]);\n";
  if (OutFile == "stderr") {
    std::cerr << out.str();
  } else if (OutFile == "stdout") {
    std::cout << out.str();
  } else {
    std::ofstream os_out(OutFile);
    os_out << out.str();
  }
  std::cerr << "ERDAHL_CutInequalities: n=" << n
            << " |cs perfect|=" << l_cs.size() << " lifted=" << n_lifted
            << " facet orbits from polytopes=" << n_facet_orbit << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_CUT_INEQUALITIES();
    if (argc != 2) {
      std::cerr << "ERDAHL_CutInequalities [file.nml]\n";
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
    std::cerr << "Normal termination of ERDAHL_CutInequalities\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ERDAHL_CutInequalities\n";
    exit(e.eVal);
  }
  runtime(time);
}
