// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  IgusaParameters<T> params{BlockDATA.get_int("NormBound"),
                            BlockDATA.get_int("NbRandomFacet"),
                            BlockDATA.get_int("NbCleanRound"),
                            BlockDATA.get_string("IlpMethod"),
                            BlockDATA.get_bool("NormClosure"),
                            BlockDATA.get_int("MaxNormEnumeration"),
                            BlockDATA.get_string("FileInstancePrefix"),
                            BlockDATA.get_bool("StopBeforeDualDescription"),
                            BlockDATA.get_string("FileDualDescription"),
                            BlockDATA.get_int("NbSampleFacet"),
                            BlockDATA.get_int("StopMinRays"),
                            BlockDATA.get_bool("OrbitClosure"),
                            BlockDATA.get_int("MaxEnumerationBound"),
                            BlockDATA.get_bool("LogNewRays")};
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  using Tdata = DataIgusaFunc<T, Tint, Tgroup>;
  Tdata data{std::move(space), os, {}, {}};
  std::string FileInitialVertex = BlockDATA.get_string("FileInitialVertex");
  if (FileInitialVertex != "unset") {
    data.InitialVertex = ReadMatrixFile<T>(FileInitialVertex);
  }
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  std::string OutFile = BlockDATA.get_string("OutFile");
  auto f_output = [&](auto f_print) -> void {
    if (OutFile == "stderr") {
      return f_print(std::cerr);
    }
    if (OutFile == "stdout") {
      return f_print(std::cout);
    }
    std::ofstream os_out(OutFile);
    f_print(os_out);
  };
  // Only the local cone of a given vertex
  std::string FileVertex = BlockDATA.get_string("FileVertex");
  if (FileVertex != "unset") {
    MyMatrix<T> A = ReadMatrixFile<T>(FileVertex);
    igusa_insert_initial_cuts(data.space);
    if (BlockDATA.get_bool("VertexTest")) {
      bool rigorous = BlockDATA.get_bool("VertexTestRigorous");
      bool is_vertex = igusa_vertex_test(data.space, A, 5, rigorous, os);
      os << "IGUSA: VertexTest result=" << GAP_logical(is_vertex) << "\n";
      return;
    }
    Tobj x = data.make_vertex(A);
    std::string FileKnownEdges = BlockDATA.get_string("FileKnownEdges");
    if (FileKnownEdges != "unset") {
      data.MapKnownEdge[A] = ReadListMatrixFile<T>(FileKnownEdges);
    }
    std::optional<std::vector<typename Tdata::TadjI>> opt;
    try {
      opt = data.f_adj(x);
    } catch (IgusaStopException const &) {
      os << "IGUSA: stopped after writing the dual description instance\n";
      return;
    }
    std::vector<typename Tdata::TadjI> l_adj =
        unfold_opt(opt, "The local cone computation");
    os << "IGUSA: |ListAdj|=" << l_adj.size() << "\n";
    if (BlockDATA.get_bool("CompareDualSide")) {
      (void)igusa_dual_side_description(x, os);
    }
    // The vertex, and one neighbor B = A + t D for each orbit of edges
    f_output([&](std::ostream &os_out) -> void {
      os_out << "return rec(Vertex:=";
      WriteEntryGAP(os_out, x);
      os_out << ",\nListNeighbor:=[";
      bool is_first = true;
      for (auto &adj : l_adj) {
        if (!is_first) {
          os_out << ",\n";
        }
        is_first = false;
        WriteMatrixGAP(os_out, adj.Gram);
      }
      os_out << "]);\n";
    });
    return;
  }
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = BlockDATA.get_int("max_runtime_second");
  std::optional<std::vector<Tout>> opt_l_tot;
  try {
    opt_l_tot = EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
        data, f_incorrect, max_runtime_second);
  } catch (IgusaStopException const &) {
    os << "IGUSA: stopped after writing the dual description instance\n";
    return;
  }
  std::vector<Tout> l_tot =
      unfold_opt(opt_l_tot, "EnumerateAndStore_Serial (Igusa vertices)");
  os << "IGUSA: number of orbits of vertices=" << l_tot.size() << "\n";
  //
  std::vector<Tobj> l_vert;
  for (auto &ent : l_tot) {
    l_vert.push_back(ent.x);
  }
  std::vector<IgusaFacetOrbit<T>> l_facet_orbit =
      igusa_facet_orbits(data, l_vert);
  os << "IGUSA: number of orbits of facets=" << l_facet_orbit.size() << "\n";
  //
  std::string OutFormat = BlockDATA.get_string("OutFormat");
  auto f_print = [&](std::ostream &os_out) -> void {
    auto write_list = [&](auto const &l_ent, auto f_write) -> void {
      os_out << "[";
      bool IsFirst = true;
      for (auto &ent : l_ent) {
        if (!IsFirst) {
          os_out << ",\n";
        }
        IsFirst = false;
        f_write(ent);
      }
      os_out << "]";
    };
    if (OutFormat == "GAP") {
      os_out << "return rec(ListVertex:=";
      write_list(l_tot,
                 [&](Tout const &ent) -> void { WriteEntryGAP(os_out, ent); });
      os_out << ",\nListFacetOrbit:=";
      write_list(l_facet_orbit, [&](IgusaFacetOrbit<T> const &orb) -> void {
        WriteEntryGAP(os_out, orb);
      });
      os_out << ");\n";
      return;
    }
    if (OutFormat == "PYTHON") {
      os_out << "{\"ListVertex\":";
      write_list(l_tot, [&](Tout const &ent) -> void {
        WriteEntryPYTHON(os_out, ent);
      });
      os_out << ",\n\"ListFacetOrbit\":";
      write_list(l_facet_orbit, [&](IgusaFacetOrbit<T> const &orb) -> void {
        WriteEntryPYTHON(os_out, orb);
      });
      os_out << "}\n";
      return;
    }
    std::cerr << "IGUSA: Unknown OutFormat=" << OutFormat << "\n";
    std::cerr << "IGUSA: Allowed formats: GAP, PYTHON\n";
    throw TerminalException{1};
  };
  f_output(f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ENUMERATE_IGUSA_TSPACE();
    if (argc != 2) {
      std::cerr << "IGUSA_EnumerateVertices [file.nml]\n";
      std::cerr << "with file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    std::string arithmetic = BlockDATA.get_string("arithmetic");
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    auto f_process = [&]() -> void {
      if (arithmetic == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        using Tgroup = permutalib::Group<Telt, mpz_class>;
        return process<T, Tint, Tgroup>(eFull, std::cerr);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arithmetic == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        using Tgroup = permutalib::Group<Telt, mpz_class>;
        return process<T, Tint, Tgroup>(eFull, std::cerr);
      }
#endif
      std::cerr << "IGUSA: Unknown arithmetic=" << arithmetic << "\n";
      std::cerr << "IGUSA: Allowed arithmetic: gmp";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << ", flint";
#endif
      std::cerr << "\n";
      throw TerminalException{1};
    };
    f_process();
    std::cerr << "Normal termination of IGUSA_EnumerateVertices\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_EnumerateVertices\n";
    exit(e.eVal);
  }
  runtime(time);
}
