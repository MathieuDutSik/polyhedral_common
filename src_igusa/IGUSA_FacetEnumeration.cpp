// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_face_vertices.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  int n = LinSpa.n;
  if (static_cast<int>(LinSpa.ListMat.size()) != n * (n + 1) / 2) {
    std::cerr << "IGUSA_FACET_ENUM: only the full space of symmetric "
                 "matrices (Classic) is supported\n";
    throw TerminalException{1};
  }
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  igusa_insert_initial_cuts(space);
  std::vector<MyMatrix<T>> ListF =
      ReadListMatrixFile<T>(BlockDATA.get_string("FileInitialFacets"));
  std::ifstream is_rhs(BlockDATA.get_string("FileInitialRhs"));
  std::vector<std::pair<MyMatrix<T>, T>> ListInit;
  for (auto &F : ListF) {
    std::string str;
    is_rhs >> str;
    ListInit.push_back({F, ParseScalar<T>(str)});
  }
  int max_orbit = BlockDATA.get_int("MaxOrbit");
  bool full_rank_only = BlockDATA.get_bool("FullRankOnly");
  std::string checkpoint_file = BlockDATA.get_string("CheckpointFile");
  if (checkpoint_file == "unset") {
    checkpoint_file = "";
  }
  int n_processed_init = BlockDATA.get_int("NumberProcessedInit");
  if (n_processed_init < 0) {
    std::cerr << "NumberProcessedInit should be non-negative\n";
    throw TerminalException{1};
  }
  std::string incidence_dir = BlockDATA.get_string("IncidenceDirectory");
  if (incidence_dir == "unset") {
    incidence_dir = "";
  }
  std::string file_vertex = BlockDATA.get_string("FileVertexRepresentatives");
  IgusaIncidenceFunction<T> incidence_function;
  std::vector<MyMatrix<T>> ListVertexRepr;
  if (file_vertex != "unset") {
    ListVertexRepr = ReadListMatrixFile<T>(file_vertex);
    // the vertices of the face, back from the basis B of W: X0 = B X B^T
    incidence_function = [&](MyMatrix<T> const &F,
                             T const &rhs) -> std::vector<MyMatrix<T>> {
      IgusaFaceVertices<T> info = igusa_face_vertices<T, Tint, Tgroup>(
          params, F, rhs, ListVertexRepr, os);
      if (info.rank != F.rows()) {
        std::cerr << "IGUSA_FACET_ENUM: FileVertexRepresentatives is for "
                     "the full rank facets\n";
        throw TerminalException{1};
      }
      MyMatrix<T> Binv = Inverse(info.B);
      std::vector<MyMatrix<T>> l_X;
      for (auto &X0 : info.l_vertex) {
        l_X.push_back(Binv * X0 * Binv.transpose());
      }
      return l_X;
    };
  }
  int complete_up_to = BlockDATA.get_int("LowerRankCompleteUpTo");
  std::vector<std::pair<MyMatrix<T>, T>> ListLowerRank;
  if (complete_up_to > 0) {
    std::vector<MyMatrix<T>> ListFL =
        ReadListMatrixFile<T>(BlockDATA.get_string("FileLowerRankFacets"));
    std::ifstream is_rhs_l(BlockDATA.get_string("FileLowerRankRhs"));
    for (auto &F : ListFL) {
      std::string str;
      is_rhs_l >> str;
      ListLowerRank.push_back({F, ParseScalar<T>(str)});
    }
  }
  IgusaFacetEnumOptions<T> opt;
  opt.full_rank_only = full_rank_only;
  opt.checkpoint_file = checkpoint_file;
  opt.n_processed_init = n_processed_init;
  opt.ListLowerRank = ListLowerRank;
  opt.complete_up_to = complete_up_to;
  opt.incidence_dir = incidence_dir;
  opt.incidence_function = incidence_function;
  {
    std::istringstream is_ranks(BlockDATA.get_string("LowerRankRidgeRanks"));
    int r;
    while (is_ranks >> r) {
      opt.lower_ranks.insert(r);
    }
  }
  opt.seed_initial_facets = BlockDATA.get_bool("SeedInitialFacets");
  std::string file_vertex_lower =
      BlockDATA.get_string("FileVertexRepresentativesLower");
  std::vector<MyMatrix<T>> ListVertexReprLower;
  if (file_vertex_lower != "unset") {
    ListVertexReprLower = ReadListMatrixFile<T>(file_vertex_lower);
    // the vertices of the face (F0, c) of P_r, in the coordinates of F0
    opt.lower_vertex_function =
        [&](MyMatrix<T> const &F0, T const &c) -> std::vector<MyMatrix<T>> {
      int r = F0.rows();
      std::vector<MyMatrix<T>> ListRepr;
      for (auto &A : ListVertexReprLower) {
        if (A.rows() == r) {
          ListRepr.push_back(A);
        }
      }
      if (ListRepr.empty()) {
        std::cerr << "IGUSA_FACET_ENUM: no vertex representative of rank "
                  << r << " in FileVertexRepresentativesLower\n";
        throw TerminalException{1};
      }
      IgusaFaceVertices<T> info =
          igusa_face_vertices<T, Tint, Tgroup>(params, F0, c, ListRepr, os);
      MyMatrix<T> Binv = Inverse(info.B);
      std::vector<MyMatrix<T>> l_Y;
      for (auto &Y0 : info.l_vertex) {
        l_Y.push_back(Binv * Y0 * Binv.transpose());
      }
      return l_Y;
    };
  }
  std::string file_ridges = BlockDATA.get_string("FileInitialRidges");
  if (file_ridges != "unset") {
    std::vector<MyMatrix<T>> l_mat = ReadListMatrixFile<T>(file_ridges);
    if (l_mat.size() % 5 != 0) {
      std::cerr << "IGUSA_FACET_ENUM: FileInitialRidges should contain G, "
                   "(cG), Q, (cQ), (done) for each ridge\n";
      throw TerminalException{1};
    }
    for (size_t i = 0; i < l_mat.size(); i += 5) {
      // done: 0 to do, 1 done, 2 postponed (processed after all the others)
      T d = l_mat[i + 4](0, 0);
      opt.ListInitRidge.push_back({l_mat[i], l_mat[i + 1](0, 0), l_mat[i + 2],
                                   l_mat[i + 3](0, 0), d == T(1),
                                   d == T(2)});
    }
  }
  IgusaFacetEnumResult<T> l_orbit = igusa_facet_enumeration<T, Tint, Tgroup>(
      space, ListInit, max_orbit, os, opt);
  os << "IGUSA_FACET_ENUM: number of orbits of full rank facets="
     << l_orbit.l_facet.size() << " of bounded ridges in facets of lower rank="
     << l_orbit.l_ridge.size() << "\n";
  std::string OutFile = BlockDATA.get_string("OutFile");
  if (OutFile == "stderr") {
    return WriteFacetOrbitsGAP(std::cerr, l_orbit);
  }
  if (OutFile == "stdout") {
    return WriteFacetOrbitsGAP(std::cout, l_orbit);
  }
  std::ofstream os_out(OutFile);
  WriteFacetOrbitsGAP(os_out, l_orbit);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_FACET_ENUMERATION();
    if (argc != 2) {
      std::cerr << "IGUSA_FacetEnumeration [file.nml]\n";
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
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process<T, Tint, Tgroup>(eFull, std::cerr);
        });
    dispatch.run(arithmetic);
    std::cerr << "Normal termination of IGUSA_FacetEnumeration\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_FacetEnumeration\n";
    exit(e.eVal);
  }
  runtime(time);
}
