// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_face_enum.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T>
std::vector<IgusaPair<T>> read_pairs(std::string const &FileMat,
                                     std::string const &FileRhs) {
  std::vector<IgusaPair<T>> l_pair;
  if (FileMat == "unset") {
    return l_pair;
  }
  std::vector<MyMatrix<T>> ListF = ReadListMatrixFile<T>(FileMat);
  std::ifstream is_rhs(FileRhs);
  for (auto &F : ListF) {
    std::string str;
    is_rhs >> str;
    l_pair.push_back({F, ParseScalar<T>(str)});
  }
  return l_pair;
}

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  std::vector<IgusaPair<T>> l_lower =
      read_pairs<T>(BlockDATA.get_string("FileLowerFacets"),
                    BlockDATA.get_string("FileLowerRhs"));
  std::vector<IgusaPair<T>> l_initial =
      read_pairs<T>(BlockDATA.get_string("FileInitialFacets"),
                    BlockDATA.get_string("FileInitialRhs"));
  if (l_lower.empty()) {
    std::cerr << "IGUSA_FACE_ENUM: the facets of lower rank are needed\n";
    throw TerminalException{1};
  }
  int n = l_lower[0].F.rows();
  IgusaFaceContext<T, Tint, Tgroup> ctx{params, n, {}, {}, os};
  for (auto &A :
       ReadListMatrixFile<T>(BlockDATA.get_string("FileVertexRepr"))) {
    ctx.map_repr[A.rows()].push_back(A);
  }
  int max_facet = BlockDATA.get_int("MaxFacet");
  IgusaFacetEnumeration<T> res =
      igusa_enumerate_facets<T, Tint, Tgroup>(ctx, l_lower, l_initial,
                                              max_facet);
  size_t n_full = 0;
  for (auto &G : res.l_facet) {
    if (G.rank == n) {
      n_full++;
    }
  }
  os << "IGUSA_FACE_ENUM: n=" << n << " orbits of facets=" << res.l_facet.size()
     << " of full rank=" << n_full << " faces built=" << ctx.n_face_built
     << " flips=" << ctx.n_flip << " minimizations=" << ctx.n_minimization
     << "\n";
  auto f_print = [&](std::ostream &os_out) -> void {
    os_out << "return [";
    for (size_t i = 0; i < res.l_facet.size(); i++) {
      IgusaFace<T> const &G = res.l_facet[i];
      if (i > 0) {
        os_out << ",\n";
      }
      std::vector<size_t> l_count(ctx.get_repr(G.rank).size(), 0);
      for (auto &t : G.l_type) {
        l_count[t]++;
      }
      os_out << "rec(rank:=" << G.rank << ", rhs:=" << G.pair.b
             << ", nbVertex:=" << G.l_vert.size() << ", VertexByType:=[";
      for (size_t u = 0; u < l_count.size(); u++) {
        os_out << (u > 0 ? "," : "") << l_count[u];
      }
      os_out << "], nbRidgeOrbit:=" << res.l_n_ridge[i]
             << ", nbFlip:=" << res.l_n_flip[i] << ", F:=";
      WriteMatrixGAP(os_out, G.pair.F);
      os_out << ")";
    }
    os_out << "];\n";
  };
  std::string OutFile = BlockDATA.get_string("OutFile");
  if (OutFile == "stderr") {
    return f_print(std::cerr);
  }
  if (OutFile == "stdout") {
    return f_print(std::cout);
  }
  std::ofstream os_out(OutFile);
  f_print(os_out);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_FACE_ENUMERATION();
    if (argc != 2) {
      std::cerr << "IGUSA_FaceEnumeration [file.nml]\n";
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
      std::cerr << "IGUSA_FACE_ENUM: Failed to find a matching entry for "
                   "arithmetic="
                << arithmetic << "\n";
      throw TerminalException{1};
    };
    f_process();
    std::cerr << "Normal termination of IGUSA_FaceEnumeration\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_FaceEnumeration\n";
    exit(e.eVal);
  }
  runtime(time);
}
