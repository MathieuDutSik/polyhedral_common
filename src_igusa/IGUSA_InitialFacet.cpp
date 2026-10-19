// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_initial_facet.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  MyMatrix<T> F = ReadMatrixFile<T>(BlockDATA.get_string("FileFacet"));
  T rhs = ParseScalar<T>(BlockDATA.get_string("FacetRhs"));
  std::vector<MyMatrix<T>> ListRepr =
      ReadListMatrixFile<T>(BlockDATA.get_string("FileListVertex"));
  IgusaInitialFacet<T, Tint> res =
      igusa_initial_facet<T, Tint, Tgroup>(params, F, rhs, ListRepr, os);
  os << "IGUSA_INITIAL_FACET: rank=" << res.rank
     << " |vertices of Phi0|=" << res.face.l_vertex.size()
     << " dim(Phi0)=" << res.face.dim_face << " w=" << StringVector(res.w)
     << " rank C(Phi0)=" << res.rank_C
     << " all points used=" << res.used_all_points
     << " codim of Psi_w in Phi=" << res.codim << "\n";
  auto f_print = [&](std::ostream &os_out) -> void {
    os_out << "return rec(rank:=" << res.rank << ", w:=";
    WriteVectorGAP(os_out, res.w);
    os_out << ", nbVertexPhi0:=" << res.face.l_vertex.size()
           << ", dimPhi0:=" << res.face.dim_face << ", rankC:=" << res.rank_C
           << ", codim:=" << res.codim << ", Certificate:=[";
    for (size_t i = 0; i < res.l_certificate.size(); i++) {
      if (i > 0) {
        os_out << ",";
      }
      WriteVectorGAP(os_out, res.l_certificate[i]);
    }
    os_out << "],\n Fnew:=";
    WriteMatrixGAP(os_out, res.Fnew);
    os_out << ", bnew:=" << res.bnew << ");\n";
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
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_INITIAL_FACET();
    if (argc != 2) {
      std::cerr << "IGUSA_InitialFacet [file.nml]\n";
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
    std::cerr << "Normal termination of IGUSA_InitialFacet\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_InitialFacet\n";
    exit(e.eVal);
  }
  runtime(time);
}
