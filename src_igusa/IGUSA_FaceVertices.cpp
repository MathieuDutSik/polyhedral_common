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
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  MyMatrix<T> F = ReadMatrixFile<T>(BlockDATA.get_string("FileFacet"));
  T rhs = ParseScalar<T>(BlockDATA.get_string("FacetRhs"));
  std::vector<MyMatrix<T>> ListRepr =
      ReadListMatrixFile<T>(BlockDATA.get_string("FileListVertex"));
  IgusaFaceVertices<T> info =
      igusa_face_vertices<T, Tint, Tgroup>(params, F, rhs, ListRepr, os);
  std::vector<size_t> l_count(ListRepr.size(), 0);
  for (auto &t : info.l_type) {
    l_count[t]++;
  }
  os << "IGUSA_FACE_VERTICES: rank=" << info.rank << " rhs=" << rhs
     << " slack=" << info.slack << " |support|=" << info.n_support
     << " n_node=" << info.n_node << " |vertices|=" << info.l_vertex.size()
     << " dim(Phi0)=" << info.dim_face
     << " by representative:";
  for (auto &c : l_count) {
    os << " " << c;
  }
  os << "\n";
  std::string str_compare = "unset";
  if (BlockDATA.get_bool("CompareIncidence")) {
    if (info.rank != F.rows()) {
      std::cerr << "IGUSA_FACE_VERTICES: CompareIncidence needs F positive "
                   "definite\n";
      throw TerminalException{1};
    }
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(info.rank);
    IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
    igusa_insert_initial_cuts(space);
    IgusaFacetIncidence<T> inc = igusa_facet_incidence(space, info.F0, rhs, os);
    std::set<MyMatrix<T>> set_inc(inc.l_incident.begin(), inc.l_incident.end());
    size_t n_missing = 0;
    for (auto &X : info.l_vertex) {
      if (set_inc.count(X) == 0) {
        n_missing++;
      }
    }
    os << "IGUSA_FACE_VERTICES: |Inc| by IGUSA_FacetIncidence="
       << set_inc.size() << " vertices not in Inc=" << n_missing
       << " incident forms not found=" << set_inc.size() + n_missing -
                                                info.l_vertex.size()
       << "\n";
    if (n_missing > 0) {
      std::cerr << "IGUSA_FACE_VERTICES: a vertex found is not incident\n";
      throw TerminalException{1};
    }
    std::ostringstream os_cmp;
    os_cmp << set_inc.size();
    str_compare = os_cmp.str();
  }
  auto f_print = [&](std::ostream &os_out) -> void {
    os_out << "return rec(rank:=" << info.rank << ", B:=";
    WriteMatrixGAP(os_out, info.B);
    os_out << ", F0:=";
    WriteMatrixGAP(os_out, info.F0);
    os_out << ", rhs:=" << rhs << ", RyshkovMin:=" << info.ryshkov_min
           << ", slack:=" << info.slack << ", dimFace:=" << info.dim_face
           << ", nbNode:=" << info.n_node
           << ", nbIncidence:=\"" << str_compare << "\", ListCount:=[";
    for (size_t i = 0; i < l_count.size(); i++) {
      os_out << (i > 0 ? "," : "") << l_count[i];
    }
    os_out << "],\n ListType:=[";
    for (size_t i = 0; i < info.l_type.size(); i++) {
      os_out << (i > 0 ? "," : "") << (info.l_type[i] + 1);
    }
    os_out << "],\n ListVertex:=[";
    for (size_t i = 0; i < info.l_vertex.size(); i++) {
      if (i > 0) {
        os_out << ",\n";
      }
      WriteMatrixGAP(os_out, info.l_vertex[i]);
    }
    os_out << "]);\n";
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
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_FACE_VERTICES();
    if (argc != 2) {
      std::cerr << "IGUSA_FaceVertices [file.nml]\n";
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
    add_rational_integer_arithmetics_by_library<false>(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process<T, Tint, Tgroup>(eFull, std::cerr);
        });
    dispatch.run(arithmetic);
    std::cerr << "Normal termination of IGUSA_FaceVertices\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_FaceVertices\n";
    exit(e.eVal);
  }
  runtime(time);
}
