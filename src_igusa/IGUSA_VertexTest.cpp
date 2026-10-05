// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "igusa_vertex_face.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(FullNamelist const &eFull, std::ostream &os) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa = ReadTspace<T, Tint, Tgroup>(BlockTSPACE, os);
  IgusaParameters<T> params{};
  params.IlpMethod = BlockDATA.get_string("IlpMethod");
  IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
  igusa_insert_initial_cuts(space);
  std::vector<MyMatrix<T>> ListA =
      ReadListMatrixFile<T>(BlockDATA.get_string("FileListForm"));
  std::string OutFile = BlockDATA.get_string("OutFile");
  std::ofstream os_out(OutFile);
  os_out << "return [";
  for (size_t i = 0; i < ListA.size(); i++) {
    MicrosecondTime time;
    IgusaVertexTestResult<T> res = igusa_vertex_test_face(space, ListA[i], os);
    os << "IGUSA_VERTEX_TEST: form " << i << " dim face=" << res.dim_face
       << " rank of the norm 1 vectors=" << res.rank_min
       << " iterations=" << res.n_iter
       << " vertex=" << GAP_logical(res.is_vertex) << " time=" << time
       << "\n";
    if (i > 0) {
      os_out << ",\n";
    }
    os_out << "rec(IsVertex:=" << GAP_logical(res.is_vertex)
           << ", DimFace:=" << res.dim_face << ", nIter:=" << res.n_iter;
    if (res.is_vertex) {
      os_out << ", Functional:=" << StringVectorGAP(res.h);
    } else {
      os_out << ", ListWeight:=[";
      for (size_t a = 0; a < res.l_weight.size(); a++) {
        if (a > 0) {
          os_out << ",";
        }
        os_out << res.l_weight[a];
      }
      os_out << "], ListPoint:=[";
      for (size_t a = 0; a < res.l_point.size(); a++) {
        if (a > 0) {
          os_out << ",";
        }
        WriteMatrixGAP(os_out, res.l_point[a]);
      }
      os_out << "], Remainder:=";
      WriteMatrixGAP(os_out, res.remainder);
    }
    os_out << ")";
  }
  os_out << "];\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_IGUSA_VERTEX_TEST();
    if (argc != 2) {
      std::cerr << "IGUSA_VertexTest [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    NAMELIST_ReadNamelistFile(argv[1], eFull);
    std::string arithmetic = eFull.get_block("DATA").get_string("arithmetic");
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
      std::cerr << "IGUSA_VERTEX_TEST: Failed to find a matching entry for "
                   "arithmetic="
                << arithmetic << "\n";
      throw TerminalException{1};
    };
    f_process();
    std::cerr << "Normal termination of IGUSA_VertexTest\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in IGUSA_VertexTest\n";
    exit(e.eVal);
  }
  runtime(time);
}
