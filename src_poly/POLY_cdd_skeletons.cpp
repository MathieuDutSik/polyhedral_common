// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "GRAPH_GraphicalFunctions.h"
#include "POLY_DualDesc_double_description.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

namespace {

template <typename Tgr>
void WriteEdgesGAP(std::ostream &os, Tgr const &gr) {
  size_t nbVert = gr.GetNbVert();
  os << "[";
  bool first = true;
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::vector<size_t> Adj = gr.Adjacency(iVert);
    for (auto &jVert : Adj) {
      if (jVert > iVert) {
        if (!first)
          os << ",";
        first = false;
        os << "[" << (iVert + 1) << "," << (jVert + 1) << "]";
      }
    }
  }
  os << "]";
}

template <typename Tgr>
void WriteEdgesCPP(std::ostream &os, Tgr const &gr) {
  size_t nbVert = gr.GetNbVert();
  std::vector<std::pair<size_t, size_t>> edges;
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::vector<size_t> Adj = gr.Adjacency(iVert);
    for (auto &jVert : Adj)
      if (jVert > iVert)
        edges.emplace_back(iVert, jVert);
  }
  os << nbVert << " " << edges.size() << "\n";
  for (auto &e : edges)
    os << e.first << " " << e.second << "\n";
}

} // namespace

template <typename T>
void process(std::string const &eFile, std::string const &OutFormat,
             std::ostream &os_out, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
  int rnk = RankMat(EXT);
  int nbCol = EXT.cols();
  if (rnk != nbCol) {
    std::cerr << "The polytope is not full dimensional\n";
    std::cerr << "rnk=" << rnk << " nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
  double_desc::DDA<T> eDDA =
      double_desc::DualDescriptionAdjacencies(EXT, true, true, os);
  if (OutFormat == "GAP") {
    os_out << "return rec(FAC:=";
    WriteMatrixGAP(os_out, eDDA.EXT);
    os_out << ",\n  nbVertSkel:=" << eDDA.SkelGraph->GetNbVert();
    os_out << ",\n  SkelEdges:=";
    WriteEdgesGAP(os_out, *eDDA.SkelGraph);
    os_out << ",\n  nbVertRidge:=" << eDDA.RidgeGraph->GetNbVert();
    os_out << ",\n  RidgeEdges:=";
    WriteEdgesGAP(os_out, *eDDA.RidgeGraph);
    os_out << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    WriteMatrix(os_out, eDDA.EXT);
    WriteEdgesCPP(os_out, *eDDA.SkelGraph);
    WriteEdgesCPP(os_out, *eDDA.RidgeGraph);
    return;
  }
  std::cerr << "Error in process, missing support for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileEXT;
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileEXT, OutFormat, *os_ptr, std::cerr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_cdd_skeletons [arith] [DATAEXT] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "   or\n";
      std::cerr << "POLY_cdd_skeletons [arith] [DATAEXT]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "\n";
      std::cerr << "     ------- DATAEXT -------\n";
      std::cerr << "\n";
      std::cerr << "The polytope vertices in lrs/cdd V-representation "
                   "format\n";
      std::cerr << "\n";
      std::cerr << "     ------- OutFormat -------\n";
      std::cerr << "\n";
      std::cerr << "GAP : The GAP record format\n";
      std::cerr << "CPP : The CPP polyhedral format\n";
      std::cerr << "\n";
      std::cerr << "     ------- OutFile -------\n";
      std::cerr << "\n";
      std::cerr << "stderr : output to std::cerr\n";
      std::cerr << "stdout : output to std::cout\n";
      std::cerr << "filename : output to filename (if different from stderr / "
                   "stdout)\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::string arith = argv[1];
    FileEXT = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }

    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of POLY_cdd_skeletons\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_cdd_skeletons\n";
    exit(e.eVal);
  }
  runtime(time1);
}
