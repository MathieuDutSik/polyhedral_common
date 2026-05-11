// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "GRAPH_GraphicalFunctions.h"
#include "POLY_cdd_graph.h"
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
  cdd::DDA<T> eDDA = cdd::DualDescriptionAdjacencies(EXT, os);
  if (OutFormat == "GAP") {
    os_out << "return rec(FAC:=";
    WriteMatrixGAP(os_out, eDDA.EXT);
    os_out << ",\n  nbVertSkel:=" << eDDA.SkelGraph.GetNbVert();
    os_out << ",\n  SkelEdges:=";
    WriteEdgesGAP(os_out, eDDA.SkelGraph);
    os_out << ",\n  nbVertRidge:=" << eDDA.RidgeGraph.GetNbVert();
    os_out << ",\n  RidgeEdges:=";
    WriteEdgesGAP(os_out, eDDA.RidgeGraph);
    os_out << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    WriteMatrix(os_out, eDDA.EXT);
    WriteEdgesCPP(os_out, eDDA.SkelGraph);
    WriteEdgesCPP(os_out, eDDA.RidgeGraph);
    return;
  }
  std::cerr << "Error in process, missing support for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_cdd_skeletons [arith] [DATAEXT] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "   or\n";
      std::cerr << "POLY_cdd_skeletons [arith] [DATAEXT]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      std::cerr << "safe_rational          : rational arithmetic based on "
                   "int64_t that fails\n";
      std::cerr << "    gracefully on overflowing\n";
      std::cerr << "cpp_rational           : rational arithmetic based on "
                   "boost header library\n";
      std::cerr << "mpq_rational           : rational arithmetic based on "
                   "boost mpq data type\n";
      std::cerr << "mpq_class              : rational arithmetic over GMP "
                   "mpq_class\n";
      std::cerr
          << "Qsqrt2                 : arithmetic over the field Q(sqrt(2))\n";
      std::cerr
          << "Qsqrt5                 : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
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
    std::string FileEXT = argv[2];
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }

    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        using T_rat = mpq_class;
        std::string const &FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(FileEXT, OutFormat, os_out, std::cerr);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    //
    std::cerr << "Normal termination of POLY_cdd_skeletons\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_cdd_skeletons\n";
    exit(e.eVal);
  }
  runtime(time1);
}
