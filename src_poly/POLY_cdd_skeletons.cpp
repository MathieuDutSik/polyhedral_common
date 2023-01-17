// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRAPH_GraphicalFunctions.h"
#include "NumberTheory.h"
#include "POLY_cdd_graph.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_cdd_skeletons [DATAEXT] [prefix]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT (in) : The polytope vertices\n";
      std::cerr << "prefix (out): where the skeletons are outputed\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
    int rnk = RankMat(EXT);
    int nbCol = EXT.cols();
    if (rnk != nbCol) {
      std::cerr << "The polytope is not full dimensional\n";
      std::cerr << "rnk=" << rnk << " nbCol=" << nbCol << "\n";
      exit(1);
    }
    //
    cdd::DDA<T> eDDA = cdd::DualDescriptionAdjacencies(EXT);
    std::string prefix_out = argv[2];
    //
    std::ofstream os1(prefix_out + "facet");
    WriteMatrix(os1, eDDA.EXT);
    //
    std::ofstream os2(prefix_out + "_skel_graph");
    GRAPH_PrintOutputGAP(os2, eDDA.SkelGraph);
    //
    std::ofstream os3(prefix_out + "_ridge_graph");
    GRAPH_PrintOutputGAP(os3, eDDA.RidgeGraph);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_cdd_skeletons\n";
    exit(e.eVal);
  }
  runtime(time1);
}
