#include "Temp_cdd.h"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  std::ifstream EXTfs;
  std::ofstream OUTfs;
  vector<vector<int> > TheReturn;
  VectVectInt TheOutput;
  if (argc != 2) {
    fprintf(stderr, "Number of argument is = %d\n", argc);
    fprintf(stderr, "This program is used as\n");
    fprintf(stderr, "DualDescriptionAdjacencies [DATAEXT]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DATAEXT: The input data of the polytope vertices\n");
    return -1;
  }
  //
  fprintf(stderr, "Reading input\n");
  std::string eStr=std::string(argv[1]);
  EXTfs.open(eStr);
  MyMatrix<mpq_class> EXT=ReadMatrix<mpq_class>(EXTfs);
  MyMatrix<mpq_class> EXTred=ColumnReduction<mpq_class>(EXT);
  EXTfs.close();
  //
  mpq_class smallVal;
  smallVal=0;
  cdd::DDA<mpq_class> TheDDA=cdd::DualDescriptionAdjacencies(EXTred, smallVal);
  std::cout << "EXT=\n";
  WriteMatrix(std::cout, TheDDA.EXT);
  std::cout << "SkelGraph=\n";
  GRAPH_PrintOutput(std::cout, TheDDA.SkelGraph);
  std::cout << "RidgeGraph=\n";
  GRAPH_PrintOutput(std::cout, TheDDA.RidgeGraph);
  //
  fprintf(stderr, "Completion of the program\n");
}
