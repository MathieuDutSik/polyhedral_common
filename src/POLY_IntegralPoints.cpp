#include "POLY_PolytopeInt.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TestStrictPositivity [FAC] [EXT] [OUT]\n";
      std::cerr << "\n";
      std::cerr << "FAC: list of facets (IN)\n";
      std::cerr << "EXT: list of vertices (IN)\n";
      std::cerr << "OUT: list of integral points (OUT)\n";
      return -1;
    }
    //
    //  std::cerr << "Reading input\n";
    //
    std::ifstream isFAC(argv[1]);
    MyMatrix<mpq_class> FAC=ReadMatrix<mpq_class>(isFAC);
    //
    std::ifstream isEXT(argv[2]);
    MyMatrix<mpq_class> EXT=ReadMatrix<mpq_class>(isEXT);
    //
    std::vector<MyVector<mpq_class>> ListPoint=GetListIntegralPoint(FAC, EXT);
    //
    int len=ListPoint.size();
    std::ofstream os(argv[3]);
    os << "return [";
    for (int i=0; i<len; i++) {
      if (i>0)
	os << ",\n";
      WriteVectorGAP(os, ListPoint[i]);
    }
    os << "];\n";
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
