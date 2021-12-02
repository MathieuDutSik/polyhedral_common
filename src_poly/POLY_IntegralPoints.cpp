#include "POLY_PolytopeInt.h"
#include "NumberTheory.h"

int main(int argc, char *argv[])
{
  try {
    using T=mpq_class;
    using Tint=int;
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IntegralPoints [FAC] [EXT] [OUT]\n";
      std::cerr << "\n";
      std::cerr << "FAC: list of facets (IN)\n";
      std::cerr << "EXT: list of vertices (IN)\n";
      std::cerr << "OUT: list of integral points (OUT)\n";
      return -1;
    }
    //
    //  std::cerr << "Reading input\n";
    //
    MyMatrix<T> FAC=ReadMatrixFile<T>(argv[1]);
    //
    MyMatrix<T> EXT=ReadMatrixFile<T>(argv[2]);
    //
    std::vector<MyVector<Tint>> ListPoint = GetListIntegralPoint<T,Tint>(FAC, EXT);
    //
    std::ofstream os(argv[3]);
    os << "return [";
    bool IsFirst = true;
    for (auto & ePt : ListPoint) {
      if (!IsFirst)
	os << ",\n";
      IsFirst = false;
      WriteVectorGAP(os, ePt);
    }
    os << "];\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
