#include "NumberTheory.h"
#include "POLY_PolytopeFct.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_FaceLatticeDirect [eFileI] [eFileO]\n";
      std::cerr << "\n";
      std::cerr << "eFileI: the input file\n";
      std::cerr << "eFileO: the output file\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T=mpq_class;
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT=ReadMatrix<T>(is);
    MyMatrix<T> FAC=ReadMatrix<T>(is);
    std::cerr << "After read matrix\n";
    //
    std::string eFileO=argv[2];
    ComputeFileFaceLatticeInfo(eFileO, EXT, FAC);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
