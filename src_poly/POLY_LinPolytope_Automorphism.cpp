#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Automorphism [EXTIN] [OutGroup]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      std::cerr << "OutGroup : The automorphism group file\n";
      return -1;
    }
    //
    using T=mpq_class;
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT=ReadMatrix<T>(is);
    int nbCol=EXT.cols();
    int nbRow=EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    TheGroupFormat GRP = LinPolytope_Automorphism(EXT);
    //
    std::ofstream os(argv[2]);
    WriteGroup(os, GRP);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
