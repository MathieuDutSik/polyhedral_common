#include "POLY_c_cddlib.h"
#include "NumberTheory.h"
#include "POLY_PolytopeFct.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_c_lcdd [DATAEXT]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT (in) : The polytope vertices\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    using T=mpq_class;
    MyMatrix<T> EXT=ReadMatrixLrsCdd<T>(is);
    //
    std::vector<Face> ListIncd=cbased_cdd::DualDescription_incd(EXT);
    //
    int nbFace=ListIncd.size();
    std::ofstream os(argv[2]);
    std::cerr << "nbFace=" << nbFace << "\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
