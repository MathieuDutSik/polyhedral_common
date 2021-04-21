#include "POLY_cddlib.h"
#include "NumberTheory.h"
#include "POLY_PolytopeFct.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_lcdd [DATAEXT] [DATAFAC]\n";
      std::cerr << "\n";
      std::cerr << "DATAEXT (in) : The polytope vertices\n";
      std::cerr << "DATAEXT (out): The polytope facets\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    using T=mpq_class;
    MyMatrix<T> EXT=ReadMatrixLrsCdd<T>(is);
    int rnk=RankMat(EXT);
    int nbCol=EXT.cols();
    if (rnk != nbCol) {
      std::cerr << "The polytope is not full dimensional\n";
      std::cerr << "rnk=" << rnk << " nbCol=" << nbCol << "\n";
      exit(1);
    }
    //
    vectface ListIncd=cdd::DualDescription_incd(EXT);
    //
    int nbFace=ListIncd.size();
    std::ofstream os(argv[2]);
    os << nbFace << " " << nbCol << "\n";
    for (auto & eFace : ListIncd) {
      MyVector<T> eV=FindFacetInequality(EXT, eFace);
      WriteVector(os, eV);
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
