#include "POLY_cddlib.h"
#include "POLY_LinearProgramming.h"

int main(int argc, char *argv[])
{
  try {
    std::vector<std::vector<int> > TheReturn;
    VectVectInt TheOutput;
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
    MyMatrix<mpq_class> EXT=ReadMatrixLrsCdd<mpq_class>(is);
    int rnk=RankMat(EXT);
    int nbCol=EXT.cols();
    if (rnk != nbCol) {
      std::cerr << "The polytope is not full dimensional\n";
      std::cerr << "rnk=" << rnk << " nbCol=" << nbCol << "\n";
      exit(1);
    }
    //
    mpq_class smallVal=0;
    std::vector<Face> ListIncd=cdd::DualDescription_incd(EXT, smallVal);
    //
    int nbFace=ListIncd.size();
    std::ofstream os(argv[2]);
    os << nbFace << " " << nbCol << "\n";
    for (auto & eFace : ListIncd) {
      MyVector<mpq_class> eV=FindFacetInequality(EXT, eFace);
      WriteVector(os, eV);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
