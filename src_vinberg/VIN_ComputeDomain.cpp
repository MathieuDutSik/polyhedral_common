#include "Permlib_specific.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "vinberg_code.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "VIN_ComputeDomain [FileI] [FileO]\n";
      throw TerminalException{1};
    }
    using T=mpq_class;
    using Tint=mpz_class;

    std::string FileI = argv[1];
    std::ifstream is(FileI);
    //
    MyMatrix<Tint> G = ReadMatrix<Tint>(is);
    MyMatrix<Tint> v0 = ReadVector<Tint>(is);
    VinbergTot<T,Tint> Vtot = GetVinbergAux<T,Tint>(G, v0);
    //
    std::vector<MyVector<Tint>> ListVect = FindRoots(Vtot);
    //
    std::ostream& os = std::cout;
    int nVect = ListVect.size();
    os << "|ListVect|=" << nVect << "\n";
    for (int i=0; i<nVect; i++) {
      WriteVector(os, ListVect[i]);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
