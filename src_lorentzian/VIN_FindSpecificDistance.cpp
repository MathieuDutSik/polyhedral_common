#include "Permlib_specific.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "vinberg_code.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "VIN_FindSpecificDistance [FileIn]\n";
      throw TerminalException{1};
    }
    using T=mpz_class;
    using Tint=mpz_class;

    std::string FileIn  = argv[1];
    std::ifstream is(FileIn);
    //
    MyMatrix<T> G = ReadMatrix<T>(is);
    MyMatrix<T> v0 = ReadVector<T>(is);
    VinbergInput<T,Tint> Vin{G, v0};
    VinbergTot<T,Tint> Vtot = GetVinbergAux(Vin);
    //
    MyVector<Tint> a = ReadVector<Tint>(is);
    T n;
    is >> n;
    std::vector<MyVector<Tint>> ListVect = Roots_decomposed_into(Vtot, a, n);
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
