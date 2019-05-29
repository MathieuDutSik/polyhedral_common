#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_canonicalize [GramI] [GramO]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "GramO (output): The gram matrix on output\n";
      return -1;
    }
    //    using T=mpq_class;
    using T=long;
    //    using T=mpz_class;
    
    //    using Tint=mpz_class;
    using Tint=long;
    //
    //    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    MyMatrix<T> eMat=ReadMatrix<T>(is);
    MyMatrix<T> eMatCan=ComputeCanonicalForm<T,Tint>(eMat).second;
    //
    std::ofstream os(argv[2]);
    WriteMatrix(os, eMatCan);
    std::cerr << "Normal termination of LATT_canonicalize\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Raised exception led to premature end of LATT_canonicalize\n";
    exit(e.eVal);
  }
}
