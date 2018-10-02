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
    using T=mpq_class;
    using Tint=mpz_class;
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    MyMatrix<T> eMat=ReadMatrix<T>(is);
    MyMatrix<T> eMatCan=ComputeCanonicalForm<T,Tint>(eMat);
    //
    std::ofstream os(argv[2]);
    WriteMatrix(os, eMatCan);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
