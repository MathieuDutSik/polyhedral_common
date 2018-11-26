#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_cannspanset [GramI] [CannSpanSetO]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "CannSpanSetO (output): The canonical spanning vector set\n";
      return -1;
    }
    using T=mpq_class;
    using Tint=mpz_class;
    //
    std::cerr << "Reading input\n";
    std::ifstream is(argv[1]);
    MyMatrix<T> eMat=ReadMatrix<T>(is);
    MyMatrix<Tint> SHVcan=ComputeCanonicalSpanningSet<T,Tint>(eMat);
    //
    std::ofstream os(argv[2]);
    WriteMatrix(os, SHVcan);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
