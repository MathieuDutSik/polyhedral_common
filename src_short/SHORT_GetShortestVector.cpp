#include "Permlib_specific.h"
#include "SHORT_ShortestConfig.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_GetShortestVector [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]   : The input file of the system\n";
      std::cerr << "[FileOut]  : The output file of the program (GAP readable)\n";
      return -1;
    }
    //
    using T=mpq_class;
    //    using Tint=mpz_class;
    using Tint=int;
    //    using Tgroup=TheGroupFormat<mpz_class>;
    std::ifstream is(argv[1]);
    MyMatrix<T> M = ReadMatrix<T>(is);
    //
    Tshortest<T,Tint> RecSHV=T_ShortestVector<T,Tint>(M);
    int nbRow=RecSHV.SHV.rows();
    std::cerr << "nbRow=" << nbRow << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
