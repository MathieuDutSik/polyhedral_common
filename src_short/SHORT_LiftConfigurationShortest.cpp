#include "SHORT_ShortestConfig.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LiftConfigurationShortest [NPROC] [DATAIN] [TheMethod] [DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "NPROC     : The number of processor used\n";
      std::cerr << "DATAIN    : The list of configuration of shortest vector as input\n";
      std::cerr << "TheMethod : can be cdd or glpk_secure\n";
      std::cerr << "DATAOUT   : The list of configuration of shortest vector as output\n";
      return -1;
    }
    //
    using T=mpq_class;
    using Tint=int;
    using Tgroup=TheGroupFormat<mpz_class>;
    int NPROC;
    sscanf(argv[1], "%d", &NPROC);
    //
    std::string eFileIN(argv[2]);
    std::string TheMethod(argv[3]);
    std::vector<MyMatrix<Tint>> ListConfIn=ReadListConfigurationShortestVector<Tint>(eFileIN);
    //
    std::vector<MyMatrix<Tint>> ListConfOut=SHORT_SimplicialEnumeration<T,Tint,Tgroup>(ListConfIn, NPROC, TheMethod);
    //
    std::string eFileOUT(argv[4]);
    WriteListConfigurationShortestVector(eFileOUT, ListConfOut);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
