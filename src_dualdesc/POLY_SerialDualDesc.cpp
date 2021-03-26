#include "Group.h"
#include "POLY_ThreadDualDescription.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_RecursiveDualDescription();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ThreadedADM [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    using T=mpq_class;
    using Tgroup=permutalib::Group<mpz_class>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionSerialDualDesc<T,Tgroup>(eFull);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
