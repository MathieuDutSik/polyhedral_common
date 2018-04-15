#include "POLY_ThreadDualDescription.h"
int main(int argc, char *argv[])
{
  try {
    Eigen::initParallel();
    FullNamelist eFull=NAMELIST_GetStandard_TEMP_THREADED_ADM();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ThreadedADM [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionComputeDualDesc<mpq_class>(eFull);
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
