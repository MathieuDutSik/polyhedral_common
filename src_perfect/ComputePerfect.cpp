#include "Permlib_specific.h"
#include "Temp_PerfectForm.h"
int main(int argc, char *argv[])
{
  try {
    Eigen::initParallel();
    FullNamelist eFull=NAMELIST_GetStandard_COMPUTE_PERFECT();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ComputePerfect [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    using T=mpq_class;
    using Tint=mpz_class;
    using Tgroup=TheGroupFormat<mpz_class>;
    //
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    TreatPerfectLatticesEntry<T,Tint,Tgroup>(eFull);
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
