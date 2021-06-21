#include "Permutation.h"
#include "Group.h"
#include "POLY_RecursiveDualDesc.h"
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
    using T = mpq_class;
    //    using Tidx = uint8_t;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    using Tidx_value = int16_t;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionSerialDualDesc<T,Tgroup,Tidx_value>(eFull);
    std::cerr << "Normal termination of the program\n";
  }
  catch (netCDF::exceptions::NcInvalidCoords & e) {
    std::cerr << "e complaint=" << e.what() << "\n";
    exit(1);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
