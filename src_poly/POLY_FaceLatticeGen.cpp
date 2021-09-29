#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
#include "POLY_Kskeletton.h"
int main(int argc, char *argv[])
{
  try {
    using T=mpq_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    FullNamelist eFull = NAMELIST_GetStandard_FaceLattice();
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_FaceLatticeGen [file.nml]\n";
      std::cerr << "with file.nml a namelist\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string FileNML=argv[1];
    NAMELIST_ReadNamelistFile(FileNML, eFull);
    //
    MainFunctionFaceLattice<T,Tgroup>(eFull);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
