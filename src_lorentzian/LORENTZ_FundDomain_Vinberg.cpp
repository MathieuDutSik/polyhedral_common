#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "vinberg.h"


int main(int argc, char* argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_VINBERG();
    if (argc != 2) {
      std::cerr << "LORENTZ_FundDomain_Vinberg [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      throw TerminalException{1};
    }
    std::string eFileName=argv[1];
    using T=mpq_class;
    using Tint=mpz_class;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionVinberg<T,Tint>(eFull);
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
