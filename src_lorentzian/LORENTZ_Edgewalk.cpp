#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"
#include "edgewalk.h"


int main(int argc, char* argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_EDGEWALK();
    if (argc != 2) {
      std::cerr << "LORENTZ_EdgeWalk [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      throw TerminalException{1};
    }
    std::string eFileName=argv[1];
    using T=mpq_class;
    using Tint=mpz_class;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionEdgewalk(eFull);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
