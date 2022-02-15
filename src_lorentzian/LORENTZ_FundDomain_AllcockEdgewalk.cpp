#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"
#include "Group.h"



#include "edgewalk.h"


int main(int argc, char* argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_EDGEWALK();
    if (argc != 2) {
      std::cerr << "LORENTZ_FundDomain_AllcockEdgeWalk [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      throw TerminalException{1};
    }
    std::string eFileName=argv[1];
    //    using T = long;
    //    using Tint = long;
    using T = mpq_class;
    using Tint = mpz_class;
    //    using Tint = mpq_class;
    //    using T = boost::multiprecision::cpp_rational;
    //    using Tint = boost::multiprecision::cpp_int;
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionEdgewalk<T,Tint,Tgroup>(eFull);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
