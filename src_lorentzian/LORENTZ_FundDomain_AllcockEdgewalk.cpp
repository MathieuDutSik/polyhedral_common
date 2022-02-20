#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"
#include "Group.h"



#include "edgewalk.h"


int main(int argc, char* argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
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
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
