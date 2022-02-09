//#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryCommon.h"
#include "Permutation.h"
#include "Group.h"
#include "POLY_RecursiveDualDesc.h"


template<typename T, typename Tidx>
void Process_eFull(FullNamelist const& eFull)
{
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint = mpz_class;
  using Tgroup = permutalib::Group<Telt,Tint>;
  //    using Tidx_value = int16_t;
  using Tidx_value = int32_t;
  MainFunctionSerialDualDesc<T,Tgroup,Tidx_value>(eFull);
}


int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull = NAMELIST_GetStandard_RecursiveDualDescription();
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
    //    using T = mpq_class;
    using T = boost::multiprecision::cpp_rational;
    MyMatrix<T> EXT = GetEXT_from_efull<T>(eFull);
    //
    auto process=[&]() -> void {
      if (size_t(EXT.rows()) < std::numeric_limits<uint8_t>::max())
        return Process_eFull<T,uint8_t>(eFull);
      if (size_t(EXT.rows()) < std::numeric_limits<uint16_t>::max())
        return Process_eFull<T,uint16_t>(eFull);
      if (size_t(EXT.rows()) < std::numeric_limits<uint32_t>::max())
        return Process_eFull<T,uint32_t>(eFull);
#if !defined __APPLE__
      if (size_t(EXT.rows()) < std::numeric_limits<uint64_t>::max())
        return Process_eFull<T,uint64_t>(eFull);
#endif
      std::cerr << "Failed to find a numeric type that matches\n";
      throw TerminalException{1};
    };
    process();
    //
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
