// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Group.h"
#include "Permutation.h"
#include "edgewalk.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_EDGEWALK();
    if (argc != 2) {
      std::cerr << "LORENTZ_FundDomain_AllcockEdgeWalk [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      throw TerminalException{1};
    }
    std::string eFileName = argv[1];
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
    using Tgroup = permutalib::Group<Telt, Tint>;
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    MainFunctionEdgewalk<T, Tint, Tgroup>(eFull, std::cerr);
    std::cerr << "Normal termination of LORENTZ_FundDomain_AllcockEdgewalk\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_FundDomain_AllcockEdgewalk\n";
    exit(e.eVal);
  }
  runtime(time);
}
