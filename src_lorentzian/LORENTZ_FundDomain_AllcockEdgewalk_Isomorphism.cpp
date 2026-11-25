// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "edgewalk.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint>
void process(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in mpz_class, whichever the arithmetic of
  // the coordinates.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MainFunctionEdgewalk_Isomorphism<T, Tint, Tgroup>(eFull, std::cerr);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_EDGEWALK_Isomorphism();
    if (argc != 2) {
      std::cerr << "LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism [FileNML]\n";
      std::cerr << "with fileNML a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      throw TerminalException{1};
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    std::string arith = eFull.get_block("PROC").get_string("Arithmetic");
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process<T, Tint>(eFull);
        });
    dispatch.run(arith);
    std::cerr << "Normal termination of "
                 "LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
