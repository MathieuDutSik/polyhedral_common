// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tidx>
void Process_eFull(FullNamelist const &eFull) {
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint>;
  //    using Tidx_value = int16_t;
  using Tidx_value = int32_t;
  MainFunctionSerialDualDesc<T, Tgroup, Tidx_value>(eFull, std::cerr);
}

template <typename T> void Process(FullNamelist const &eFull) {
  MyMatrix<T> EXT = GetEXT_from_efull<T>(eFull);
  //
  if (size_t(EXT.rows()) < std::numeric_limits<uint8_t>::max())
    return Process_eFull<T, uint8_t>(eFull);
  if (size_t(EXT.rows()) < std::numeric_limits<uint16_t>::max())
    return Process_eFull<T, uint16_t>(eFull);
  if (size_t(EXT.rows()) < std::numeric_limits<uint32_t>::max())
    return Process_eFull<T, uint32_t>(eFull);
#if !defined __APPLE__
  if (size_t(EXT.rows()) < std::numeric_limits<uint64_t>::max())
    return Process_eFull<T, uint64_t>(eFull);
#endif
  std::cerr << "Failed to find a numeric type that matches\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_RecursiveDualDescription();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SerialDualDesc [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    std::string NumericalType = GetNumericalType(eFull);
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void { Process<T>(eFull); };
    add_rational_arithmetics(dispatch, f_arith);
    add_integer_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
#ifdef QUADRATIC_FIELDS
    dispatch.add("Qsqrt3", "arithmetic over the field Q(sqrt(3))",
                 [&]() -> void {
                   f_arith.template operator()<QuadField<Trat_quad_field, 3>>();
                 });
#endif
    dispatch.run(NumericalType);
    //
    std::cerr << "Normal termination of POLY_SerialDualDesc\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SerialDualDesc\n";
    exit(e.eVal);
  } catch (RuntimeException const &e) {
    std::cerr << "The maximum runtime has been reached, exiting "
                 "POLY_SerialDualDesc\n";
    // exit(e.eVal);
  } catch (SafeIntException const &e) {
    std::cerr << "The int64_t arithmetic overflowed: rerun with an "
                 "arbitrary precision NumericalType (mpz_class, mpq_class)\n";
    std::cerr << "Error in POLY_SerialDualDesc\n";
    exit(e.eVal);
  }
  runtime(time1);
}
