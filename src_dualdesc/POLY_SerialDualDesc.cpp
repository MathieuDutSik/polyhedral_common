// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryRealField.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
#include "QuadField.h"
// clang-format on

template <typename T, typename Tidx>
void Process_eFull(FullNamelist const &eFull) {
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint>;
  //    using Tidx_value = int16_t;
  using Tidx_value = int32_t;
  MainFunctionSerialDualDesc<T, Tgroup, Tidx_value>(eFull);
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
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_RecursiveDualDescription();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SerialDualDesc [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    std::string NumericalType = GetNumericalType(eFull);
    /*
    if (NumericalType == "integer") {
      using T = mpz_class;
      Process<T>(eFull);
    }
    */
    if (NumericalType == "rational") {
      using T = mpq_class;
      //    using T = boost::multiprecision::cpp_rational;
      //    using T = boost::multiprecision::mpq_rational;
      Process<T>(eFull);
    }
    if (NumericalType == "Qsqrt5") {
      using Trat = mpq_class;
      using T = QuadField<Trat, 5>;
      Process<T>(eFull);
    }
    if (NumericalType == "Qsqrt2") {
      using Trat = mpq_class;
      using T = QuadField<Trat, 2>;
      Process<T>(eFull);
    }
    if (NumericalType == "RealAlgebraic") {
      using T_rat = mpq_class;
      SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
      std::string FileAlgebraicField =
          BlockDATA.ListStringValues.at("FileAlgebraicField");
      if (!IsExistingFile(FileAlgebraicField)) {
        std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                  << " is missing\n";
        throw TerminalException{1};
      }
      HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
      int const idx_real_algebraic_field = 1;
      insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
      using T = RealField<idx_real_algebraic_field>;
      Process<T>(eFull);
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  } catch (RuntimeException const &e) {
    std::cerr << "The maximum runtime has been reached, exiting the program\n";
    exit(e.eVal);
  }
  runtime(time1);
}
