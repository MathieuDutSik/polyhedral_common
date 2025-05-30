// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
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
    auto process = [&]() -> void {
      /*
        if (NumericalType == "integer") {
        using T = mpz_class;
        Process<T>(eFull);
        }
      */
      if (NumericalType == "safe_rational") {
        using T = Rational<SafeInt64>;
        return Process<T>(eFull);
      }
      if (NumericalType == "rational") {
        using T = mpq_class;
        return Process<T>(eFull);
      }
      if (NumericalType == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return Process<T>(eFull);
      }
      if (NumericalType == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return Process<T>(eFull);
      }
      if (NumericalType == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return Process<T>(eFull);
      }
      if (NumericalType == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return Process<T>(eFull);
      }
      if (NumericalType == "RealAlgebraic") {
        using T_rat = mpq_class;
        SingleBlock const& BlockDATA = eFull.get_block("DATA");
        std::string FileAlgebraicField = BlockDATA.get_string("FileAlgebraicField");
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return Process<T>(eFull);
      }
      std::cerr << "Failed to find a matching type entry\n";
      throw TerminalException{1};
    };
    process();
    //
    std::cerr << "Normal termination of POLY_SerialDualDesc\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SerialDualDesc\n";
    exit(e.eVal);
  } catch (RuntimeException const &e) {
    std::cerr << "The maximum runtime has been reached, exiting POLY_SerialDualDesc\n";
    // exit(e.eVal);
  }
  runtime(time1);
}
