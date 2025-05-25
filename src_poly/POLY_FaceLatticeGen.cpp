// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "Group.h"
#include "POLY_Kskeletton.h"
#include "Permutation.h"
// clang-format on

template <typename Tgroup>
void MainFunctionFaceLattice(FullNamelist const &eFull) {
  SingleBlock const& BlockPROC = eFull.get_block("PROC");
  std::string const& arith = BlockPROC.get_string("Arithmetic");
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  }
  if (arith == "rational") {
    using T = mpq_class;
    return MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  }
  std::optional<std::string> opt_realalgebraic =
      get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  }
  std::cerr << "Failed to find a matching arithmetic for arith=" << arith << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    FullNamelist eFull = NAMELIST_GetStandard_FaceLattice();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_FaceLatticeGen [file.nml]\n";
      std::cerr << "with file.nml a namelist\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string FileNML = argv[1];
    NAMELIST_ReadNamelistFile(FileNML, eFull);
    //
    MainFunctionFaceLattice<Tgroup>(eFull);
    std::cerr << "Normal termination of POLY_FaceLatticeGen\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_FaceLatticeGen\n";
    exit(e.eVal);
  }
  runtime(time1);
}
