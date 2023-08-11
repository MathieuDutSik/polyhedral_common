// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "poincare_polyhedron.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on


template <typename Tgroup>
void Process_rec_option(RecOption const &rec_option) {
  std::string arith = rec_option.Arithmetic;
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T,Tgroup>(rec_option);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_type<T,Tgroup>(rec_option);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_type<T,Tgroup>(rec_option);
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
    return full_process_type<T,Tgroup>(rec_option);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}


int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetPoincareInput();
    if (argc != 2) {
      std::cerr << "POINCARE_Initial_Computation [FileNML]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      throw TerminalException{1};
    }
    using Tidx = int32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    RecOption rec_option = ReadInitialData(eFull);
    Process_rec_option<Tgroup>(rec_option);
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POINCARE_Initial_Computation time=" << time << "\n";
    exit(e.eVal);
  }
}
