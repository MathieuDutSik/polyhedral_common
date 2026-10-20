// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "poincare_polyhedron.h"
#include "Group.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename Tgroup>
void Process_rec_option(RecOption const &rec_option, std::ostream &os) {
  std::string arith = rec_option.Arithmetic;
  ArithmeticDispatch dispatch;
  auto f_arith = [&]<typename T>() -> void {
    full_process_type<T, Tgroup>(rec_option, os);
  };
  add_rational_arithmetics(dispatch, f_arith);
  add_algebraic_field_arithmetics(dispatch, f_arith);
  dispatch.run(arith);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetPoincareInput();
    if (argc != 2) {
      std::cerr << "POINCARE_Initial_Computation [FileNML]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      throw TerminalException{1};
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    RecOption rec_option = ReadInitialData(eFull);
    Process_rec_option<Tgroup>(rec_option, std::cerr);
    std::cerr << "Normal termination of POINCARE_Initial_Computation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POINCARE_Initial_Computation\n";
    exit(e.eVal);
  }
}
