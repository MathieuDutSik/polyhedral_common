// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_LinearProgramming.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::ostream &os_out,
             std::ostream &os) {
  MyMatrix<T> ListVect = ReadMatrixFile<T>(eFileI);
  PosRelRes<T> prr = SearchPositiveRelationSimple_Direct(ListVect, os);
  if (prr.eTestExist) {
    os_out << "relation found\n";
    MyVector<T> const &V = *prr.TheRelat;
    WriteVector(os_out, V);
  } else {
    os_out << "no relation found\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string eFileI;
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileI, *os_ptr, std::cerr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IsPointedCone arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    eFileI = argv[2];
    auto compute_pointedness = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    compute_pointedness(std::cerr);
    std::cerr << "Normal termination of POLY_IsPointedCone\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_IsPointedCone\n";
    exit(e.eVal);
  }
  runtime(time1);
}
