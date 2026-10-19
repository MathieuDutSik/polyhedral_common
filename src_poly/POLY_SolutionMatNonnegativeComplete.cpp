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
void process(std::string const &eFileFAC, std::string const &eFileINEQ,
             std::ostream &os) {
  MyMatrix<T> ListVect = ReadMatrixFile<T>(eFileFAC);
  MyVector<T> eVect = ReadVectorFile<T>(eFileINEQ);
  SolutionMatNonnegativeComplete<T> eSol =
      GetSolutionMatNonnegativeComplete(ListVect, eVect, std::cerr);
  if (eSol.ExtremeRay) {
    os << "eEXT=" << StringVector(*eSol.ExtremeRay) << "\n";
  } else {
    os << "No extreme ray found\n";
  }
  //
  if (eSol.SolNonnegative) {
    os << "eSol=" << StringVector(*eSol.SolNonnegative) << "\n";
  } else {
    os << "No nonnegative solution found\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string eFileFAC;
    std::string eFileINEQ;
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileFAC, eFileINEQ, *os_ptr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_SolutionMatNonnegativeComplete arith [FileFAC] [FileINEQ]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    eFileFAC = argv[2];
    eFileINEQ = argv[3];
    auto compute_solution = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    compute_solution(std::cerr);
    std::cerr << "Normal termination of POLY_SolutionMatNonnegativeComplete\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SolutionMatNonnegativeComplete\n";
    exit(e.eVal);
  }
  runtime(time1);
}
