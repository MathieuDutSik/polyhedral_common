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
void full_process_type(std::string const &FileFAC,
                       std::string const &FileINEQ) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
  std::cerr << "We have FAC, |FAC|=" << FAC.rows() << " / " << FAC.cols()
            << "\n";
  MyVector<T> Ineq = ReadVectorFile<T>(FileINEQ);
  //
  std::optional<MyVector<T>> opt = SolutionMatNonnegative(FAC, Ineq, std::cerr);
  //
  if (opt) {
    MyVector<T> V = *opt;
    std::cerr << "Found a solution V=" << StringVector(V) << "\n";
  } else {
    std::cerr << "No solution found\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileFAC, FileINEQ;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      full_process_type<T>(FileFAC, FileINEQ);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SolutionMatNonnegative Arith [FAC] [INEQ]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "\n";
      std::cerr << "FAC : The list of defining inequalities\n";
      std::cerr << "INEQ: A single inequality\n";
      return -1;
    }
    std::string arith = argv[1];
    FileFAC = argv[2];
    FileINEQ = argv[3];
    dispatch.run(arith);
    //
    std::cerr << "Normal termination of POLY_SolutionMatNonnegative\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SolutionMatNonnegative\n";
    exit(e.eVal);
  }
  runtime(time1);
}
