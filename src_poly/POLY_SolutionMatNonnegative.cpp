// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_LinearProgramming.h"
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
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SolutionMatNonnegative Arith [FAC] [INEQ]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully if overflowing\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "FAC : The list of defining inequalities\n";
      std::cerr << "INEQ: A single inequality\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileFAC = argv[2];
    std::string FileINEQ = argv[3];
    auto evaluation = [&]() -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return full_process_type<T>(FileFAC, FileINEQ);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return full_process_type<T>(FileFAC, FileINEQ);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return full_process_type<T>(FileFAC, FileINEQ);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return full_process_type<T>(FileFAC, FileINEQ);
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
        return full_process_type<T>(FileFAC, FileINEQ);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    evaluation();
    //
    std::cerr << "Normal termination of POLY_SolutionMatNonnegative\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SolutionMatNonnegative\n";
    exit(e.eVal);
  }
  runtime(time1);
}
