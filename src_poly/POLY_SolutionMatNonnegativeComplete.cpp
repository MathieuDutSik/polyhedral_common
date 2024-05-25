// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_LinearProgramming.h"
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
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SolutionMatNonnegativeComplete arith [FileFAC] [FileINEQ]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflowing\n";
      std::cerr << "rational : rational arithmetic on gmp\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileFAC = argv[2];
    std::string eFileINEQ = argv[3];
    auto compute_solution = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileFAC, eFileINEQ, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileFAC, eFileINEQ, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileFAC, eFileINEQ, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileFAC, eFileINEQ, os);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        using T_rat = mpq_class;
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFileFAC, eFileINEQ, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    compute_solution(std::cerr);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SolutionMatNonnegativeComplete\n";
    exit(e.eVal);
  }
  runtime(time1);
}
