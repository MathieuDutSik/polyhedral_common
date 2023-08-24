// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_LinearProgramming.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::ostream &os) {
  MyMatrix<T> ListVect = ReadMatrixFile<T>(eFileI);
  PosRelRes<T> prr = SearchPositiveRelationSimple_Direct(ListVect);
  if (prr.eTestExist) {
    os << "relation found\n";
    MyVector<T> const& V = *prr.TheRelat;
    WriteVector(os, V);
  } else {
    os << "no relation found\n";
  }
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IsPointedCone arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational : rational based on int64_t that fails gracefully of overflowing\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a "
                   "field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileI = argv[2];
    auto compute_pointedness = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileI, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, os);
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
        return process<T>(eFileI, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    compute_pointedness(std::cerr);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs\n";
    exit(e.eVal);
  }
  runtime(time1);
}
