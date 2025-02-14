// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "POLY_lrslib.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  T volume = lrs::Kernel_VolumePolytope(EXT);
  if (OutFormat == "Direct") {
    os << "volume=" << volume << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return " << volume << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching output\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_lrs_volume arith [DATAIN] [OutFormat] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "POLY_lrs_volume arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith   : the chosen arithmetic\n";
      std::cerr << "DATAIN  : The polytope vertices\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr
          << "safe_integer  : integer arithmetic based on int64_t that fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "integer  : integer arithmetic on input\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a ";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileI = argv[2];
    std::string OutFormat = "Direct";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    auto call_lrs = [&](std::ostream &os) -> void {
      if (arith == "safe_integer") {
        using T = SafeInt64;
        return process<T>(eFileI, OutFormat, os);
      }
      /*
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, OutFormat, os);
      }
      */
      if (arith == "integer") {
        using T = mpz_class;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "double") {
        using T = double;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, OutFormat, os);
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
        return process<T>(eFileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      call_lrs(std::cerr);
    } else {
      if (FileO == "stdout") {
        call_lrs(std::cout);
      } else {
        std::ofstream os(FileO);
        call_lrs(os);
      }
    }
    std::cerr << "Normal termination of POLY_lrs_volume\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs_volume\n";
    exit(e.eVal);
  }
  runtime(time1);
}
