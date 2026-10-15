// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
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
// clang-format on

template <typename T>
void process(std::string const &FileEXT1, std::string const &FileEXT2,
             std::string const &OutFormat, std::ostream &os_out,
             std::ostream &os) {
  MyMatrix<T> EXT1 = ReadMatrixFile<T>(FileEXT1);
  MyMatrix<T> EXT2 = ReadMatrixFile<T>(FileEXT2);
  bool test = has_empty_intersection(EXT1, EXT2, os);
  if (OutFormat == "GAP") {
    os_out << "return " << GAP_logical(test) << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    os_out << test << "\n";
    return;
  }
  std::cerr << "Error in process, missing support for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_check_empty_intersection arith [EXT1] [EXT2] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "   or\n";
      std::cerr << "POLY_check_empty_intersection arith [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "Whether the convex hulls of the rows of EXT1 and of EXT2 "
                   "(points with first\n";
      std::cerr << "coordinate 1) are disjoint.\n";
      std::cerr << "\n";
      std::cerr << "with\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "cpp_rational           : rational arithmetic based on "
                   "boost header library\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "mpq_rational           : rational arithmetic based on "
                   "boost mpq data type\n";
#endif
      std::cerr << "mpq_class              : rational arithmetic over GMP "
                   "mpq_class\n";
#ifdef QUADRATIC_FIELDS
      std::cerr
          << "Qsqrt2                 : arithmetic over the field Q(sqrt(2))\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr
          << "Qsqrt5                 : arithmetic over the field Q(sqrt(5))\n";
#endif
      std::cerr
          << "RealAlgebraic=FileDesc : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "     ------- EXT1, EXT2 -------\n";
      std::cerr << "\n";
      std::cerr << "The files containing the two lists of points\n";
      std::cerr << "\n";
      std::cerr << "     ------- OutFormat -------\n";
      std::cerr << "\n";
      std::cerr << "GAP : return true; or return false;\n";
      std::cerr << "CPP : 1 or 0\n";
      std::cerr << "\n";
      std::cerr << "     ------- OutFile -------\n";
      std::cerr << "\n";
      std::cerr << "stderr : output to std::cerr\n";
      std::cerr << "stdout : output to std::cout\n";
      std::cerr << "filename : output to filename (if different from stderr / "
                   "stdout)\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileEXT1 = argv[2];
    std::string FileEXT2 = argv[3];
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }

    auto f = [&](std::ostream &os_out) -> void {
#ifdef ENABLE_BOOST_TYPES
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "fmpq_class") {
        using T = fmpq_class;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
#endif
#ifdef ENABLE_BOOST_TYPES
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
#endif
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
#ifdef QUADRATIC_FIELDS
      if (arith == "Qsqrt5") {
        using Trat = Trat_quad_field;
        using T = QuadField<Trat, 5>;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
      if (arith == "Qsqrt2") {
        using Trat = Trat_quad_field;
        using T = QuadField<Trat, 2>;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
#endif
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!FILE_IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat_real_field> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(FileEXT1, FileEXT2, OutFormat, os_out,
                          std::cerr);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    //
    std::cerr << "Normal termination of POLY_check_empty_intersection\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_check_empty_intersection\n";
    exit(e.eVal);
  }
  runtime(time1);
}
