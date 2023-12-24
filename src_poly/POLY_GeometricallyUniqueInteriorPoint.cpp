// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "POLY_LinearProgramming.h"
// clang-format on

template <typename T>
void process(std::string const &eFile, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(eFile);
  //
  MyVector<T> eSol = GetGeometricallyUniqueInteriorPoint(FAC);
  if (OutFormat == "GAP") {
    os << "return ";
    WriteVectorGAP(os, eSol);
    os << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    return WriteVector(os, eSol);
  }
  std::cerr << "Error in process, missing support for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_GeometricallyUniqueInteriorPoint arith [DATAFAC] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "   or\n";
      std::cerr << "POLY_GeometricallyUniqueInteriorPoint arith [DATAFAC]\n";
      std::cerr << "\n";
      std::cerr << "with\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      std::cerr << "safe_rational          : rational arithmetic based on "
                   "int64_t that fails\n";
      std::cerr << "    gracefully on overflowing\n";
      std::cerr << "cpp_rational           : rational arithmetic based on "
                   "boost header library\n";
      std::cerr << "mpq_rational           : rational arithmetic based on "
                   "boost mpq data type\n";
      std::cerr << "rational               : rational arithmetic on input\n";
      std::cerr << "Qsqrt2                 : arithmetic over the Q(sqrt(2))\n";
      std::cerr << "Qsqrt5                 : arithmetic over the Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc : For the real algebraic case\n";
      std::cerr << "    of a field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "     ------- DATAFAC -------\n";
      std::cerr << "\n";
      std::cerr << "The file containing the matrix of the inequalities\n";
      std::cerr << "\n";
      std::cerr << "     ------- OutFormat -------\n";
      std::cerr << "\n";
      std::cerr << "GAP : The GAP file\n";
      std::cerr << "CPP : The CPP polyhedral format\n";
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
    std::cerr << "Reading input\n";
    std::string arith = argv[1];
    std::string FileFAC = argv[2];
    std::string OutFormat = "CPP";
    std::string FileSPA = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileSPA = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileFAC, OutFormat, os);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileFAC, OutFormat, os);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileFAC, OutFormat, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(FileFAC, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(FileFAC, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(FileFAC, OutFormat, os);
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
        return process<T>(FileFAC, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    if (FileSPA == "stderr") {
      f(std::cerr);
    } else {
      if (FileSPA == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(FileSPA);
        f(os);
      }
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_LinearDetermineByInequalities\n";
    exit(e.eVal);
  }
  runtime(time1);
}
