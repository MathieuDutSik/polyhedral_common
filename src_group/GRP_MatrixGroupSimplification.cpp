// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "MatrixGroupSimplification.h"
// clang-format on

template <typename T>
void process(std::string const &FileMatrGroup, std::string const &OutFormat,
             std::ostream &os_out) {
  std::vector<MyMatrix<T>> ListM = ReadListMatrixFile<T>(FileMatrGroup);
  std::vector<MyMatrix<T>> ListMred =
      ExhaustiveReductionComplexityGroupMatrix<T>(ListM, std::cerr);
  //
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, ListMred);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    WriteListMatrix(os_out, ListMred);
    return;
  }
  if (OutFormat == "Stat") {
    os_out << "INPUT : " << compute_complexity_listmat(ListM) << "\n";
    os_out << "OUTPUT: " << compute_complexity_listmat(ListMred) << "\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_MatrixGroupSimplification [Arith] [FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_MatrixGroupSimplification [Arith] [FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "Arith         : mpq_class / mpz_class / mpq_rational / cpp_rational / mpz_int / cpp_int\n";
      std::cerr << "FileMatrGroup : The file containing the list of matrices\n";
      std::cerr << "OutFormat     : Optional parameter GAP / CPP\n";
      std::cerr
          << "FileOut       : File where to write the simplified matrices\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileMatrGroup = argv[2];
    std::string OutFormat = "Stat";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "mpz_class") {
        using T = mpz_class;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "safe_integer") {
        using T = SafeInt64;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "mpz_int") {
        using T = boost::multiprecision::mpz_int;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "cpp_int") {
        using T = boost::multiprecision::cpp_int;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      std::cerr << "Failed to find a matching arith. Allowed is mpq_class / "
                   "mpz_class / mpq_rational / cpp_rational / mpz_int / cpp_int\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of GRP_MatrixGroupSimplification\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_MatrixGroupSimplification\n";
    exit(e.eVal);
  }
  runtime(time);
}
