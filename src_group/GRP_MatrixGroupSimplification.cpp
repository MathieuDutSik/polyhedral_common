// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "MatrixGroupSimplification.h"
#include "Arithmetic_dispatch_types.h"
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
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileMatrGroup;
    std::string OutFormat = "Stat";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileMatrGroup, OutFormat, *os_ptr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_integer_arithmetics(dispatch, f_arith);
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_MatrixGroupSimplification [Arith] [FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_MatrixGroupSimplification [Arith] [FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "FileMatrGroup : The file containing the list of matrices\n";
      std::cerr << "OutFormat     : Optional parameter GAP / CPP / Stat\n";
      std::cerr
          << "FileOut       : File where to write the simplified matrices\n";
      std::cerr << "\n";
      std::cerr << "     ------- Arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileMatrGroup = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of GRP_MatrixGroupSimplification\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_MatrixGroupSimplification\n";
    exit(e.eVal);
  }
  runtime(time);
}
