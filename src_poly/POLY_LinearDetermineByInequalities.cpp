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
void process(std::string const &eFile, std::string const &OutFormat,
             std::ostream &os_out, std::ostream &os) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(eFile);
  //
  MyMatrix<T> LinSpace = LinearDeterminedByInequalities(FAC, os);
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteMatrixGAP(os_out, LinSpace);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    return WriteMatrix(os_out, LinSpace);
  }
  std::cerr << "Error in process, missing support for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string FileFAC;
    std::string OutFormat = "CPP";
    std::string FileSPA = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileFAC, OutFormat, *os_ptr, std::cerr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinearDetermineByInequalities arith [DATAFAC] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "   or\n";
      std::cerr << "POLY_LinearDetermineByInequalities arith [DATAFAC]\n";
      std::cerr << "\n";
      std::cerr << "with\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
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
    FileFAC = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileSPA = argv[4];
    }

    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileSPA, f);
    //
    std::cerr << "Normal termination of POLY_LinearDetermineByInequalities\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_LinearDetermineByInequalities\n";
    exit(e.eVal);
  }
  runtime(time1);
}
