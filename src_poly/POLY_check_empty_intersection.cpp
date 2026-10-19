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
#include "Arithmetic_dispatch_types.h"
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
    std::string FileEXT1;
    std::string FileEXT2;
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileEXT1, FileEXT2, OutFormat, *os_ptr,
                          std::cerr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
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
      dispatch.print_descriptions(std::cerr);
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
    FileEXT1 = argv[2];
    FileEXT2 = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }

    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
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
