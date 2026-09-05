// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "zero_one_solution.h"
// clang-format on

template <typename T>
void write_solutions(std::vector<Face> const &ListSol, int n_col,
                     std::string const &OutFormat, std::ostream &os) {
  if (OutFormat == "ZeroOne") {
    // One line of n_col characters 0 or 1 per solution
    for (auto &sol : ListSol) {
      for (int j = 0; j < n_col; j++)
        os << static_cast<int>(sol[j]);
      os << "\n";
    }
    return;
  }
  if (OutFormat == "matrix") {
    MyMatrix<T> M = ZeroOneSolutionsAsMatrix<T>(ListSol, n_col);
    WriteMatrix(os, M);
    return;
  }
  if (OutFormat == "GAP") {
    os << "return [";
    bool is_first = true;
    for (auto &sol : ListSol) {
      if (!is_first)
        os << ",\n";
      is_first = false;
      os << "[";
      for (int j = 0; j < n_col; j++) {
        if (j > 0)
          os << ",";
        os << static_cast<int>(sol[j]);
      }
      os << "]";
    }
    os << "];\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat=" << OutFormat << "\n";
  std::cerr << "Allowed values are ZeroOne, matrix and GAP\n";
  throw TerminalException{1};
}

template <typename T>
void process_zero_one_solutions(std::string const &FileMatrix,
                                std::string const &FileRHS,
                                std::string const &OutFormat,
                                std::string const &FileOut,
                                ZeroOneOptions const &options) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileMatrix);
  MyVector<T> b = ReadVectorFile<T>(FileRHS);
  if (A.rows() != b.size()) {
    std::cerr << "The matrix has " << A.rows() << " rows but the right hand "
              << "side has " << b.size() << " entries\n";
    throw TerminalException{1};
  }
  int n_col = A.cols();
  std::cerr << "n_row=" << A.rows() << " n_col=" << n_col << "\n";
  std::pair<ZeroOneResult, std::vector<Face>> pair =
      GetAllZeroOneSolutions(A, b, options, std::cerr);
  ZeroOneResult const &result = pair.first;
  std::cerr << "n_node=" << result.n_node
            << " n_solution=" << result.n_solution << "\n";
  if (!result.resolved) {
    std::cerr << "The enumeration is UNRESOLVED: the node budget of "
              << options.max_node << " was exhausted\n";
  } else {
    std::cerr << "The enumeration is complete\n";
  }
  if (FileOut == "stderr") {
    write_solutions<T>(pair.second, n_col, OutFormat, std::cerr);
  } else {
    if (FileOut == "stdout") {
      write_solutions<T>(pair.second, n_col, OutFormat, std::cout);
    } else {
      std::ofstream os(FileOut);
      write_solutions<T>(pair.second, n_col, OutFormat, os);
    }
  }
}

void process_arithmetic(std::string const &arithmetic,
                        std::string const &FileMatrix,
                        std::string const &FileRHS,
                        std::string const &OutFormat,
                        std::string const &FileOut,
                        ZeroOneOptions const &options) {
  if (arithmetic == "gmp") {
    using T = mpz_class;
    return process_zero_one_solutions<T>(FileMatrix, FileRHS, OutFormat,
                                         FileOut, options);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpz_int;
    return process_zero_one_solutions<T>(FileMatrix, FileRHS, OutFormat,
                                         FileOut, options);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_int;
    return process_zero_one_solutions<T>(FileMatrix, FileRHS, OutFormat,
                                         FileOut, options);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc < 6 || argc > 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ZeroOneSolutions [arithmetic] [FileMatrix] [FileRHS] "
                   "[OutFormat] [FileOut] <max_node> <lp_max_depth>\n";
      std::cerr << "\n";
      std::cerr << "It enumerates all the x in {0,1}^n with A x = b\n";
      std::cerr << "\n";
      std::cerr << "arithmetic  : gmp, gmp_boost, multi_boost\n";
      std::cerr << "FileMatrix  : the matrix A, in the format of "
                   "ReadMatrixFile\n";
      std::cerr << "FileRHS     : the vector b, in the format of "
                   "ReadVectorFile\n";
      std::cerr << "OutFormat   : ZeroOne, matrix or GAP\n";
      std::cerr << "FileOut     : the output file, or stdout or stderr\n";
      std::cerr << "max_node    : the node budget of the search, default is "
                   "unbounded\n";
      std::cerr << "lp_max_depth: the depth down to which the linear "
                   "programming\n";
      std::cerr << "              relaxation is used for pruning, default 0 "
                   "(the root\n";
      std::cerr << "              only), a negative value disabling it\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string FileMatrix = argv[2];
    std::string FileRHS = argv[3];
    std::string OutFormat = argv[4];
    std::string FileOut = argv[5];
    ZeroOneOptions options;
    if (argc >= 7) {
      std::string max_node_str = argv[6];
      options.max_node = ParseScalar<size_t>(max_node_str);
    }
    if (argc >= 8) {
      std::string lp_max_depth_str = argv[7];
      options.lp_max_depth = ParseScalar<int>(lp_max_depth_str);
    }
    process_arithmetic(arithmetic, FileMatrix, FileRHS, OutFormat, FileOut,
                       options);
    std::cerr << "Normal termination of LATT_ZeroOneSolutions\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ZeroOneSolutions\n";
    exit(e.eVal);
  }
  runtime(time);
}
