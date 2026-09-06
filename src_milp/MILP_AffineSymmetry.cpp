// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "affine_symmetry.h"
// clang-format on

template <typename T, typename Tgroup>
void process_affine_symmetry(std::string const &FileMatrix,
                             std::string const &FileRHS,
                             std::string const &OutFormat,
                             std::string const &FileOut) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileMatrix);
  MyVector<T> b = ReadVectorFile<T>(FileRHS);
  if (A.rows() != b.size()) {
    std::cerr << "The matrix has " << A.rows() << " rows but the right hand "
              << "side has " << b.size() << " entries\n";
    throw TerminalException{1};
  }
  std::cerr << "n_row=" << A.rows() << " n_col=" << A.cols() << "\n";
  if (!SolutionMat(TransposedMat(A), b)) {
    std::cerr << "The system has no solution at all, its affine subspace is "
              << "empty and its stabilizer is meaningless\n";
    throw TerminalException{1};
  }
  Tgroup GRP = ComputeAffineSymmetry<T, Tgroup>(A, b, std::cerr);
  std::cerr << "|G_aff|=" << GRP.size() << "\n";
  std::cerr << "n_generator=" << GRP.GeneratorsOfGroup().size() << "\n";
  auto f_print = [&](std::ostream &os) -> void {
    if (OutFormat == "GAP") {
      os << "return " << GRP.GapString() << ";\n";
      return;
    }
    if (OutFormat == "order") {
      os << GRP.size() << "\n";
      return;
    }
    std::cerr << "Failed to find a matching OutFormat=" << OutFormat << "\n";
    std::cerr << "Allowed values are GAP and order\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(FileOut, f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "MILP_AffineSymmetry [FileMatrix] [FileRHS] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "        or\n";
      std::cerr << "MILP_AffineSymmetry [FileMatrix] [FileRHS]\n";
      std::cerr << "\n";
      std::cerr << "It computes G_aff, the group of the permutations of the "
                   "coordinates\n";
      std::cerr << "that preserve the affine subspace {x : A x = b}. It "
                   "contains G_mat.\n";
      std::cerr << "\n";
      std::cerr << "FileMatrix    : the matrix A, in the format of "
                   "ReadMatrixFile\n";
      std::cerr << "FileRHS       : the vector b, in the format of "
                   "ReadVectorFile\n";
      std::cerr << "OutFormat     : GAP or order\n";
      std::cerr << "FileOut       : the output file, or stdout or stderr\n";
      return -1;
    }
    using T = mpq_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string FileMatrix = argv[1];
    std::string FileRHS = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    process_affine_symmetry<T, Tgroup>(FileMatrix, FileRHS, OutFormat,
                                       FileOut);
    std::cerr << "Normal termination of MILP_AffineSymmetry\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in MILP_AffineSymmetry\n";
    exit(e.eVal);
  }
  runtime(time);
}
