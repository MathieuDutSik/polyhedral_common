// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "system_symmetry.h"
// clang-format on

template <typename T, typename Tgroup>
void process_system_symmetry(std::string const &FileMatrix,
                             std::string const &FileRHS,
                             std::string const &FileColColour,
                             std::string const &OutFormat,
                             std::string const &FileOut) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileMatrix);
  MyVector<T> b = ReadVectorFile<T>(FileRHS);
  if (A.rows() != b.size()) {
    std::cerr << "The matrix has " << A.rows() << " rows but the right hand "
              << "side has " << b.size() << " entries\n";
    throw TerminalException{1};
  }
  std::vector<T> ColColour;
  if (FileColColour != "none") {
    MyVector<T> V = ReadVectorFile<T>(FileColColour);
    if (V.size() != A.cols()) {
      std::cerr << "The colour vector has " << V.size() << " entries but the "
                << "matrix has " << A.cols() << " columns\n";
      throw TerminalException{1};
    }
    for (int j = 0; j < V.size(); j++)
      ColColour.push_back(V(j));
  }
  std::cerr << "n_row=" << A.rows() << " n_col=" << A.cols() << "\n";
  Tgroup GRP = ComputeSystemSymmetry<T, Tgroup>(A, b, ColColour, std::cerr);
  std::cerr << "|G_mat|=" << GRP.size() << "\n";
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
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "MILP_SystemSymmetry [FileMatrix] [FileRHS] "
                   "[FileColColour] [OutFormat] [FileOut]\n";
      std::cerr << "        or\n";
      std::cerr << "MILP_SystemSymmetry [FileMatrix] [FileRHS] "
                   "[FileColColour]\n";
      std::cerr << "\n";
      std::cerr << "It computes G_mat, the group of the permutations of the "
                   "columns of A\n";
      std::cerr << "for which some permutation of the rows leaves A and b "
                   "unchanged.\n";
      std::cerr << "\n";
      std::cerr << "FileMatrix    : the matrix A, in the format of "
                   "ReadMatrixFile\n";
      std::cerr << "FileRHS       : the vector b, in the format of "
                   "ReadVectorFile\n";
      std::cerr << "FileColColour : a colour per column, in the format of\n";
      std::cerr << "                ReadVectorFile, or none. Colouring the "
                   "columns\n";
      std::cerr << "                asks for the subgroup preserving that "
                   "colouring,\n";
      std::cerr << "                a partial assignment for instance\n";
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
    std::string FileColColour = argv[3];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    process_system_symmetry<T, Tgroup>(FileMatrix, FileRHS, FileColColour,
                                       OutFormat, FileOut);
    std::cerr << "Normal termination of MILP_SystemSymmetry\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in MILP_SystemSymmetry\n";
    exit(e.eVal);
  }
  runtime(time);
}
