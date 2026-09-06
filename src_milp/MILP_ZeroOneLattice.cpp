// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "zero_one_lattice.h"
// clang-format on

template <typename T, typename Tint>
void process_zero_one_lattice(std::string const &FileMatrix,
                              std::string const &FileRHS,
                              std::string const &choice,
                              std::string const &FileOut) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileMatrix);
  MyVector<T> b = ReadVectorFile<T>(FileRHS);
  if (A.rows() != b.size()) {
    std::cerr << "The matrix has " << A.rows() << " rows but the right hand "
              << "side has " << b.size() << " entries\n";
    throw TerminalException{1};
  }
  std::cerr << "n_row=" << A.rows() << " n_col=" << A.cols() << "\n";
  ZeroOneLattice<T, Tint> lattice = BuildZeroOneLattice<T, Tint>(A, b, std::cerr);
  std::cerr << "dim=" << lattice.dim << " bound=" << lattice.bound << "\n";
  ZeroOneLatticeEstimate est = EstimateEnumerationCost(lattice);
  std::cerr << "gs_min=" << est.gs_min << " gs_max=" << est.gs_max
            << " log10_det=" << est.log10_det << "\n";
  std::cerr << "log10_max_node=" << est.log10_max_node
            << " at level " << est.level << " of " << lattice.dim << "\n";
  std::cerr << "log10_max_node_flat=" << est.log10_max_node_flat
            << " at level " << est.level_flat << "\n";
  if (choice == "profile") {
    auto f_print = [&](std::ostream &os) -> void {
      os << "return rec(dim:=" << lattice.dim << ", gs_min:=" << est.gs_min
         << ", gs_max:=" << est.gs_max
         << ", log10_max_node:=" << est.log10_max_node
         << ", log10_max_node_flat:=" << est.log10_max_node_flat << ");\n";
    };
    FILE_PrintStderrStdoutFile(FileOut, f_print);
    return;
  }
  if (choice == "enumerate") {
    std::vector<Face> ListSol =
        EnumerateZeroOneByLattice<T, Tint>(lattice, std::cerr);
    std::cerr << "n_solution=" << ListSol.size() << "\n";
    int n_col = A.cols();
    auto f_print = [&](std::ostream &os) -> void {
      for (auto &sol : ListSol) {
        for (int j = 0; j < n_col; j++)
          os << static_cast<int>(sol[j]);
        os << "\n";
      }
    };
    FILE_PrintStderrStdoutFile(FileOut, f_print);
    return;
  }
  std::cerr << "Failed to find a matching choice=" << choice << "\n";
  std::cerr << "Allowed values are profile and enumerate\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "MILP_ZeroOneLattice [FileMatrix] [FileRHS] [choice] "
                   "[FileOut]\n";
      std::cerr << "\n";
      std::cerr << "It builds the lattice ker([A | -b]) with the form\n";
      std::cerr << "  Q(z,s) = ||2z - s 1||^2 + s^2\n";
      std::cerr << "for which the 0/1 solutions of A x = b are the vectors "
                   "of norm\n";
      std::cerr << "n+1 with s = +-1, and either reports what enumerating it "
                   "would\n";
      std::cerr << "cost, or enumerates.\n";
      std::cerr << "\n";
      std::cerr << "FileMatrix : the matrix A, in the format of "
                   "ReadMatrixFile\n";
      std::cerr << "FileRHS    : the vector b, in the format of "
                   "ReadVectorFile\n";
      std::cerr << "choice     : profile or enumerate\n";
      std::cerr << "             profile stops after the reduction and "
                   "reports the\n";
      std::cerr << "             Gram-Schmidt profile and the estimated "
                   "number of\n";
      std::cerr << "             nodes of the widest level of the search "
                   "tree, both\n";
      std::cerr << "             for the basis at hand and for a perfectly "
                   "flat one,\n";
      std::cerr << "             which no reduction can beat\n";
      std::cerr << "             enumerate runs the enumeration, which is "
                   "only worth\n";
      std::cerr << "             trying when the profile says so\n";
      std::cerr << "FileOut    : the output file, or stdout or stderr\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    std::string FileMatrix = argv[1];
    std::string FileRHS = argv[2];
    std::string choice = argv[3];
    std::string FileOut = argv[4];
    process_zero_one_lattice<T, Tint>(FileMatrix, FileRHS, choice, FileOut);
    std::cerr << "Normal termination of MILP_ZeroOneLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in MILP_ZeroOneLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}
