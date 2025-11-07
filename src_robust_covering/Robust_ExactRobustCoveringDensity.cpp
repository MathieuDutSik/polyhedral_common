// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "Group.h"
#include "NumberTheory.h"
#include "enum_robust_covering.h"
// clang-format on

template <typename T, typename Tint>
void process_B(std::string const& MatFile, std::string const& OutFormat, std::string const& OutFile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  //
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
    AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  DataLattice<T, Tint, Tgroup> eData =
    GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);

  T sqr_dist = compute_square_robust_covering_radius(eData);

  T det = DeterminantMat(GramMat);
  int dim = GramMat.rows();
  ResultCov<T> rc = ComputeCoveringDensityFromDimDetCov(dim, det, sqr_dist);
  auto f_print=[&](std::ostream& osf) -> void {
    if (OutFormat == "GAP") {
      osf << "return ";
      osf << to_stringGAP(rc);
      osf << ";\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: GAP\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

void process_A(std::string const &arithmetic,
               std::string MatFile, std::string OutFormat, std::string OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T,Tint>(MatFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_ExactRobustCoveringDensity arithmetic [MaFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_ExactRobustCoveringDensity arithmetic [MaFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string MatFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    process_A(arithmetic, MatFile, OutFormat, OutFile);
    std::cerr << "Normal termination of Robust_ExactRobustCoveringDensity\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_ExactRobustCoveringDensity\n";
    exit(e.eVal);
  }
  runtime(time);
}
