// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "enum_robust_covering.h"
// clang-format on

template <typename T, typename Tint>
void process_B(size_t const& n_iter, std::string const& MatFile, std::string const& OutFormat, std::string const& OutFile) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);

  T sqr_dist = random_estimation_robust_covering<T,Tint>(GramMat, n_iter, std::cerr);
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
               size_t const& n_iter,
               std::string MatFile, std::string OutFormat, std::string OutFile) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T,Tint>(n_iter, MatFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Robust_RandomEstimation arithmetic n_iter [MaFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "Robust_RandomEstimation arithmetic n_iter [MaFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic: gmp\n";
      std::cerr << "n_iter: 1, 10, or whatever you want\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    size_t n_iter = ParseScalar<size_t>(argv[2]);
    std::string MatFile = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    process_A(arithmetic, n_iter, MatFile, OutFormat, OutFile);
    std::cerr << "Normal termination of Robust_RandomEstimation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Robust_RandomEstimation\n";
    exit(e.eVal);
  }
  runtime(time);
}
