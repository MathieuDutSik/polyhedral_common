// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "lorentzian_perfect.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename Tint>
void WriteGenerators(std::vector<MyMatrix<Tint>> const &l_gen,
                     std::string const &OutFormat, std::ostream &os_out) {
  if (OutFormat == "CPP") {
    return WriteListMatrix(os_out, l_gen);
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, l_gen);
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching entry in WriteGenerators\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
void process_C(std::string const &FileMatrix, std::string const &OutFormat,
               std::string const &FileOut) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;

  MyMatrix<T> LorMat = ReadMatrixFile<T>(FileMatrix);
  std::vector<MyMatrix<Tint>> l_gen =
      LORENTZ_GetGeneratorsAutom<T, Tint, Tgroup>(LorMat, std::cerr);
  auto f = [&](std::ostream &os_out) -> void {
    WriteGenerators(l_gen, OutFormat, os_out);
  };
  print_stderr_stdout_file(FileOut, f);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    Eigen::initParallel();
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "LORENTZ_PERF_Automorphism [fileMatrix] [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "LORENTZ_PERF_Automorphism [fileMatrix]\n";
      std::cerr << "With file.nml a namelist file\n";
      return -1;
    }
    //
    std::string FileMatrix = argv[1];
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileOut = argv[3];
    }
    using T = mpq_class;
    using Tint = mpz_class;
    process_C<T, Tint>(FileMatrix, OutFormat, FileOut);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_PERF_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
