// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "lorentzian_perfect.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template<typename Tint>
void WriteEquivalence(std::optional<MyMatrix<Tint>> const& opt, std::string const& OutFormat, std::ostream & os_out) {
  if (OutFormat == "CPP") {
    if (opt) {
      WriteMatrix(os_out, *opt);
    } else {
      os_out << "-1\n";
    }
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (opt) {
      WriteMatrixGAP(os_out, *opt);
    } else {
      os_out << "fail";
    }
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

template<typename T, typename Tint>
void process_C(std::string const& FileMatrix1, std::string const& FileMatrix2, std::string const& OutFormat, std::string const& FileOut) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;

  MyMatrix<T> LorMat1 = ReadMatrixFile<T>(FileMatrix1);
  MyMatrix<T> LorMat2 = ReadMatrixFile<T>(FileMatrix2);
  std::optional<MyMatrix<Tint>> opt = LORENTZ_TestEquivalenceMatrices<T,Tint,Tgroup>(LorMat1, LorMat2, std::cerr);
  auto f=[&](std::ostream& os_out) -> void {
    WriteEquivalence(opt, OutFormat, os_out);
  };
  print_stderr_stdout_file(FileOut, f);
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    Eigen::initParallel();
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_PERF_Isomorphism [fileMatrix1] [FileMatrix2] [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "LORENTZ_PERF_Isomorphism [fileMatrix1] [FileMatrix2]\n";
      std::cerr << "With file.nml a namelist file\n";
      return -1;
    }
    //
    std::string FileMatrix1 = argv[1];
    std::string FileMatrix2 = argv[2];
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    using T = mpq_class;
    using Tint = mpz_class;
    process_C<T,Tint>(FileMatrix1, FileMatrix2, OutFormat, FileOut);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_PERF_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
