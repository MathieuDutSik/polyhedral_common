// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_General.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T>
void write_result(std::optional<MyMatrix<T>> const &opt,
                  std::string const &OutFormat, std::ostream &os) {
  if (OutFormat == "simple") {
    if (opt) {
      os << "It is equivalence\n";
    } else {
      os << "It is not equivalence\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    if (opt) {
      os << "return " << StringMatrixGAP_line(*opt) << ";\n";
    } else {
      os << "return fail;\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    FullNamelist eFull = NAMELIST_GetOneTSPACE();
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Equivalence [FileLinSpa] [FileMat1] [FileMat2]\n";
      std::cerr << "or\n";
      std::cerr << "TSPACE_Equivalence [FileLinSpa] [FileMat1] [FileMat2] "
                   "[OutFormat] [FileOut]\n";
      return -1;
    }
    std::string FileTspace = argv[1];
    std::string FileGram1 = argv[2];
    std::string FileGram2 = argv[3];
    std::string OutFormat = "simple";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    //
    LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace, std::cerr);
    MyMatrix<T> eMat1 = ReadMatrixFile<T>(FileGram1);
    MyMatrix<T> eMat2 = ReadMatrixFile<T>(FileGram2);
    std::optional<MyMatrix<Tint>> opt =
        LINSPA_TestEquivalenceGramMatrix<T, Tint, Tgroup>(LinSpa, eMat1, eMat2,
                                                          std::cerr);
    auto f = [&](std::ostream &os_out) -> void {
      write_result(opt, OutFormat, os_out);
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of TSPACE_Equivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Equivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}
